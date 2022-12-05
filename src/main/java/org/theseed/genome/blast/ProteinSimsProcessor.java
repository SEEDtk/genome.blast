package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.ProteinCompare;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;
import org.theseed.utils.StringPair;

/**
 * This command compares all the proteins in a set of protein FASTA files in a directory.  All the proteins
 * in each file are compared with all the proteins in each other file.  The FASTA files can be meaningful or
 * they can simply be a giant list divided into chunks for batching performance.
 *
 * The positional parameter is the name of the source directory for the FASTA files, which should have
 * file name extensions of ".fa", ".faa", or ".fasta".
 *
 * It is presumed this will be run with the Diamond protein comparison tool installed and comparison method
 * BLAST.  It will run correctly with NCBI BLAST as the underlying tool, but it is optimized for Diamond,
 * which is heavily parallelized.  As a result, this program is single-threaded, though design decisions have been
 * made to allow for parallelization in the future (in case different methods or tools become prominent).
 *
 * Extra information is written to the report to allow for resume processing.  Each output line contains a file
 * pair name, a sequence number, and a total number of protein pairs for the file pairing.  When we resume, we
 * load all the protein pairs for a file pairing, and then if we found all of them, we write them out and denote
 * the file pairing is already processed.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 * -K	kmer size for kmer-based comparisons (default 8)
 *
 * --method		comparison method to use (default BLAST)
 * --resume		if specified, a file from a previous run; this can be used to restart a failed job
 *
 * @author Bruce Parrello
 *
 */
public class ProteinSimsProcessor extends BaseReportProcessor implements ProteinCompare.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinSimsProcessor.class);
    /** comparison method engine */
    private ProteinCompare method;
    /** set of file-pair requests to process */
    private Set<Request> filePairs;
    /** counter of comparisons output */
    private int batchCounter;
    /** counter of comparisons required */
    private int batchTotal;
    /** file name filter for protein FASTA files */
    private static FilenameFilter FASTA_FILTER = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            String namelc = name.toLowerCase();
            return namelc.endsWith(".faa") || namelc.endsWith(".fa") || namelc.endsWith(".fasta");
        }
    };

    // COMMAND-LINE OPTIONS

    /** protein kmer size */
    @Option(name = "--kmer", aliases = { "-K" }, metaVar = "9", usage = "protein kmer size")
    private int kmerSize;

    /** comparison method */
    @Option(name = "--method", usage = "protein comparison method")
    private ProteinCompare.Type methodType;

    /** name of a file containing output from a prior, if we are resuming */
    @Option(name = "--resume", metaVar = "oldOutput.tbl", usage = "if specified, the name of a file from a prior run that was interrupted")
    private File oldOutputFile;

    /** input directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory containing FASTA files", required = true)
    private File inDir;

    /**
     * This object tracks a pair of FASTA files that need to be processed.  It is keyed by the name portions of the
     * files separated by a slash.
     */
    protected class Request implements Comparable<Request> {

        /** name of the first file */
        private File file1;
        /** name of the second file */
        private File file2;

        /**
         * Create a request to process a file pair.
         *
         * @param f1	query file
         * @param f2	subject file
         */
        protected Request(File f1, File f2) {
            this.file1 = f1;
            this.file2 = f2;
        }

        /**
         * Compute the file-pair request for a specified key string.  The
         * key string is produced by our "toString" method.
         *
         * @param key	key string for the request
         */
        public Request(String key) throws IOException {
            String[] pair = StringUtils.split(key, " / ");
            if (pair.length != 2)
                throw new IOException("Invalid file-pair descriptor \"" + key + "\".");
            this.file1 = new File(ProteinSimsProcessor.this.inDir, pair[0]);
            this.file2 = new File(ProteinSimsProcessor.this.inDir, pair[1]);
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.file1 == null) ? 0 : this.file1.hashCode());
            result = prime * result + ((this.file2 == null) ? 0 : this.file2.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Request)) {
                return false;
            }
            Request other = (Request) obj;
            if (this.file1 == null) {
                if (other.file1 != null) {
                    return false;
                }
            } else if (!this.file1.equals(other.file1)) {
                return false;
            }
            if (this.file2 == null) {
                if (other.file2 != null) {
                    return false;
                }
            } else if (!this.file2.equals(other.file2)) {
                return false;
            }
            return true;
        }

        @Override
        public String toString() {
            return this.file1.getName() + " / " + this.file2.getName();
        }

        @Override
        public int compareTo(Request o) {
            int retVal = this.file1.compareTo(o.file1);
            if (retVal == 0)
                retVal = this.file2.compareTo(o.file2);
            return retVal;
        }

        /**
         * @return the name of the subject file
         */
        protected File getFile1() {
            return this.file1;
        }

        /**
         * @return the name of the query file
         */
        protected File getFile2() {
            return this.file2;
        }

    }

    @Override
    public int getKmerSize() {
        return this.kmerSize;
    }

    @Override
    protected void setReporterDefaults() {
        this.kmerSize = 8;
        this.methodType = ProteinCompare.Type.BLAST;
        this.oldOutputFile = null;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        //Insure the input directory exists and has FASTA files.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        File[] inFiles = this.inDir.listFiles(FASTA_FILTER);
        log.info("{} files found in input directory.", inFiles.length);
        // Sort the files so that the left-right order doesn't change between runs.
        Arrays.sort(inFiles);
        this.filePairs = new HashSet<Request>(inFiles.length * inFiles.length * 2 / 3 + 1);
        // Form the list of pair requests.
        for (int i = 0; i < inFiles.length; i++) {
            for (int j = i; j < inFiles.length; j++) {
                this.filePairs.add(this.new Request(inFiles[i], inFiles[j]));
            }
        }
        log.info("{} file pairs need to be processed.", this.filePairs.size());
        // Validate the resume file.
        if (this.oldOutputFile != null && ! this.oldOutputFile.canRead())
            throw new FileNotFoundException("Prior output file " + this.oldOutputFile + " is not found or unreadable.");
        // Validate the kmer size.
        if (this.kmerSize <= 0)
            throw new ParseFailureException("Kmer size must be at least 1.");
        // Create the method.
        log.info("Initializing {} comparison engine.", this.methodType);
        this.method = this.methodType.create(this);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Prepare the output file.
        writer.println("prot1\tprot2\tsimilarity\tfile_pair\tfp_total");
        // Set up our progress counter.
        this.batchCounter = 0;
        this.batchTotal = this.filePairs.size();
        // Determine which file pairs have already been processed.
        if (this.oldOutputFile !=  null) try (var oldStream = new TabbedLineReader(this.oldOutputFile)) {
            log.info("Reading {} for prior-run results.", this.oldOutputFile);
            // Check the headers to make sure this is a good file.
            int p1Idx = oldStream.findField("prot1");
            int fpIdx = oldStream.findField("file_pair");
            int totIdx = oldStream.findField("fp_total");
            if (p1Idx != 0 || fpIdx != 3)
                throw new IOException(this.oldOutputFile + " does not appear to be in the proper format for a prior-run output file.");
            // Denote we have no current request.  We set the total to 1 so the null request appears incomplete and is not output.
            Request curr = null;
            String currName = "";
            int currTotal = 1;
            Map<StringPair, Double> sims = new HashMap<StringPair, Double>(1000);
            for (var line : oldStream) {
                String reqName = line.get(fpIdx);
                if (! reqName.contentEquals(currName)) {
                    // Here we have a new batch.  Check the old batch to see if it needs to be output.
                    this.checkResumeBatch(curr, writer, sims, currTotal);
                    // Set up for the next batch.
                    currName = reqName;
                    curr = new Request(currName);
                    currTotal = line.getInt(totIdx);
                    sims.clear();
                }
                // Add this record to the sims map.
                StringPair pair = new StringPair(line.get(0), line.get(1));
                double sim = line.getDouble(2);
                sims.put(pair, sim);
            }
            // Output the residual.  This is the one most likely to be incomplete.
            this.checkResumeBatch(curr, writer, sims, currTotal);
            log.info("{} file pairs remaining after resume processing.", this.filePairs.size());
        }
        // Loop through the file pairs.
        this.filePairs.stream().forEach(x -> this.runCompare(x, writer));
    }

    /**
     * If the specified resume batch is complete, write it to the output and remove its file pairing
     * from the request set.
     *
     * @param curr			request being checked
     * @param writer		output writer for the report
     * @param sims			map of protein-id pairs to similarity scores
     * @param currTotal		expected number of pairs
     *
     * @throws IOException
     */
    private void checkResumeBatch(Request curr, PrintWriter writer, Map<StringPair, Double> sims, int currTotal) throws IOException {
        if (sims.size() < currTotal) {
            // Here the batch is incomplete.
            if (log.isInfoEnabled() && curr != null)
                log.info("Incomplete pairing {} not output.", curr);
        } else {
            // Here we can write the batch.
            this.writeSims(curr, writer, sims);
            // Remove it from the set of unprocessed requests.
            this.filePairs.remove(curr);
        }
    }

    /**
     * Run a comparison between two files and write the output.
     *
     * @param req		file-processing request
     * @param writer	output writer for report
     */
    private void runCompare(Request req, PrintWriter writer) {
        try {
            long start = System.currentTimeMillis();
            File f1 = req.getFile1();
            File f2 = req.getFile2();
            // Compute the results for this file pairing.
            Map<StringPair, Double> sims = this.method.computeSim(f1, f2);
            // Write all the results for this file pairing.
            this.writeSims(req, writer, sims);
            if (log.isInfoEnabled()) {
                Duration d = Duration.ofMillis(System.currentTimeMillis() - start);
                log.info("{} required to process {}.", d, req);
            }
        } catch (IOException e) {
            // Percolate the exception.  We can't have an unchecked one in a parallel stream.
            throw new RuntimeException(e);
        }
    }

    /**
     * Write out all the protein similarities for a file pairing.
     *
     * @param req		processing request for the file pairing
     * @param writer	output writer for the similarity results
     * @param sims		map of protein ID pairs to similarity scores
     */
    private void writeSims(Request req, PrintWriter writer, Map<StringPair, Double> sims) {
        // Compute the file pairing ID and the total to be used for the resume processing.
        // For each line written, we output the file pair name and the total number of protein
        // pairings output for that file pair.  This information allows us to determine if the file
        // pair needs to be rerun after a resume.  Currently, we do not parallelize, but if we did, the
        // fact that access to the output is single-threaded guarantees that each file pairing's
        // records are grouped, simplifying the resume logic.
        String fpName = req.toString();
        int fpTotal = sims.size();
        synchronized (writer) {
            log.info("{} pairs returned from {}.", fpTotal, fpName);
            for (var mapEntry : sims.entrySet()) {
                StringPair pair = mapEntry.getKey();
                double sim = mapEntry.getValue();
                writer.format("%s\t%s\t%8.6f\t%s\t%d%n", pair.getString1(), pair.getString2(), sim,
                        fpName, fpTotal);
            }
            // Flush the output to increase the likelihood of seeing whole output batches.
            writer.flush();
            // Count this output.
            this.batchCounter++;
            log.info("{} of {} comparisons output.", this.batchCounter, this.batchTotal);
        }
    }

}
