package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.stream.IntStream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
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
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 * -K	kmer size for kmer-based comparisons (default 8)
 *
 * --method		comparison method to use (default KMERS)
 *
 * @author Bruce Parrello
 *
 */
public class ProteinDistanceProcessor extends BaseReportProcessor implements ProteinCompare.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinDistanceProcessor.class);
    /** comparison method engine */
    private ProteinCompare method;
    /** list of input FASTA files */
    private File[] inFiles;
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

    /** input directory */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory containing FASTA files")
    private File inDir;

    @Override
    public int getKmerSize() {
        return this.kmerSize;
    }

    @Override
    protected void setReporterDefaults() {
        this.kmerSize = 8;
        this.methodType = ProteinCompare.Type.KMERS;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        //Insure the input directory exists and has FASTA files.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        this.inFiles = this.inDir.listFiles(FASTA_FILTER);
        log.info("{} files found in input directory.", this.inFiles.length);
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
        writer.println("prot1\tprot2\tdistance");
        // Loop through the file pairs.
        for (int i = 0; i < this.inFiles.length; i++) {
            File file1 = this.inFiles[i];
            log.info("Processing file {} of {}: {}.",  i+1, this.inFiles.length, file1);
            IntStream.range(i, this.inFiles.length).parallel()
                    .forEach(j -> this.runCompare(file1, this.inFiles[j], writer));
        }
    }

    /**
     * Run a comparison between two files and write the output.
     *
     * @param file1		first FASTA file
     * @param file2		second FASTA file
     * @param writer	output writer for report
     */
    private void runCompare(File file1, File file2, PrintWriter writer) {
        try {
            Map<StringPair, Double> distances = this.method.computeDistance(file1, file2);
            synchronized (writer) {
                log.info("{} pairs returned from {} and {}.", distances.size(), file1, file2);
                for (var mapEntry : distances.entrySet()) {
                    StringPair pair = mapEntry.getKey();
                    double dist = mapEntry.getValue();
                    if (dist < 1.0)
                        writer.println(pair.getString1() + "\t" + pair.getString2() + "\t" + Double.toString(dist));
                }
            }
        } catch (IOException e) {
            // Percolate the exception.  We can't have an unchecked one in a parallel stream.
            throw new RuntimeException(e);
        }
    }

}
