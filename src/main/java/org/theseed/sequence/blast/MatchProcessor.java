/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.proteins.DnaTranslator;
import org.theseed.reports.MatchReporter;
import org.theseed.reports.NaturalSort;
import org.theseed.sequence.DnaDataStream;
import org.theseed.sequence.DnaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceDataStream;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.utils.BaseProcessor;

/**
 * This command takes DNA from the input and BLASTS it against the contigs of a genome, then parses proteins
 * out of the DNA that has good hits.
 *
 * The positional parameter is the name of the genome file.  The DNA FASTA file should come in via the standard
 * input.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -i	name of the input file for the DNA (the default is STDIN)
 * -b	number of input sequences to process at a time (the default is 10)
 * -x	distance to extend the genome hit on either side (the default is 50)
 *
 * --maxE		maximum permissible e-value (the default is 1e-10)
 * --minPct		minimum percent coverage of an incoming DNA sequence for a match (the default is 95)
 * --tempDir	temporary directory for BLAST databases; the default is "Temp" in the current directory
 * --maxGap		maximum gap between adjacent hits when they are to be joined.  The default is 500.
 * --filter		file name of a protein BLAST database to use for filtering out the input sequences. The
 * 				default is "uniRoles.faa" in the current directory.
 * --upstream	upstream region to include with output proteins
 * --outFormat	output format
 *
 * @author Bruce Parrello
 *
 */
public class MatchProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static final Logger log = LoggerFactory.getLogger(MatchProcessor.class);
    /** DNA input stream */
    private DnaInputStream inStream;
    /** input genome */
    private Genome genome;
    /** filtering BLAST datanase */
    private BlastDB filterDb;
    /** BLAST parameters for the main query */
    private BlastParms mainParms;
    /** output report */
    private MatchReporter reporter;

    // COMMAND-LINE OPTION

    /** input file */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "rna.fasta", usage = "input file (if not STDIN)")
    private File inFile;

    /** number of DNA sequences to submit in each BLAST call */
    @Option(name = "-b", aliases = { "--batchSize", "--batch" }, metaVar = "1",
            usage = "number of input sequences to submit to each BLAST call")
    private int batchSize;

    /** distance to extend the genome hit on either side */
    @Option(name = "-x", aliases = { "--extend" }, metaVar = "20",
            usage = "distance to extend the genome hit on either side")
    private int extend;

    /** maximum permissible e-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-20", usage = "maximum permissible e-value for a match")
    private double eValue;

    /** temporary directory for BLAST database */
    @Option(name = "--tempDir", metaVar = "Tmp", usage = "temporary directory for BLAST databases")
    private File tempDir;

    /** maximum gap between adjacent sequences */
    @Option(name = "--maxGap", metaVar = "100", usage = "maximum gap between hits to join")
    private int maxGap;

    /** protein filter file */
    @Option(name = "--filter", metaVar = "proteins.faa", usage = "protein database for filtering input sequences")
    private File filterFile;

    /** upstream region size */
    @Option(name = "--upstream", metaVar = "50", usage = "upstream DNA to include in FASTA output")
    private int upstream;

    /** output format */
    @Option(name = "--format", aliases = { "--outFormat", "--outFmt" }, usage = "output format")
    private MatchReporter.Type outFormat;

    /** file containing the target genome */
    @Argument(index = 0, metaVar = "genome.gto", usage = "target genome file containing the proteins",
            required = true)
    private File genomeFile;

    @Override
    protected void setDefaults() {
        this.inFile = null;
        this.batchSize = 10;
        this.eValue = 1e-100;
        this.maxGap = 500;
        this.extend = 50;
        this.tempDir = new File(System.getProperty("user.dir"), "Temp");
        this.filterFile = new File(System.getProperty("user.dir"), "uniRoles.faa");
        this.upstream = 40;
        this.outFormat = MatchReporter.Type.TRAINING;
        this.reporter = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Read the genome.
        log.info("Loading genome from {}.", this.genomeFile);
        this.genome = new Genome(this.genomeFile);
        // Connect to the DNA input stream.
        if (this.inFile == null) {
            this.inStream = new DnaInputStream(System.in, this.genome.getGeneticCode());
            log.info("DNA sequences will be read from standard input.");
        } else {
            this.inStream = new DnaInputStream(this.inFile, this.genome.getGeneticCode());
            log.info("DNA sequences will be read from {}.", this.inFile);
        }
        if (! this.tempDir.isDirectory()) {
            log.info("Creating temporary file directory {}.", this.tempDir);
            FileUtils.forceMkdir(this.tempDir);
        }
        // Validate the filter file.
        if (! this.filterFile.canRead())
            throw new FileNotFoundException("Protein FASTA " +
                    this.filterFile + " not found or unreadable.  Use --filter parameter to specify.");
        else {
            log.info("Input sequences will be filtered by BLAST database in {}.", this.filterFile);
            try {
                this.filterDb = ProteinBlastDB.createOrLoad(this.filterFile);
            } catch (InterruptedException e) {
                throw new IOException("Interruption in BLAST database creation: " + e.getMessage());
            }
        }
        // Validate the parameters.
        if (this.eValue >= 1.0)
            throw new IllegalArgumentException("Invalid eValue specified.  Must be less than 1.");
        if (this.batchSize <= 0)
            throw new IllegalArgumentException("Batch size must be 1 or more.");
        if (this.maxGap < 0)
            throw new IllegalArgumentException("Maximum gap must be 0 or more.");
        if (this.upstream < 0)
            throw new IllegalArgumentException("Upstream length must be 0 or more.");
        if (this.extend < 0)
            throw new IllegalArgumentException("Extension length must be 0 or more.");
        if (! this.genomeFile.canRead())
            throw new FileNotFoundException("Input file" + this.genomeFile + " is not found or unreadable.");
        return true;
    }

    @Override
    public void run() {
        try {
            // Create the BLAST database.
            BlastDB blastDB = this.createBlastDb();
            // Create the BLAST parameters.
            this.mainParms = new BlastParms().maxE(this.eValue);
            // Create the output report.
            this.reporter = this.outFormat.create(System.out, this.genome);
            reporter.initialize(this);
            // Now we loop through the input FASTA stream, building batches.
            int batchCount = 0;
            int seqCount = 0;
            Iterator<SequenceDataStream> batcher = this.inStream.batchIterator(batchSize);
            while (batcher.hasNext()) {
                DnaDataStream batch = (DnaDataStream) batcher.next();
                batchCount++;
                seqCount += batch.size();
                log.info("Processing input batch {} with {} sequences.", batchCount, batch.size());
                this.processBatch(batch, blastDB);
            }
            // Finish the output report.
            reporter.finish();
            log.info("All done. {} sequences in {} batches processed.", seqCount, batchCount);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        } finally {
            this.inStream.close();
            if (this.reporter != null)
                this.reporter.close();
        }
    }

    /**
     * Process a batch of input sequences and write their output.
     *
     * @param batch		batch of DNA sequences to process
     * @param blastDB	blast database to blast against
     *
     * @throws InterruptedException
     * @throws IOException
     */
    private void processBatch(DnaDataStream batch, BlastDB blastDB) throws IOException, InterruptedException {
        // This will be the filtered query stream.
        DnaDataStream queryStream;
        // Here we have to filter using the filter database.
        Map<String, List<BlastHit>> hitMap = BlastHit.sort(batch.blast(this.filterDb, this.mainParms));
        queryStream = new DnaDataStream(this.batchSize, batch.getGeneticCode());
        // For each RNA sequence we need to pull out the hits that are close to each other and create a
        // set of the ones that are hit backward.
        Set<String> minusSeqs = new HashSet<String>(batch.size());
        for (Sequence seq : batch) {
            if (hitMap.containsKey(seq.getLabel())) {
                queryStream.add(seq);
                // Compute the direction of this sequence's hits.
                int dirCount = 0;
                for (BlastHit hit : hitMap.get(seq.getLabel())) {
                    if (hit.getQueryLoc().getDir() == '+')
                        dirCount++;
                    else
                        dirCount--;
                }
                if (dirCount < 0) minusSeqs.add(seq.getLabel());
            }
        }
        log.info("{} sequences left in batch after filtering.", queryStream.size());
        // Perform the BLAST.
        hitMap = BlastHit.sort(queryStream.blast(blastDB, this.mainParms));
        // Get a sorted list of the query sequence IDs that produced results.
        ArrayList<String> queryList = new ArrayList<String>(hitMap.keySet());
        queryList.sort(new NaturalSort());
        // Get a map of query sequence IDs to sequences.
        Map<String, String> qMap = queryStream.stream().collect(Collectors.toMap(Sequence::getLabel, Sequence::getSequence));
        // Get the DNA translator for the query sequences.
        DnaTranslator xlator = new DnaTranslator(this.genome.getGeneticCode());
        // Loop through the query sequences.
        for (String seqId : queryList) {
            // Now we need to get the proteins.  First, we need to orient the sequence properly.
            String targetDna = qMap.getOrDefault(seqId, "").toLowerCase();
            if (minusSeqs.contains(seqId)) targetDna = Contig.reverse(targetDna);
            List<Sequence> prots = xlator.opTranslate(targetDna, 1, targetDna.length(), seqId, this.upstream);
            // Produce the output.
            if (prots.size() > 0)
                this.reporter.processSequence(seqId, targetDna, prots);
        }
    }

    /**
     * Read in the genome and create the blast database.
     *
     * @throws IOException
     * @throws InterruptedException
     */
    private BlastDB createBlastDb() throws IOException, InterruptedException {
        File tempFile = File.createTempFile("blast", ".fasta", this.tempDir);
        tempFile.deleteOnExit();
        BlastDB retVal = DnaBlastDB.create(tempFile, genome);
        return retVal;
    }

}
