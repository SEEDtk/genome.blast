/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;
import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.reports.BlastReporter;
import org.theseed.sequence.SequenceInputStream;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command performs a BLAST between two files.  Each file can be a full BLAST database, a FASTA file,
 * or a genome GTO file.
 *
 * The positional parameters are the query type (db, dna, prot, contigs, pegs, features), the query file
 * name, the database type (db, dna, prot, contigs, pegs, features), and the database file name.
 *
 * The types are as follows
 *
 * db			existing BLAST database
 * dna			DNA FASTA file
 * prot			protein FASTA file
 * contigs		genome GTO file, use the DNA contigs
 * pegs			genome GTO file, use the protein translations of protein-encoding features
 * features		genome GTO file, use the DNA sequences of the features
 *
 * The following command-line options are supported
 *
 * -h		display command-line usage
 * -v		show more detailed progress messages
 * -d		working directory for storing created files; the default is the current directory
 * -t		number of threads to use; the default is 1
 * -b		number of sequences to submit per batch; the default is 20;
 *
 * --gc			genetic code for type "dna" sequence files; the default is 11
 * --maxE		maximum permissible e-value; the default is 1e-10
 * --minQuery	minimum percent coverage for query sequences; the default is 0
 * --minSubject	minimum percent coverage for subject sequences; the default is 0
 * --max		maximum number of results to return per query sequence; the default is 100
 * --format		output format; the default is HTML
 * --sort		sort order of output (QUERY or SUBJECT); the default is QUERY
 * --keep		if specified, the BLAST database will be kept (ignored if database type is "db")
 * --color		color computation scheme for HTML reports; the default is "sim"
 * --minIdent	minimum percent identity for hits; the default is 0
 *
 * @author Bruce Parrello
 *
 */
public class BlastProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BlastProcessor.class);
    /** BLAST database */
    private BlastDB subject;
    /** query sequence stream */
    private SequenceInputStream query;
    /** output reporter */
    private BlastReporter reporter;

    // COMMAND-LINE OPTIONS
    /** working directory for created files */
    @Option(name = "-d", aliases = { "--dir", "--workDir" }, metaVar = "DB", usage = "working directory for created files")
    private File workDir;

    /** number of threads to use while BLASTing */
    @Option(name = "-t", aliases = { "--threads" }, metaVar = "8", usage = "number of threads to use")
    private int numThreads;

    /** query batch size */
    @Option(name = "-b", aliases = { "--batch", "--batchSize" }, metaVar = "50", usage = "batch size for query submission")
    private int batchSize;

    /** genetic code for DNA sequences */
    @Option(name = "--gc", aliases = { "--genCode", "--geneticCode" }, metaVar = "4",
            usage = "genetic code for DNA sequences not in genomes")
    private int geneticCode;

    /** maximum E-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-5",
            usage = "maximum permissible e-value for a query")
    private double eValue;

    /** minimum percent query coverage for a legitimate hit */
    @Option(name = "--minQuery", aliases = { "--minQ" }, metaVar = "75",
            usage  = "minimum percent of query sequence that must be hit")
    private double minPctQuery;

    /** minimum percent subject coverage for a legitimate hit */
    @Option(name = "--minSubject", aliases = { "--minS" }, metaVar = "50",
            usage  = "minimum percent of subject sequence that must be hit")
    private double minPctSubject;

    /** minimum percent indentity */
    @Option(name = "--minIdent", aliases = { "--percIdentity", "--minI" }, metaVar = "75",
            usage = "minimum percent identity for a hit")
    private double minPctIdentity;

    /** maximum number of results to return per query */
    @Option(name = "--max", aliases = { "--maxHSP", "--maxPerQuery" }, metaVar = "10",
            usage = "maximum number of hits to return for each query")
    private int maxPerQuery;

    /** output format */
    @Option(name = "--format", aliases = { "--outFormat", "--outFmt" }, usage = "output report format")
    private BlastReporter.Type format;

    /** sort type of output format */
    @Option(name = "--sort", usage = "type of sequence to sort on in output report")
    private BlastDB.SortType sortType;

    /** TRUE to keep created files for BLAST databases, else FALSE */
    @Option(name = "--keep", usage = "keep BLAST database if one is created")
    private boolean keepDB;

    /** color computation scheme for HTML reports */
    @Option(name = "--color", usage = "color computation scheme for hits")
    private BlastDB.ColorType colorType;

    /** type of query file */
    @Argument(index = 0, usage = "type of query input file", required = true)
    private Source queryType;

    /** name of query file */
    @Argument(index = 1, usage = "file containing query sequences", required = true)
    private File queryFile;

    /** type of subject file */
    @Argument(index = 2, usage = "type of database input file", required = true)
    private Source subjectType;

    /** name of subject file */
    @Argument(index = 3, usage = "file containing subject sequences", required = true)
    private File subjectFile;


    @Override
    protected void setDefaults() {
        this.workDir = new File(System.getProperty("user.dir"));
        this.eValue = 1e-10;
        this.format = BlastReporter.Type.HTML;
        this.geneticCode = 11;
        this.keepDB = false;
        this.maxPerQuery = 100;
        this.minPctQuery = 0.0;
        this.minPctSubject = 0.0;
        this.minPctIdentity = 0.0;
        this.numThreads = 1;
        this.sortType = BlastDB.SortType.QUERY;
        this.colorType = BlastDB.ColorType.ident;
        this.batchSize = 20;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the work directory exists.
        if (! this.workDir.isDirectory()) {
            log.info("Creating working file directory {}.", this.workDir);
            FileUtils.forceMkdir(this.workDir);
        }
        // Validate the BLAST parameters.
        if (this.eValue < 0.0)
            throw new ParseFailureException("Maximum e-value cannot be negative.");
        if (this.maxPerQuery < 1)
            throw new ParseFailureException("Cannot keep less than one result per query.");
        if (this.numThreads < 1)
            throw new ParseFailureException("At least one thread is required.");
        // Validate the filters.
        if (this.minPctQuery < 0.0 || this.minPctQuery > 100.0)
            throw new ParseFailureException("Minimum query coverage must be between 0 and 100.");
        if (this.minPctSubject < 0.0 || this.minPctSubject > 100.0)
            throw new ParseFailureException("Minimum subject coverage must be between 0 and 100.");
        if (this.minPctIdentity < 0.0 || this.minPctIdentity > 100.0)
            throw new ParseFailureException("Minimum percent identity must be between 0 and 100.");
        // Validate the tuning parameters.
        if (this.batchSize < 1)
            throw new ParseFailureException("Query batch size must be at least 1");
        // Create the query sequence stream.
        log.info("Opening query sequence stream of type {} in {}.", this.queryType, this.queryFile);
        this.query = this.queryType.query(this.workDir, this.queryFile, this.geneticCode);
        // Create the subject database.
        try {
            log.info("Connecting to subject database of type {} in {}.", this.subjectType, this.subjectFile);
            this.subject = this.subjectType.subject(this.workDir, this.subjectFile, this.geneticCode, keepDB);
        } catch (InterruptedException e) {
            throw new RuntimeException("Error creating subject database: " + e.getMessage());
        }
        // Create the reporting object.
        this.reporter = this.format.create(System.out, this.sortType);
        log.info("Report format is {} sorted by {}.", this.format, this.sortType);
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        try {
            // Create the BLAST parameters.
            BlastParms parms = new BlastParms().maxE(this.eValue).num_threads(this.numThreads)
                    .maxPerQuery(this.maxPerQuery).pctLenOfQuery(this.minPctQuery)
                    .pctLenOfSubject(this.minPctSubject).minPercent(this.minPctIdentity);
            this.subject.runBlast(this.query, this.batchSize, parms, this.reporter);
        } finally {
            this.query.close();
            this.reporter.close();
        }
    }

}
