/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.reports.BlastHtmlReporter;
import org.theseed.reports.BlastReporter;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.BlastProcessor;
import org.theseed.sequence.blast.ProteinProfiles;
import org.theseed.sequence.blast.Source;
import org.theseed.utils.BaseProcessor;

/**
 * This command performs a profile search against a BLAST database.  The profiles must be stored in a
 * profile directory, such as one made by "ProfileMakeProcessor".  The positional parameters are the
 * profile directory, the database type (db, dna, prot, contigs, pegs, features), and the database file name.
 *
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
 * --minSubject	minimum percent coverage for subject sequences; the default is 0
 * --format		output format; the default is HTML
 * --sort		sort order of output (QUERY or SUBJECT); the default is QUERY
 * --keep		if specified, the BLAST database will be kept (ignored if database type is "db")
 * --color		color computation scheme for HTML reports; the default is "sim"
 * --minIdent	minimum percent identity for hits; the default is 0
 *
 * @author Bruce Parrello
 *
 */
public class ProfileProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BlastProcessor.class);
    /** BLAST database */
    private BlastDB subject;
    /** profile directory */
    private ProteinProfiles profiler;
    /** output reporter */
    private BlastReporter reporter;

    // COMMAND-LINE OPTIONS
    /** working directory for created files */
    @Option(name = "-d", aliases = { "--dir", "--workDir" }, metaVar = "DB", usage = "working directory for created files")
    private File workDir;

    /** number of threads to use while BLASTing */
    @Option(name = "-t", aliases = { "--threads" }, metaVar = "8", usage = "number of threads to use")
    private int numThreads;

    /** genetic code for DNA sequences */
    @Option(name = "--gc", aliases = { "--genCode", "--geneticCode" }, metaVar = "4",
            usage = "genetic code for DNA sequences not in genomes")
    private int geneticCode;

    /** maximum E-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-5",
            usage = "maximum permissible e-value for a query")
    private double eValue;

    /** minimum percent subject coverage for a legitimate hit */
    @Option(name = "--minSubject", aliases = { "--minS" }, metaVar = "50",
            usage  = "minimum percent of subject sequence that must be hit")
    private double minPctSubject;

    /** minimum percent identity */
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
    private BlastReporter.SortType sortType;

    /** TRUE to keep created files for the BLAST database, else FALSE */
    @Option(name = "--keep", usage = "keep BLAST database if one is created")
    private boolean keepDB;

    /** color computation scheme for HTML reports */
    @Option(name = "--color", usage = "color computation scheme for hits")
    private BlastHtmlReporter.ColorType colorType;

    /** profile direcory */
    @Argument(index = 0, usage = "protein profile directory", required = true)
    private File protFile;

    /** type of subject file */
    @Argument(index = 1, usage = "type of database input file", required = true)
    private Source subjectType;

    /** name of subject file */
    @Argument(index = 2, usage = "file containing subject sequences", required = true)
    private File subjectFile;


    @Override
    protected void setDefaults() {
        this.workDir = new File(System.getProperty("user.dir"));
        this.eValue = 1e-10;
        this.format = BlastReporter.Type.HTML;
        this.geneticCode = 11;
        this.keepDB = false;
        this.maxPerQuery = 100;
        this.minPctSubject = 0.0;
        this.minPctIdentity = 0.0;
        this.numThreads = 1;
        this.sortType = BlastReporter.SortType.SUBJECT;
        this.colorType = BlastHtmlReporter.ColorType.ident;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Insure the work directory exists.
        if (! this.workDir.isDirectory()) {
            log.info("Creating working file directory {}.", this.workDir);
            FileUtils.forceMkdir(this.workDir);
        }
        // Validate the BLAST parameters.
        if (this.eValue < 0.0)
            throw new IllegalArgumentException("Maximum e-value cannot be negative.");
        if (this.maxPerQuery < 1)
            throw new IllegalArgumentException("Cannot keep less than one result per query.");
        if (this.numThreads < 1)
            throw new IllegalArgumentException("At least one thread is required.");
        // Validate the filters.
        if (this.minPctSubject < 0.0 || this.minPctSubject > 100.0)
            throw new IllegalArgumentException("Minimum subject coverage must be between 0 and 100.");
        if (this.minPctIdentity < 0.0 || this.minPctIdentity > 100.0)
            throw new IllegalArgumentException("Minimum percent identity must be between 0 and 100.");
        // Create the profiler.
        log.info("Opening profile directory {}.", this.protFile);
        this.profiler = new ProteinProfiles(this.protFile);
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
    public void run() {
        try {
            // Create the BLAST parameters.
            BlastParms parms = new BlastParms().maxE(this.eValue).num_threads(this.numThreads)
                    .maxPerQuery(this.maxPerQuery).pctLenOfSubject(this.minPctSubject)
                    .minPercent(this.minPctIdentity);
            // Get the hits for all the profiles.
            Map<String, List<BlastHit>> profileHits = this.profiler.profile(subject, parms);
            log.info("{} profiles found hits.", profileHits.size());
            int hitCount = 0;
            for (Map.Entry<String, List<BlastHit>> profileEntry : profileHits.entrySet()) {
                List<BlastHit> list = profileEntry.getValue();
                log.debug("Recording {} hits for {}.", list.size(), profileEntry.getKey());
                hitCount += list.size();
                for (BlastHit hit : list)
                    this.reporter.recordHit(hit);
            }
            log.info("Writing report. {} total hits recorded.", hitCount);
            int seqCount = this.profiler.size();
            BlastReporter.Info blastInfo = new BlastReporter.Info(subject.getBlastParms(), seqCount,
                    seqCount - profileHits.size(), hitCount);
            this.reporter.writeReport(subject.getBlastType().toUpperCase()
                    + " run against " + this.subjectFile.getName(), blastInfo);
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

}
