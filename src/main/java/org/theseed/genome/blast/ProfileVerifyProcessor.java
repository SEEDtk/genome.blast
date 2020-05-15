/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.ProteinProfiles;
import org.theseed.utils.BaseProcessor;

/**
 * This command verifies a profile against a genome.  The profile will be applied to the genome's contigs.
 * We will output the features containing profiled roles that are not hit, and the hits that do not
 * correspond to protein locations.
 *
 * The positional parameters are the name of the profile directory and the name of a genome file.
 *
 * The following command-line options are supported
 *
 * -h		display command-line usage
 * -v		show more detailed progress messages
 * -d		working directory for storing created files; the default is the current directory
 * -t		number of threads to use; the default is 1
 *
 * --gc			genetic code for type "dna" sequence files; the default is 11
 * --maxE		maximum permissible e-value; the default is 1e-10
 * --minIdent	minimum percent identity for hits; the default is 0
 * --minQIdent	minimum query identity fraction for hits; the default is 0.75
 * --minQbsc	minimum query-scaled bit score for hits; the default is 0.35
 *
 * @author Bruce Parrello
 *
 */
public class ProfileVerifyProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProfileVerifyProcessor.class);
    /** profile directory manager */
    private ProteinProfiles profiler;
    /** target genome */
    private Genome genome;

    // COMMAND-LINE OPTIONS

    /** working directory for created files */
    @Option(name = "-d", aliases = { "--dir", "--workDir" }, metaVar = "DB", usage = "working directory for created files")
    private File workDir;

    /** number of threads to use while BLASTing */
    @Option(name = "-t", aliases = { "--threads" }, metaVar = "8", usage = "number of threads to use")
    private int numThreads;

    /** maximum E-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-5",
            usage = "maximum permissible e-value for a query")
    private double eValue;

    /** minimum percent identity */
    @Option(name = "--minIdent", aliases = { "--percIdentity", "--minI" }, metaVar = "75",
            usage = "minimum percent identity for a hit")
    private double minPctIdentity;

    /** minimum query-scaled bit score */
    @Option(name = "--minQbsc", metaVar = "1.2", usage = "minimum acceptable query-scaled bit score")
    private double minQbsc;

    /** minimum query identity fraction */
    @Option(name = "--minQIdent", metaVar = "0.5", usage = "minimum acceptable query identity fraction")
    private double minQIdent;

    /** profile direcory */
    @Argument(index = 0, metaVar = "profileDir", usage = "protein profile directory", required = true)
    private File protFile;

    /** target genome files */
    @Argument(index = 1, metaVar = "genome.gto", usage = "target genome file")
    private List<File> genomeFiles;

    @Override
    protected void setDefaults() {
        this.workDir = new File(System.getProperty("user.dir"));
        this.eValue = 1e-10;
        this.minPctIdentity = 0.0;
        this.minQbsc = 0.0;
        this.minQIdent = 0.0;
        this.numThreads = 1;
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
        if (this.numThreads < 1)
            throw new IllegalArgumentException("At least one thread is required.");
        // Validate the filters.
        if (this.minPctIdentity < 0.0 || this.minPctIdentity > 100.0)
            throw new IllegalArgumentException("Minimum percent identity must be between 0 and 100.");
        if (this.minQbsc < 0.0 || this.minQbsc > 10.0)
            throw new IllegalArgumentException("Minimum query-scaled bit score must be between 0 and 10.");
        if (this.minQIdent < 0.0 || this.minQIdent > 1.0)
            throw new IllegalArgumentException("Minimum query identity fraction must be between 0 and 1");
        // Create the profiler.
        log.info("Opening profile directory {}.", this.protFile);
        this.profiler = new ProteinProfiles(this.protFile);
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        throw new RuntimeException("Testing");
        // TODO Auto-generated method stub
    }

}
