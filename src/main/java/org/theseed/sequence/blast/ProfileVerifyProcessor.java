/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.FeatureList;
import org.theseed.genome.Genome;
import org.theseed.io.GtoFilter;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;
import org.theseed.proteins.Role;
import org.theseed.reports.ProfileVerifyReporter;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command verifies a profile against a genome.  The profile will be applied to the genome's contigs.
 * We will output the hits that do not correspond to protein locations and count the missed features
 * with roles that should have been found.
 *
 * The positional parameters are the name of the profile directory and the name of one or more genome files.
 * If a directory is specified, all of the GTO files in the directory will be processed.
 *
 * The following command-line options are supported
 *
 * -h		display command-line usage
 * -v		show more detailed progress messages
 * -d		working directory for storing created files; the default is the current directory
 * -t		number of threads to use; the default is 1
 *
 * --maxE			maximum permissible e-value; the default is 1e-10
 * --minIdent		minimum percent identity for hits; the default is 0
 * --minQIdent		minimum query identity fraction for hits; the default is 0.50
 * --minQbsc		minimum query-scaled bit score for hits; the default is 0.0
 * --minQuery		minimum percent query match; the default is 0.0
 * --roleFilter		if specified, a file containing roles to use in its first column (tab-delimited with headers); use this
 * 					to restrict processing to a subset of the profile roles
 * --outFormat		output format (default LIST)
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
    /** count of missed features */
    private int missCount;
    /** count of eligible features */
    private int goodCount;
    /** count of bad hits */
    private int badCount;
    /** temporary file for blast databases */
    private File blastFile;
    /** empty list of blast hits */
    private static final List<BlastHit> NO_HITS = new ArrayList<BlastHit>();
    /** list of genome files to process */
    private List<File> genomeFiles;
    /** output report processor */
    private ProfileVerifyReporter reporter;

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
    @Option(name = "--minQbsc", metaVar = "1.1", usage = "minimum acceptable query-scaled bit score")
    private double minQbsc;

    /** minimum query identity fraction */
    @Option(name = "--minQIdent", metaVar = "0.5", usage = "minimum acceptable query identity fraction")
    private double minQIdent;

    /** minimum percent query coverage for a legitimate hit */
    @Option(name = "--minQuery", aliases = { "--minQ" }, metaVar = "75",
            usage  = "minimum percent of query sequence that must be hit")
    private double minPctQuery;

    /** include both good and bad hits in the output */
    @Option(name = "--all", usage = "include good hits as well as bad hits in the output")
    private boolean showAll;

    /** filtering file to restrict roles used */
    @Option(name = "--roleFilter", metaVar = "sours.tbl", usage = "file containing roles to use")
    private File roleFilter;

    /** output report format */
    @Option(name = "--outFormat", usage = "output report format")
    private ProfileVerifyReporter.Type reportType;

    /** profile directory */
    @Argument(index = 0, metaVar = "profileDir", usage = "protein profile directory", required = true)
    private File protFile;

    /** target genome files and directories */
    @Argument(index = 1, metaVar = "genome.gto ...", usage = "target genome file/directory", multiValued = true)
    private List<File> sourceFiles;

    @Override
    protected void setDefaults() {
        this.workDir = new File(System.getProperty("user.dir"));
        this.eValue = 1e-10;
        this.minPctIdentity = 0.0;
        this.minQbsc = 0.0;
        this.minQIdent = 0.5;
        this.numThreads = 1;
        this.showAll = false;
        this.roleFilter = null;
        this.reportType = ProfileVerifyReporter.Type.LIST;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the work directory exists.
        if (! this.workDir.isDirectory()) {
            log.info("Creating working file directory {}.", this.workDir);
            FileUtils.forceMkdir(this.workDir);
        }
        // Insure the role-filter file exists.
        Set<String> filterRoles = null;
        if (this.roleFilter != null) {
            log.info("Reading list of roles to use from {}.", this.roleFilter);
            filterRoles = TabbedLineReader.readSet(this.roleFilter, "1");
            log.info("{} roles found in filter file.", filterRoles.size());
        }
        // Create the temporary file for the BLAST databases.
        this.blastFile = File.createTempFile("blast", ".fa", this.workDir);
        // Validate the BLAST parameters.
        if (this.eValue < 0.0)
            throw new ParseFailureException("Maximum e-value cannot be negative.");
        if (this.numThreads < 1)
            throw new ParseFailureException("At least one thread is required.");
        // Validate the filters.
        if (this.minPctIdentity < 0.0 || this.minPctIdentity > 100.0)
            throw new ParseFailureException("Minimum percent identity must be between 0 and 100.");
        if (this.minQbsc < 0.0 || this.minQbsc > 10.0)
            throw new ParseFailureException("Minimum query-scaled bit score must be between 0 and 10.");
        if (this.minQIdent < 0.0 || this.minQIdent > 1.0)
            throw new ParseFailureException("Minimum query identity fraction must be between 0 and 1");
        if (this.minPctQuery < 0.0 || this.minPctQuery > 100.0)
            throw new ParseFailureException("Minimum query percentation must be between 0 and 100.");
        // Create the profiler.
        log.info("Opening profile directory {}.", this.protFile);
        this.profiler = new ProteinProfiles(this.protFile, filterRoles);
        // Create the report writer.
        this.reporter = this.reportType.create(System.out);
        // Now we must convert the incoming directories to file lists.
        this.genomeFiles = new ArrayList<File>(this.sourceFiles.size() * 10);
        for (File gFile : this.sourceFiles) {
            if (gFile.isDirectory()) {
                File[] gFiles = GtoFilter.getAll(gFile);
                this.genomeFiles.addAll(Arrays.asList(gFiles));
            } else if (! gFile.canRead())
                throw new FileNotFoundException("Input file " + gFile + " not found or unreadable.");
            else
                this.genomeFiles.add(gFile);
        }
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        // Initialize the result counters.
        this.missCount = 0;
        this.goodCount = 0;
        this.badCount = 0;
        // Start the report.
        this.reporter.openReport();
        // Create the BLAST parameters.
        BlastParms parms = new BlastParms().maxE(this.eValue).num_threads(this.numThreads).minPercent(this.minPctIdentity)
                .minQueryBitScore(this.minQbsc).minQueryIdentity(this.minQIdent).pctLenOfQuery(this.minPctQuery);
        // Loop through the genomes.
        for (File gFile : this.genomeFiles) {
            Genome genome = new Genome(gFile);
            log.info("Testing genome {}.", genome);
            this.reporter.openGenome(genome.getId());
            // Create the blast database for this genome.
            DnaBlastDB blastDB = DnaBlastDB.create(this.blastFile, genome);
            blastDB.deleteOnExit();
            // Profile the genome.
            Map<String, List<BlastHit>> hitMap = this.profiler.profile(blastDB, parms);
            // We need to process each contig in the genome.
            for (Contig contig : genome.getContigs()) {
                // Get the contig's hits.
                String contigId = contig.getId();
                List<BlastHit> hits = hitMap.getOrDefault(contigId, NO_HITS);
                // Get the contig's feature list.
                FeatureList contigFeatures = genome.getContigFeatures(contigId);
                // Fill the missing-set with all the features that have useful roles.  Any features that were hit will
                // be removed.
                Set<String> missedFeatures = new HashSet<String>(contigFeatures.size());
                for (Feature feat : contigFeatures) {
                    List<Role> useful = feat.getUsefulRoles(this.profiler.roleMap());
                    if (useful.size() > 0)
                        missedFeatures.add(feat.getId());
                }
                this.goodCount += missedFeatures.size();
                log.debug("Processing contig {} with {} hits against {} features.", contigId, hits.size(), missedFeatures.size());
                // Find the features hit and remove them from the good-features set.  Save the bad hits.
                for (BlastHit hit : hits) {
                    List<Role> profileRoles = Feature.usefulRoles(this.profiler.roleMap(), hit.getQueryDef());
                    Location hitLoc = hit.getSubjectLoc();
                    Collection<Feature> overlap = contigFeatures.inRegion(hitLoc.getLeft(), hitLoc.getRight());
                    // Each feature hit must be removed from the good-features set.  If the roles of the feature don't
                    // overlap the roles of the profile, then it is a bad hit.
                    boolean badHit = true;
                    for (Feature feat : overlap) {
                        if (feat.getLocation().getDir() == hitLoc.getDir() &&
                                this.goodHit(feat, profileRoles)) {
                            // Here the hit is good.
                            missedFeatures.remove(feat.getId());
                            badHit = false;
                        }
                    }
                    // Report the hit.
                    String type = "good";
                    if (badHit) {
                        type = "bad";
                        this.badCount++;
                    }
                    this.reporter.processHit(type, genome.getId(), hit.getQueryId(), profileRoles.get(0).getId(), hit, overlap);
                }
                // Count the missed features.
                int misses = missedFeatures.size();
                this.missCount += misses;
                // Add them to the report.
                for (String fid : missedFeatures) {
                    Feature feat = genome.getFeature(fid);
                    List<Role> roles = feat.getUsefulRoles(this.profiler.roleMap());
                    this.reporter.processMiss(roles.get(0).getId(), feat);
                }
                log.info("{} features not found in contig {}.", misses, contigId);
            }
            this.reporter.closeGenome();
        }
        log.info("{} total features missed out of {} possible. {} bad hits.", this.missCount,
                this.goodCount, this.badCount);
    }

    /**
     * @return TRUE if the roles in the feature include a role in the profile-role list
     *
     * @param feat			feature of interest
     * @param profileRoles	list of roles in the current profile
     */
    private boolean goodHit(Feature feat, List<Role> profileRoles) {
        List<Role> fRoles = feat.getUsefulRoles(this.profiler.roleMap());
        boolean retVal = false;
        for (int i = 0; i < fRoles.size() && ! retVal; i++)
            retVal = profileRoles.contains(fRoles.get(i));
        return retVal;
    }

}
