/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
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
import org.theseed.genome.Feature;
import org.theseed.genome.FeatureList;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Role;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.sequence.blast.ProteinProfiles;
import org.theseed.utils.BaseProcessor;

/**
 * This command verifies a profile against a genome.  The profile will be applied to the genome's contigs.
 * We will output the hits that do not correspond to protein locations and count the missed features
 * with roles that should have been found.
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

    /**
     * Representation of a bad blast hit
     */
    private static class Miss {
        public BlastHit hit;
        public String neighborhood;

        public Miss(BlastHit hit, Collection<Feature> neighbors) {
            this.hit = hit;
            this.neighborhood = neighbors.stream().map(x -> x.getId()).collect(Collectors.joining(", "));
        }

    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProfileVerifyProcessor.class);
    /** profile directory manager */
    private ProteinProfiles profiler;
    /** list of bad hits */
    private List<Miss> badHits;
    /** count of missed features */
    private int missCount;
    /** count of eligible features */
    private int goodCount;
    /** temporary file for blast databases */
    private File blastFile;
    /** empty list of blast hits */
    private static final List<BlastHit> NO_HITS = new ArrayList<BlastHit>();

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
        // Create the temporary file for the BLAST databases.
        this.blastFile = File.createTempFile("blast", "fa", this.workDir);
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
        // Initialize the result buffers.
        this.badHits = new ArrayList<Miss>(100);
        this.missCount = 0;
        this.goodCount = 0;
        // Create the BLAST parameters.
        BlastParms parms = new BlastParms().maxE(this.eValue).num_threads(this.numThreads).minPercent(this.minPctIdentity)
                .minQueryBitScore(this.minQbsc).minQueryIdentity(this.minQIdent);
        // Loop through the genomes.
        for (File gFile : this.genomeFiles) {
            Genome genome = new Genome(gFile);
            log.info("Testing genome {}.", genome);
            // Create the blast database for this genome.
            DnaBlastDB blastDB = DnaBlastDB.create(this.blastFile, genome);
            // Profile the genome.
            Map<String, List<BlastHit>> hitMap = this.profiler.profile(blastDB, parms);
            // We need to process each contig in the genome.
            for (Contig contig : genome.getContigs()) {
                // Get the contig's hits.
                String contigId = contig.getId();
                List<BlastHit> hits = hitMap.getOrDefault(contigId, NO_HITS);
                // Get the contig's feature list.
                FeatureList contigFeatures = genome.getContigFeatures(contigId);
                // Fill this set with all the features that have useful roles.
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
                    // If there is nothing being hit, save the hit as bad.
                    if (badHit)
                        this.badHits.add(new Miss(hit, overlap));
                }
                // Count the missed features.
                int misses = missedFeatures.size();
                this.missCount += misses;
                log.info("{} features not found in contig {}.", misses, contigId);
            }
        }
        log.info("{} total features missed out of {} possible. {} bad hits.", this.missCount,
                this.goodCount, this.badHits.size());
        // Now we produce the report.
        System.out.println("profile\tq_len\thit_loc\te_value\tp_ident\tbit_score\tq_bit_score\tq_ident");
        for (Miss miss : this.badHits) {
            BlastHit hit = miss.hit;
            System.out.format("%s\t%d\t%s\t%4.2e\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%s%n", hit.getQueryId(), hit.getQueryLen(),
                    hit.getSubjectLoc(), hit.getEvalue(), hit.getPercentIdentity(), hit.getBitScore(),
                    hit.getQueryBitScore(), hit.getQueryIdentity(), miss.neighborhood);
        }
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
