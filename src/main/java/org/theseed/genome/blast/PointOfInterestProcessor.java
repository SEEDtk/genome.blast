/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableSet;
import java.util.TreeSet;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.FeatureList;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.locations.Location;
import org.theseed.sequence.DnaDataStream;
import org.theseed.sequence.DnaStream;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.utils.BaseProcessor;

/**
 * This command adds points of interest to a genome.  A point of interest is a feature belonging to a class defined by a percent
 * identity with a known reference sequence.  Each point of interest is created as a new feature in the genome. It is assigned a family
 * ID from the sequence ID and a function from the sequence comment.
 *
 * The positional parameters are the names of the definition FASTA file and the names of the input and output directories.  All the genomes
 * from the input directory will be processed and updated versions written to the output directory.  If a genome already exists in the
 * output directory, it will be skipped.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more detailed log messages
 * -m	minimum percent identity for a match (default 80)
 * -l	minimum percent query coverage for a match (default 65)
 * -b	number of sequences to submit per blast invocation (default 20)
 *
 * --type		feature type to use (default "poi")
 * --clear		erase output directory before processing
 * --maxE		maximum permissible e-value (default 1e-10)
 * --maxGap		maximum permissible gap between related hits (default 150)
 * --reset		delete existing features of the specified type before processing a new genome
 *
 * @author Bruce Parrello
 *
 */
public class PointOfInterestProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(PointOfInterestProcessor.class);
    /** blast database of input sequences */
    private BlastDB referenceDB;
    /** parameters for the BLAST */
    private BlastParms parms;
    /** map of contig IDs to feature lists */
    private Map<String, FeatureList> contigMaps;
    /** sorted set of new features for the current genome */
    private NavigableSet<Feature> newFeats;
    /** map of minimum feature lengths for each family */
    private Map<String, Integer> familyLengths;

    // COMMAND-LINE OPTIONS

    /** minimum percent identity */
    @Option(name = "-m", aliases = { "--minIdent" }, metaVar = "95.0", usage = "minimum percent identity for a match")
    private double minIdentity;

    /** minimum percent match length */
    @Option(name = "-l", aliases = { "--minLen", "--minCoverage" }, metaVar = "80.0", usage = "minimum percent length coverage")
    private double minCoverage;

    /** maximum E-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-5",
            usage = "maximum permissible e-value for a query")
    private double eValue;

    /** type of feature to create */
    @Option(name = "--type", metaVar = "prna", usage = "feature type to create")
    private String featureType;

    /** TRUE if the output directory should be cleared before processing */
    @Option(name = "--clear", usage = "clear output directory before processing")
    private boolean clearFlag;

    /** TRUE if we are rebuilding the genomes, and old points of interest should be deleted */
    @Option(name = "--reset", usage = "remove old points of interest before processing a genome")
    private boolean resetFlag;

    /** maximum gap between adjacent hits */
    @Option(name = "--maxGap", metaVar = "200", usage = "maximum gap between related hits")
    private int maxGap;

    /** input FASTA file with family reference sequences */
    @Argument(index = 0, metaVar = "referenceSeqs.fna", usage = "FASTA file containing reference sequences, family IDs, and functions")
    private File seqFile;

    /** input genome directory */
    @Argument(index = 1, metaVar = "inDir", usage = "genome input directory")
    private File inDir;

    /** output genome directory */
    @Argument(index = 2, metaVar = "outDir", usage = "genome output directory")
    private File outDir;


    @Override
    protected void setDefaults() {
        this.minIdentity = 80.0;
        this.minCoverage = 65.0;
        this.clearFlag = false;
        this.resetFlag = false;
        this.featureType = "poi";
        this.eValue = 1e-10;
        this.maxGap = 150;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Insure the feature type is only letters and digits.
        if (! StringUtils.isAlpha(this.featureType))
            throw new IllegalArgumentException("Feature type must contain only letters.");
        // Validate the blast parameters.
        if (this.minIdentity > 100.0)
            throw new IllegalArgumentException("MinIdent cannot be greater than 100%.");
        if (this.minCoverage > 100.0)
            throw new IllegalArgumentException("MinCoverage cannot be greater than 100%.");
        if (this.eValue <= 0.0)
            throw new IllegalArgumentException("MaxE cannot be less than zero.");
        // Store them in the blast parameters.
        this.parms = new BlastParms().minPercent(this.minIdentity).maxE(this.eValue);
        // Verify the input file.
        if (! this.seqFile.canRead())
            throw new FileNotFoundException("Input FASTA file " + this.seqFile + " not found or unreadable.");
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " not found or invalid.");
        // Set up the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else {
            log.info("Output will be to directory {}.", this.outDir);
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Initialize the stat counters.
        int newFeatures = 0;
        int genomeCount = 0;
        int skipCount = 0;
        // Create the contig hash map and the new-feature set.  They will be cleared and re-used for each genome.
        this.contigMaps = new HashMap<String, FeatureList>(100);
        this.newFeats = new TreeSet<Feature>(new Feature.StrandComparator());
        // Create the family-length map.
        this.familyLengths = new HashMap<String, Integer>(1000);
        // Create the blast database from the reference sequences.  Note that the genetic code does not matter here,
        // since we are looking at RNA.
        log.info("Connecting to BLAST database from reference sequences in {}.", this.seqFile);
        this.referenceDB = DnaBlastDB.createOrLoad(this.seqFile, 11);
        // Get all the input genomes.
        log.info("Scanning input directory {}.", this.inDir);
        GenomeDirectory genomes = new GenomeDirectory(this.inDir);
        for (Genome genome : genomes) {
            String genomeId = genome.getId();
            // Verify that this genome is not already in the output directory.
            File outFile = new File(this.outDir, genomeId + ".gto");
            if (outFile.exists())
                log.info("{} already exists in {}-- skipping.", genomeId, this.outDir);
            else {
                log.info("Processing genome {}.", genome);
                if (this.resetFlag)
                    this.clearGenome(genome);
                // Create a DNA sequence stream from this genome's contigs.
                DnaStream contigs = new DnaDataStream(genome);
                // Run the blast to find hits.
                List<BlastHit> hits = this.referenceDB.blast(contigs, this.parms);
                // Loop through the hits, creating new features.  We need the next available feature index
                // for this feature type.
                int nextIdx = genome.getNextIdNum(this.featureType);
                for (BlastHit hit : hits) {
                    // The subject sequence is the protein family.  The query sequence is the contig.
                    String contigId = hit.getQueryId();
                    FeatureList contigData = this.contigMaps.computeIfAbsent(contigId, k -> genome.getContigFeatures(k));
                    // Use the feature list to figure out if we hit an intergenic region.
                    Location newLoc = hit.getQueryLoc();
                    if (! contigData.isOccupied(newLoc)) {
                        log.debug("Possible POI at {}, coverage {}, family ID {}.", newLoc, hit.getSubjectPercentMatch(), hit.getSubjectId());
                        // Here we have a possible point of interest.  Create a new feature for it.
                        String fid = "fig|" + genomeId + "." + this.featureType + "." + Integer.toString(nextIdx++);
                        Feature feat = new Feature(fid, hit.getSubjectDef(), newLoc);
                        String pgFam = hit.getSubjectId();
                        feat.setPgfam(pgFam);
                        // Verify we have a known minimum length for this family.
                        if (! this.familyLengths.containsKey(pgFam))
                            this.familyLengths.put(pgFam, (int) Math.ceil(hit.getSubjectLen() * this.minCoverage / 100));
                        // Check to see if the new feature should be combined with its neighbors.
                        checkNeighbor(this.newFeats.floor(feat), newLoc, pgFam);
                        checkNeighbor(this.newFeats.ceiling(feat), newLoc, pgFam);
                        // Add it to the feature list.
                        newFeats.add(feat);
                        log.debug("Provisional {} added at {} for family {}.", feat, feat.getLocation(), feat.getPgfam());
                    }
                }
                // Having processed all the hits, keep only the features that are long enough.
                int count = 0;
                for (Feature feat : this.newFeats) {
                    Location loc = feat.getLocation();
                    String pgFam = feat.getPgfam();
                    if (loc.getLength() >= this.familyLengths.get(pgFam)) {
                        genome.addFeature(feat);
                        log.info("{} at {} added to {} for family {}.", feat, loc, genome, pgFam);
                        count++;
                    }
                }
                // Now we are done with the genome.
                log.info("Updating {} to {} with {} new features.", genome, outFile, count);
                genome.update(outFile);
                this.contigMaps.clear();
                this.newFeats.clear();
                newFeatures += count;
                if (count > 0) genomeCount++; else skipCount++;
                log.info("{} features added to {} genomes to this point, {} genomes unchanged.", newFeatures, genomeCount, skipCount);
            }
        }
        log.info("All done. {} features added to {} of {} genomes.", newFeatures, genomeCount, genomes.size());
    }

    /**
     * Delete all of the features of our type from the current genome.
     *
     * @param genome	genome whose features are to be deleted
     */
    private void clearGenome(Genome genome) {
        Collection<Feature> deletes = new ArrayList<Feature>(1000);
        log.info("Deleting {} features from {}.", this.featureType, genome);
        String prefix = "fig|" + genome.getId() + "." + this.featureType;
        for (Feature feat : genome.getFeatures()) {
            if (feat.getId().startsWith(prefix))
                deletes.add(feat);
        }
        for (Feature feat : deletes)
            genome.deleteFeature(feat);
        log.info("{} features deleted.", deletes.size());
    }

    /**
     * Check to determine if an existing feature should be combined with the current feature because they are
     * neighbors.  If the existing feature should be combined, its location is merged and it is removed from
     * the existing-feature structure.
     *
     * @param neighbor	possible neighboring feature
     * @param newLoc	location of the current feature
     * @param pgFam		family ID of the current feature
     */
    private void checkNeighbor(Feature neighbor, Location newLoc, String pgFam) {
        if (neighbor != null && neighbor.getPgfam().contentEquals(pgFam) &&
                neighbor.getLocation().strandDistance(newLoc) < this.maxGap) {
            this.newFeats.remove(neighbor);
            newLoc.merge(neighbor.getLocation());
        }
    }

}
