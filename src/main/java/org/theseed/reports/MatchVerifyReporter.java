/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.security.NoSuchAlgorithmException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.similarity.LevenshteinDetailedDistance;
import org.apache.commons.text.similarity.LevenshteinResults;
import org.theseed.counters.CountMap;
import org.theseed.counters.EnumCounter;
import org.theseed.genome.Feature;
import org.theseed.genome.FeatureList;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.MD5Hex;
import org.theseed.sequence.ProteinKmers;

/**
 * This report displays each protein sequence found and compares it to the closest known protein in the
 * target region of the genome.  For genomes where the proteins are well-annotated, this is a good
 * indication of how accurate the output proteins are.
 *
 * @author Bruce Parrello
 *
 */
public class MatchVerifyReporter extends MatchReporter {

    public enum ErrorType {
        EXACT, TOO_SHORT, TOO_LONG, CHANGED, NOT_FOUND;
    }
    // FIELDS
    /** utility object for counting changes */
    private LevenshteinDetailedDistance computer;
    /** number of too-long proteins */
    private EnumCounter<ErrorType> counters;
    /** total number of proteins */
    private int protCount;
    /** total number of records */
    private int gtiCount;
    /** protein ID computer */
    private MD5Hex idFactory;
    /** duplicate-peg counter */
    private CountMap<String> pegCounts;

    /**
     * Create a new verification report.
     *
     * @param output		output stream
     * @param genome		target genome that should contain the proteins found
     */
    public MatchVerifyReporter(OutputStream output) {
        super(output);
        this.computer = new LevenshteinDetailedDistance();
        this.counters = new EnumCounter<ErrorType>(ErrorType.class);
        try {
            this.idFactory = new MD5Hex();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("Error initializing MD5 engine: " + e.getMessage());
        }
    }

    @Override
    public void initialize() throws IOException {
        this.println("sample\tprot_id\trna_id\tlocation\tbest_peg\tdistance\tnotes\tORF");
        this.pegCounts = new CountMap<String>();
        clearCounters();
    }

    /**
     * Reset all the counters to 0.
     */
    public void clearCounters() {
        this.counters.clear();
        this.protCount = 0;
        this.gtiCount = 0;
    }

    /**
     * @return the error counters.
     */
    public EnumCounter<ErrorType> getCounters() {
        return this.counters;
    }

    /**
     * @return the number of GTI records.
     */
    public int getRecordCount() {
        return this.gtiCount;
    }

    /**
     * @return the nubmer of proteins.
     */
    public int getProteinCount() {
        return this.protCount;
    }


    @Override
    public void processSequence(String id, Location loc, String dna, List<String> prots)
            throws IOException, InterruptedException {
        this.gtiCount++;
        Genome genome = this.getGenome();
        // We process each protein individually.  First, however, we need to get the features in
        // the genome that overlap the target region.
        Map<String, ProteinKmers> protMap = getProteinMap(loc);
        // Now loop through the proteins from the input.
        for (String prot : prots) {
            this.protCount++;
            // Get a proteinkmers object for this protein.
            ProteinKmers protKmers = new ProteinKmers(prot);
            // Find the closest protein in our list.  We default to not finding anything.
            double distance = 1.0;
            String fid = "";
            String comment = "No match found.";
            String orfLoc = "";
            ErrorType error = ErrorType.NOT_FOUND;
            for (Map.Entry<String, ProteinKmers> featEntry : protMap.entrySet()) {
                double newDist = protKmers.distance(featEntry.getValue());
                if (newDist < distance) {
                    distance = newDist;
                    fid = featEntry.getKey();
                    if (distance == 0.0) {
                        comment = "";
                        error = ErrorType.EXACT;
                        orfLoc = "";
                    } else {
                        /// Here we need to analyze the nature of the disagreement.  The default is simply a change.
                        String protSeq = protKmers.getProtein();
                        String featSeq = featEntry.getValue().getProtein();
                        LevenshteinResults comparison = this.computer.apply(featSeq, protSeq);
                        comment = String.format("Found sequence has %d insertions, %d deletions, and %d substitutions.",
                                comparison.getInsertCount(), comparison.getDeleteCount(), comparison.getSubstituteCount());
                        error = ErrorType.CHANGED;
                        // We want to highlight cases where the protein is simply longer or shorter.  To do this, we
                        // have to strip off the starts, since the start translates differently than the rest of the
                        // protein.
                        String protCodons = protSeq.substring(1);
                        String featCodons = featSeq.substring(1);
                        if (protSeq.length() > featSeq.length()) {
                            int longer = protSeq.length() - featSeq.length();
                            if (StringUtils.endsWith(protCodons, featCodons)) {
                                comment = String.format("Found sequence has %d extra codons.", longer);
                                error = ErrorType.TOO_LONG;
                            }
                        } else {
                            int shorter = featSeq.length() - protSeq.length();
                            if (StringUtils.endsWith(featCodons, protCodons)) {
                                comment = String.format("Found sequence has %d fewer codons.", shorter);
                                error = ErrorType.TOO_SHORT;
                            }
                        }
                        orfLoc = genome.getFeature(fid).getLocation().extendToOrf(genome).toString();
                    }
                }
            }
            // Count the match type.
            this.counters.count(error);
            // Count the peg.
            this.pegCounts.count(fid);
            // Write this sequence.
            this.print("%s\t%s\t%s\t%s\t%s\t%4.4f\t%s\t%s", this.getSampleId(), this.idFactory.sequenceMD5(prot),
                    id, loc.toString(), fid, distance, comment, orfLoc);
        }
    }

    /**
     * Return protein kmer objects for each protein in the specified region of the genome.
     *
     * @param loc	location of interest
     *
     * @return a map from feature IDs to protein kmer objects
     */
    private Map<String, ProteinKmers> getProteinMap(Location loc) {
        FeatureList contigFeatures = this.getGenome().getContigFeatures(loc.getContigId());
        Collection<Feature> found = contigFeatures.inRegion(loc.getLeft(), loc.getRight());
        // For each feature in the region, we need a proteinkmers object so we can compute
        // the distance.
        Map<String, ProteinKmers> retVal = new HashMap<String, ProteinKmers>(contigFeatures.size());
        for (Feature feat : found) {
            String protein = feat.getProteinTranslation();
            if (protein != null)
                retVal.put(feat.getId(), new ProteinKmers(protein));
        }
        return retVal;
    }

    @Override
    public void finish() { }

    /**
     * @return the map of peg counts
     */
    public CountMap<String> getPegCounts() {
        return this.pegCounts;
    }

    @Override
    protected void beginSection() { }

    @Override
    public void endSection() { }

}
