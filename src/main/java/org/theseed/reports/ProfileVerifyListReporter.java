/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

import org.theseed.genome.Feature;
import org.theseed.sequence.blast.BlastHit;

/**
 * This produces a simple list report of the profile hits and misses.  The output is sorted by role within genome.
 *
 * @author Bruce Parrello
 *
 */
public class ProfileVerifyListReporter extends ProfileVerifyReporter {

    // FIELDS
    /** list of lines for each role */
    private Map<String, List<String>> roleMap;

    public ProfileVerifyListReporter(OutputStream output) {
        super(output);
        this.roleMap = new TreeMap<String, List<String>>();
    }

    @Override
    public void openReport() {
        this.println("genome\ttype\trole\tprofile\tq_len\thit_loc\te_value\tp_ident\tbit_score\tq_bit_score\tq_ident\tq_percent\tneighborhood");
    }

    @Override
    public void openGenome(String genomeId) {
        this.roleMap.clear();
    }

    @Override
    public void processHit(String type, String genome, String profile, String role, BlastHit hit, Collection<Feature> feats) {
        String neighborhood = feats.stream().map(x -> x.getId()).collect(Collectors.joining(" ,"));
        this.saveLine(role, String.format("%s\t%s\t%s\t%s\t%d\t%s\t%4.2e\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%4.2f\t%s", genome, type, role, hit.getQueryId(),
                hit.getQueryLen(), hit.getSubjectLoc(), hit.getEvalue(), hit.getPercentIdentity(), hit.getBitScore(),
                hit.getQueryBitScore(), hit.getQueryIdentity(), hit.getQueryPercentMatch(), neighborhood));
    }

    @Override
    public void processMiss(String role, Feature feat) {
        this.saveLine(role, String.format("%s\t%s\t%s\t \t \t%s\t \t \t \t \t  \t \t%s", feat.getParent().getId(), "miss", role,
                feat.getLocation(), feat.getId()));
    }

    private void saveLine(String role, String line) {
        List<String> roleLines = this.roleMap.computeIfAbsent(role, x -> new ArrayList<String>());
        roleLines.add(line);
    }

    @Override
    public void closeGenome() {
        for (List<String> lineList : this.roleMap.values()) {
            for (String line : lineList)
                this.println(line);
        }
    }

    @Override
    public void closeReport() {
    }


}
