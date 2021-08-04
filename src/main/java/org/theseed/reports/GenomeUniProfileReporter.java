/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastHit;

/**
 * This class produces a summary of the universal role occurrences in each genome.  For each genome we want to
 * know the number of features annotated with the role and the number found to contain the role using profiles.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeUniProfileReporter extends UniProfileReporter {

    // FIELDS
    /** map of roles to annotated features */
    private Map<String, Set<String>> annotationMap;
    /** map of roles to profile-identified features */
    private Map<String, Set<String>> profileMap;

    public GenomeUniProfileReporter(File outDir) throws FileNotFoundException {
        super(new File(outDir, "genomeSourReport.tbl"));
    }

    @Override
    protected void startReport() {
        this.println("genome_id\tgenome_name\trole_id\trole_name\tannotated\tprofiled\ttotal");
    }

    @Override
    protected void startGenome(Genome genome) {
        Set<String> roles = this.getRoles();
        // For a new genome, we re-create the role maps.
        this.annotationMap = new HashMap<String, Set<String>>(roles.size());
        this.profileMap = new HashMap<String, Set<String>>(roles.size());
        // Give each role an empty set.  This simplifies a lot of logic.
        for (String role : roles) {
            this.annotationMap.put(role, new TreeSet<String>());
            this.profileMap.put(role, new TreeSet<String>());
        }
    }

    @Override
    public void recordAnnotation(Feature feat, String role) {
        Set<String> featSet = this.annotationMap.get(role);
        featSet.add(feat.getId());
    }

    @Override
    public void recordProfileHit(String role, BlastHit hit) {
        Set<String> featSet = this.profileMap.get(role);
        featSet.add(hit.getSubjectId());
    }

    @Override
    protected void finishGenome() {
        // Here we can write out all the counts for the current genome.
        String genomeId = this.getGenomeId();
        String genomeName = this.getGenome().getName();
        // Loop through the roles.
        SortedSet<String> roles = this.getRoles();
        for (String role : roles) {
            String roleName = this.getRoleName(role);
            int annotated = this.annotationMap.get(role).size();
            int profiled = this.profileMap.get(role).size();
            if (annotated != 1 || profiled > 1)
                this.print("%s\t%s\t%s\t%s\t%d\t%d\t%d", genomeId, genomeName, role, roleName, annotated, profiled, annotated + profiled);
        }
    }

    @Override
    protected void finishReport() {
    }

}
