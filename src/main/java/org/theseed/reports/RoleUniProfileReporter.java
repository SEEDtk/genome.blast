/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashSet;
import java.util.Set;

import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastHit;

/**
 * This report shows for each role how many genomes contain a single annotated copy, how many contain at least one profiled copy,
 * how many contain multiple annotated copies, and how many contain no copies.
 *
 * @author Bruce Parrello
 *
 */
public class RoleUniProfileReporter extends UniProfileReporter {

    // FIELDS
    /** count of genomes where the role occurs once */
    private CountMap<String> goodCounts;
    /** count of genomes where the role occurs multiple times */
    private CountMap<String> multiCounts;
    /** count of genomes where the role was found using profiles */
    private CountMap<String> profiledCounts;
    /** count of genomes where the role was missing */
    private CountMap<String> missingCounts;
    /** number of annotations for each role in the current genome */
    private CountMap<String> genomeAnnotationCounts;
    /** set of roles profiled in the current genome */
    private Set<String> genomeProfiledRoles;

    public RoleUniProfileReporter(File outDir) throws FileNotFoundException {
        super(new File(outDir, "roleSourReport.tbl"));
        // Initialize the global count maps.
        this.goodCounts = new CountMap<String>();
        this.multiCounts = new CountMap<String>();
        this.profiledCounts = new CountMap<String>();
        this.missingCounts = new CountMap<String>();
    }

    @Override
    protected void startReport() {
        this.println("role_id\trole_name\tgood\tmultiple\tprofiled\tmissing");
    }

    @Override
    protected void startGenome(Genome genome) {
        // Set up the counters for this genome.
        this.genomeAnnotationCounts = new CountMap<String>();
        this.genomeProfiledRoles = new HashSet<String>(this.getRoles().size());
    }

    @Override
    public void recordAnnotation(Feature feat, String role) {
        this.genomeAnnotationCounts.count(role);
    }

    @Override
    public void recordProfileHit(String role, BlastHit hit) {
        this.genomeProfiledRoles.add(role);
    }

    @Override
    protected void finishGenome() {
        // Now we must analyze what we found in this genome.
        Set<String> roles = this.getRoles();
        for (String role : roles) {
            int roleCount = this.genomeAnnotationCounts.getCount(role);
            if (roleCount > 1)
                this.multiCounts.count(role);
            else if (roleCount == 1)
                this.goodCounts.count(role);
            else if (this.genomeProfiledRoles.contains(role))
                this.profiledCounts.count(role);
            else
                this.missingCounts.count(role);
        }
    }

    @Override
    protected void finishReport() {
    }

}
