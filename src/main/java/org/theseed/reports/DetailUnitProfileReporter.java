/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.FileNotFoundException;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastHit;

/**
 * This is the detail report for profile hits.  For each such hit, it contains the role id and name, the cluster, the hit strength, and
 * the feature ID and function.
 *
 * @author Bruce Parrello
 *
 */
public class DetailUnitProfileReporter extends UniProfileReporter {

    public DetailUnitProfileReporter(File outDir) throws FileNotFoundException {
        super(new File(outDir, "detailSourReport.tbl"));
    }

    @Override
    protected void startReport() {
        this.println("genome_id\trole_id\trole_name\tcluster_id\te_value\tbit_score\tfeature_coverage\tfeature_id\tfeature_function");
    }

    @Override
    protected void startGenome(Genome genome) {
    }

    @Override
    public void recordAnnotation(Feature feat, String role) {
    }

    @Override
    public void recordProfileHit(String role, BlastHit hit) {
        this.print("%s\t%s\t%s\t%s\t%2.6e\t%8.2f\t%8.4f\t%s\t%s", this.getGenomeId(), role, this.getRoleName(role), hit.getQueryId(),
                hit.getEvalue(), hit.getBitScore(), hit.getSubjectPercentMatch(), hit.getSubjectId(), hit.getSubjectDef());
    }

    @Override
    protected void finishGenome() {
    }

    @Override
    protected void finishReport() {
    }

}
