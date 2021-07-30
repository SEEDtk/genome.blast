/**
 *
 */
package org.theseed.reports;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.blast.BlastHit;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This is the base class for a universal-role profiling report.
 *
 * @author Bruce Parrello
 *
 */
public abstract class UniProfileReporter extends BaseReporter {

    // FIELDS
    /** current genome */
    private Genome genome;
    /** role definition map */
    private RoleMap roleMap;
    /** sorted set of roles of interest, created when a subclass needs one */
    private SortedSet<String> roleList;

    /**
     * THis enumeration describes the different report types.
     */
    public static enum Type {
        DETAIL {
            @Override
            public UniProfileReporter create(File outDir) throws FileNotFoundException {
                return new DetailUnitProfileReporter(outDir);
            }
        }, ROLES {
            @Override
            public UniProfileReporter create(File outDir) throws FileNotFoundException {
                return new RoleUniProfileReporter(outDir);
            }
        }, GENOMES {
            @Override
            public UniProfileReporter create(File outDir) throws FileNotFoundException {
                return new GenomeUniProfileReporter(outDir);
            }
        };

        /**
         * @return a processor for this report type
         *
         * @param outDir		output directory to contain the report file
         */
        public abstract UniProfileReporter create(File outDir) throws FileNotFoundException;
    }

    public UniProfileReporter(File outFile) throws FileNotFoundException {
        super(new FileOutputStream(outFile));
        this.roleList = null;
    }

    /**
     * Initialize the report.
     *
     * @param roles		role definition map
     */
    public void openReport(RoleMap roles) {
        this.roleMap = roles;
        this.startReport();
    }

    /**
     * Do header processing for the report.
     */
    protected abstract void startReport();

    /**
     * Begin processing a genome.
     *
     * @param genome	current genome
     */
    public void openGenome(Genome genome) {
        this.genome = genome;
        this.startGenome(genome);
    }

    /**
     * Do header processing for a new genome.
     */
    protected abstract void startGenome(Genome genome);

    /**
     * Record an annotated role.
     *
     * @param feat		feature containing the role
     * @param role		ID of the role of interest
     */
    public abstract void recordAnnotation(Feature feat, String role);

    /**
     * Record a found role.  Note that the ID of the feature hit is the subject ID and the function of the feature hit is the
     * subject comment.
     *
     * @param role		ID of the role of interest
     * @param hit		BLAST hit for the profile that found the role
     */
    public abstract void recordProfileHit(String role, BlastHit hit);

    /**
     * Finish processing a genome.
     */
    public void closeGenome() {
        this.finishGenome();
    }

    /**
     * Do trailer processing for a genome.
     */
    protected abstract void finishGenome();

    /**
     * Close this report.
     */
    public void closeReport() {
        this.finishReport();
    }

    /**
     * Do trailer processing for the whole report.
     */
    protected abstract void finishReport();

    /**
     * @return the current genome
     */
    protected Genome getGenome() {
        return this.genome;
    }

    /**
     * @return the ID of the current genome
     */
    protected String getGenomeId() {
        return this.genome.getId();
    }

    /**
     * @return the name of a role
     *
     * @param role	role ID
     */
    public String getRoleName(String role) {
        return this.roleMap.getName(role);
    }

    /**
     * @return the list of roles of interest
     */
    protected SortedSet<String> getRoles() {
        SortedSet<String> retVal = this.roleList;
        if (retVal == null)
            retVal = new TreeSet<String>(this.roleMap.keySet());
        return retVal;
    }



}
