/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.Collection;

import org.theseed.genome.Feature;
import org.theseed.sequence.blast.BlastHit;

/**
 * This is the base class for reports from the profile verification processor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ProfileVerifyReporter extends BaseReporter {

    /**
     * This enum specifies the type of output from the profile verification command.
     */
    public static enum Type {
        LIST {
            @Override
            public ProfileVerifyReporter create(OutputStream outStream) {
                return new ProfileVerifyListReporter(outStream);
            }
        };

        public abstract ProfileVerifyReporter create(OutputStream outStream);
    }

    /**
     * Construct a reporting class.
     *
     * @param output	output stream for the report
     */
    public ProfileVerifyReporter(OutputStream output) {
        super(output);
    }

    /**
     * Start the report.
     */
    public abstract void openReport();

    /**
     * Start processing a genome.
     */
    public abstract void openGenome(String genomeId);

    /**
     * Report on a blast hit.
     *
     * @param type		type of hit (good or bad)
     * @param genome	ID of genome hit
     * @param profile	name of the profile used
     * @param role		role ID for the profile making the hit
     * @param hit		descriptor for the blast hit
     * @param feats		collection of neighboring features
     */
    public abstract void processHit(String type, String genome, String profile, String role, BlastHit hit, Collection<Feature> feats);

    /**
     * Report on a missed feature.
     *
     * @param role		role ID for the feature's expected role
     * @param feat		feature missed
     */
    public abstract void processMiss(String role, Feature feat);

    /**
     * Finish processing a genome.
     */
    public abstract void closeGenome();

    /**
     * Finish the report.
     */
    public abstract void closeReport();
}
