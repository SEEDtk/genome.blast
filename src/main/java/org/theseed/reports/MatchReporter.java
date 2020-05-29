/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;

/**
 * This is the base reporting class for MatchProcessor.  The main heading is an input sequence.  The
 * details are the individual proteins found in that sequence.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MatchReporter extends BaseReporter {

    public enum Type {
        GTI, VERIFY;

        /**
         * Construct an RNA Sequence Matching report of this type.
         *
         * @param output	output stream
         */
        public MatchReporter create(OutputStream output) {
            MatchReporter retVal = null;
            switch (this) {
            case GTI :
                retVal = new MatchGtiReporter(output);
                break;
            case VERIFY :
                retVal = new MatchVerifyReporter(output);
                break;
            }
            return retVal;
        }
    }

    // FIELDS
    /** genome of interest */
    private Genome genome;
    /** ID of current sample */
    private String sampleId;

    /**
     * Construct an RNA Sequence Matching report.
     *
     * @param output	output stream
     * @param genome	input genome
     */
    public MatchReporter(OutputStream output) {
        super(output);
        this.genome = null;
        this.sampleId = null;
    }

    /**
     * Initialize the report.
     *
     * @throws IOException
     */
    public abstract void initialize() throws IOException;

    /**
     * Start a section.
     *
     * @param genome	target genome
     * @param sample	source sample ID
     */
    public void startSection(Genome genome, String sampleId) {
        this.genome = genome;
        this.sampleId = sampleId;
    }
    /**
     * Process the output from a sequence.
     *
     * @param id		input sequence ID
     * @param loc 		genome location
     * @param dna		genome nucleotides
     * @param prots		list of the operon proteins in the RNA (must be nonempty)
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract void processSequence(String id, Location loc, String dna, List<String> prots) throws IOException, InterruptedException;

    /**
     * Finish the report.
     */
    public abstract void finish();

    /**
     * @return the genome
     */
    public Genome getGenome() {
        return genome;
    }

    /**
     * @return the sample ID
     */
    public String getSampleId() {
        return sampleId;
    }

}
