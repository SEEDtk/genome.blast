/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.MatchBaseProcessor;

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
         * @param genome	input genome
         */
        public MatchReporter create(OutputStream output, Genome genome) {
            MatchReporter retVal = null;
            switch (this) {
            case GTI :
                retVal = new MatchGtiReporter(output, genome);
                break;
            case VERIFY :
                retVal = new MatchVerifyReporter(output, genome);
                break;
            }
            return retVal;
        }
    }

    // FIELDS
    /** genome of interest */
    private Genome genome;

    /**
     * Construct an RNA Sequence Matching report.
     *
     * @param output	output stream
     * @param genome	input genome
     */
    public MatchReporter(OutputStream output, Genome genome) {
        super(output);
        this.genome = genome;
    }

    /**
     * Initialize the report.
     *
     * @throws IOException
     */
    public abstract void initialize(MatchBaseProcessor base) throws IOException;

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
    public abstract void processSequence(String id, Location loc, String dna, List<Sequence> prots) throws IOException, InterruptedException;

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

}
