/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.MatchProcessor;

/**
 * This is the base reporting class for MatchProcessor.  The main heading is an input sequence.  The
 * details are the individual proteins found in that sequence.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MatchReporter extends BaseReporter {

    public enum Type {
        TRAINING, GTI;

        /**
         * Construct an RNA Sequence Matching report of this type.
         *
         * @param output	output stream
         * @param genome	input genome
         */
        public MatchReporter create(OutputStream output, Genome genome) {
            MatchReporter retVal = null;
            switch (this) {
            case TRAINING :
                retVal = new MatchTrainingReporter(output, genome);
                break;
            case GTI :
                retVal = new MatchGtiReporter(output, genome);
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
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract void initialize(MatchProcessor base) throws IOException, InterruptedException;

    /**
     * Process the output from a sequence.
     *
     * @param id		input sequence ID
     * @param dna		RNA sequence nucleotides
     * @param prots		list of the operon proteins in the RNA (must be nonempty)
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract void processSequence(String id, String dna, List<Sequence> prots) throws IOException, InterruptedException;

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
