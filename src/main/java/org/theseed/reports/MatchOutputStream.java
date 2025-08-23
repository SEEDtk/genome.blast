/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastHit;

/**
 * This is the base class for an object used to output data regarding RNA/genome matches to a file.  A MatchBaseProcessor
 * client can have many of these attached, and each produces a different file for a given sample.  This is distinct from
 * the MatchReporter objects, that tend to produce a single file for an entire run, and can be used by both a
 * MatchBaseProcessor and the MatchVerifyProcessor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MatchOutputStream implements AutoCloseable {

    // FIELDS
    /** genome of interest */
    private final Genome genome;
    /** ID of current sample */
    private final String sampleId;

    /**
     * Type of output stream.
     */
    public static enum Type {
        GTI, FASTA;

        /**
         * Create an output stream of the specified type.
         * @throws IOException
         */
        public MatchOutputStream create(File outputFile, Genome genome, String sampleId) throws IOException {
            MatchOutputStream retVal;
            switch (this) {
            case GTI -> retVal = new GtiMatchOutputStream(outputFile, genome, sampleId);
            case FASTA -> retVal = new FastaMatchOutputStream(outputFile, genome, sampleId);
            default -> throw new IOException("Unknown output stream type: " + this);
            }
            retVal.initialize();
            return retVal;
        }
    }

    /**
     * Construct an output stream from a file.
     *
     * @param outFile		output file
     * @param genome		genome of interest
     * @param sampleId		ID of the current RNA sample
     *
     * @throws IOException
     */
    public MatchOutputStream(File outFile, Genome genome, String sampleId) throws IOException {
        this.genome = genome;
        this.sampleId = sampleId;
    }

    /**
     * Initialize the report.
     */
    protected abstract void initialize() throws IOException;

    /* Process the output from a sequence.
    *
    * @param id		input sequence fragment ID, in the form "r.XXXXXXXX.NNNN", where XXXXXXXX isa the RNA contig
    * 				ID and NNNN is the fragment sequence number
    * @param loc 	genome location
    * @param dna	genome nucleotides
    * @param prots	list of the operon proteins in the RNA in the form of sequence objects (must be nonempty);
    * 				the ID of each protein must be unique, and the comment should be its location in the RNA
    * @param hit	RNA blast hit to the genome
    *
    * @throws InterruptedException
    * @throws IOException
    */
   public abstract void processSequence(String id, Location loc, String dna, List<Sequence> prots, BlastHit hit)
           throws IOException, InterruptedException;

    /**
     * Finish the report.
     */
    protected abstract void finish();

    @Override
    public void close() {
        this.finish();
    }

    /**
     * @return the genome
     */
    protected Genome getGenome() {
        return genome;
    }

    /**
     * @return the sample ID
     */
    protected String getSampleId() {
        return sampleId;
    }

}
