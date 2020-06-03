/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastHit;

/**
 * This output stream produces a simple FASTA file of the proteins.
 *
 * @author Bruce Parrello
 *
 */
public class FastaMatchOutputStream extends MatchOutputStream {

    // FIELDS
    /** stream for writing sequences in FASTA format */
    private FastaOutputStream outStream;

    /**
     * Construct a match output stream that produces a FASTA file.
     *
     * @param outFile	output file
     * @param genome	genome of interest
     * @param sampleId	ID of the RNA sample
     *
     * @throws IOException
     */
    public FastaMatchOutputStream(File outFile, Genome genome, String sampleId) throws IOException {
        super(outFile, genome, sampleId);
        this.outStream = new FastaOutputStream(outFile);
    }

    @Override
    protected void initialize() throws IOException { }

    @Override
    public void processSequence(String id, Location loc, String dna, List<Sequence> prots, BlastHit hit)
            throws IOException, InterruptedException {
        this.outStream.write(prots);
    }

    @Override
    protected void finish() {
        this.outStream.close();
    }

}
