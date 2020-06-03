/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.stream.Collectors;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastHit;

/**
 * This output stream produces raw GTIs.  Each consists of some identifying information, a DNA sequence from the
 * genome, and a list of protein sequences derived from the DNA.
 *
 * @author Bruce Parrello
 */
public class GtiMatchOutputStream extends MatchOutputStream {

    // FIELDS
    /** output writer */
    private PrintWriter writer;

    /**
     * Construct a GTI output stream.
     *
     * @param outFile	file to receive the output
     * @param genome	genome from which the DNA is derived
     * @param sampleId	ID of the RNA sample
     *
     * @throws IOException
     */
    public GtiMatchOutputStream(File outFile, Genome genome, String sampleId) throws IOException {
        super(outFile, genome, sampleId);
        this.writer = new PrintWriter(outFile);
    }

    @Override
    protected void initialize() throws IOException { }

    @Override
    public void processSequence(String id, Location loc, String dna, List<Sequence> prots, BlastHit hit)
            throws IOException, InterruptedException {
        String proteinList = prots.stream().map(x -> x.getSequence()).collect(Collectors.joining(","));
        this.writer.format("%s\t%s\t%s\t%s\t%s%n", this.getSampleId(), id, loc.toString(), dna, proteinList);

    }

    @Override
    protected void finish() {
        this.writer.close();
    }

}
