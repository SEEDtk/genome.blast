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
 * This outputs simple ground truth instances.  Each consists of a DNA sequence representing a
 * portion of the genome that matches the RNA sequence and the verified proteins found therein.
 *
 * @author Bruce Parrello
 *
 */
public class MatchGtiReporter extends MatchReporter {

    /**
     * @param output
     * @param genome
     */
    public MatchGtiReporter(OutputStream output, Genome genome) {
        super(output, genome);
        // TODO Auto-generated constructor stub
    }

    @Override
    public void initialize(MatchProcessor base) throws IOException, InterruptedException {
        // TODO Auto-generated method stub

    }

    @Override
    public void processSequence(String id, String dna, List<Sequence> prots) throws IOException, InterruptedException {
        // TODO Auto-generated method stub

    }

    @Override
    public void finish() {
        // TODO Auto-generated method stub

    }

}
