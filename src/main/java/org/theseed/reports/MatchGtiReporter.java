/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.stream.Collectors;

import org.theseed.genome.Genome;
import org.theseed.locations.Location;
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
     * Create a GTI report.
     *
     * @param output	output stream
     * @param genome	controlling genome
     */
    public MatchGtiReporter(OutputStream output, Genome genome) {
        super(output, genome);
    }

    @Override
    public void initialize(MatchProcessor base) throws IOException, InterruptedException {
    }

    @Override
    public void processSequence(String id, Location loc, String dna, List<Sequence> prots) throws IOException, InterruptedException {
        String proteinList = prots.stream().map(x -> x.getSequence()).collect(Collectors.joining("\t"));
        this.print("%s\t%s\t%s\t%s", id, loc.toString(), dna, proteinList);
    }

    @Override
    public void finish() {
    }

}
