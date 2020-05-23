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
import org.theseed.sequence.blast.MatchBaseProcessor;

/**
 * This outputs simple ground truth instances.  Each consists of a DNA sequence representing a
 * portion of the genome that matches the RNA sequence and the verified proteins found therein.
 *
 * @author Bruce Parrello
 *
 */
public class MatchGtiReporter extends MatchReporter {

    // FIELDS
    /** base processor containing sample ID */
    private MatchBaseProcessor base;

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
    public void initialize(MatchBaseProcessor base) throws IOException {
        this.base = base;
    }

    @Override
    public void processSequence(String id, Location loc, String dna, List<Sequence> prots) throws IOException, InterruptedException {
        String proteinList = prots.stream().map(x -> x.getSequence()).collect(Collectors.joining(","));
        this.print("%s\t%s\t%s\t%s\t%s", this.base.getSampleID(), id, loc.toString(), dna, proteinList);
    }

    @Override
    public void finish() {
    }

}
