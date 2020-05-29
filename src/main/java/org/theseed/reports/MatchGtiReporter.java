/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.stream.Collectors;

import org.theseed.locations.Location;

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
    public MatchGtiReporter(OutputStream output) {
        super(output);
    }

    @Override
    public void initialize() throws IOException { }

    @Override
    public void processSequence(String id, Location loc, String dna, List<String> prots) throws IOException, InterruptedException {
        String proteinList = prots.stream().collect(Collectors.joining(","));
        this.print("%s\t%s\t%s\t%s\t%s", this.getSampleId(), id, loc.toString(), dna, proteinList);
    }

    @Override
    public void finish() { }

}
