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
 * This is a simple summary report that just counts the inputs and outputs.
 *
 * @author Bruce Parrello
 *
 */
public class MatchSummaryReporter extends MatchReporter {

    // FIELDS
    /** number of proteins in the current section */
    private int protCount;
    /** number of RNA fragments in the current section */
    private int rnaCount;
    /** total number of proteins */
    private int protTotal;
    /** total number of RNA fragments */
    private int rnaTotal;
    /** total number of samples */
    private int sampleCount;

    /**
     * Construct a summary report.
     *
     * @param output	output stream
     */
    public MatchSummaryReporter(OutputStream output) {
        super(output);
    }

    @Override
    public void initialize() throws IOException {
        this.protTotal = 0;
        this.rnaTotal = 0;
        this.sampleCount = 0;
        this.print("sample_id\tgenome_id\tgenome_name\tgen_code\trna_count\tprot_count");
    }

    @Override
    protected void beginSection() { }

    @Override
    public void processSequence(String id, Location loc, String dna, List<String> prots)
            throws IOException, InterruptedException {
        this.protCount += prots.size();
        this.rnaCount++;
    }

    @Override
    public void endSection() {
        Genome genome = this.getGenome();
        this.print("%s\t%s\t%s\t%d\t%d\t%d", this.getSampleId(), genome.getId(), genome.getName(), genome.getGeneticCode(),
                this.rnaCount, this.protCount);
        this.rnaTotal += this.rnaCount;
        this.protTotal += this.protCount;
        this.sampleCount++;
    }

    @Override
    public void finish() {
        this.println();
        this.print("TOTAL\t%d\t \t \t \t%d\t%d", this.sampleCount, this.rnaTotal, this.protTotal);
    }


}
