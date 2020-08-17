/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;

import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

/**
 * This is a tabular report that shows the blast results in alignment format.  The first column gives the
 * eValue, and the remaining columns show data from the anchor sequence atop the same data from the
 * target sequence.
 *
 * @author Bruce Parrello
 *
 */
public class BlastAlignReporter extends BlastReporter {

    /**
     * Construct an alignment report.
     *
     * @param output	output stream to receive the report
     * @param sort		sort type (query or subject)
     */
    public BlastAlignReporter(OutputStream output, BlastDB.SortType sort) {
        super(output, sort);
    }

    @Override
    protected void openReport(String title) {
        // Display the headings.
        this.println("e_value\tid\tdescription\tstart\tend\tmatch_len\ttot_len\tcoverage\tsequence");
    }

    @Override
    protected void openSection(SeqData data) {
    }

    @Override
    protected void processHit(SeqData target, SeqData anchor, BlastHit hit) {
        // Space before this group.
        this.println();
        // Write the e-value and the anchor sequence.
        this.showData(String.format("%4.2e", hit.getEvalue()),
                anchor, hit);
        // Fix up the target sequence so that the snips are visible.
        this.findSnips(anchor, target);
        // Write the target sequence with the e-column blank.
        this.showData("", target, hit);
    }

    /**
     * Display an output line.
     *
     * @param eColumn	text for the e-value column
     * @param data		sequence data for this line
     * @param hit		blast hit being displayed
     */
    private void showData(String eColumn, SeqData data, BlastHit hit) {
        Location loc = data.getLoc();
        double coverage = hit.getNumSimilar() * 100.0 / data.getLen();
        this.print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%4.2f\t%s", eColumn, data.getId(), data.getDef(),
                loc.getBegin(), loc.getEnd(), loc.getLength(), data.getLen(),
                coverage, data.getAlignment());
    }

    @Override
    protected void closeSection() {
    }

    /**
     * Change the target sequence so that only the SNPs are visible.
     *
     * @param anchor	anchor data
     * @param target	target data
     */
    private void findSnips(SeqData anchor, SeqData target) {
        StringBuilder snips = new StringBuilder(target.getAlignment());
        String anchorAlignment = anchor.getAlignment();
        String targetAlignment = target.getAlignment();
        for (int i = 0; i < anchorAlignment.length(); i++) {
            if (anchorAlignment.charAt(i) == targetAlignment.charAt(i))
                snips.setCharAt(i, '.');
        }
        target.setAlignment(snips.toString());
    }

    @Override
    protected void closeReport() {
    }

}
