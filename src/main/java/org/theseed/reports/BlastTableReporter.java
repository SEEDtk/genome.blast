/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;

import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

/**
 * This reporter produces a tabular flat file report of the BLAST hits.
 *
 * @author Bruce Parrello
 *
 */
public class BlastTableReporter extends BlastReporter {

    /**
     * Construct a new tabular reporting stream.
     *
     * @param output	output stream to contain the report
     * @param sort		type of sequence (subject, query) to sort the report on
     */
    public BlastTableReporter(OutputStream output, BlastDB.SortType sort) {
        super(output, sort);
    }

    @Override
    protected void openReport(String title) {
        // We don't use the title; instead, we simply display the header line.
        this.println(BlastHit.PRINT_HEADER);
    }

    @Override
    protected void openSection(SeqData data) {
        // Space before each new section.
        this.println();
    }

    @Override
    protected void processHit(SeqData target, SeqData anchor, BlastHit hit) {
        // Display this hit.
        this.println(hit.getPrintLine());
    }

    @Override
    protected void closeSection() {
    }

    @Override
    protected void closeReport() {
    }

}
