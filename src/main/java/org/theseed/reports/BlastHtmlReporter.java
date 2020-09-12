/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

import static j2html.TagCreator.*;

/**
 * @author Bruce Parrello
 *
 */
public class BlastHtmlReporter extends BlastReporter {

    // FIELDS
    /** type of color display */
    private BlastDB.ColorType colorType;
    /** current sequence object */
    private HtmlFullSequence container;
    /** link object to use */
    private LinkObject linker;

    /**
     * Construct a report stream for HTML reports.
     *
     * @param output	output stream to contain report
     * @param sort		sequence type (query, subject) to sort on
     */
    public BlastHtmlReporter(OutputStream output, BlastDB.SortType sort) {
        super(output, sort);
        this.colorType = BlastDB.ColorType.sim;
        this.linker = new LinkObject.None();
    }

    @Override
    protected void openReport(String title) {
        this.println("<html>");
        this.println(head(title(title)).render());
        this.println("<body>");
        this.println(h1(title).render());
        this.println(BlastHtmlUtilities.showColorInfo(this.colorType).renderFormatted());
    }

    @Override
    protected void showSubtitle(BlastInfo blastInfo) {
        this.println(BlastHtmlUtilities.showBlastInfo(blastInfo, this.getSequencesHit(), this.getSortType(),
                this.getRedundant()).render());
    }

    @Override
    protected void openSection(SeqData data) {
        this.container = new HtmlFullSequence(1, data.getLen(), data.getId() + " " + data.getDef());
    }

    @Override
    protected void processHit(SeqData target, SeqData anchor, BlastHit hit) {
        Color color = this.colorType.computeColor(target, hit);
        HtmlHitSequence hitDescriptor = BlastHtmlUtilities.getHitArrow(target, anchor, hit, color);
        this.container.add(hitDescriptor);
    }

    @Override
    protected void closeSection() {
        String text = this.container.draw(this.linker).render();
        this.println(text);
    }

    @Override
    protected void closeReport() {
        this.println("</body></html>");
    }

    /**
     * Specify the scheme for computing the hit color.
     *
     * @param colorType the color computation scheme to use
     */
    public void setColorType(BlastDB.ColorType colorType) {
        this.colorType = colorType;
    }

}
