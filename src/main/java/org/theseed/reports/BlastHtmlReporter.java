/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

import static j2html.TagCreator.*;

/**
 * @author Bruce Parrello
 *
 */
public class BlastHtmlReporter extends BlastReporter {

    /**
     * This enum indicates how the color of a hit is determined.
     */
    public enum ColorType {
        /** percent identity */
        ident("percent identity"),
        /** percent similarity */
        sim("percent similarity"),
        /** percent coverage of target sequence */
        covg("percent coverage");

        // FIELDS
        private String name;

        private ColorType(String name) {
            this.name = name;
        }

        /**
         * @return the color for a blast hit by the specified target sequence.
         *
         * @param target	target sequence forming the denominator
         * @param hit		blast hit whose color is needed
         */
        public Color computeColor(BlastHit.SeqData target, BlastHit hit) {
            double fraction = 1.0;
            switch (this) {
            case ident:
                fraction = ((double) hit.getNumIdentical()) / hit.getAlignLen();
                break;
            case sim:
                fraction = ((double) hit.getNumSimilar()) / hit.getAlignLen();
                break;
            case covg:
                fraction = ((double) hit.getNumSimilar()) / target.getLen();
                break;
            }
            Color retVal;
            if (fraction >= 1.0)
                retVal = Color.BLUE;
            else if (fraction >= 0.9)
                retVal = Color.DARK_GREEN.brighten((1.0 - fraction)*5);
            else if (fraction >= 0.7)
                retVal = Color.ORANGE.brighten((0.9 - fraction)*2.5);
            else if (fraction >= 0.5)
                retVal = Color.RED.brighten((0.7 - fraction)*2.5);
            else
                retVal = Color.DARK_GRAY.brighten(0.5 - fraction);
            return retVal;
        }

        public String description() {
            return this.name;
        }

    }

    // FIELDS
    /** type of color display */
    private ColorType colorType;
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
    public BlastHtmlReporter(OutputStream output, SortType sort) {
        super(output, sort);
        this.colorType = ColorType.sim;
        this.linker = new LinkObject.None();
    }

    @Override
    protected void openReport(String title) {
        this.println("<html>");
        this.println(header(title(title)).render());
        this.println("<body>");
        this.println(h1(title).render());
        this.println(p("Color is determined by " + this.colorType.description() + ".")
                .with(table(tr(
                        th("100%").withStyle("background-color: " + Color.BLUE.html() + "; color: white"),
                        th("90% to 99%").withStyle("background-color: " + Color.DARK_GREEN.html() + ";"),
                        th("70% to 89%").withStyle("background-color: " + Color.ORANGE.html() + ";"),
                        th("50% to 69%").withStyle("background-color: " + Color.RED.html() + ";"),
                        th("0% to 49%").withStyle("background-color: " + Color.DARK_GRAY.html() + "; color: white")
                        ))).renderFormatted());
    }

    @Override
    protected void showSubtitle(BlastReporter.Info blastInfo) {
        this.println(ul(
                li(blastInfo.getParms()),
                li(String.format("%d queries produced %d hits.", blastInfo.getQueriesIn(),
                        blastInfo.getHitCount())),
                li(String.format("%d queries had no hits.", blastInfo.getMissCount())
                        )).render());
    }

    @Override
    protected void openSection(SeqData data) {
        this.container = new HtmlFullSequence(1, data.getLen(), data.getId() + " " + data.getDef());
    }

    @Override
    protected void processHit(SeqData target, SeqData anchor, BlastHit hit) {
        Color color = this.colorType.computeColor(target, hit);
        Location hitLoc = target.getLoc();
        char dir = (hitLoc.getDir() != anchor.getLoc().getDir() ? '-' : '+');
        String label = String.format("[e=%4.2e, ident=%4.1f%%, gap=%d, loc=(%d,%d)/%d] %s",
                hit.getEvalue(), hit.getPercentIdentity(), hit.getNumGap(), hitLoc.getBegin(),
                hitLoc.getEnd(), target.getLen(), target.getDef());
        HtmlHitSequence hitDescriptor = new HtmlHitSequence(target.getId(), label,
                anchor.getLoc(), dir, color);
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
    public void setColorType(ColorType colorType) {
        this.colorType = colorType;
    }

}
