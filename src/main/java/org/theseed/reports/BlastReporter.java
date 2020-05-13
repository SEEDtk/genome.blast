/**
 *
 */
package org.theseed.reports;

import java.io.Closeable;
import java.io.OutputStream;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.Shuffler;
import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

/**
 * This is the base class for BLAST output reports.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BlastReporter extends BaseReporter implements Closeable, AutoCloseable {

    /**
     * This enum indicates the high-level sort for output.
     */
    public enum SortType {
        /** sort by query, list subjects within query */
        QUERY(BlastHit.QUERY, "queries"),
        /** sort by subject, list queries within subject */
        SUBJECT(BlastHit.SUBJECT, "subject sequences");

        // FIELDS
        private int sortIdx;
        private int otherIdx;
        private String plural;

        private SortType(int idx, String plural) {
            this.sortIdx = idx;
            this.otherIdx = 1 - idx;
            this.plural = plural;
        }

        /**
         * @return the ID of the sequence being sorted on
         *
         * @param hit	blast hit whose sort ID is desired
         */
        public String idOf(BlastHit hit) {
            return hit.getData(sortIdx).getId();
        }

        /**
         * @return the ID of the sequence being hit
         *
         * @param hit	blast hit whose sequence ID is desired
         */
        public String targetOf(BlastHit hit) {
            return hit.getData(otherIdx).getId();
        }

        /**
         * @return TRUE if the hit has the specified target ID
         *
         * @param hit		blast hit to check
         * @param otherId	ID of the desired target
         */
        public boolean targetEquals(BlastHit hit, String otherId) {
            return otherId.contentEquals(this.targetOf(hit));
        }

        /**
         * @return the sorting sequence data for the specified hit
         */
        public BlastHit.SeqData data(BlastHit hit) {
            return hit.getData(this.sortIdx);
        }

        /**
         * @return the target sequence data for the specified hit
         */
        public BlastHit.SeqData target(BlastHit hit) {
            return hit.getData(this.otherIdx);
        }

        /**
         * @return the plural phrase for the things being sorted
         */
        public String getPlural() {
            return this.plural;
        }

    }

    /**
     * This enum indicates the type of BLAST report.
     */
    public enum Type {
        TABLE, HTML, ALIGN;

        /**
         * Create a BLAST reporter of this type sorting by the specified sequence.
         *
         * @param stream	output stream to contain the report
         * @param type		type of sequence (SUBJECT, QUERY) by which to sort the report
         */
        public BlastReporter create(OutputStream stream, SortType type) {
            BlastReporter retVal = null;
            switch (this) {
            case TABLE:
                retVal = new BlastTableReporter(stream, type);
                break;
            case HTML:
                retVal = new BlastHtmlReporter(stream, type);
                break;
            case ALIGN:
                retVal = new BlastAlignReporter(stream, type);
                break;
            }
            return retVal;
        }
    }

    /**
     * This simple object describes the BLAST info.
     */
    public static class Info {

        // FIELDS
        /** BLAST parameters */
        private String parms;
        /** query count */
        private int queriesIn;
        /** number of misses */
        private int missCount;
        /** number of hits */
        private int hitCount;

        /**
         * Construct the information object.
         *
         * @param parms		parameters for the BLAST
         * @param queries	number of incoming queries
         * @param misses	number of queries without a hit
         * @param hits		total number of hits
         */
        public Info(String parms, int queries, int misses, int hits) {
            this.parms = parms;
            this.queriesIn = queries;
            this.missCount = misses;
            this.hitCount = hits;
        }

        /**
         * @return the BLAST parameters
         */
        public String getParms() {
            return parms;
        }

        /**
         * @return the input query count
         */
        public int getQueriesIn() {
            return queriesIn;
        }

        /**
         * @return the query miss count
         */
        public int getMissCount() {
            return missCount;
        }

        /**
         * @return the hit count
         */
        public int getHitCount() {
            return hitCount;
        }

    }

    // FIELDS
    /** message log */
    protected Logger log = LoggerFactory.getLogger(BlastReporter.class);
    /** target sort type (QUERY or SUBJECT) */
    private SortType sortType;
    /** map of sort sequences to hits */
    private SortedMap<String, List<BlastHit>> hitMap;
    /** rejection count */
    private int rejected;

    /**
     * Construct a new BLAST reporting facility.
     *
     * @param output	output stream to receive report
     * @param sort		type of sequence (subject or query) to sort on
     */
    public BlastReporter(OutputStream output, SortType sort) {
        super(output);
        this.sortType = sort;
        this.hitMap = new TreeMap<String, List<BlastHit>>(new NaturalSort());
        this.rejected = 0;

    }

    /**
     * Record a BLAST hit.  Hits wholly contained inside others are discarded.
     *
     * @param hit	BLAST hit to record
     */
    public void recordHit(BlastHit hit) {
        String id = this.sortType.idOf(hit);
        List<BlastHit> hitList = hitMap.get(id);
        if (hitList == null) {
            hitList = new Shuffler<BlastHit>(5).add1(hit);
            this.hitMap.put(id, hitList);
        } else {
            int i = 0;
            String otherId = this.sortType.targetOf(hit);
            Location otherLoc = this.sortType.target(hit).getLoc();
            // Find the hit for this query/subject pair.
            boolean keep = true;
            while (i < hitList.size() && keep) {
                if (! this.sortType.targetEquals(hitList.get(i), otherId)) {
                    // Different target sequence, keep looking.
                    i++;
                } else {
                    Location hitLoc = this.sortType.target(hitList.get(i)).getLoc();
                    if (hitLoc.contains(otherLoc)) {
                        // Old location subsumes this one, so we're done.
                        keep = false;
                        this.rejected++;
                    } else if (otherLoc.contains(hitLoc)) {
                        // This location subsumes old one, so remove it.
                        hitList.remove(i);
                        this.rejected++;
                    } else {
                        // Locations are different, so keep moving.
                        i++;
                    }
                }
            }
            if (keep)
                hitList.add(hit);
        }
    }

    /**
     * Output the report.
     *
     * @param title		title to put on report
     * @param runInfo 	statistics on the run
     */
    public void writeReport(String title, Info runInfo) {
        this.openReport(title);
        this.showSubtitle(runInfo);
        for (String id : this.hitMap.keySet()) {
            // Get the list of hits for this sequence.
            List<BlastHit> hitList = this.hitMap.get(id);
            this.openSection(this.sortType.data(hitList.get(0)));
            // Process each hit.
            for (BlastHit hit : hitList)
                this.processHit(this.sortType.target(hit),
                        this.sortType.data(hit), hit);
            // Finish up this section.
            this.closeSection();
        }
        // Finish the report.
        this.closeReport();
        log.info("{} hits removed due to overlap.", this.rejected);
    }

    /**
     * Display the subtitle.  (This method is optional.)
     *
     * @param runInfo	subtitle containing the BLAST parameters
     */
    protected void showSubtitle(Info runInfo) { }

    /**
     * Start the report with the specified title.
     *
     * @param title		title string for report
     */
    protected abstract void openReport(String title);

    /**
     * Begin a new section for the specified sort sequence.
     *
     * @param data	data describing the sequence anchoring this section
     */
    protected abstract void openSection(SeqData data);

    /**
     * Process a single hit.
     *
     * @param target	target sequence being hit by the sort sequence
     * @param anchor	anchor sequence hit by the sort sequence
     * @param hit		BLAST hit
     */
    protected abstract void processHit(SeqData target, SeqData anchor, BlastHit hit);

    /**
     * Finish the section for the current sort sequence.
     */
    protected abstract void closeSection();

    /**
     * Finish the entire report.
     */
    protected abstract void closeReport();

    /**
     * @return the sort type
     */
    protected SortType getSortType() {
        return sortType;
    }

    /**
     * @return the number of sort sequences hit
     */
    protected int getSequencesHit() {
        return this.hitMap.size();
    }


}
