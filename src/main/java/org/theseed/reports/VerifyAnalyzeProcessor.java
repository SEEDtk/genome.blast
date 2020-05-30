/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.FloatList;

/**
 * This command reads a tab-delimited file into memory with multiple named numeric columns and a single column that says "good" or "bad".
 * For each numeric column, a value will be found where a certain percentage of good results are above the value (sensitivity), and a value
 * where no more than a certain percentage of bad results are above the value (specificity).
 *
 * The positional parameters are the names of the numeric columns to which the analysis should be applied.  The input file should come from
 * the standard input.
 *
 * The command-line options are as follows.
 *
 * -h	display usage
 * -v	display more detailed progress messages
 * -i	input file (if not STDIN)
 * -c	column label for good/bad indicator (default "type")
 * -g	type value for "good" records (default "good")
 *
 * --thresholds		comma-delimited list of thresholds to display (default: 99,95,90,80)
 *
 * @author Bruce Parrello
 *
 */
public class VerifyAnalyzeProcessor extends BaseProcessor {

    // FIELDS
    /** logging object */
    protected static Logger log = LoggerFactory.getLogger(VerifyAnalyzeProcessor.class);
    /** main list of records */
    private AnalysisList records;
    /** list of thresholds */
    private FloatList thresholds;

    // COMMAND-LINE OPTIONS

    /** input file (if not STDIN) */
    @Option(name = "-i", aliases = { "--input" }, metaVar = "input.tbl", usage = "input file (if not STDIN)")
    private File inFile;

    /** column label for type column */
    @Option(name = "-c", aliases = { "--typeCol" }, metaVar = "qual", usage = "index (1-based) or label of type column")
    private String typeCol;

    /** type value for positive (good) records */
    @Option(name = "-g", aliases = { "--goodValue" }, metaVar = "1",
            usage = "value of a positive (good) record in the type column")
    private String goodValue;

    /** list of thresholds to target */
    @Option(name = "--thresholds", aliases = { "--levels" }, metaVar="80,50,40",
            usage = "comma-delimited list of thresholds to target")
    private void setThresholds(String tList) {
        this.thresholds = new FloatList(tList);
    }

    /** list of numeric column headings */
    @Argument(index = 0, metaVar = "evalue percent ...", usage = "labels of numeric columns to analyze")
    private List<String> numLabels;

    @Override
    protected void setDefaults() {
        this.thresholds = new FloatList(99, 95, 90, 80);
        this.inFile = null;
        this.typeCol = "type";
        this.goodValue = "good";
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Verify that all the thresholds are valid.
        for (double thresh : this.thresholds) {
            if (thresh <= 0.0 || thresh >= 100.0)
                throw new IllegalArgumentException("Invalid threshold " + Double.toString(thresh) +
                        ": must be between 0 and 100.");
        }
        // Read in the records.
        TabbedLineReader inStream = null;
        if (this.inFile == null) {
            log.info("Records to analyze will be read from standard input.");
            inStream = new TabbedLineReader(System.in);
        } else {
            log.info("Records to analyze will be read from {}.", this.inFile);
            inStream = new TabbedLineReader(this.inFile);
        }
        try {
            this.records = new AnalysisList(inStream, this.typeCol, this.numLabels, this.goodValue);
            log.info("{} records read: {} good, {} bad.", this.records.size(), this.records.getGoodCount(),
                    this.records.getBadCount());
        } finally {
            inStream.close();
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Write the output header.
        System.out.println("label\tcutoff\tspecificity\tsensitivity\taccuracy");
        // We process each number value separately.  We compute the accuracy cutoff, then the specificity and
        // sensitivity cutoffs for each threshold.  The cutoffs are output in order with all three values.
        for (int idx = 0; idx < this.numLabels.size(); idx++) {
            String label = this.numLabels.get(idx);
            log.info("Processing {} values.", label);
            // Sort by this label's values.
            this.records.sort(idx);
            // This will hold the cutoffs.
            SortedSet<Double> cutoffs = new TreeSet<Double>();
            // Compute the accuracy cutoff.
            double cutoff = this.records.getAccurate();
            cutoffs.add(cutoff);
            // Loop through the thresholds.
            for (double thresh : this.thresholds) {
                cutoff = this.records.getSensitive(thresh);
                cutoffs.add(cutoff);
                cutoff = this.records.getSpecific(thresh);
                if (Double.isNaN(cutoff))
                    log.debug("Could not find a cutoff for specificity threshold {} of {}.", thresh, label);
                else
                    cutoffs.add(cutoff);
            }
            // Output the data for these cutoffs.
            for (double cutoffVal : cutoffs) {
                System.out.format("%s\t%8.6f\t%4.1f\t%4.1f\t%4.1f%n", label, cutoffVal,
                        this.records.getSpecificity(cutoffVal, idx),
                        this.records.getSensitivity(cutoffVal, idx),
                        this.records.getAccuracy(cutoffVal, idx));
            }
        }

    }

}
