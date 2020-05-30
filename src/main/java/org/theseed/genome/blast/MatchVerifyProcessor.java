/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.EnumCounter;
import org.theseed.genome.Genome;
import org.theseed.io.GtoFilter;
import org.theseed.reports.MatchReporter;
import org.theseed.reports.MatchVerifyReporter;
import org.theseed.reports.MatchVerifyReporter.ErrorType;
import org.theseed.sequence.blast.GtiFile;
import org.theseed.utils.BaseProcessor;

/**
 * This command runs through the directory produced by the "mrun" command and creates a verification report
 * for all the GTI files found.  In addition, a "results.txt" file is produced that summarizes each sample.
 *
 * The profiles used in the match run are prokaryotic, so only prokaryotic genomes are fully verified.  Other
 * genomes will have basic counts displayed in results.txt, but no error reporting.
 *
 * The positional parameter is the name of the match-run directory.  Each GTI file must have a GTO file
 * with a corresponding name.  The first part of the file number must be the sample ID.
 *
 * @author Bruce Parrello
 *
 */
public class MatchVerifyProcessor extends BaseProcessor {

    // FIELDS
    /** log facility */
    protected Logger log = LoggerFactory.getLogger(MatchVerifyProcessor.class);
    /** error totals */
    private EnumCounter<ErrorType> totals;
    /** protein count */
    private int protTotal;
    /** gti count */
    private int gtiTotal;
    /** current sample ID */
    private String sampleId;
    /** current genome */
    private Genome genome;


    // COMMAND-LINE OPTIONS

    /** directory containing match-run output */
    @Argument(index = 0, metaVar = "runDir", usage = "directory containing GTI and GTO files")
    private File runDir;

    @Override
    protected void setDefaults() {	}

    @Override
    protected boolean validateParms() throws IOException {
        // Verify the run directory.
        if (! this.runDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.runDir + " not found or invalid.");
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // This will be the result summary file.
        File resultFile = new File(this.runDir, "results.txt");
        // Initialize the totals.
        this.totals = new EnumCounter<ErrorType>(ErrorType.class);
        this.gtiTotal = 0;
        this.protTotal = 0;
        // Open the files and produce the reports.
        try (MatchVerifyReporter reporter = (MatchVerifyReporter) MatchReporter.Type.VERIFY.create(System.out);
                PrintWriter resultStream = new PrintWriter(resultFile))  {
            // Output the report headers.
            reporter.initialize();
            resultStream.println("sample\tgenome_id\tgenome_name\tgc\trecords\tproteins\tgood\tshorter\tlonger\tchanged\tinvalid");
            // Now we loop through each of the genomes in the input directory.  We can't use a genome directory
            // here because the file name is the sample ID.
            File[] gFiles = GtoFilter.getAll(this.runDir);
            log.info("{} genomes found in {}.", gFiles.length, this.runDir);
            for (File gFile : gFiles) {
                // Get the GTI file and the genome.
                this.sampleId = StringUtils.removeEnd(gFile.getName(), ".gto");
                File gtiFile = new File(this.runDir, sampleId + ".gti");
                log.info("Loading genome from {}.", gFile);
                this.genome = new Genome(gFile);
                log.info("Genome is {}.", genome);
                switch (genome.getDomain()) {
                case "Archaea" :
                case "Bacteria" :
                    this.verifySample(gtiFile, reporter, resultStream);
                    break;
                default :
                    this.countSample(gtiFile, resultStream);
                }
            }
            // Write the totals.
            resultStream.println();
            resultStream.format("TOTALS\t \t \t \t%d\t%d\t%d\t%d\t%d\t%d\t%d%n", this.gtiTotal, this.protTotal,
                    this.totals.getCount(ErrorType.EXACT), this.totals.getCount(ErrorType.TOO_SHORT),
                            this.totals.getCount(ErrorType.TOO_LONG), this.totals.getCount(ErrorType.CHANGED),
                            this.totals.getCount(ErrorType.NOT_FOUND));
            log.info("All done.");
        }
    }

    /**
     * Process a sample with unreliable proteins.  In this case, we simply count the records in the GTI file and
     * produce a line in the results.txt file.
     *
     * @param gtiFile		input GTI file
     * @param resultStream	output stream for results.txt
     *
     * @throws IOException
     */
    private void countSample(File gtiFile, PrintWriter resultStream) throws IOException {
        log.info("Counting proteins in {}.", gtiFile);
        try (GtiFile gtiStream = new GtiFile(gtiFile)) {
            int gtiCount = 0;
            int protCount = 0;
            for (GtiFile.Record record : gtiStream) {
                gtiCount++;
                protCount += record.getProts().size();
            }
            showResult(resultStream, gtiCount, protCount, null);
        }
    }

    /**
     * Process a sample with reliable proteins.  Here we produce a full report on all the records and proteins.
     *
     * @param gtiFile		input file containing GTIs
     * @param reporter		output reporting object
     * @param resultStream	output stream for results.txt
     *
     * @throws IOException
     * @throws InterruptedException
     */
    private void verifySample(File gtiFile, MatchVerifyReporter reporter, PrintWriter resultStream)
            throws IOException, InterruptedException {
        log.info("Verifying proteins in {}.", gtiFile);
        try (GtiFile gtiStream = new GtiFile(gtiFile)) {
            // Initialize the current section of the report.
            reporter.clearCounters();
            reporter.startSection(this.genome, this.sampleId);
            // We report on each record, then we write the totals to the result file.
            for (GtiFile.Record record : gtiStream) {
                reporter.processSequence(record.getRnaId(), record.getDnaLoc(), record.getDna(), record.getProts());
            }
            showResult(resultStream, reporter.getRecordCount(), reporter.getProteinCount(), reporter.getCounters());
        }
    }

    /**
     * Display a report line and update the totals.
     *
     * @param recordCount
     * @param proteinCount
     * @param counters
     */
    private void showResult(PrintWriter resultStream, int recordCount, int proteinCount, EnumCounter<ErrorType> counters) {
        // Process the counters.
        String counterString;
        if (counters == null) {
            counterString = " \t \t \t \t \t ";
        } else {
            // Create the counter output text.
            counterString = String.format("%d\t%d\t%d\t%d\t%d", counters.getCount(ErrorType.EXACT),
                    counters.getCount(ErrorType.TOO_SHORT), counters.getCount(ErrorType.TOO_LONG),
                    counters.getCount(ErrorType.CHANGED), counters.getCount(ErrorType.NOT_FOUND));
            // Accumulate the counters in the totals.
            this.totals.sum(counters);
        }
        this.gtiTotal += recordCount;
        this.protTotal += proteinCount;
        resultStream.format("%s\t%s\t%s\t%d\t%d\t%d\t%s%n", this.sampleId, this.genome.getId(),
                this.genome.getName(), this.genome.getGeneticCode(), recordCount, proteinCount, counterString);
    }

}
