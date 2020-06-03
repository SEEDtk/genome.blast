/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.io.GtoFilter;
import org.theseed.reports.MatchReporter;
import org.theseed.sequence.blast.BlastDB;

/**
 * This method runs the MatchProcessor against all rna/genome pairs in a directory.  It has
 * most of the sample command-line options (the exception being "--sample".  The two positional
 * parameters are the name of the profile directory and the name of the directory containing the
 * samples.  Each sample must consist of two files-- XXXXXX.assembled.fasta, which is the RNA sequence
 * FASTA file, and XXXXXX.gto, which is the target genome.  XXXXXX should be the ID of the sample.
 * The output files are XXXXXX.gti, which contains DNA from the genome and the associated proteins found
 * from the RNA sequences, and XXXXXX.faa, which contains the proteins and the location in each RNA they
 * came from.
 *
 * The command-line options are as follows.
 *
 * -v	show more detail on the log
 * -h	display usage information
 * -b	number of input sequences to process at a time (the default is 10)
 * -x	distance to extend the genome hit on either side (the default is 50)
 *
 * --maxE			maximum permissible e-value for a profile hit; the default is 1e-10
 * --minPct			minimum percent coverage of an incoming RNA sequence for genome hits; the default is 95.0
 * --minPctIdent	minimum percent identity for an RNA-to-genome hit; the default is 90.0
 * --tempDir		temporary directory for BLAST databases; the default is "Temp" in the current directory
 * --maxGap			maximum gap between adjacent hits when they are to be joined; the default is 500
 * --minIdent		minimum percent identity for profile hits; the default is 90.0
 * --minQIdent		minimum query identity fraction for profile hits; the default is 0.0
 * --minQbsc		minimum query-scaled bit score for profile hits; the default is 1.1
 * --minQuery		minimum percent query match for profile hits; the default is 65.0
 * --starts			algorithm for finding starts; the default is NEAREST
 * --verify			if specified, a verification report will be produced in the file "report.txt"
 * --missing		if specified, only samples without an existing GTI file will be processed
 *
 * @author Bruce Parrello
 *
 */
public class MatchRunProcessor extends MatchBaseProcessor {

    // COMMAND-LINE OPTIONS

    /** if specified, a verification report will be produced in the directory after the run */
    @Option(name = "--verify", usage = "if specified, a verification report will be produced")
    private boolean verify;

    /** if specified, only samples that do not already have a GTI file will be processed */
    @Option(name = "--missing", usage = "if specified, only new samples will be processed")
    private boolean missingMode;

    /** directory containing samples */
    @Argument(index = 1, usage = "directory containing sample FASTA files and genomes")
    private File inDir;

    @Override
    protected void setDefaults() {
        this.setupDefaults();
        this.verify = false;
        this.missingMode = false;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " not found or invalid.");
        this.validateCommonParms(MatchReporter.Type.SUMMARY, System.out);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        int sampCount = 0;
        // Insure the BLAST messages aren't too numerous.
        BlastDB.configureLogging("INFO");
        // Loop through the input directory.
        File[] gtoFiles = GtoFilter.getAll(this.inDir);
        log.info("{} GTO files found in {}.", gtoFiles.length, this.inDir);
        for (File gtoFile : gtoFiles) {
            // Find the corresponding RNA file.
            String sampleID = StringUtils.removeEnd(gtoFile.getName(), ".gto");
            File rnaFile = new File(this.inDir, sampleID + ".assembled.fasta");
            File outFile = new File(this.inDir, sampleID + ".gti");
            if (! rnaFile.canRead()) {
                log.warn("RNA file {} for {} not found or unreadable:  skipping.", rnaFile, sampleID);
            } else if (this.missingMode && outFile.exists()) {
                log.info("Sample {} skipped-- already processed.", sampleID);
            } else {
                this.runSample(this.inDir, gtoFile, sampleID, rnaFile);
                sampCount++;
            }
        }
        finish();
        log.info("{} samples processed.", sampCount);
        if (this.verify) {
            log.info("Producing verification reports.");
            MatchVerifyProcessor.verifyDirectory(this.inDir);
        }
    }

}
