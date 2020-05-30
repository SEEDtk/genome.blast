/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.theseed.io.GtoFilter;

/**
 * This method runs the MatchProcessor against all rna/genome pairs in a directory.  It has
 * most of the sample command-line options (the exception being "--sample".  The two positional
 * parameters are the name of the profile directory and the name of the directory containing the
 * samples.  Each sample must consist of two files-- XXXXXX.assembled.fasta, which is the RNA sequence
 * FASTA file, and XXXXXX.gto, which is the target genome.  XXXXXX should be the ID of the sample.
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
 * --format			output format
 * --starts			algorithm for finding starts; the default is NEAREST
 *
 * @author Bruce Parrello
 *
 */
public class MatchRunProcessor extends MatchBaseProcessor {

    // COMMAND-LINE OPTIONS

    /** directory containing samples */
    @Argument(index = 1, usage = "directory containing sample FASTA files and genomes")
    private File inDir;

    @Override
    protected void setDefaults() {
        this.setupDefaults();
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " not found or invalid.");
        this.validateCommonParms();
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        int sampCount = 0;
        // Loop through the input directory.
        File[] gtoFiles = GtoFilter.getAll(this.inDir);
        log.info("{} GTO files found in {}.", gtoFiles.length, this.inDir);
        for (File gtoFile : gtoFiles) {
            // Find the corresponding RNA file.
            String sampleID = StringUtils.removeEnd(gtoFile.getName(), ".gto");
            File rnaFile = new File(this.inDir, sampleID + ".assembled.fasta");
            if (! rnaFile.canRead()) {
                log.warn("RNA file {} for {} not found or unreadable:  skipping.", rnaFile, sampleID);
            } else {
                // Compute the output file.
                File outFile = new File(this.inDir, sampleID + ".gti");
                OutputStream outStream = new FileOutputStream(outFile);
                // Set up and run the sample.
                long start = System.currentTimeMillis();
                this.setup(gtoFile, rnaFile, outStream);
                log.info("Processing RNA file {} for genome {}.", rnaFile, super.getGenome());
                this.runGenome(sampleID, rnaFile);
                sampCount++;
                log.info("{} samples processed.  {} took {} seconds.", sampCount, sampleID,
                        (System.currentTimeMillis() - start + 500) / 1000);
            }
        }
        finish();
        log.info("All done.  {} samples processed.", sampCount);
    }

}
