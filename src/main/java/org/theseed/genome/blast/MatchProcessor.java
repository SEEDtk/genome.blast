/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import org.kohsuke.args4j.Argument;
import org.theseed.basic.ParseFailureException;
import org.theseed.reports.MatchReporter;
import org.theseed.sequence.blast.BlastDB;

/**
 * This is designed to find ground-truth proteins from RNA sequences that correspond to a known genome.
 *
 * The RNA sequences are profiled to isolate protein locations.  Each protein location is extended to a
 * start and a stop and translated to produce actual proteins.  Adjacent proteins are grouped together,
 * and then the RNA that produced them is blasted against the genome contigs to find the source location.
 * The source location DNA is then output along with the protein sequences.
 *
 * The positional parameters are the name of the profile directory, the name of the input directory,
 * and the names of one or more samples to process.  Each sample XXXXXXXX must have a GTO file with
 * the name XXXXXXXX.gto and an RNA FASTA file with the name XXXXXXXX.assembled.fasta.
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
 * --minIdent		minimum percent identity for genome hits; the default is 90.0
 * --minQIdent		minimum query identity fraction for profile hits; the default is 0.0
 * --minQbsc		minimum query-scaled bit score for profile hits; the default is 1.1
 * --minQuery		minimum percent query match for profile hits; the default is 65.0
 * --starts			algorithm for finding starts; the default is NEAREST
 *
 * @author Bruce Parrello
 *
 */
public class MatchProcessor extends MatchBaseProcessor {

    /**
     * This is a dinky little object that contains the names of the files for a sample.
     */
    protected class SampleFiles {
        public File rnaFile;
        public File gtoFile;

        public SampleFiles(String sample) {
            this.rnaFile = new File(inDir, sample + ".assembled.fasta");
            this.gtoFile = new File(inDir, sample + ".gto");
        }
    }

    // FIELDS
    /** list of samples to process */
    private SortedMap<String, SampleFiles> sampleMap;

    // COMMAND-LINE OPTION

    /** RNA sequence input file */
    @Argument(index = 1, metaVar = "inDir", usage = "input directory containing samples", required = true)
    private File inDir;

    /** IDs of samples to process */
    @Argument(index = 2, metaVar = "sample1 sample2 ...", usage = "IDs of samples to process",
            multiValued = true)
    private List<String> samples;


    @Override
    protected void setDefaults() {
        this.setupDefaults();
        this.sampleMap = new TreeMap<>();
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        validateCommonParms(MatchReporter.Type.SUMMARY, System.out);
        // Verify all the input samples.
        for (String sample : samples) {
            SampleFiles files = this.new SampleFiles(sample);
            if (! files.gtoFile.canRead())
                throw new FileNotFoundException("GTO file " + files.gtoFile + " not found or unreadable.");
            if (! files.rnaFile.canRead())
                throw new FileNotFoundException("RNA file " + files.rnaFile + " not found or unreadable.");
            this.sampleMap.put(sample, files);
        }
    }

    @Override
    public void runCommand() throws Exception {
        // Turn off BLAST details.
        BlastDB.configureLogging("INFO");
        // Run all the samples.
        for (Map.Entry<String, SampleFiles> sampleEntry : this.sampleMap.entrySet()) {
            SampleFiles files = sampleEntry.getValue();
            String sampleId = sampleEntry.getKey();
            this.runSample(this.inDir, files.gtoFile, sampleId, files.rnaFile);
        }
        finish();
    }


}
