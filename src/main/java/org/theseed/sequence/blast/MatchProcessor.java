/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;

/**
 * This is designed to find ground-truth proteins from RNA sequences that correspond to a known genome.
 *
 * The RNA sequences are profiled to isolate protein locations.  Each protein location is extended to a
 * start and a stop and translated to produce actual proteins.  Adjacent proteins are grouped together,
 * and then the RNA that produced them is blasted against the genome contigs to find the source location.
 * The source location DNA is then output along with the protein sequences.
 *
 * The positional parameters are the name of the profile directory, the name of the RNA sequence file,
 * and the name of the genome file.
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
 * --format			output format
 * --starts			algorithm for finding starts; the default is NEAREST
 * --sample			ID of the sample that produced these RNAs
 *
 * @author Bruce Parrello
 *
 */
public class MatchProcessor extends MatchBaseProcessor {




    // COMMAND-LINE OPTION

    /** ID of the sample that produced the RNAs */
    @Option(name = "--sample", metaVar="SRR134798", usage = "ID of source sample, to be displayed in report")
    private String sampleID;

    /** RNA sequence input file */
    @Argument(index = 1, metaVar = "rna.fasta", usage = "RNA sequence input file", required = true)
    private File inFile;

    /** file containing the target genome */
    @Argument(index = 2, metaVar = "genome.gto", usage = "target genome file containing the proteins",
            required = true)
    private File genomeFile;

    @Override
    protected void setDefaults() {
        this.sampleID = "N/K";
    }

    @Override
    protected boolean validateParms() throws IOException {
        setup(this.genomeFile, this.inFile, System.out);
        validateCommonParms();
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        runGenome(this.sampleID, this.inFile);
    }


}
