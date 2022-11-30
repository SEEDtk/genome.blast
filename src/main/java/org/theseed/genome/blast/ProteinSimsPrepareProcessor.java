/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.p3api.P3Genome;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command will prepare the genomes from a genome source so that they can be processed by
 * the ProteinSimsProcessor.  All the proteins will be collected and output into FASTA files
 * with the MD5 as the sequence label.  A command-line option determines the number of sequences
 * put in each file.
 *
 * The positional parameters are the name of the genome source (directory or file) and the name
 * of the output directory.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -b	number of protein sequences to put in each FASTA file (default 1200000)
 *
 * --source		type of genome source (default DIR)
 * --clear		erase the target directory before beginning
 *
 * @author Bruce Parrello
 *
 */
public class ProteinSimsPrepareProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinSimsPrepareProcessor.class);
    /** incoming genome source */
    private GenomeSource genomes;
    /** set of proteins already output */
    private Set<String> prots;
    /** number of protein files output */
    private int fileIdx;

    // COMMAND-LINE OPTIONS

    /** type of input */
    @Option(name = "--source", usage = "type of genome input (master directory, normal, PATRIC ID file)")
    private GenomeSource.Type inType;

    /** optimal number of sequences per FASTA file */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "4000", usage = "maximum number of protein sequences per file")
    private int batchSize;

    /** TRUE to clear the output directory before beginning */
    @Option(name = "--clear", usage = "if specified, the output directory will be cleared before processing")
    private boolean clearFlag;

    /** genome source directory */
    @Argument(index = 0, metaVar = "inDir", usage = "genome source input file or directory")
    private File inDir;

    /** output directory */
    @Argument(index = 1, metaVar = "outDir", usage = "output directory for FASTA files")
    private File outDir;

    @Override
    protected void setDefaults() {
        this.inType = GenomeSource.Type.DIR;
        this.batchSize = 1200000;
        this.clearFlag = false;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Verify the batch size.
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be at least 1 (but >= 1000 is better).");
        // Validate the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        } else
            log.info("Output will be written to {}.", this.outDir);
        // Validate the input.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " does not exist.");
        log.info("Connecting to genome source {}.", this.inDir);
        this.genomes = this.inType.create(this.inDir);
        this.genomes.setDetailLevel(P3Genome.Details.PROTEINS);
        // Create the protein set.  This prevents us from writing redundant proteins.
        this.prots = new HashSet<String>(5500 * this.genomes.size());
        // Denote no files have been output.
        this.fileIdx = 0;
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // We loop through the genome source, creating FASTA files. Generate the first output file.
        FastaOutputStream outStream = null;
        try {
            outStream = this.getNextFile();
            // Track the number of proteins in this file.
            int protCount = 0;
            // We need counters for the genomes, the missing proteins, and the proteins found.
            int gCount = 0;
            int missing = 0;
            int found = 0;
            int skipped = 0;
            // Here is the main genome loop.
            for (Genome genome : this.genomes) {
                gCount++;
                log.info("Processing genome {} of {}: {}", gCount, this.genomes.size(), genome);
                for (Feature feat : genome.getPegs()) {
                    String protMd5 = feat.getMD5();
                    if (protMd5 == null)
                        missing++;
                    else if (this.prots.contains(protMd5))
                        skipped++;
                    else {
                        // Here we have a new protein to output.  Verify there is room in this file.
                        if (protCount >= this.batchSize) {
                            outStream.close();
                            outStream = this.getNextFile();
                            protCount = 0;
                        }
                        // Write out the protein.
                        Sequence seq = new Sequence(protMd5, "", feat.getProteinTranslation());
                        outStream.write(seq);
                        found++;
                        protCount++;
                        // Insure we don't write it again.
                        this.prots.add(protMd5);
                    }
                }
                log.info("{} genomes processed.  {} proteins output, {} redundant, {} missing.",
                        gCount, found, skipped, missing);
            }
        } finally {
            // Close the last output file.
            if (outStream != null)
                outStream.close();
        }

    }

    /**
     * @return an output stream for a new FASTA file
     *
     * @throws FileNotFoundException
     */
    private FastaOutputStream getNextFile() throws FileNotFoundException {
        this.fileIdx++;
        String fileName = String.format("proteins%04d.faa", this.fileIdx);
        FastaOutputStream retVal = new FastaOutputStream(new File(this.outDir, fileName));
        return retVal;
    }

}
