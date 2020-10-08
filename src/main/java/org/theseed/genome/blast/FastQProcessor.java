/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.LineReader;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseProcessor;

/**
 * This is a simple command that converts one or more FASTQ files to a single FASTA file.  This is a crude process that simply outputs
 * the raw sequences with a computed label and a quality estimate in the comments.
 *
 * The positional parameters are the names of the incoming FASTQ files.  The FASTA file will be produced on the standard output.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	show more detailed progress messages
 * -m	minimum mean quality for a read to be acceptable; the default is 20
 * -k	minimum number of DNA kmers in common required for a sequence to be kept; the default is 25
 * -s	if specified, the name of a FASTA file containing a DNA sequence; only incoming sequences that have kmers in common will be kept
 *
 * @author Bruce Parrello
 *
 */
public class FastQProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FastQProcessor.class);
    /** current sequence number, for generating IDs */
    private int seqNum;
    /** kmer object for sequence filtering (if any) */
    private DnaKmers kmerFilter;

    // COMMAND-LINE OPTIONS

    /** minimum mean quality */
    @Option(name = "-m", aliases = { "--min" }, metaVar = "20.0", usage = "minimum mean quality number for acceptable reads")
    private double minQual;

    /** minimum kmers in common when sequence-filtering */
    @Option(name = "-k", aliases = { "--ksim" }, metaVar = "20", usage = "minimum number of kmers-in-common required during sequence filtering")
    private int minSim;

    /** file containing optional sequence for sequence filtering */
    @Option(name = "-s", aliases = { "--seqFilter" }, metaVar = "seq.fa", usage = "optional FASTA containing a sequence for filtering")
    private File filterFile;

    /** input files */
    @Argument(index = 0, metaVar = "file1.fq file2.fq ...", usage = "input FASTQ files", multiValued = true)
    private List<File> inFiles;

    @Override
    protected void setDefaults() {
        this.minQual = 20.0;
        this.minSim = 25;
        this.filterFile = null;
        this.kmerFilter = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
        if (this.minQual < 0.0 || this.minQual > 40.0)
            throw new IllegalArgumentException("Minimum quality must be between 0.0 and 4.0 inclusive.");
        // Set up the sequence filter.
        if (this.filterFile != null) {
            if (! this.filterFile.canRead())
                throw new FileNotFoundException("Filter file " + this.filterFile + " is not found or unreadable.");
            try (FastaInputStream filterStream = new FastaInputStream(this.filterFile)) {
                Sequence seq = filterStream.next();
                this.kmerFilter = new DnaKmers(seq.getSequence());
                log.info("Filter sequence read from {}.", this.filterFile);
            }
        }
        // Verify the input files.
        for (File inFile : this.inFiles) {
            if (! inFile.canRead())
                throw new FileNotFoundException("Input file " + inFile + " is not found or unreadable.");
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Set up the output.
        try (FastaOutputStream outStream = new FastaOutputStream(System.out)) {
            // Loop through the input files.
            for (File inFile : this.inFiles) {
                log.info("Processing input file {}.", inFile);
                // These counters will be used for progress messages.
                int lineOut = 0;
                int lineSkipped = 0;
                int lineRead = 0;
                // Loop through the file itself.
                try (LineReader inStream = new LineReader(inFile)) {
                    Iterator<String> iter = inStream.iterator();
                    // Each output sequence is four input lines.  We skip the marker line.
                    while (iter.hasNext()) {
                        iter.next();
                        String sequence = iter.next();
                        iter.next();
                        double qual = this.parseQuality(iter.next());
                        boolean good = true;
                        if (qual < this.minQual)
                            good = false;
                        else if (this.kmerFilter != null) {
                            // Here we have to check the kmer distance.
                            DnaKmers kmers = new DnaKmers(sequence);
                            good = (this.kmerFilter.similarity(kmers) >= this.minSim);
                         }
                          if (! good)
                              lineSkipped++;
                          else {
                            String label = String.format("Sequence.%06d", ++this.seqNum);
                            Sequence seq = new Sequence(label, String.format("mean_q=%1.4f", qual), sequence);
                            outStream.write(seq);
                            lineOut++;
                        }
                        lineRead++;
                        if (log.isInfoEnabled() && lineRead % 50000 == 0)
                            log.info("{} sequences written, {} skipped.", lineOut, lineSkipped);
                    }
                }
            }
        }
    }

    /**
     * @return the mean quality from the current quality line
     *
     * @param qline		quality line to parse
     */
    private double parseQuality(String qline) {
        int sum = qline.chars().sum();
        double retVal = (sum - (qline.length() * 33)) / (double) qline.length();
        return retVal;
    }

}
