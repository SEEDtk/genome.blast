/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPInputStream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.LineReader;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

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
 * -p	is specified, the fraction of sequences to be written to the output; the default is 1.0, so all sequences will be kept;
 * 		a randomizer is used, so this will only be an estimate
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
    /** random number generator */
    private Random rand;

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

    /** fraction of sequences to pick */
    @Option(name = "-p", aliases = { "--fraction", "--pickFraction" }, metaVar = "0.10", usage = "fraction of sequences to output")
    private double pickFraction;

    /** input files */
    @Argument(index = 0, metaVar = "file1.fq file2.fq ...", usage = "input FASTQ files", multiValued = true)
    private List<File> inFiles;

    @Override
    protected void setDefaults() {
        this.minQual = 20.0;
        this.minSim = 25;
        this.filterFile = null;
        this.kmerFilter = null;
        this.pickFraction = 1.0;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (this.minQual < 0.0 || this.minQual > 40.0)
            throw new ParseFailureException("Minimum quality must be between 0.0 and 4.0 inclusive.");
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
        // Validate the pick fraction.
        if (this.pickFraction <= 0.0 || this.pickFraction > 1.0)
            throw new ParseFailureException("Pick fraction must be greater than 0 and not greater than 1.");
        this.rand = new Random();
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
                // Determine the file type.  If it's .gz we need to decompress.
                try (InputStream byteStream = new FileInputStream(inFile)) {
                    InputStream dataStream;
                    if (inFile.getName().endsWith(".gz"))
                        dataStream = new GZIPInputStream(byteStream);
                    else
                        dataStream = byteStream;
                    try (LineReader inStream = new LineReader(dataStream)) {
                        // Loop through the input lines.
                        Iterator<String> iter = inStream.iterator();
                        // Each output sequence is four input lines.
                        while (iter.hasNext()) {
                            // Skip first marker label.
                            iter.next();
                            // Get the sequence.
                            String sequence = iter.next();
                            // Skip second marker label.
                            iter.next();
                            // Get the quality.
                            double qual = this.parseQuality(iter.next());
                            // Now we filter for quality and kmers.
                            boolean good = true;
                            if (qual < this.minQual)
                                good = false;
                            else if (this.kmerFilter != null) {
                                // Here we have to check the kmer distance.
                                DnaKmers kmers = new DnaKmers(sequence);
                                good = (this.kmerFilter.similarity(kmers) >= this.minSim);
                            }
                            if (! good || this.rand.nextDouble() > this.pickFraction)
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
