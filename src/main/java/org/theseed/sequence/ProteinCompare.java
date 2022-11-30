package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.utils.StringPair;

/**
 * This is a base class for protein comparisons.  The client presents two protein FASTA files, and
 * we come up with a similarity for each pair between the left and right.  The results are returned
 * in a list.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ProteinCompare {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinCompare.class);


    /**
     * This defines the parameters that must be supplied by any client doing protein comparisons.
     */
    public interface IParms {

        /**
         * @return the size of each kmer to use
         */
        int getKmerSize();

    }

    /**
     * This enum contains the supported protein comparison types.
     */
    public enum Type {
        /** compare using Jaccard distance and kmers */
        KMERS {
            @Override
            public ProteinCompare create(IParms processor) {
                return new KmerProteinCompare(processor);
            }
        },
        /** compare using BLAST scores */
        BLAST {
            @Override
            public ProteinCompare create(IParms processor) {
                return new BlastProteinCompare(processor);
            }
        };

        public abstract ProteinCompare create(IParms processor);
    }

    /**
     * Construct a protein comparison method.
     *
     * @param processor		controlling command processor
     */
    public ProteinCompare(IParms processor) {
    }

    /**
     * Compute the similarities between proteins in two FASTA files.
     *
     * @param file1		first file to compare
     * @param file2		second file to compare
     *
     * @return a map from ID pairs to distances
     *
     * @throws IOException
     */
    public Map<StringPair, Double> computeSim(File file1, File file2) throws IOException {
        log.info("Comparing proteins in {} with {}.", file1, file2);
        return this.computeSimInternal(file1, file2);
    }

    /**
     * Compute the similarities for all the protein pairs in the two specified sets.
     *
     * @param file1		first file to compare
     * @param file2		second file to compare
     *
     * @return a map from ID pairs to non-trivial similarities
     *
     * @throws IOException
     */
    protected abstract Map<StringPair, Double> computeSimInternal(File file1, File file2) throws IOException;

}
