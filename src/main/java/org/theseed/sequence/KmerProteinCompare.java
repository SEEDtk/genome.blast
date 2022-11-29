/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.theseed.utils.StringPair;

/**
 * This method compares proteins using Kmer distance.  The kmer size is a parameter.
 *
 * @author Bruce Parrello
 *
 */
public class KmerProteinCompare extends ProteinCompare {

    // FIELDS
    /** kmer size */
    private int kmerSize;

    public KmerProteinCompare(IParms processor) {
        super(processor);
        this.kmerSize = processor.getKmerSize();
    }

    @Override
    protected Map<StringPair, Double> computeDistanceInternal(File file1, File file2) throws IOException {
        // Read in the first set of proteins.
        List<Sequence> prots1 = FastaInputStream.readAll(file1);
        log.info("{} sequences read from {}.", prots1.size(), file1);
        List<Sequence> prots2;
        int numCompares;
        // We do a little optimization for recursive comparisons, which occur a lot in the first client
        // that used this facility.
        boolean recursive = (file1.equals(file2));
        if (recursive) {
            prots2 = prots1;
            numCompares = prots1.size() * prots2.size() / 2;
        } else {
            prots2 = FastaInputStream.readAll(file2);
            log.info("{} sequences read from {}.", prots2.size(), file2);
            numCompares = prots1.size() * prots2.size();
        }
        // Create kmers for all the second-list proteins.
        Map<String, ProteinKmers> kmerMap2 = new HashMap<String, ProteinKmers>(prots2.size() * 4 / 3);
        prots2.stream().forEach(x -> kmerMap2.put(x.getLabel(), new ProteinKmers(x.getSequence(), this.kmerSize)));
        // Create the result hash.
        var retVal = new HashMap<StringPair, Double>(numCompares * 4 / 3 + 1);
        // Run through all the comparisons.
        int count = 0;
        long start = System.currentTimeMillis();
        for (Sequence seq1 : prots1) {
            // Note we re-use the kmers in the map if we can find then.  This is another optimization
            // for the recursive case.
            String id1 = seq1.getLabel();
            ProteinKmers kmers1 = kmerMap2.getOrDefault(id1, new ProteinKmers(seq1.getSequence(), this.kmerSize));
            for (var entry2 : kmerMap2.entrySet()) {
                String id2 = entry2.getKey();
                StringPair pair = new StringPair(id1, id2);
                // Only proceed if we haven't done this comparison yet.
                if (! id1.equals(id2) && ! retVal.containsKey(pair)) {
                    ProteinKmers kmers2 = entry2.getValue();
                    double dist = kmers1.distance(kmers2);
                    retVal.put(pair, dist);
                    count++;
                }
                if (log.isInfoEnabled() && count % 10000 == 0  && count > 0) {
                    Duration d = Duration.ofMillis(System.currentTimeMillis() - start);
                    log.info("{} pairs processed out of {} for {} and {} in {}.", count, numCompares, file1, file2, d);
                }
            }
        }
        return retVal;
    }

}
