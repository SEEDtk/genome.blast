/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.util.HashMap;
import java.util.Map;

import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.ProteinBlastDB;
import org.theseed.utils.StringPair;

/**
 * This performs protein comparisons using BLAST.   We blast the first FASTA against the second, and
 * keep the top result for each sequence pair.  The bit score over the length is returned as a similarity
 * score.
 *
 * @author Bruce Parrello
 *
 */
public class BlastProteinCompare extends ProteinCompare {

    // FIELDS
    /** blast parameters to use */
    private BlastParms parms;

    public BlastProteinCompare(IParms processor) {
        super(processor);
        this.parms = new BlastParms();
        this.parms.maxPerPair(1).maxE(1e-1).pctLenOfQuery(50);
    }

    @Override
    protected Map<StringPair, Double> computeSimInternal(File file1, File file2) throws IOException {
        // Insure there is a database for the second FASTA.
        BlastDB db2;
        synchronized (this) {
            log.info("Creating blast database for {}.", file2);
            try {
                db2 = ProteinBlastDB.createOrLoad(file2);
            } catch (InterruptedException e) {
                throw new RuntimeException("Error creating BLAST DB: " +e.toString());
            }
        }
        ProteinStream inStream = new ProteinInputStream(file1);
        log.info("Performing blast between {} and {}.", file1, file2);
        var hits = db2.blast(inStream, this.parms);
        final int nHits = hits.size();
        log.info("{} results from blast between {} and {}.", nHits, file1, file2);
        // We will store the hits in here.
        int count = 0;
        long start = System.currentTimeMillis();
        long timer = start;
        Map<StringPair, Double> retVal = new HashMap<StringPair, Double>(hits.size() * 4 / 3);
        for (var hit : hits) {
            String id1 = hit.getQueryId();
            String id2 = hit.getSubjectId();
            count++;
            if (! id1.equals(id2)) {
                StringPair pair = new StringPair(id1, id2);
                double maxLen = Math.max(hit.getQueryLen(), hit.getSubjectLen());
                double score = hit.getBitScore() / maxLen;
                double oldScore = retVal.getOrDefault(pair, 0.0);
                if (score > oldScore) {
                    retVal.put(pair, score);
                    if (log.isInfoEnabled() && System.currentTimeMillis() - timer >= 5000) {
                        timer = System.currentTimeMillis();
                        Duration d = Duration.ofMillis(timer - start);
                        log.info("{} of {} hits processed for {} and {} in {}.", count, nHits, file1, file2, d);
                    }
                }
            }
        }
        return retVal;
    }

}
