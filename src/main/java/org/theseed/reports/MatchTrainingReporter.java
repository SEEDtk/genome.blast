/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Genome;
import org.theseed.sequence.ProteinDataStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.MatchProcessor;
import org.theseed.sequence.blast.ProteinBlastDB;

/**
 * This report produces a training set from the RNA sequence matching.  The training set can be used to
 * identify RNA starts.  The RNA proteins are blasted against the protein sequences in the genome.
 * For each RNA protein that matches at least one genome protein, we mark it as "coding".  For each
 * other protein, we mark it as "other".  The upstream DNA is used as the training set data.
 *
 * @author Bruce Parrello
 *
 */
public class MatchTrainingReporter extends MatchReporter {

    // FIELDS
    /** BLAST parameters */
    private BlastParms parms;
    /** genome protein database */
    private BlastDB protDB;
    /** TRUE if headers have been output */
    private boolean headerFlag;

    /**
     * Construct the report.
     *
     * @param output	output stream
     * @param genome	base genome
     */
    public MatchTrainingReporter(OutputStream output, Genome genome) {
        super(output, genome);
    }

    @Override
    public void initialize(MatchProcessor base) throws IOException, InterruptedException {
        // Create the BLAST database for verifying the proteins.
        File tempFile = File.createTempFile("pegs", "fa");
        this.protDB = ProteinBlastDB.create(tempFile, this.getGenome());
        this.protDB.deleteOnExit();
        // Create the BLAST parameters.
        this.parms = new BlastParms().maxE(1e-30).minPercent(95.0).maxPerQuery(1);
        // Denote we have not written the column headers.
        this.headerFlag = false;

    }

    @Override
    public void processSequence(String id, String dna, List<Sequence> prots) throws IOException, InterruptedException {
        // Insure we have a header.
        if (! this.headerFlag) {
            String upstream = prots.get(0).getComment();
            String header = "label\thit\t" + IntStream.rangeClosed(1, upstream.length())
                    .mapToObj(i -> String.format("p%d", i)).collect(Collectors.joining("\t"));
            this.println(header);
            this.headerFlag = true;
        }
        // Now we must verify the proteins to determine the good ones.
        ProteinDataStream protStream = new ProteinDataStream(prots);
        Map<String, List<BlastHit>> hitMap = BlastHit.sort(protStream.blast(this.protDB, this.parms));
        // Loop through the proteins, generating the training rows.
        for (Sequence prot : prots) {
            String label = prot.getLabel();
            String type = (hitMap.containsKey(label) ? "coding" : "other");
            String line = label + "\t" + type + "\t" + StringUtils.join(prot.getComment().toCharArray(), '\t');
            this.println(line);
        }
    }

    @Override
    public void finish() {
    }

}
