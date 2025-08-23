/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.contains;
import static org.hamcrest.Matchers.equalTo;
import org.theseed.locations.Location;
import org.theseed.sequence.blast.GtiFile;

import junit.framework.TestCase;

/**
 * Test the GTI file.
 *
 * @author Bruce Parrello
 *
 */
public class TestGtiFile extends TestCase {

    /**
     * test small GTI file
     *
     * @throws IOException
     */
    public void testFile() throws IOException {
        File inFile = new File("data", "test.gti");
        try (GtiFile inStream = new GtiFile(inFile)) {
            Iterator<GtiFile.Record> iter = inStream.iterator();
            GtiFile.Record record = iter.next();
            assertTrue(iter.hasNext());
            assertThat(record.getSampleId(), equalTo("SRR10210920"));
            assertThat(record.getRnaId(), equalTo("r.TRINITY_GG_131_c0_g1_i1.0001"));
            assertThat(record.getDna(), equalTo("acaacattacgtttgatatcgtccgatatcgatttactatttccatttagacgttgttcaagctcgtgacaattcccaagtttttggtgttgctagaatttacgcttctttcaacgatactttcgttcatgttaccgatttatctggtaaggaaaccatcgccagagttactggtggtatgaaggttaaggctgacagagatgaatcttctccatacgctgctatgttggctgcccaagatgttgccgctaagtgtaaggaagtcggtatcactgccgttcacgttaagatcagagctaccggtggtactagaaccaagactccaggtccaggtggtcaagctgctttgagagctttggccagatctggtttgagaattggccgtatcgaagatgttaccccagttccatctgactccaccagaaagaagggtggtagaagaggtagaagattatgagttatgcatgtattgtacttgtattgccgtattattttttacagttaaaa"));
            Location loc = record.getDnaLoc();
            Location loc2 = Location.create("559292.28.con.0003", "-", 177450, 177956);
            assertThat(loc, equalTo(loc2));
            List<String> prots = record.getProts();
            assertThat(prots.size(), equalTo(1));
            assertThat(prots.get(0), equalTo("MSNVVQARDNSQVFGVARIYASFNDTFVHVTDLSGKETIARVTGGMKVKADRDESSPYAAMLAAQDVAAKCKEVGITAVHVKIRATGGTRTKTPGPGGQAALRALARSGLRIGRIEDVTPVPSDSTRKKGGRRGRRL"));
            record = iter.next();
            assertTrue(iter.hasNext());
            assertThat(record.getSampleId(), equalTo("SRR10210920"));
            record = iter.next();
            assertFalse(iter.hasNext());
            assertThat(record.getSampleId(), equalTo("SRR10210920"));
            assertThat(record.getProts(), contains("MDSAKIINIILSLFLPPVAVFLARGWGTDCIVDIILTILAWFPGMLYALYIVLQDIAHSINESNYQMK",
                    "MDSAKIINIILSLFLPPVAVFLARGWGTDCIVDIILTILAWFPGMLYALYIVLQD"));
        }
    }

}
