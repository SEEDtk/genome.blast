package org.theseed.genome.blast;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.theseed.io.Shuffler;
import org.theseed.io.TabbedLineReader;
import org.theseed.reports.AnalysisList;

/**
 * Unit test for simple App.
 */
public class AnalyzeTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AnalyzeTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AnalyzeTest.class );
    }

    /**
     * test analysis list
     *
     * @throws IOException
     */
    public void testAnalyzer() throws IOException {
        File inFile = new File("src/test", "analysis.txt");
        TabbedLineReader inStream = new TabbedLineReader(inFile);
        List<String> cols = new Shuffler<String>(2).add1("eval").add1("percent");
        AnalysisList aList = new AnalysisList(inStream, "type", cols, "good");
        assertThat(aList.getTypeName(), equalTo("type"));
        assertThat(aList.getGoodCount(), equalTo(10));
        assertThat(aList.getBadCount(), equalTo(10));
        assertThat(aList.size(), equalTo(20));
        assertThat(aList.getSpecificity(65.0, 1), closeTo(72.727, 0.001));
        assertThat(aList.getSpecificity(68.0, 1), closeTo(66.666, 0.001));
        assertThat(aList.getSensitivity(65.0, 1), closeTo(80.000, 0.001));
        assertThat(aList.getSensitivity(68.0, 1), closeTo(60.000, 0.001));
        assertThat(aList.getGoodValue(), equalTo("good"));
        assertThat(aList.getLabel(0), equalTo("eval"));
        assertThat(aList.getLabel(1), equalTo("percent"));
        for (int idx = 0; idx < 2; idx++) {
            aList.sort(idx);
            for (double pct = 5; pct < 100; pct += 5) {
                double cutoff = aList.getSensitive(pct);
                assertThat(String.format("%s at %4.2f", aList.getLabel(idx), pct),
                        aList.getSensitivity(cutoff, idx), greaterThanOrEqualTo(pct));
                cutoff = aList.getSpecific(pct);
                if (! Double.isNaN(cutoff))
                    assertThat(String.format("%s at %4.2f", aList.getLabel(idx), pct),
                            aList.getSpecificity(cutoff, idx), greaterThanOrEqualTo(pct));
            }

        }
        double cutoff = aList.getAccurate();
        assertThat(cutoff, equalTo(55.0));
        assertThat(aList.getAccuracy(55.0, 1), closeTo(80.000, 0.001));
        assertThat(aList.getAccuracy(80.0, 1), closeTo(55.000, 0.001));
        assertThat(aList.getAccuracy(68.0, 1), closeTo(65.000, 0.001));
    }

}
