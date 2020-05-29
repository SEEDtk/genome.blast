/**
 *
 */
package org.theseed.sequence.blast;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;

/**
 * @author Bruce Parrello
 *
 */
public class GtiFile implements Iterable<GtiFile.Record>, AutoCloseable, Closeable, Iterator<GtiFile.Record> {

    // FIELDS
    /** input stream */
    private TabbedLineReader inStream;

    /**
     * This class describes a GTI record.
     */
    public class Record {

        // FIELDS
        /** ID of the source sample */
        private String sampleId;
        /** ID of the RNA sequence fragment */
        private String rnaId;
        /** location of the DNA in the genome */
        private Location dnaLoc;
        /** DNA string */
        private String dna;
        /** list of proteins */
        private List<String> prots;

        /**
         * Construct a GTI record from an input line.
         *
         * @param line		input line containing the GTI data
         */
        protected Record(TabbedLineReader.Line line) {
            this.sampleId = line.get(0);
            this.rnaId = line.get(1);
            this.dnaLoc = Location.fromString(line.get(2));
            this.dna = line.get(3);
            this.prots = Arrays.asList(StringUtils.split(line.get(4), ','));
        }

        /**
         * @return the sample ID
         */
        public String getSampleId() {
            return sampleId;
        }

        /**
         * @return the rna fragment ID
         */
        public String getRnaId() {
            return rnaId;
        }

        /**
         * @return the dna location
         */
        public Location getDnaLoc() {
            return dnaLoc;
        }

        /**
         * @return the dna sequence
         */
        public String getDna() {
            return dna;
        }

        /**
         * @return the protein sequences
         */
        public List<String> getProts() {
            return prots;
        }
    }


    /**
     * Construct a GTI file object for reading from a specified file.
     *
     * @param inFile	GTI file to iterate through
     *
     * @throws IOException
     */
    public GtiFile(File inFile) throws IOException {
        // Note we use a headerless tab-delimited file.
        this.inStream = new TabbedLineReader(inFile, 5);
    }

    @Override
    public void close() {
        // Close the input stream.
        this.inStream.close();
    }

    @Override
    public Iterator<Record> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return this.inStream.hasNext();
    }

    @Override
    public Record next() {
        Record retVal = null;
        TabbedLineReader.Line line = this.inStream.next();
        if (line != null)
            retVal = new Record(line);
        return retVal;
    }

}
