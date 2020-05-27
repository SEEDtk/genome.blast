/**
 *
 */
package org.theseed.reports;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.theseed.io.TabbedLineReader;

/**
 * This is a utility class used to create analysis records.  To create the class, the client specifies
 * the names of the numeric fields, the name of the type field, and the good/bad values.
 *
 * An analysis record is either "good" or "bad".  Our goal is to find filtering values that maximize the
 * number of good records and minimize the number of bad records.
 *
 * Analysis records are designed to be sorted on any one of the numeric fields.  Within a single value, they
 * are sorted in input order.  Methods are provided to do the sorts and to compute metrics for each of the
 * number sets.
 *
 * @author Bruce Parrello
 *
 */
public class AnalysisList {

    // FIELDS
    /** array of names of numeric fields */
    private String[] fieldName;
    /** array of indices of numeric fields */
    private int[] fieldIdx;
    /** name of type field */
    private String typeName;
    /** sortNum of type field */
    private int typeIdx;
    /** list of records */
    private List<Record> records;
    /** value for a good record */
    private String goodValue;
    /** index of the current sort for the list */
    private int sortIdx;
    /** total good values in list */
    private int totalGood;

    /**
     * This object represents a single input record for the analysis process.  It contains a fixed number of numeric fields
     * and a good/bad flag.  It can be sorted on any of the fields.
     */
    public class Record {

        // FIELDS
        /** numeric fields for sorting */
        private double[] numbers;
        /** TRUE if good */
        private boolean good;
        /** sortNum number to resolve ties */
        private int sortNum;

        /**
         * Create an analysis record.
         *
         * @param goodFlag	TRUE if this is a good record, else FALSE
         * @param sortNum	sort number to assign to the record
         * @param numbers	array of numbers to use (will not be cloned)
         */
        protected Record(boolean goodFlag, int sortNum, double[] numbers) {
            this.numbers = numbers;
            this.good = goodFlag;
            this.sortNum = sortNum;
        }

        /**
         * @return TRUE if this is a good record, else FALSE
         */
        public boolean isGood() {
            return good;
        }

        /**
         * @return the indicated number value for this record
         *
         * @param idx	sortNum of the number value desired
         */
        public double getNumber(int idx) {
            return this.numbers[idx];
        }

        /**
         * @return the input number of this record
         */
        public int getSortNumber() {
            return this.sortNum;
        }

    }

    /**
     * This class generates a comparator based on a given number field.  The comparator sorts high
     * values to the beginning.
     */
    public class Compare implements Comparator<Record> {

        // FIELDS
        /** sortNum of comparison field */
        private int index;

        /**
         * Construct a comparator for sorting on a specific number field.
         *
         * @param sortNum		sortNum of field to sort
         */
        public Compare(int index) {
            this.index = index;

        }

        @Override
        public int compare(Record o1, Record o2) {
            int retVal = Double.compare(o2.getNumber(this.index), o1.getNumber(this.index));
            if (retVal == 0)
                retVal = o1.getSortNumber() - o2.getSortNumber();
            return retVal;
        }

    }


    /**
     * Construct an analysis list from an input file.
     *
     * @param inStream		input stream containing the records to analyze
     * @param typeCol		label of the column containing the good/bad type
     * @param numLabels		labels of the columns containing the number values
     * @param goodValue		type value for a good column
     * @throws IOException
     */
    public AnalysisList(TabbedLineReader inStream, String typeCol, List<String> numLabels, String goodValue) throws IOException {
        // Determine the number of numeric fields.
        int n = numLabels.size();
        this.fieldIdx = new int[n];
        this.fieldName = new String[n];
        this.fieldName = numLabels.toArray(this.fieldName);
        // Fill in the column indices for the number fields.
        for (int i = 0; i < n; i++)
            this.fieldIdx[i] = inStream.findField(this.fieldName[i]);
        // Get the column index for the type column.
        this.typeName = typeCol;
        this.typeIdx = inStream.findField(typeCol);
        // Save the "good" label.
        this.goodValue = goodValue;
        // Now we need to read in the records.
        this.totalGood = 0;
        this.records = new ArrayList<Record>(5000);
        for (TabbedLineReader.Line line : inStream) {
            boolean goodFlag = (line.get(this.typeIdx).contentEquals(goodValue));
            if (goodFlag) this.totalGood++;
            double[] numbers = new double[n];
            for (int i = 0; i < n; i++)
                numbers[i] = line.getDouble(this.fieldIdx[i]);
            Record newRecord = new Record(goodFlag, this.records.size(), numbers);
            this.records.add(newRecord);
        }
        // Denote the list is unsorted.
        this.sortIdx = -1;
        // Abort if the input is empty.
        if (this.records.size() == 0)
            throw new IOException("Input file cannot be empty.");
    }

    /**
     * @return the number of records in the list
     */
    public int size() {
        return this.records.size();
    }

    /**
     * Sort the list by the specified numeric field.
     *
     * @param idx		index of the number field by which to sort
     */
    public void sort(int idx) {
        Comparator<Record> sorter = this.new Compare(idx);
        this.records.sort(sorter);
        this.sortIdx = idx;
    }

    /**
     * @return the number of good records
     */
    public int getGoodCount() {
        return this.totalGood;
    }

    /**
     * @return the number of bad records
     */
    public int getBadCount() {
        return (this.records.size() - this.totalGood);
    }

    /**
     * @return the label for the specified number column
     *
     * @param idx	index of desired number
     */
    public String getLabel(int idx) {
        return this.fieldName[idx];
    }

    /**
     * Compute the cutoff for the specified sensitivity using the field currently sorted.
     *
     * @param sensitivePct	the sensitivity percentage
     *
     * @return the value of the specified field above which the specified percentage of good values can be found
     */
    public double getSensitive(double sensitivePct) {
        // Compute the minimum number of good values needed.
        int minGoods = (int) Math.ceil(sensitivePct * this.getGoodCount() / 100);
        // Loop through the input, counting good records until we reach the limit.  This is tricky, as we need
        // a copy of the last record found available at the end.  The loop below works only because we have
        // guaranteed there is at least one record in the file.
        int goodFound = 0;
        Iterator<Record> iter = this.records.iterator();
        Record record = null;
        while (goodFound < minGoods && iter.hasNext()) {
            record = iter.next();
            if (record.isGood()) goodFound++;
        }
        // If the user specified a cutoff of 0, we return the value in the first record.
        if (record == null) record = iter.next();
        return record.getNumber(this.sortIdx);
    }

    /**
     * Compute the cutoff for the specified specificity using the field currently sorted.
     *
     * @param specificPct	the specificity percentage
     *
     * @return the value of the specified field above which the specified percentage of records are good,
     * 			or NaN if we can't find such a value
     *
     */
    public double getSpecific(double specificPct) {
        // Loop through the input, counting good and bad records until we reach the specificity.  We remember
        // the last value for which the specificity is larger than our minimum.  The value only counts if
        // the next record has a different sort value than the current one, so we only save when we are
        // one record further than the target.
        double retVal = Double.NaN;
        int goodCount = 0;
        Iterator<Record> iter = this.records.iterator();
        Record record = iter.next();
        if (record.isGood()) goodCount++;
        int totCount = 1;
        double oldNumber = record.getNumber(this.sortIdx);
        while (iter.hasNext()) {
            record = iter.next();
            double newNumber = record.getNumber(this.sortIdx);
            // If the number of interest has changed, check to see if we need to save the cutoff.
            if (newNumber < oldNumber) {
                double newSpecificity = specificity(goodCount, totCount);
                if (newSpecificity >= specificPct)
                    retVal = oldNumber;
            }
            // Process the new record.
            oldNumber = newNumber;
            if (record.isGood()) goodCount++;
            totCount++;
        }
        // Make a final check here for the last record.
        if (specificity(goodCount, totCount) >= specificPct)
            retVal = oldNumber;
        return retVal;

    }

    /**
     * @return the specificity percent for a given good count and total count
     *
     * @param goodCount		number of good records
     * @param totCount		total number of records
     */
    protected double specificity(int goodCount, int totCount) {
        return goodCount * 100.0 / totCount;
    }

    /**
     * Compute the sensitivity percentage of the specified cutoff value using the specified field.
     *
     * @param cutoff	the cutoff value
     * @param idx		index of the relevant field
     *
     * @return the percentage of good values above the cutoff
     */
    public double getSensitivity(double cutoff, int idx) {
        int goodCount = 0;
        for (Record record : this.records) {
            if (record.getNumber(idx) >= cutoff && record.isGood())
                goodCount++;
        }
        return (goodCount * 100.0) / this.totalGood;
    }

    /**
     * @return the accuracy associated with a given cutoff for the specified index
     *
     * @param cutoff	the cutoff value
     * @param idx		index of the relevant field
     *
     * @return the percentage of correctly-called values for the cutoff
     */
    public double getAccuracy(double cutoff, int idx) {
        int falsePositive = 0;
        int falseNegative = 0;
        for (Record record : this.records) {
            if (record.getNumber(idx) >= cutoff) {
                // Here we think the record is good.
                if (! record.isGood()) falsePositive++;
            } else {
                // Here we think the record is bad.
                if (record.isGood()) falseNegative++;
            }
        }
        return this.accuracy(falseNegative, falsePositive);
    }


    /**
     * Compute the cutoff the produces the highest accuracy rate using the field currently sorted.
     */
    public double getAccurate() {
        // Get the first record.  We only record an accuracy when the value of the sorted number changes.
        Iterator<Record> iter = this.records.iterator();
        Record record = iter.next();
        double oldNumber = record.getNumber(this.sortIdx);
        // Start above the first record.  Remember the false positive and false negative counts.
        // Note the initial value is slightly above the highest in the list.
        double retVal = oldNumber;
        if (retVal == 0) retVal = 0.1;
        else if (retVal > 0) retVal = 1.01 * retVal;
        else retVal = 0.99 * retVal;
        int falseNegative = this.totalGood;
        int falsePositive = 0;
        double bestAccuracy = this.accuracy(falseNegative, falsePositive);
        // Now loop through the records, accumulating the best accuracy.
        while (iter.hasNext()) {
            double newNumber = record.getNumber(this.sortIdx);
            if (newNumber < oldNumber) {
                double newAccuracy = this.accuracy(falseNegative, falsePositive);
                if (newAccuracy > bestAccuracy) {
                    retVal = oldNumber;
                    bestAccuracy = newAccuracy;
                }
            }
            // Process this record.
            if (record.isGood())
                falseNegative--;
            else
                falsePositive++;
            oldNumber = newNumber;
            record = iter.next();
        }
        // Check the last record's value.
        if (this.accuracy(falseNegative, falsePositive) > bestAccuracy)
            retVal = oldNumber;
        return retVal;
    }

    /**
     * @return the accuracy fraction given the false negative and positive counts
     *
     * @param falseNegative		number of false negatives
     * @param falsePositive		number of false positives
     */
    private double accuracy(int falseNegative, int falsePositive) {
        return (this.records.size() - falseNegative - falsePositive) * 100.0 / this.records.size();
    }

    /**
     * Compute the specificity percentage of the specified cutoff value using the specified field.
     *
     * @param cutoff	the cutoff value
     * @param idx		index of the relevant field
     *
     * @return the percentage of values above the cutoff that are good
     */
    public double getSpecificity(double cutoff, int idx) {
        int goodCount = 0;
        int badCount = 0;
        for (Record record : this.records) {
            if (record.getNumber(idx) >= cutoff) {
                if (record.isGood())
                    goodCount++;
                else
                    badCount++;
            }
        }
        return (goodCount * 100.0) / (goodCount + badCount);
    }

    /**
     * @return the name of the type column
     */
    public String getTypeName() {
        return typeName;
    }

    /**
     * @return the type value for a good record
     */
    public String getGoodValue() {
        return goodValue;
    }



}
