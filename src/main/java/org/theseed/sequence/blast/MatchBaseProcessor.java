package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.locations.BLocation;
import org.theseed.locations.FLocation;
import org.theseed.locations.Location;
import org.theseed.proteins.LocationFixer;
import org.theseed.reports.MatchReporter;
import org.theseed.sequence.DnaDataStream;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseProcessor;

/**
 * This is a base class that runs a single genome against a single RNA sequence file.  It is used by
 * both MatchProcessor and MatchRunProcessor.
 *
 * @author Bruce Parrello
 *
 */
public abstract class MatchBaseProcessor extends BaseProcessor {

    // FIELDS

    /** logging facility */
    protected static final Logger log = LoggerFactory.getLogger(MatchProcessor.class);
    /** RNA sequence BLAST database */
    private DnaBlastDB rnaDB;
    /** genome contig BLAST database */
    private DnaBlastDB genomeDB;
    /** input genome */
    private Genome genome;
    /** filtering protein profiles */
    private ProteinProfiles profiler;
    /** BLAST parameters for the main query */
    private BlastParms mainParms;
    /** output report */
    private MatchReporter reporter;
    /** RNA sequences in the current batch */
    private DnaDataStream rnaStream;
    /** map of sequence IDs to proteins for the current batch */
    private Map<String, List<Sequence>> protMap;
    /** DNA translator for the genome's genetic code */
    private LocationFixer xlator;
    /** number of batches processed */
    private int batchNum;
    /** current protein number */
    private int protNum;
    /** current sample ID */
    private String sample;

    // COMMAND-LINE OPTIONS
    /** number of DNA sequences to submit in each BLAST call */
    @Option(name = "-b", aliases = { "--batchSize",
            "--batch" }, metaVar = "1", usage = "number of input sequences to submit to each BLAST call")
    private int batchSize;
    /** distance to extend the genome hit on either side */
    @Option(name = "-x", aliases = {
            "--extend" }, metaVar = "20", usage = "distance to extend the genome hit on either side")
    private int extend;
    /** minimum percent of query coverage for an RNA-to-genome match */
    @Option(name = "--minPct", metaVar = "95", usage = "minimum percent of RNA query that must match genome DNA")
    private double minPct;
    /** maximum permissible e-value */
    @Option(name = "--maxE", aliases = {
            "--evalue" }, metaVar = "1e-20", usage = "maximum permissible e-value for a match")
    private double eValue;
    /** temporary directory for BLAST database */
    @Option(name = "--tempDir", metaVar = "Tmp", usage = "temporary directory for BLAST databases")
    private File tempDir;
    /** maximum gap between adjacent sequences */
    @Option(name = "--maxGap", metaVar = "100", usage = "maximum gap between proteins to join into an operon")
    private int maxGap;
    /** output format */
    @Option(name = "--format", aliases = { "--outFormat", "--outFmt" }, usage = "output format")
    private MatchReporter.Type outFormat;
    /** minimum percent identity */
    @Option(name = "--minIdent", aliases = { "--percIdentity",
            "--minI" }, metaVar = "75", usage = "minimum percent identity for a genome hit")
    private double minPctIdentity;
    /** minimum query-scaled bit score */
    @Option(name = "--minQbsc", metaVar = "1.1", usage = "minimum acceptable query-scaled bit score for profile hits")
    private double minQbsc;
    /** minimum query identity fraction */
    @Option(name = "--minQIdent", metaVar = "0.5", usage = "minimum acceptable query identity fraction for profile hits")
    private double minQIdent;
    /** minimum percent query coverage for a legitimate hit */
    @Option(name = "--minQuery", aliases = {
            "--minQ" }, metaVar = "75", usage = "minimum percent of query sequence that must match in a profile hit")
    private double minPctQuery;
    /** algorithm for finding starts */
    @Option(name = "--starts", usage = "algorithm for finding start codons from profile hits")
    private LocationFixer.Type algorithm;

    /** profile direcory */
    @Argument(index = 0, metaVar = "profileDir", usage = "protein profile directory", required = true)
    private File profileDir;

    public MatchBaseProcessor() {
        super();
    }

    /**
     * Set the default values for common parameters.
     */
    protected void setupDefaults() {
        this.batchSize = 10;
        this.eValue = 1e-10;
        this.minPct = 95.0;
        this.maxGap = 500;
        this.extend = 50;
        this.minPctIdentity = 90.0;
        this.minPctQuery = 65.0;
        this.minQbsc = 1.1;
        this.minQIdent = 0.0;
        this.tempDir = new File(System.getProperty("user.dir"), "Temp");
        this.outFormat = MatchReporter.Type.GTI;
        this.reporter = null;
        this.protNum = 1;
        this.algorithm = LocationFixer.Type.NEAREST;
    }

    /**
     * Validate the common parameters for a match run.
     *
     * @throws FileNotFoundException
     * @throws IOException
     */
    protected void validateCommonParms() throws FileNotFoundException, IOException {
        // Validate the profile directory.
        if (! this.profileDir.isDirectory())
            throw new FileNotFoundException("Profile directory " + this.profileDir +
                    " not found or invalid.");
        else {
            log.info("Input sequences will be profiled from {}.", this.profileDir);
            this.profiler = new ProteinProfiles(this.profileDir);
        }
        // Set up the temp directory.
        if (! this.tempDir.isDirectory()) {
            log.info("Creating temporary file directory {}.", this.tempDir);
            FileUtils.forceMkdir(this.tempDir);
        }
        // Validate the parameters.
        if (this.eValue >= 1.0)
            throw new IllegalArgumentException("Invalid eValue specified.  Must be less than 1.");
        if (this.batchSize <= 0)
            throw new IllegalArgumentException("Batch size must be 1 or more.");
        if (this.maxGap < 0)
            throw new IllegalArgumentException("Maximum gap must be 0 or more.");
        if (this.extend < 0)
            throw new IllegalArgumentException("Extension length must be 0 or more.");
        if (this.minPct < 0.0 || this.minPct > 100.0)
            throw new IllegalArgumentException("Minimum RNA match percent must be between 0 and 100.");
        if (this.minPctIdentity < 0.0 || this.minPctIdentity > 100.0)
            throw new IllegalArgumentException("Minimum percent identity must be between 0 and 100.");
        if (this.minQbsc < 0.0 || this.minQbsc > 10.0)
            throw new IllegalArgumentException("Minimum query-scaled bit score must be between 0 and 10.");
        if (this.minQIdent < 0.0 || this.minQIdent > 1.0)
            throw new IllegalArgumentException("Minimum query identity fraction must be between 0 and 1");
        if (this.minPctQuery < 0.0 || this.minPctQuery > 100.0)
            throw new IllegalArgumentException("Minimum query percentation must be between 0 and 100.");
    }

    /**
     * Set up the match processor for a single sample's run.
     *
     * @param sampleId		ID of the current sample
     * @param gFile			file containing the genome
     * @param rnaFile		file containing the RNA input
     * @param outStream 	output stream for report
     *
     * @throws IOException
     * @throws FileNotFoundException
     * @throws InterruptedException
     */
    protected void setup(File gFile, File rnaFile, OutputStream outStream) throws IOException, FileNotFoundException {
        // Read the genome.
        log.info("Loading genome from {}.", gFile);
        this.genome = new Genome(gFile);
        // Get the appropriate translator.
        this.xlator = this.algorithm.create(this.getGenome().getGeneticCode());
        // Connect to the RNA input stream.
        if (! rnaFile.canRead())
            throw new FileNotFoundException("RNA input file " + rnaFile + " not foudn or unreadable.");
        else
            try {
                log.info("RNA sequences will be read from {}.", rnaFile);
                this.rnaDB = DnaBlastDB.create(rnaFile, this.getGenome().getGeneticCode());
                this.rnaDB.cleanOnExit();
            } catch (InterruptedException e) {
                throw new IOException("Interruption during BLAST DB creation: " + e.getMessage());
            }
        // Create the output report.
        this.reporter = this.outFormat.create(outStream, this.getGenome());
        reporter.initialize(this);
    }

    /**
     * Process the RNA for a single genome.
     *
     * @throws IOException
     * @throws InterruptedException
     * @throws FileNotFoundException
     */
    protected void runGenome(String sample, File rnaFile)
            throws IOException, InterruptedException, FileNotFoundException {
                try {
                    // Save the sample ID.
                    this.sample = sample;
                    // Create the BLAST database from the genome.
                    File tempFile = File.createTempFile("blast", ".fasta", this.tempDir);
                    this.genomeDB = DnaBlastDB.create(tempFile, getGenome());
                    this.genomeDB.deleteOnExit();
                    // Create the BLAST parameters.
                    this.mainParms = new BlastParms().maxE(this.eValue).pctLenOfQuery(this.minPct).minPercent(this.minPctIdentity);
                    // First, we must find the protein regions in the DNA.  We use the profiler to get a list
                    // of blast hits for each RNA sequence.
                    log.info("Profiling the RNA sequences for protein regions.");
                    BlastParms parms = new BlastParms().maxE(this.eValue).pctLenOfQuery(this.minPctQuery)
                            .minQueryBitScore(this.minQbsc).minQueryIdentity(this.minQIdent);
                    Map<String, List<BlastHit>> hitMap = this.profiler.profile(this.rnaDB, parms);
                    log.info("{} contigs contained proteins.", hitMap.size());
                    // This DNA stream will hold the current batch of RNA sequence fragments found.
                    this.rnaStream = new DnaDataStream(this.batchSize, this.getGenome().getGeneticCode());
                    // For each sequence, this will contain the proteins found.
                    this.protMap = new HashMap<String, List<Sequence>>();
                    // Denote no batches have been processed yet.
                    this.batchNum = 0;
                    // Now we loop through all the RNA sequences, processing the ones that have hits.
                    try (FastaInputStream rnaSequences = new FastaInputStream(rnaFile)) {
                        for (Sequence rnaSeq : rnaSequences) {
                            // Get the hits for this sequence.
                            List<BlastHit> hits = hitMap.get(rnaSeq.getLabel());
                            if (hits != null) {
                                // Insure there is room in this batch for a new sequence.
                                if (rnaStream.size() >= this.batchSize) {
                                    this.processBatch();
                                    rnaStream.clear();
                                    protMap.clear();
                                }
                                // Find all the proteins in the contig and compute their DNA.
                                this.processSequence(rnaSeq, hits);
                            }
                        }
                    }
                    this.processBatch();
                    // Finish the output report.
                    reporter.finish();
                } finally {
                    if (this.reporter != null)
                        this.reporter.close();
                }
            }

    /**
     * Extract protein regions from the RNA.  Each hit is extended to a full protein.  RNA regions that are close
     * together are combined and the combined sequence is saved as an RNA fragment.  The fragment is then associated
     * with the protein sequences extracted.
     *
     * It is worth noting that the number of hits is not going to be large.  Most of the time we get one hit, and
     * a large RNA sequence will have only ten.  This means sequential searches are ok.
     *
     * Finally, the RNA sequence is the subject of the hit, not the query.
     *
     * @param rnaSeq	incoming full RNA sequence
     * @param hits		list of protein hits against the sequence
     */
    private void processSequence(Sequence rnaSeq, List<BlastHit> hits) {
        String rnaSequence;
        // These variables are used to build protein IDs.  The first part is an RNA seq ID and the second part is
        // a counter value.
        String rnaLabel = rnaSeq.getLabel();
        log.info("Processing RNA sequence {}.", rnaLabel);
        // First we need to pick a direction.  We will keep only the locations in that direction.  Most
        // of the time, all of the hits will be in the same direction.
        int dirCount = hits.stream().mapToInt(x -> (x.getSubjectLoc().getDir() == '+' ? 1 : -1)).sum();
        int seqLen = rnaSeq.length();
        List<ProteinHit> kept = null;
        if (dirCount >= 0) {
            // Here we are working in the plus direction.
            kept = hits.stream().filter(x -> x.getSubjectLoc() instanceof FLocation)
                    .map(x -> new ProteinHit(x.getQueryId(), x.getSubjectLoc(), seqLen)).collect(Collectors.toList());
            rnaSequence = rnaSeq.getSequence().toLowerCase();
            log.debug("{} forward hits found in {}.", kept.size(), rnaLabel);
        } else {
            // Here we are working in the minus direction.  We have to flip the sequence and convert the locations.
            rnaSequence = Contig.reverse(rnaSeq.getSequence().toLowerCase());
            kept = hits.stream().filter(x -> x.getSubjectLoc() instanceof BLocation)
                    .map(x -> new ProteinHit(x.getQueryId(), x.getSubjectLoc(), seqLen)).collect(Collectors.toList());
            log.debug("{} backward hits found in {}.", kept.size(), rnaLabel);
        }
        // At this point we have a bunch of forward locations known to correspond to proteins. We need to expand each
        // one to a start and a stop.  If we can't find a start and a stop, we toss the hit.  If the resulting
        // location matches another hit, we throw away the shorter one.
        List<ProteinHit> processed = new ArrayList<ProteinHit>(kept.size());
        for (ProteinHit protHit : kept) {
            boolean ok = this.xlator.fix(protHit.getLoc(), rnaSequence);
            // Only proceed if we successfully extended to a start and stop.
            if (ok) {
                // Find out if there is already a location with the same stop.
                final int stopLoc = protHit.getLoc().getRight();
                Optional<ProteinHit> rival = processed.stream().filter(x -> x.getLoc().getRight() == stopLoc).findFirst();
                if (rival.isPresent()) {
                    ProteinHit otherHit = rival.get();
                    // Here we have another location with the same stop.  There can be only one.  If it is longer than this
                    // one, throw this one away.  Otherwise, delete it and add this one.
                    if (otherHit.getLoc().getLength() < protHit.getLoc().getLength()) {
                        processed.remove(otherHit);
                        processed.add(protHit);
                    }
                } else {
                    // This is the only location with this stop.
                    processed.add(protHit);
                }
            }
        }
        // Now "processed" contains each location with a protein.  The next step is to create the proteins.  Note we strip
        // off the stop codon in the translation.  The protein is encapsulated by a Sequence object (in case we want to
        // put it in a FASTA file), with the location as the comment and a manufactured label.  Each protein is associated
        // with a containing location in the RNA sequence.  Close locations are combined and their sequence lists merged.
        log.debug("{} proteins found in {}.", processed.size(), rnaLabel);
        Map<Location, List<Sequence>> proteinMap = new HashMap<Location, List<Sequence>>(processed.size());
        for (ProteinHit protHit : processed) {
            Sequence protein = new Sequence(String.format("%s.%04d", protHit.getProfile(), this.protNum++),
                    protHit.getLoc().toString(), this.xlator.pegTranslate(rnaSequence, protHit.getLoc().getLeft(),
                    protHit.getLoc().getLength() - 3));
            // Now we do the location-combining.  If we find a close location in the location list, we
            // merge it with this one.  We need a list to hold all the proteins for the merged locations.
            List<Sequence> proteins = new ArrayList<Sequence>(10);
            proteins.add(protein);
            // Since we are deleting while looping, we have to use an iterator.
            Iterator<Map.Entry<Location, List<Sequence>>> mapIter = proteinMap.entrySet().iterator();
            while (mapIter.hasNext()) {
                Map.Entry<Location, List<Sequence>> protEntry = mapIter.next();
                Location other = protEntry.getKey();
                List<Sequence> otherSeqs = protEntry.getValue();
                if (other.distance(protHit.getLoc()) < this.maxGap) {
                    // Combine the two locations.
                    protHit.getLoc().merge(other);
                    // Remove the other location from the map and add its proteins to ours.
                    mapIter.remove();
                    proteins.addAll(otherSeqs);
                }
            }
            // Now we add the combined location and its proteins to the map.
            proteinMap.put(protHit.getLoc(), proteins);
        }
        log.debug("{} operons found in {}.", proteinMap.size(), rnaLabel);
        // We are almost done.  Now each location gets turned into an RNA fragment and gets added to the main protein map
        // along with its associated proteins.
       int rnaNum = 1;
        for (Map.Entry<Location, List<Sequence>> protEntry : proteinMap.entrySet()) {
            Location loc = protEntry.getKey();
            String fragment = loc.getDna(rnaSequence);
            String fragmentId = String.format("r.%s.%04d", rnaLabel, rnaNum++);
            Sequence rnaFragment = new Sequence(fragmentId, loc.toString(), fragment);
            this.rnaStream.add(rnaFragment);
            this.protMap.put(fragmentId, protEntry.getValue());
        }
    }

    /**
     * BLAST the RNA sequence fragments against the genome to find the corresponding DNA, then output the
     * DNA with the associated proteins.  The RNA sequence fragments are in this.rnaStream.  The proteins
     * are in this.protMap.  The comment for each protein is its location in the original RNA sequence.
     *
     * @throws InterruptedException
     * @throws IOException
     */
    private void processBatch() throws IOException, InterruptedException {
        // First we need the DNA from the genome.  We use the BLAST to do this.
        this.batchNum++;
        log.debug("Processing RNA batch {} with {} fragments.", this.batchNum, this.rnaStream.size());
        Map<String, List<BlastHit>> hitMap = BlastHit.sort(this.genomeDB.blast(this.rnaStream, this.mainParms));
        // Now for each RNA sequence fragment we have a list of the blast hits associated with it.  We take the
        // best hit and use its DNA.
        int outputCount = 0;
        BlastHit.Longest getBest = new BlastHit.Longest();
        for (Map.Entry<String, List<BlastHit>> hitEntry : hitMap.entrySet()) {
            // Get the longest hit.
            Optional<BlastHit> bestHit = hitEntry.getValue().stream().max(getBest);
            if (bestHit.isPresent()) {
                BlastHit hit = bestHit.get();
                // We have our longest hit.  Get its DNA location in the genome.
                Location loc = hit.getSubjectLoc();
                int contigLen = this.getGenome().getContig(loc.getContigId()).length();
                // Add the requested padding on each side.
                loc.expand(this.extend, this.extend, contigLen);
                // Compute the DNA.
                String dna = this.getGenome().getDna(loc);
                // Output this result.
                String rnaId = hit.getQueryId();
                List<Sequence> prots = this.protMap.get(rnaId);
                this.reporter.processSequence(rnaId, loc, dna, prots);
                outputCount++;
            }
        }
        log.info("Batch {} contained {} fragments and {} were output.", this.batchNum, this.rnaStream.size(), outputCount);
    }

    /**
     * @return the ID of the source RNA sample (if any)
     */
    public String getSampleID() {
        return this.sample;
    }

	/**
	 * @return the genome
	 */
	public Genome getGenome() {
		return genome;
	}

}
