/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.reports.UniProfileReporter;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.ProteinBlastDB;
import org.theseed.sequence.blast.ProteinProfiles;

/**
 * This class checks for singly-occurring universal roles in all genomes from a genome source.  If a role is missing, a profile blast will be
 * used to try to find it in existing proteins.  If the profile hits, we will display the annotation of the feature hit to determine if the
 * annotation is faulty.
 *
 * The positional parameters are the name of the input genome source, the name of the profile directory, the file containing the list
 * of universal roles, and the name of the output directory.  The universal role file should be tab-delimited with the universal roles
 * in the first column.  The output directory will contain various reports about the findings.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -d		working directory for storing created files; the default is the current directory
 *
 * --clear			erase the output directory before starting
 * --source			type of genome source (default DIR)
 * --maxE			maximum permissible e-value; the default is 1e-10
 * --minIdent		minimum percent identity for hits; the default is 0
 * --minQIdent		minimum query identity fraction for hits; the default is 0.50
 * --minQbsc		minimum query-scaled bit score for hits; the default is 0.0
 * --minQuery		minimum percent query match; the default is 50.0
 *
 * @author Bruce Parrello
 *
 */
public class UniProfileProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UniProfileProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** protein profile manager */
    private ProteinProfiles profiles;
    /** role definition map */
    private RoleMap roleMap;
    /** list of reporting objects */
    private List<UniProfileReporter> reporters;
    /** empty set of blast hits */
    private static final Set<BlastHit> NO_HITS = Collections.emptySet();

    // COMMAND-LINE OPTIONS

    /** if TRUE, the output directory will be erased before processing */
    @Option(name = "--clear", usage = "if specified, the output directory will be erased before processing")
    private boolean clearFlag;

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** maximum E-value */
    @Option(name = "--maxE", aliases = { "--evalue" }, metaVar = "1e-5",
            usage = "maximum permissible e-value for a query")
    private double eValue;

    /** minimum percent identity */
    @Option(name = "--minIdent", aliases = { "--percIdentity", "--minI" }, metaVar = "75",
            usage = "minimum percent identity for a hit")
    private double minPctIdentity;

    /** minimum query-scaled bit score */
    @Option(name = "--minQbsc", metaVar = "1.1", usage = "minimum acceptable query-scaled bit score")
    private double minQbsc;

    /** minimum query identity fraction */
    @Option(name = "--minQIdent", metaVar = "0.5", usage = "minimum acceptable query identity fraction")
    private double minQIdent;

    /** minimum percent query coverage for a legitimate hit */
    @Option(name = "--minQuery", aliases = { "--minQ" }, metaVar = "75",
            usage  = "minimum percent of query sequence that must be hit")
    private double minPctQuery;

    /** input genome directory */
    @Argument(index = 0, metaVar = "genomeDir", usage = "file or directory containing input genomes")
    private File inDir;

    /** profile directory */
    @Argument(index = 1, metaVar = "profileDir", usage = "directory containing role profiles")
    private File profileDir;

    /** SOUR role file */
    @Argument(index = 2, metaVar = "sours.tbl", usage = "file of universal roles to check")
    private File sourFile;

    /** output directory */
    @Argument(index = 3, metaVar = "outDir", usage = "output directory")
    private File outDir;


    @Override
    protected void setDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
        this.eValue = 1e-10;
        this.minPctIdentity = 0.50;
        this.minQbsc = 0.0;
        this.minQIdent = 0.0;
        this.minPctQuery = 0.0;
        this.clearFlag = false;
    }

    @Override
    protected void validateParms() throws IOException, ParseFailureException {
        // Verify the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else {
            if (this.clearFlag) {
                log.info("Erasing output directory {}.", this.outDir);
                FileUtils.cleanDirectory(this.outDir);
            } else
                log.info("Output will be placed in directory {}.", this.outDir);
        }
        // Validate the BLAST parameters.
        if (this.eValue < 0.0)
            throw new ParseFailureException("Maximum e-value cannot be negative.");
        // Validate the filters.
        if (this.minPctIdentity < 0.0 || this.minPctIdentity > 100.0)
            throw new ParseFailureException("Minimum percent identity must be between 0 and 100.");
        if (this.minQbsc < 0.0 || this.minQbsc > 10.0)
            throw new ParseFailureException("Minimum query-scaled bit score must be between 0 and 10.");
        if (this.minQIdent < 0.0 || this.minQIdent > 1.0)
            throw new ParseFailureException("Minimum query identity fraction must be between 0 and 1");
        if (this.minPctQuery < 0.0 || this.minPctQuery > 100.0)
            throw new ParseFailureException("Minimum query percentation must be between 0 and 100.");
        // Get the list of SOURS.
        log.info("Reading universal roles from {}.", this.sourFile);
        Set<String> sourRoles = TabbedLineReader.readSet(this.sourFile, "1");
        log.info("{} universal roles found.", sourRoles.size());
        // Load the input genomes.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " not found in file system.");
        log.info("Loading genomes from {} source {}.", this.sourceType.toString(), this.inDir);
        this.genomes = this.sourceType.create(this.inDir);
        // Create the reporting objects.
        this.reporters = new ArrayList<>(UniProfileReporter.Type.values().length);
        for (UniProfileReporter.Type reportType : UniProfileReporter.Type.values()) {
            UniProfileReporter reporter = reportType.create(this.outDir);
            this.reporters.add(reporter);
        }
        // Finally, load the profiles.
        log.info("Loading protein profiles from {}.", this.profileDir);
        this.profiles = new ProteinProfiles(this.profileDir, sourRoles);
        this.roleMap = this.profiles.roleMap();
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the temporary FASTA file.
        File fastaFile = File.createTempFile("pegs", ".faa");
        try {
            // Initialize the reports.
            log.info("Initializing reports.");
            this.reporters.stream().forEach(x -> x.openReport(this.roleMap));
            // Turn off BLAST logging.
            BlastDB.setQuiet();
            // Create the BLAST parameters.
            BlastParms parms = new BlastParms().maxE(this.eValue).minPercent(this.minPctIdentity)
                    .minQueryBitScore(this.minQbsc).minQueryIdentity(this.minQIdent).pctLenOfQuery(this.minPctQuery);
            // This set will track the roles found in each genome.
            Set<String> foundRoles = new HashSet<>(this.roleMap.size());
            // Loop through the genomes.
            int gCount = 0;
            int gTotal = this.genomes.size();
            long start = System.currentTimeMillis();
            for (Genome genome : this.genomes) {
                gCount++;
                log.info("Processing genome {} of {}: {}.", gCount, gTotal, genome);
                this.reporters.stream().forEach(x -> x.openGenome(genome));
                // We want to loop through the genome, finding roles in pegs.
                // We will keep a set of the roles found in here.
                foundRoles.clear();
                // We will create a FASTA of the pegs without found roles in here.
                try (FastaOutputStream proteinFastaStream = new FastaOutputStream(fastaFile)) {
                    for (Feature feat : genome.getPegs()) {
                        boolean found = false;
                        for (Role role : feat.getUsefulRoles(this.roleMap)) {
                            String roleId = role.getId();
                            this.reporters.stream().forEach(x -> x.recordAnnotation(feat, roleId));
                            foundRoles.add(roleId);
                            found = true;
                        }
                        if (! found) {
                            // Here this protein does not have a desired role in it, so save it in the FASTA file.
                            Sequence seq = new Sequence(feat.getId(), feat.getPegFunction(), feat.getProteinTranslation());
                            proteinFastaStream.write(seq);
                        }
                    }
                }
                // Get a list of the missing roles.
                Set<String> missingRoles = this.roleMap.keySet().stream().filter(x -> ! foundRoles.contains(x)).collect(Collectors.toSet());
                if (missingRoles.isEmpty())
                    log.info("NO MISSING ROLES FOUND.");
                else {
                    log.info("Missing {} roles:  {} found. Applying profiles.", missingRoles.size(), foundRoles.size());
                    // Create a BLAST database from the protein FASTA.  This particular method erases any existing support files.
                    ProteinBlastDB blastDB = ProteinBlastDB.create(fastaFile);
                    Map<String, Set<BlastHit>> profileHits = this.profiles.profile(missingRoles, blastDB, parms);
                    // Now record the profile hits.
                    for (String role : missingRoles) {
                        Set<BlastHit> hits = profileHits.getOrDefault(role, NO_HITS);
                        for (BlastHit hit : hits)
                            this.reporters.stream().forEach(x -> x.recordProfileHit(role, hit));
                    }
                }
                // Allow the reports to finish processing this genome.
                this.reporters.stream().forEach(x -> x.closeGenome());
                if (log.isInfoEnabled()) {
                    double seconds = (System.currentTimeMillis() - start) / (1000.0 * gCount);
                    log.info("{} genomes remaining, {} seconds / genome.", gTotal - gCount, seconds);
                }
            }
            // Finish up the reports.
            this.reporters.stream().forEach(x -> x.closeReport());
        } finally {
            // Close the reporters to flush out the final output.
            this.reporters.stream().forEach(x -> x.close());
            // Delete the temporary files.
            log.info("Cleaning up blast database for {}.", fastaFile.getAbsolutePath());
            FileUtils.forceDelete(fastaFile);
            String fastaPath = fastaFile.getPath();
            for (String suffix : ProteinBlastDB.SUFFIXES)
                FileUtils.forceDelete(new File(fastaPath + suffix));
        }
    }

}
