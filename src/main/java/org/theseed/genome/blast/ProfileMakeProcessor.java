/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;

/**
 * This command builds a new profile directory from a role map and raw profiles.  The profile directory
 * will be suitable for use in the "profile" and "match" commands of this application.  Only roles specified
 * in the role map will be loaded into the new profile directory.
 *
 * The positional parameters are the names of the role file, the input directory and the output directory,
 * respectively. The role file should be a standard role definition file (3-column, tab-delimited,
 * role ID in the first column and role name in the third).
 *
 * The following command-line options are supported.
 *
 * -h	display usage
 * -v	display more detailed progress messages
 * -m	minimum plurality count for an acceptable role (the default is 10)
 * -p	minimum percentage of the total count for an acceptable role's plurality count (the default is 90)
 *
 * @author Bruce Parrello
 *
 */
public class ProfileMakeProcessor extends BaseProcessor {

    // FIELDS
    /** suffix for assigned-count files */
    private static final String ASSIGN_SUFFIX = ".assign_cnt";
    /** logging facility */
    protected Logger log = LoggerFactory.getLogger(ProfileMakeProcessor.class);
    /** role map */
    private RoleMap roleMap;
    /** map file writer */
    private PrintWriter mapOutput;
    /** output role map */
    private RoleMap rolesOut;
    /** pattern for finding local ID from file */
    private static final Pattern CLUSTER_ID = Pattern.compile(".+(nr\\d+([/\\\\])cl\\d+)\\" + ASSIGN_SUFFIX);

    // COMMAND-LINE OPTIONS

    @Option(name = "-m", aliases = { "--minPlurality", "--minP" }, metaVar = "12",
            usage = "minimum acceptable plurality count for a profile")
    private int minPlurality;

    @Option(name = "-p", aliases = { "--minPercent", "--minPct" }, metaVar = "95",
            usage = "minimum acceptable percentage of plurality occurrences for a profile")
    private double minPercent;

    /** role map input file */
    @Argument(index = 0, metaVar = "roleFile.tbl", usage = "role map input file", required = true)
    private File roleFile;

    /** input profile directory */
    @Argument(index = 1, metaVar = "inDir", usage = "raw profiles input directory", required = true)
    private File inDir;

    /** output profile directory */
    @Argument(index = 2, metaVar = "outDir", usage = "refined profiles output directory", required = true)
    private File outDir;

    @Override
    protected void setDefaults() {
        this.minPercent = 90.0;
        this.minPlurality = 10;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Read in the role map.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role input file " + this.roleFile + " not found or invalid.");
        log.info("Reading roles from {}.", this.roleFile);
        this.roleMap = RoleMap.load(this.roleFile);
        // Verify the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Profile input directory " + this.inDir + " not found or invalid.");
        log.info("Profiles will be read from {}.", this.inDir);
        // Verify the output directory.
        if (! this.outDir.isDirectory()) {
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else {
            log.info("Erasing output directory {}.", this.outDir);
            FileUtils.cleanDirectory(this.outDir);
        }
        // Create the map file.
        this.mapOutput = new PrintWriter(new File(this.outDir, "_map.tbl"));
        // Create the processed-roles set.
        this.rolesOut = new RoleMap();
        return true;
    }

    @Override
    public void runCommand() throws Exception {
        try {
            // Create the map file and walk through all the smp files in the input.
            try (Stream<Path> paths = Files.walk(Paths.get(this.inDir.getPath()), 2)) {
                paths.map(p -> p.toString()).filter(f -> f.endsWith(ASSIGN_SUFFIX)).forEach(f -> this.process(f));
                this.mapOutput.flush();
            }
            // Write out the role map.
            this.rolesOut.save(new File(this.outDir, "_roles.tbl"));
        } finally {
            this.mapOutput.close();
        }
    }

    /**
     * Process a single profile.  The incoming file name is the assigned-count file. This is a 2-column
     * headerless tab-delimited file.  Each line contains a count followed by a protein function.
     * The first line is the primary function.  If the primary function's count is minPlurality or more, and if
     * it is minPercent of the sum of all the counts, the profile is considered good.  The ".smp" file with the
     * same name will be modified to contain the role ID and description and then copied to the
     * compiled profile directory. If more than one profile specifies the same role, the one with
     * the highest primary function count will be the winner.
     *
     * @param profileFile	name of the assigned-count file to read
     */
    private void process(String profileFile) {
        // Get the cluster identifier for this profile.
        Matcher m = CLUSTER_ID.matcher(profileFile);
        if (! m.matches())
            throw new RuntimeException("Invalid profile file name " + profileFile + ".");
        else {
            // To create the cluster ID, we take the full identifier and replace the path separator found
            // with a period, a complicated operation.
            String localId = StringUtils.replace(m.group(1), m.group(2), ".");
            // Now read the assigned-count file to determine the quality of this profile.
            try (TabbedLineReader countReader = new TabbedLineReader(new File(profileFile), 2)) {
                log.debug("Processing profile {}.", localId);
                Iterator<TabbedLineReader.Line> lineIter = countReader.iterator();
                // Get the plurality line.
                TabbedLineReader.Line line = lineIter.next();
                if (line == null)
                    log.warn("Assignment file for {} is empty.", localId);
                else {
                    String function = line.get(1);
                    int plurality = line.getInt(0);
                    int total = plurality;
                    // Parse the function for a good role.
                    List<Role> roles = Feature.usefulRoles(this.roleMap, function);
                    if (roles.size() > 0 && plurality >= this.minPlurality) {
                        // Here the plurality function contains useful roles.  Compute the total count.
                        while (lineIter.hasNext())
                            total += lineIter.next().getInt(0);
                        log.debug("Plurality for {} is {} out of {}.", localId, plurality, total);
                        if (plurality * 100.0 / total >= this.minPercent) {
                            // We can output this profile. We give it a new file name based on the cluster ID, but
                            // first we need its current name.
                            String profileSmp = StringUtils.substringBeforeLast(profileFile, ASSIGN_SUFFIX) + ".smp";
                            File oldProfile = new File(profileSmp);
                            File newProfile = new File(this.outDir, localId + ".smp");
                            log.info("Copying profile from {} to {}.", oldProfile, newProfile);
                            FileUtils.copyFile(oldProfile, newProfile);
                            // Update the maps.
                            mapOutput.format("%s\t%s%n", localId, function);
                            for (Role role : roles)
                                if (! this.rolesOut.containsKey(role.getId()))
                                    this.rolesOut.put(role);
                        }
                    }
                }
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
        }
    }

}
