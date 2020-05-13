/**
 *
 */
package org.theseed.genome.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
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
 *
 * @author Bruce Parrello
 *
 */
public class ProfileMakeProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected Logger log = LoggerFactory.getLogger(ProfileMakeProcessor.class);
    /** role map */
    private RoleMap roleMap;
    /** map file writer */
    private PrintWriter mapOutput;
    /** pattern for finding local ID line */
    private static final Pattern ID_LINE = Pattern.compile("\\s+local str \"([^\"]+)\"");
    /** pattern for finding title line */
    private static final Pattern TITLE_LINE = Pattern.compile("\\s+title \"\\S+ ([^\"]+)\"");
    /** set of roles already processed */
    private Set<String> processedRoles;

    // COMMAND-LINE OPTIONS

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
        // Create the process set.
        this.processedRoles = new HashSet<String>(this.roleMap.size());
        return true;
    }

    @Override
    public void run() {
        try {
            // Create the map file and walk through all the smp files in the input.
            try (Stream<Path> paths = Files.walk(Paths.get(this.inDir.getPath()), 2)) {
                paths.map(p -> p.toString()).filter(f -> f.endsWith(".smp")).forEach(f -> this.process(f));
                this.mapOutput.flush();
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
        } finally {
            this.mapOutput.close();
        }
    }

    /**
     * Process a single profile file.  The file is a simple text file, and it contains among other things
     * the description of the role.  This information is used to update the map file, find the role ID,
     * and modify the file so that the role ID is returned as the sequence ID.
     *
     * @param profileFile	name of the profile file to read
     *
     * @throws IOException
     */
    private void process(String profileFile) {
        try {
            // Read all the lines from the file.
            List<String> lines = Files.readAllLines(Paths.get(profileFile));
            log.debug("Processing file {}.", profileFile);
            // We need to search for the title and ID lines, and we need to memorize the local ID.
            String description = null;
            String localId = null;
            int localIdLine = 0;
            for (int i = 0; description == null; i++) {
                if (i >= lines.size())
                    throw new IndexOutOfBoundsException("Invalid profile:  local ID or title missing in " + profileFile + ".");
                // Get the current line.
                String line = lines.get(i);
                if (localId == null) {
                    // Here we have not yet found the local ID for this profile.
                    Matcher lineMatcher = ID_LINE.matcher(line);
                    if (lineMatcher.matches()) {
                        localIdLine = i;
                        localId = lineMatcher.group(1);
                    }
                } else {
                    // Here we are looking for the title line.
                    Matcher lineMatcher = TITLE_LINE.matcher(line);
                    if (lineMatcher.matches()) {
                        description = lineMatcher.group(1);
                        // Use the description to find the role.
                        Role role = roleMap.getByName(description);
                        if (role != null && ! this.processedRoles.contains(role.getId())) {
                            // Here the role is one we care about.  Replace the local ID with the role
                            // ID in the file and write it to the output directory.
                            this.fixLine(lines, localId, i, role);
                            this.fixLine(lines, localId, localIdLine, role);
                            File outFile = new File(this.outDir, role.getId() + ".smp");
                            log.info("Creating {} from {} for: {}.", outFile, profileFile, description);
                            try (FileWriter newProfileOut = new FileWriter(outFile)) {
                                for (String oLine : lines)
                                    newProfileOut.write(oLine);
                            }
                            // Store the mapping in the mapping file.
                            this.mapOutput.format("%s\t%s%n", role.getId(), description);
                            // Remember this role.
                            this.processedRoles.add(role.getId());
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Replace the local ID with the role ID in the specified line.
     *
     * @param lines		line list
     * @param localId	local ID to replace
     * @param i			index of line to update
     * @param role		role whose ID is to be replaced
     */
    private void fixLine(List<String> lines, String localId, int i, Role role) {
        lines.set(i, StringUtils.replaceOnce(lines.get(i), localId, role.getId()));
    }

}
