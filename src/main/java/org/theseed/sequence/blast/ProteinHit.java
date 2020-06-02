/**
 *
 */
package org.theseed.sequence.blast;

import org.theseed.locations.FLocation;
import org.theseed.locations.Location;

/**
 * This class represents a protein hit in an RNA match run/
 *
 * @author Bruce Parrello
 */
public class ProteinHit {

    // FIELDS
    private String profile;
    private FLocation loc;

    /**
     * Construct a new protein hit.
     *
     * @param profile		name of profile that hit the protein
     * @param loc			location of the protein
     * @param seqLen		length of the sequence containing the location
     */
    public ProteinHit(String profile, Location loc, int seqLen) {
        if (loc instanceof FLocation) {
            this.loc = (FLocation) loc;
        } else {
            this.loc = (FLocation) loc.converse(seqLen);
        }
        this.profile = profile;
    }

    /**
     * @return the ID of the profile that made the hit
     */
    public String getProfile() {
        return this.profile;
    }

    /**
     * @return the location hit
     */
    public FLocation getLoc() {
        return this.loc;
    }

    @Override
    public String toString() {
        return (profile != null ? profile : "null") + "->" + (loc != null ? loc.toString() : "null");
    }
}
