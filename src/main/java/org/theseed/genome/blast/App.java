package org.theseed.genome.blast;

import java.util.Arrays;

import org.theseed.sequence.blast.BlastProcessor;
import org.theseed.sequence.blast.MatchProcessor;
import org.theseed.sequence.blast.MatchRunProcessor;
import org.theseed.utils.ICommand;

/**
 * These are various commands related to genome BLAST and matching.
 *
 * match		match RNA sequences to proteins in a genome
 * blast		BLAST a sequence source against a BLAST database
 * profile		run profiles against a DNA BLAST database
 * makepdb		build a profile directory from the raw profiles
 * pverify		verify a profile against a genome
 * mrun			match RNA sequences to proteins in a genome for a whole directory
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "match" :
            processor = new MatchProcessor();
            break;
        case "blast" :
            processor = new BlastProcessor();
            break;
        case "profile" :
            processor = new ProfileProcessor();
            break;
        case "pverify" :
            processor = new ProfileVerifyProcessor();
            break;
        case "makepdb" :
            processor = new ProfileMakeProcessor();
            break;
        case "mrun" :
            processor = new MatchRunProcessor();
            break;
        case "vanalyze" :
            processor = new VerifyAnalyzeProcessor();
            break;
        default :
            throw new IllegalArgumentException("Invalid command " + command);
        }

        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
