package org.theseed.genome.blast;

import java.util.Arrays;

import org.theseed.basic.BaseProcessor;
import org.theseed.sequence.blast.BlastProcessor;
import org.theseed.sequence.blast.ProfileMakeProcessor;
import org.theseed.sequence.blast.ProfileProcessor;
import org.theseed.sequence.blast.ProfileVerifyProcessor;

/**
 * These are various commands related to genome BLAST and matching.
 *
 * match		match RNA sequences to proteins in a genome
 * blast		BLAST a sequence source against a BLAST database
 * profile		run profiles against a DNA BLAST database
 * makepdb		build a profile directory from the raw profiles
 * pverify		verify a profile against a genome
 * mrun			match RNA sequences to proteins in a genome for a whole directory
 * mverify		verify GTI files from "mrun"
 * vanalyze		determine a cutoff value for regression results
 * poi			add points of interest to a genome
 * splice		merge DNA into a reference sequence
 * fastq		hack a pair of FASTQ files into a FASTA file
 * uniProfile	use profiles to verify universal roles
 * proteins		compute pairwise distances between proteins in FASTA files
 * protPrep		package a genome source into protein FASTA files
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        // Parse the parameters.
        switch (command) {
        case "match" -> processor = new MatchProcessor();
        case "blast" -> processor = new BlastProcessor();
        case "profile" -> processor = new ProfileProcessor();
        case "pverify" -> processor = new ProfileVerifyProcessor();
        case "makepdb" -> processor = new ProfileMakeProcessor();
        case "mrun" -> processor = new MatchRunProcessor();
        case "mverify" -> processor = new MatchVerifyProcessor();
        case "vanalyze" -> processor = new VerifyAnalyzeProcessor();
        case "poi" -> processor = new PointOfInterestProcessor();
        case "fastq" -> processor = new FastQProcessor();
        case "uniProfile" -> processor = new UniProfileProcessor();
        case "proteins" -> processor = new ProteinSimsProcessor();
        case "protPrep" -> processor = new ProteinSimsPrepareProcessor();
        default -> throw new IllegalArgumentException("Invalid command " + command);
        }

        processor.parseCommand(newArgs);
        processor.run();
    }
}
