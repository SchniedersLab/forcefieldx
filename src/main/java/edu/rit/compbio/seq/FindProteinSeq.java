//******************************************************************************
//
// File:    FindProteinSeq.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.FindProteinSeq
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package edu.rit.compbio.seq;

import edu.rit.pj.Comm;

import java.io.File;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Class FindProteinSeq is a sequential program that finds matches for a protein
 * in a protein sequence database. The query sequence is stored in a file in
 * FASTA format; for further information, see class {@linkplain
 * ProteinSequence}. The protein sequence database is stored in two files, a
 * database file in FASTA format and an index file; for further information, see
 * class {@linkplain ProteinDatabase}. The program uses the Smith-Waterman
 * algorithm to compute a local alignment between the query sequence and each
 * subject sequence in the database. The program uses the BLOSUM-62 protein
 * substitution matrix. The program uses affine gap penalties with a gap
 * existence penalty of &minus;11 and a gap extension penalty of &minus;1. The
 * program prints on the standard output the resulting alignments from highest
 * to lowest score. The program only prints alignments with an <I>E</I>-value
 * below the given threshold; if not specified, the default is 10.
 * <P>
 * The program does all the alignments sequentially, in a single thread. The
 * program uses class {@linkplain ProteinLocalAlignmentSeq} to do the
 * alignments.
 * <P>
 * Usage: java edu.rit.compbio.seq.FindProteinSeq <I>queryfile</I>
 * <I>databasefile</I> <I>indexfile</I> [ <I>expect</I> ]
 * <BR><I>queryfile</I> = Query sequence file
 * <BR><I>databasefile</I> = Protein sequence database file
 * <BR><I>indexfile</I> = Protein sequence index file
 * <BR><I>expect</I> = <I>E</I>-value threshold (default: 10)
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class FindProteinSeq
	{

// Prevent construction.

	private FindProteinSeq()
		{
		}

// Global variables.

	// Command line arguments.
	static File queryfile;
	static File databasefile;
	static File indexfile;
	static double expect;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long t1 = System.currentTimeMillis();

		Comm.init (args);

		// Parse command line arguments.
		if (3 > args.length || args.length > 4) usage();
		queryfile = new File (args[0]);
		databasefile = new File (args[1]);
		indexfile = new File (args[2]);
		expect = 10.0;
		if (args.length == 4) expect = Double.parseDouble (args[3]);

		// Set up query sequence.
		ProteinSequence query = new ProteinSequence (queryfile);

		// Set up protein sequence database.
		ProteinDatabase database =
			new ProteinDatabase (databasefile, indexfile);

		// Set up object to compute alignment statistics.
		AlignmentStats stats =
			new DefaultAlignmentStats (database.getDatabaseLength());

		// Set up list to hold alignments.
		List<Alignment> alignments = new ArrayList<Alignment>();

		// Set up object to perform alignments.
		ProteinLocalAlignment aligner = new ProteinLocalAlignmentSeq();
		aligner.setQuerySequence (query, 0);

		long t2 = System.currentTimeMillis();

		// Align query sequence against every subject sequence.
		for (long id = 0; id < database.getProteinCount(); ++ id)
			{
			ProteinSequence subject = database.getProteinSequence (id);
			aligner.setSubjectSequence (subject, id);
			Alignment a = aligner.align();
			if (stats.eValue (a) <= expect)
				{
				alignments.add (a);
				}
			}

		long t3 = System.currentTimeMillis();

		// Sort alignments into descending order of score.
		Collections.sort (alignments);

		// Set up alignment printer.
		AlignmentPrinter printer = new AlignmentPrinter (System.out, stats);

		// Print query sequence.
		System.out.println ("Query Description:");
		System.out.println (query.description());
		System.out.println ("Length = "+query.length());
		System.out.println();

		// Print summary of each alignment.
		System.out.println ("                                                                Bit  E-");
		System.out.println ("Subject Description                                           Score  Value");
		for (Alignment a : alignments)
			{
			printer.printSummary
				(a, database.getProteinSequence (a.getSubjectId()));
			}
		System.out.println();

		// Print details of each alignment.
		for (Alignment a : alignments)
			{
			printer.printDetails
				(a, query, database.getProteinSequence (a.getSubjectId()));
			}

		// Print various information about the alignment procedure.
		System.out.println ("Query file: "+queryfile);
		System.out.println ("Database file: "+databasefile);
		System.out.println ("Database index file: "+indexfile);
		System.out.println ("Number of sequences: "+database.getProteinCount());
		System.out.println ("Number of matches: "+alignments.size());
		System.out.println ("Query length: "+query.length());
		System.out.println ("Database length: "+database.getDatabaseLength());
		stats.print (System.out);
		System.out.println();

		// All done.
		database.close();

		long t4 = System.currentTimeMillis();
		System.out.println ((t2-t1)+" msec pre");
		System.out.println ((t3-t2)+" msec calc");
		System.out.println ((t4-t3)+" msec post");
		System.out.println ((t4-t1)+" msec total");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.seq.FindProteinSeq <queryfile> <databasefile> <indexfile> [<expect>]");
		System.err.println ("<queryfile> = Query sequence file");
		System.err.println ("<databasefile> = Protein sequence database file");
		System.err.println ("<indexfile> = Protein sequence index file");
		System.err.println ("<expect> = E-value threshold (default: 10)");
		System.exit (1);
		}

	}
