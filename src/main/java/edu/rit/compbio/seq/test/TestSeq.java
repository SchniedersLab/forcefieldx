//******************************************************************************
//
// File:    TestSeq.java
// Package: edu.rit.compbio.seq.test
// Unit:    Class edu.rit.compbio.seq.test.TestSeq
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

package edu.rit.compbio.seq.test;

import edu.rit.compbio.seq.Alignment;
import edu.rit.compbio.seq.AlignmentPrinter;
import edu.rit.compbio.seq.AlignmentStats;
import edu.rit.compbio.seq.DefaultAlignmentStats;
import edu.rit.compbio.seq.ProteinLocalAlignment;
import edu.rit.compbio.seq.ProteinLocalAlignmentSeq;
import edu.rit.compbio.seq.ProteinSequence;

import java.io.File;

/**
 * Class TestSeq is a unit test main program for classes {@linkplain
 * edu.rit.compbio.seq.ProteinSequence} and {@linkplain
 * edu.rit.compbio.seq.ProteinLocalAlignmentSeq}. The program reads in a query
 * sequence and a subject sequence, aligns them, and prints the alignment.
 * <P>
 * Usage: java edu.rit.compbio.seq.test.TestSeq <I>queryfile</I>
 * <I>subjectfile</I>
 * <BR><I>queryfile</I> = Query sequence file
 * <BR><I>subjectfile</I> = Subject sequence file
 *
 * @author  Alan Kaminsky
 * @version 02-Jul-2008
 */
public class TestSeq
	{

// Prevent construction.

	private TestSeq()
		{
		}

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 2) usage();
		File queryfile = new File (args[0]);
		File subjectfile = new File (args[1]);

		// Set up query sequence and subject sequence.
		ProteinSequence query = new ProteinSequence (queryfile);
		ProteinSequence subject = new ProteinSequence (subjectfile);

		// Set up object to compute alignment statistics.
		AlignmentStats stats = new DefaultAlignmentStats (subject.length());

		// Perform alignment.
		ProteinLocalAlignment aligner = new ProteinLocalAlignmentSeq();
		aligner.setQuerySequence (query, 0);
		aligner.setSubjectSequence (subject, 0);
		Alignment a = aligner.align();

		// Set up alignment printer.
		AlignmentPrinter printer = new AlignmentPrinter (System.out, stats);

		// Print query sequence.
		System.out.println();
		System.out.println ("Query Description:");
		System.out.println (query.description());
		System.out.println ("Length = "+query.length());
		System.out.println();

		// Print details of alignment.
		printer.printDetails (a, query, subject);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.seq.test.TestSeq <queryfile> <subjectfile>");
		System.err.println ("<queryfile> = Query sequence file");
		System.err.println ("<subjectfile> = Subject sequence file");
		System.exit (1);
		}

	}
