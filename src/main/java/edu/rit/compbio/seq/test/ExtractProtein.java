//******************************************************************************
//
// File:    ExtractProtein.java
// Package: edu.rit.compbio.seq.test
// Unit:    Class edu.rit.compbio.seq.test.ExtractProtein
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

import edu.rit.compbio.seq.ProteinDatabase;
import edu.rit.compbio.seq.ProteinSequence;

import java.io.File;

/**
 * Class ExtractProtein is a unit test main program for classes {@linkplain
 * edu.rit.compbio.seq.ProteinDatabase} and {@linkplain
 * edu.rit.compbio.seq.ProteinSequence}. The program extracts the protein
 * sequence at a given index from the protein sequence database and prints it on
 * the standard output in FASTA format.
 * <P>
 * Usage: java edu.rit.compbio.seq.test.ExtractProtein <I>databasefile</I>
 * <I>indexfile</I> <I>index</I>
 * <BR><I>databasefile</I> = Protein sequence database file
 * <BR><I>indexfile</I> = Protein sequence index file
 * <BR><I>index</I> = Protein sequence index
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class ExtractProtein
	{

// Prevent construction.

	private ExtractProtein()
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
		if (args.length != 3) usage();
		File databasefile = new File (args[0]);
		File indexfile = new File (args[1]);
		long index = Long.parseLong (args[2]);

		// Set up protein sequence database.
		ProteinDatabase db = new ProteinDatabase (databasefile, indexfile);

		// Get and print protein sequence.
		ProteinSequence seq = db.getProteinSequence (index);
		System.out.print (seq.description());
		String s = seq.elementsToString();
		int n = s.length();
		int i = 0;
		while (i < n)
			{
			if (i % 60 == 0) System.out.println();
			System.out.print (s.charAt(i));
			++ i;
			}
		System.out.println();

		// All done.
		db.close();
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.seq.test.ExtractProtein <databasefile> <indexfile> <index>");
		System.err.println ("<databasefile> = Protein sequence database file");
		System.err.println ("<indexfile> = Protein sequence index file");
		System.err.println ("<index> = Protein sequence index");
		System.exit (1);
		}

	}
