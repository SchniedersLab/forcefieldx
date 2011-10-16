//******************************************************************************
//
// File:    ListDatabase.java
// Package: edu.rit.compbio.seq.test
// Unit:    Class edu.rit.compbio.seq.test.ListDatabase
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
 * Class ListDatabase is a unit test main program for classes {@linkplain
 * edu.rit.compbio.seq.ProteinDatabase} and {@linkplain
 * edu.rit.compbio.seq.ProteinSequence}. The program prints the indexes and
 * descriptions of all the protein sequences in the protein sequence database.
 * <P>
 * Usage: java edu.rit.compbio.seq.test.ListDatabase <I>databasefile</I>
 * <I>indexfile</I>
 * <BR><I>databasefile</I> = Protein sequence database file
 * <BR><I>indexfile</I> = Protein sequence index file
 *
 * @author  Alan Kaminsky
 * @version 03-Jul-2008
 */
public class ListDatabase
	{

// Prevent construction.

	private ListDatabase()
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
		File databasefile = new File (args[0]);
		File indexfile = new File (args[1]);

		// Set up protein sequence database.
		ProteinDatabase db = new ProteinDatabase (databasefile, indexfile);

		// Get and print all protein sequences.
		for (long i = 0; i < db.getProteinCount(); ++ i)
			{
			ProteinSequence seq = db.getProteinSequence (i);
			System.out.print (i);
			System.out.print (' ');
			System.out.print (seq.description());
			System.out.println();
			}

		// All done.
		db.close();
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.seq.test.ListDatabase <databasefile> <indexfile>");
		System.err.println ("<databasefile> = Protein sequence database file");
		System.err.println ("<indexfile> = Protein sequence index file");
		System.exit (1);
		}

	}
