//******************************************************************************
//
// File:    SummarizeDatabase.java
// Package: edu.rit.compbio.seq.test
// Unit:    Class edu.rit.compbio.seq.test.SummarizeDatabase
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

import java.util.HashMap;
import java.util.Map;

/**
 * Class SummarizeDatabase is a unit test main program for classes {@linkplain
 * edu.rit.compbio.seq.ProteinDatabase} and {@linkplain
 * edu.rit.compbio.seq.ProteinSequence}. The program prints summary statistics
 * for the protein sequence database, including number of sequences, total
 * sequence length, smallest sequence length, largest sequence length, and a
 * histogram of the sequence lengths.
 * <P>
 * Usage: java edu.rit.compbio.seq.test.SummarizeDatabase <I>databasefile</I>
 * <I>indexfile</I> <I>binsize</I>
 * <BR><I>databasefile</I> = Protein sequence database file
 * <BR><I>indexfile</I> = Protein sequence index file
 * <BR><I>binsize</I> = Sequence length histogram bin size
 *
 * @author  Alan Kaminsky
 * @version 09-Jul-2008
 */
public class SummarizeDatabase
	{

// Prevent construction.

	private SummarizeDatabase()
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
		int binsize = Integer.parseInt (args[2]);

		// Set up protein sequence database.
		ProteinDatabase db = new ProteinDatabase (databasefile, indexfile);

		// Histogram: Mapping from bin number (Integer) to count (Long).
		Map<Integer,Long> histogram = new HashMap<Integer,Long>();

		// Sequence length statistics.
		int minLength = Integer.MAX_VALUE;
		int maxLength = 0;

		// Get all protein sequences and accumulate statistics.
		for (long i = 0; i < db.getProteinCount(); ++ i)
			{
			ProteinSequence seq = db.getProteinSequence (i);
			int len = seq.length();
			minLength = Math.min (minLength, len);
			maxLength = Math.max (maxLength, len);
			int bin = len/binsize;
			long count = histogram.containsKey (bin) ? histogram.get (bin) : 0L;
			++ count;
			histogram.put (bin, count);
			}

		// Print results.
		System.out.println ("Database file: "+databasefile);
		System.out.println ("Database index file: "+indexfile);
		System.out.println ("Number of sequences: "+db.getProteinCount());
		System.out.println ("Total sequence length: "+db.getDatabaseLength());
		System.out.println ("Minimum sequence length: "+minLength);
		System.out.println ("Maximum sequence length: "+maxLength);
		System.out.println ("Length\tCount");
		int bin = 0;
		int binlb = bin*binsize;
		int binub = binlb + binsize - 1;
		while (binlb <= maxLength)
			{
			long count = histogram.containsKey (bin) ? histogram.get (bin) : 0L;
			System.out.println (binlb+"-"+binub+"\t"+count);
			bin += 1;
			binlb += binsize;
			binub += binsize;
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
		System.err.println ("Usage: java edu.rit.compbio.seq.test.SummarizeDatabase <databasefile> <indexfile> <database>");
		System.err.println ("<databasefile> = Protein sequence database file");
		System.err.println ("<indexfile> = Protein sequence index file");
		System.err.println ("<binsize> = Sequence length histogram bin size");
		System.exit (1);
		}

	}
