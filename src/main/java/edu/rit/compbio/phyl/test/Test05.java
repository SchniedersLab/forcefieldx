//******************************************************************************
//
// File:    Test05.java
// Package: edu.rit.compbio.phyl.test
// Unit:    Class edu.rit.compbio.phyl.test.Test05
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

package edu.rit.compbio.phyl.test;

import edu.rit.compbio.phyl.DnaSequence;
import edu.rit.compbio.phyl.DnaSequenceList;
import edu.rit.compbio.phyl.DnaSequenceTree;
import edu.rit.compbio.phyl.FitchParsimony;
import edu.rit.compbio.phyl.TreeDrawing;
import edu.rit.compbio.phyl.Upgma;

import edu.rit.draw.Drawing;

import java.io.File;

/**
 * Class Test05 is a unit test program for class {@linkplain
 * edu.rit.compbio.phyl.DnaSequenceList}. The program reads a {@linkplain
 * DnaSequenceList} from the input file in interleaved PHYLIP format, then
 * prints each sequence's name and sites on a single line.
 * <P>
 * Usage: java edu.rit.compbio.phyl.test.Test05 <I>infile</I>
 * <BR><I>infile</I> = Input DNA sequence list file name
 *
 * @author  Alan Kaminsky
 * @version 17-Jul-2008
 */
public class Test05
	{

// Prevent construction.

	private Test05()
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
		if (args.length != 1) usage();
		File infile = new File (args[0]);

		// Read DNA sequence list.
		DnaSequenceList seqlist = DnaSequenceList.read (infile);

		// Print each sequence's name and sites.
		for (DnaSequence seq : seqlist)
			{
			System.out.format ("%-10s%s%n", seq.name(), seq);
			}
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.phyl.test.Test05 <infile>");
		System.err.println ("<infile> = Input DNA sequence list file name");
		System.exit (1);
		}

	}
