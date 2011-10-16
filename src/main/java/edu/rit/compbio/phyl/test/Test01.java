//*****************************************************************************
//
// File:    Test01.java
// Package: edu.rit.compbio.phyl.test
// Unit:    Class edu.rit.compbio.phyl.test.Test01
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

import edu.rit.compbio.phyl.DnaSequenceList;

import java.io.File;

/**
 * Class Test01 is a unit test program for class {@linkplain
 * edu.rit.compbio.phyl.DnaSequenceList}. It reads DNA sequences from a
 * PHYLIP-formatted input file, then prints the DNA sequences to a
 * PHYLIP-formatted output file.
 * <P>
 * Usage: java edu.rit.compbio.phyl.test.Test01 <I>infile</I> <I>outfile</I>
 *
 * @author  Alan Kaminsky
 * @version 14-Jul-2008
 */
public class Test01
	{

// Prevent construction.

	private Test01()
		{
		}

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		if (args.length != 2) usage();
		File infile = new File (args[0]);
		File outfile = new File (args[1]);
		DnaSequenceList list = DnaSequenceList.read (infile);
		list.write (outfile, 70, true, false);
		}

// Hidden operations.

	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.phyl.test.Test01 <infile> <outfile>");
		System.exit (1);
		}

	}
