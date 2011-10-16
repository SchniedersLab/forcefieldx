//******************************************************************************
//
// File:    Test04.java
// Package: edu.rit.compbio.phyl.test
// Unit:    Class edu.rit.compbio.phyl.test.Test04
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
import edu.rit.compbio.phyl.JukesCantorDistance;
import edu.rit.compbio.phyl.TreeDrawing;
import edu.rit.compbio.phyl.Upgma;

import edu.rit.draw.Drawing;

import java.io.File;

/**
 * Class Test04 is a unit test program for classes {@linkplain
 * edu.rit.compbio.phyl.DnaSequenceList}, {@linkplain
 * edu.rit.compbio.phyl.DnaSequenceTree}, {@linkplain
 * edu.rit.compbio.phyl.FitchParsimony}, and {@linkplain
 * edu.rit.compbio.phyl.Upgma}. The program:
 * <OL TYPE=1>
 * <LI>
 * Reads a {@linkplain DnaSequenceList} from the input file in interleaved
 * PHYLIP format.
 * <LI>
 * Constructs a phylogenetic tree using the UPGMA algorithm.
 * <LI>
 * Computes and prints the tree's Fitch parsimony score.
 * <LI>
 * Reorders the list of DNA sequences into descending order of branch length and
 * stores the list in an output file.
 * <LI>
 * Stores a picture of the tree in an output file as a serialized {@linkplain
 * edu.rit.draw.Drawing} object.
 * </OL>
 * <P>
 * The {@linkplain View} program can be used to view the drawing file and
 * save it in an image file in several formats.
 * <P>
 * Usage: java edu.rit.compbio.phyl.test.Test04 <I>infile</I> <I>outfile</I>
 * <I>drawfile</I> [ <I>format</I> ]
 * <BR><I>infile</I> = Input DNA sequence list file name
 * <BR><I>outfile</I> = Output DNA sequence list file name
 * <BR><I>drawfile</I> = Output drawing file name
 * <BR><I>format</I> = DecimalFormat string for branch lengths (default:
 * <TT>"0.00"</TT>)
 *
 * @author  Alan Kaminsky
 * @version 25-Jul-2008
 */
public class Test04
	{

// Prevent construction.

	private Test04()
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
		if (args.length < 3 || args.length > 4) usage();
		File infile = new File (args[0]);
		File outfile = new File (args[1]);
		File drawfile = new File (args[2]);
		TreeDrawing artist = new TreeDrawing();
		if (args.length >= 4) artist.setBranchLengthFormat (args[3]);

		// Read DNA sequence list.
		DnaSequenceList seqlist = DnaSequenceList.read (infile);

		// Construct phylogenetic tree.
		DnaSequenceTree tree =
			Upgma.buildTree (seqlist, new JukesCantorDistance());

		// Compute and print Fitch parsimony score.
		System.out.println ("Score = "+FitchParsimony.computeScore (tree));

		// Reorder DNA sequence list and store in output file.
		tree.toList().write (outfile, 70, true, false);

		// Draw the tree and store in drawing file.
		artist.draw (tree);
		Drawing.write (drawfile);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.phyl.test.Test04 <infile> <outfile> <drawfile> [<format>]");
		System.err.println ("<infile> = Input DNA sequence list file name");
		System.err.println ("<outfile> = Output DNA sequence list file name");
		System.err.println ("<drawfile> = Drawing file name");
		System.err.println ("<format> = DecimalFormat string for branch lengths (default: \"0.00\")");
		System.exit (1);
		}

	}
