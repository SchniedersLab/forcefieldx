//*****************************************************************************
//
// File:    PhylogenyParsBnbSeq.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.PhylogenyParsBnbSeq
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

package edu.rit.compbio.phyl;

import edu.rit.pj.Comm;

import java.io.File;

/**
 * Class PhylogenyParsBnbSeq is a sequential program for maximum parsimony
 * phylogenetic tree construction using branch-and-bound search. The program
 * reads a {@linkplain DnaSequenceList} from a file in interleaved PHYLIP
 * format, constructs a list of one or more maximum parsimony phylogenetic trees
 * using branch-and-bound search, and stores the results in an output directory.
 * If the third command line argument <I>N</I> is specified, only the first
 * <I>N</I> DNA sequences in the file are used; if <I>N</I> is not specified,
 * all DNA sequences in the file are used. If the fourth command line argument
 * <I>T</I> is specified, the program will only report the first <I>T</I>
 * maximum parsimony phylogenetic trees it finds; if <I>T</I> is not specified,
 * the default is <I>T</I> = 100.
 * <P>
 * To examine the results, use a web browser to look at the
 * <TT>"index.html"</TT> file in the output directory. For further information,
 * see class {@linkplain Results}.
 * <P>
 * Usage: java edu.rit.compbio.phyl.PhylogenyParsBnbSeq <I>infile</I>
 * <I>outdir</I> [ <I>N</I> [ <I>T</I> ] ]
 * <BR><I>infile</I> = Input DNA sequence list file name
 * <BR><I>outdir</I> = Output directory name
 * <BR><I>N</I> = Number of DNA sequences to use (default: all)
 * <BR><I>T</I> = Number of trees to report (default: 100)
 *
 * @author  Alan Kaminsky
 * @version 21-Nov-2008
 */
public class PhylogenyParsBnbSeq
	{

// Prevent construction.

	private PhylogenyParsBnbSeq()
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
		// Start timing.
		long t1 = System.currentTimeMillis();

		Comm.init (args);

		// Parse command line arguments.
		if (args.length < 2 || args.length > 4) usage();
		File infile = new File (args[0]);
		File outdir = new File (args[1]);
		int T = 100;
		if (args.length >= 4) T = Integer.parseInt (args[3]);

		// Read DNA sequence list from file and truncate to N sequences if
		// necessary.
		DnaSequenceList seqList = DnaSequenceList.read (infile);
		int N = seqList.length();
		if (args.length >= 3) N = Integer.parseInt (args[2]);
		seqList.truncate (N);

		// Run the UPGMA algorithm to get an approximate solution. Calculate its
		// parsimony score.
		DnaSequenceTree upgmaTree =
			Upgma.buildTree (seqList, new JukesCantorDistance());
		int upgmaScore = FitchParsimony.computeScore (upgmaTree);

		// Put the DNA sequence list in descending order of tip node branch
		// length in the UPGMA tree.
		DnaSequenceList sortedList = upgmaTree.toList();

		// Excise uninformative sites.
		DnaSequenceList excisedList = new DnaSequenceList (sortedList);
		int uninformativeScore = excisedList.exciseUninformativeSites();

		long t2 = System.currentTimeMillis();

		// Run the branch-and-bound search. Use the UPGMA parsimony score
		// (reduced by score from uninformative sites) as the initial bound.
		MaximumParsimonyResults results = new MaximumParsimonyResults (T);
		MaximumParsimonyBnbSeq searcher =
			new MaximumParsimonyBnbSeq
				(excisedList, upgmaScore - uninformativeScore, results);
		searcher.findTrees();
		results.score (results.score() + uninformativeScore);

		long t3 = System.currentTimeMillis();

		// Report results.
		Results.report
			(/*directory      */ outdir,
			 /*programName    */ "edu.rit.compbio.phyl.PhylogenyParsBnbSeq",
			 /*hostName       */ Comm.world().host(),
			 /*K              */ 1,
			 /*infile         */ infile,
			 /*originalSeqList*/ seqList,
			 /*sortedSeqList  */ sortedList,
			 /*initialBound   */ upgmaScore,
			 /*treeStoreLimit */ T,
			 /*results        */ results,
			 /*t1             */ t1,
			 /*t2             */ t2,
			 /*t3             */ t3);

		// Stop timing.
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
		System.err.println ("Usage: java edu.rit.compbio.phyl.PhylogenyParsBnbSeq <infile> <outdir> [<N> [<T>]]");
		System.err.println ("<infile> = Input DNA sequence list file name");
		System.err.println ("<outdir> = Output directory name");
		System.err.println ("<N> = Number of DNA sequences to use (default: all)");
		System.err.println ("<T> = Number of trees to report (default: 100)");
		System.exit (1);
		}

	}
