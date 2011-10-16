//*****************************************************************************
//
// File:    PhylogenyParsBnbSmp.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.PhylogenyParsBnbSmp
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
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.pj.reduction.SharedInteger;

import java.io.File;

/**
 * Class PhylogenyParsBnbSmp is an SMP parallel program for maximum parsimony
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
 * Usage: java [ -Dpj.nt=<I>Kt</I> ] [ -Dpj.schedule=<I>schedule</I> ]
 * edu.rit.compbio.phyl.PhylogenyParsBnbSmp <I>infile</I> <I>outdir</I>
 * [ <I>N</I> [ <I>T</I> ] ]
 * <BR><I>Kt</I> = Number of parallel threads (default: number of CPUs)
 * <BR><I>schedule</I> = Load balancing schedule (default: dynamic(1))
 * <BR><I>infile</I> = Input DNA sequence list file name
 * <BR><I>outdir</I> = Output directory name
 * <BR><I>N</I> = Number of DNA sequences to use (default: all)
 * <BR><I>T</I> = Number of trees to report (default: 100)
 *
 * @author  Alan Kaminsky
 * @version 21-Nov-2008
 */
public class PhylogenyParsBnbSmp
	{

// Prevent construction.

	private PhylogenyParsBnbSmp()
		{
		}

// Hidden constants.

	// Maximum level of the search graph at which to partition the search.
	private static final int MAX_START_LEVEL = 6;

// Global variables.

	// Command line arguments.
	private static File infile;
	private static File outdir;
	private static int N;
	private static int T;

	// Original DNA sequence list.
	private static DnaSequenceList seqList;

	// DNA sequence list sorted into descending order of distance.
	private static DnaSequenceList sortedList;

	// Sorted DNA sequence list with uninformative sites removed.
	private static DnaSequenceList excisedList;

	// Shared <bound> variable.
	private static SharedInteger bound;

	// Maximum parsimony search results.
	private static MaximumParsimonyResults globalResults;

	// Search graph starting level and number of vertices at that level.
	private static int startLevel;
	private static int vertexCount;

	// Number of parallel team threads.
	private static int K;

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
		infile = new File (args[0]);
		outdir = new File (args[1]);
		T = 100;
		if (args.length >= 4) T = Integer.parseInt (args[3]);

		// Read DNA sequence list from file and truncate to N sequences if
		// necessary.
		seqList = DnaSequenceList.read (infile);
		N = seqList.length();
		if (args.length >= 3) N = Integer.parseInt (args[2]);
		seqList.truncate (N);

		// Run the UPGMA algorithm to get an approximate solution. Calculate its
		// parsimony score.
		DnaSequenceTree upgmaTree =
			Upgma.buildTree (seqList, new JukesCantorDistance());
		int upgmaScore = FitchParsimony.computeScore (upgmaTree);

		// Put the DNA sequence list in descending order of tip node branch
		// length in the UPGMA tree.
		sortedList = upgmaTree.toList();

		// Excise uninformative sites.
		excisedList = new DnaSequenceList (sortedList);
		int uninformativeScore = excisedList.exciseUninformativeSites();

		// Set up shared <bound> variable. Initial bound is the UPGMA parsimony
		// score, reduced by the score from the uninformative sites.
		bound =
			MaximumParsimonyBnbSmp.createBoundVariable
				(upgmaScore - uninformativeScore);

		// Set up maximum parsimony results object.
		globalResults = new MaximumParsimonyResults (T);

		// Determine search graph starting level and number of vertices at that
		// level.
		startLevel = Math.min (MAX_START_LEVEL, N - 1);
		vertexCount = 1;
		for (int i = 2*startLevel - 1; i > 1; i -= 2) vertexCount *= i;

		long t2 = System.currentTimeMillis();

		// Run the branch-and-bound search in parallel.
		K = ParallelTeam.getDefaultThreadCount();
		new ParallelTeam(K).execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				execute (0, vertexCount - 1, new IntegerForLoop()
					{
					// Per-thread variables.
					MaximumParsimonyResults results;
					MaximumParsimonyBnbSmp searcher;

					// Extra padding to avert cache interference.
					long p0, p1, p2, p3, p4, p5, p6, p7;
					long p8, p9, pa, pb, pc, pd, pe, pf;

					// Use a default dynamic(1) schedule.
					public IntegerSchedule schedule()
						{
						return IntegerSchedule.runtime
							(IntegerSchedule.dynamic (1));
						}

					// Set up results object and searcher object.
					public void start()
						{
						results = new MaximumParsimonyResults (T);
						searcher = new MaximumParsimonyBnbSmp
							(excisedList, bound, results);
						}

					// Perform search.
					public void run (int first, int last)
						{
						searcher.findTrees (startLevel, first, last);
						}

					// Reduce per-thread results into global results.
					public void finish()
						{
						globalResults.addAll (results);
						}
					});
				}
			});

		// Add the score from the uninformative sites back in.
		globalResults.score (globalResults.score() + uninformativeScore);

		long t3 = System.currentTimeMillis();

		// Report results.
		Results.report
			(/*directory      */ outdir,
			 /*programName    */ "edu.rit.compbio.phyl.PhylogenyParsBnbSmp",
			 /*hostName       */ Comm.world().host(),
			 /*K              */ K,
			 /*infile         */ infile,
			 /*originalSeqList*/ seqList,
			 /*sortedSeqList  */ sortedList,
			 /*initialBound   */ upgmaScore,
			 /*treeStoreLimit */ T,
			 /*results        */ globalResults,
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

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java [-Dpj.nt=<Kt>] [-Dpj.schedule=<schedule>] edu.rit.compbio.phyl.PhylogenyParsBnbSmp <infile> <outdir> [<N> [<T>]]");
		System.err.println ("<Kt> = Number of parallel threads (default: number of CPUs)");
		System.err.println ("<schedule> = Load balancing schedule (default: dynamic(1))");
		System.err.println ("<infile> = Input DNA sequence list file name");
		System.err.println ("<outdir> = Output directory name");
		System.err.println ("<N> = Number of DNA sequences to use (default: all)");
		System.err.println ("<T> = Number of trees to report (default: 100)");
		System.exit (1);
		}

	}
