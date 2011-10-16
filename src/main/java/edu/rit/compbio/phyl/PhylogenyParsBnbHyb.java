//*****************************************************************************
//
// File:    PhylogenyParsBnbHyb.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.PhylogenyParsBnbHyb
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

import edu.rit.mp.IntegerBuf;
import edu.rit.mp.ObjectBuf;

import edu.rit.mp.buf.IntegerItemBuf;
import edu.rit.mp.buf.ObjectItemBuf;

import edu.rit.pj.Comm;
import edu.rit.pj.CommStatus;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.pj.reduction.IntegerOp;

import edu.rit.pj.replica.ReplicatedInteger;

import edu.rit.util.Range;

import java.io.File;
import java.io.IOException;

/**
 * Class PhylogenyParsBnbHyb is a hybrid parallel program for maximum parsimony
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
 * Usage: java [ -Dpj.np=<I>Kp</I> ] [ -Dpj.nt=<I>Kt</I> ] [
 * -Dpj.schedule=<I>schedule</I> ] edu.rit.compbio.phyl.PhylogenyParsBnbHyb
 * <I>infile</I> <I>outdir</I> [ <I>N</I> [ <I>T</I> ] ]
 * <BR><I>Kp</I> = Number of parallel processes (default: 1)
 * <BR><I>Kt</I> = Number of parallel threads per process (default: number of
 * CPUs)
 * <BR><I>schedule</I> = Load balancing schedule (default: dynamic(1))
 * <BR><I>infile</I> = Input DNA sequence list file name
 * <BR><I>outdir</I> = Output directory name
 * <BR><I>N</I> = Number of DNA sequences to use (default: all)
 * <BR><I>T</I> = Number of trees to report (default: 100)
 *
 * @author  Alan Kaminsky
 * @version 21-Nov-2008
 */
public class PhylogenyParsBnbHyb
	{

// Prevent construction.

	private PhylogenyParsBnbHyb()
		{
		}

// Hidden constants.

	// Maximum level of the search graph at which to partition the search.
	private static final int MAX_START_LEVEL = 6;

// Global variables.

	// World communicator.
	private static Comm world;
	private static int size;
	private static int rank;

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

	// Shared, replicated <bound> variable.
	private static ReplicatedInteger bound;

	// Maximum parsimony search results.
	private static MaximumParsimonyResults globalResults;

	// Search graph starting level and number of vertices at that level.
	private static int startLevel;
	private static int vertexCount;

	// Number of worker threads in this process.
	private static int Kt;

	// Number of worker threads in all lower-ranked processes.
	private static int Klower;

	// Number of worker threads in all processes.
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

		// Initialize world communicator.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

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

		// Set up shared, replicated <bound> variable. Initial bound is the
		// UPGMA parsimony score, reduced by the score from the uninformative
		// sites.
		bound =
			MaximumParsimonyBnbHyb.createBoundVariable
				(upgmaScore - uninformativeScore,
				 world,
				 Integer.MAX_VALUE);

		// Set up maximum parsimony results object.
		globalResults = new MaximumParsimonyResults (T);

		// Determine search graph starting level and number of vertices at that
		// level.
		startLevel = Math.min (MAX_START_LEVEL, N - 1);
		vertexCount = 1;
		for (int i = 2*startLevel - 1; i > 1; i -= 2) vertexCount *= i;

		// Determine number of worker threads in this process.
		Kt = ParallelTeam.getDefaultThreadCount();

		// Determine number of worker threads in all lower-ranked processes.
		IntegerItemBuf Kbuf = IntegerBuf.buffer (Kt);
		world.exclusiveScan (Kbuf, IntegerOp.SUM, null);
		Klower = Kbuf.item;

		// Determine number of worker threads in all processes. (Only process 0
		// determines this.)
		Kbuf.item = Kt;
		world.reduce (0, Kbuf, IntegerOp.SUM);
		K = Kbuf.item;

		long t2 = System.currentTimeMillis();

		// Run the branch-and-bound search using the master-worker pattern for
		// load balancing.
		new ParallelTeam (rank == 0 ? Kt+1 : Kt) .execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				int thr = getThreadIndex();
				if (thr == Kt)
					{
					master();
					}
				else
					{
					worker (Klower + thr);
					}
				}
			});

		// Add the score from the uninformative sites back in.
		globalResults.score (globalResults.score() + uninformativeScore);

		long t3 = System.currentTimeMillis();

		// Master process reports results.
		if (rank == 0)
			{
			Results.report
				(/*directory      */ outdir,
				 /*programName    */ "edu.rit.compbio.phyl.PhylogenyParsBnbHyb",
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
			}

		// Stop timing.
		long t4 = System.currentTimeMillis();
		System.out.println ((t2-t1)+" msec pre "+rank);
		System.out.println ((t3-t2)+" msec calc "+rank);
		System.out.println ((t4-t3)+" msec post "+rank);
		System.out.println ((t4-t1)+" msec total "+rank);
		}

// Hidden operations.

	/**
	 * Perform the master section. Only process 0 does this.
	 *
	 * @param  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private static void master()
		throws IOException
		{
		// Set up schedule for partitioning the search graph into tasks.
		IntegerSchedule schedule =
			IntegerSchedule.runtime
				(IntegerSchedule.dynamic (1));
		schedule.start (K, new Range (0, vertexCount - 1));

		// Set up buffer for receiving results from worker.
		ObjectItemBuf<MaximumParsimonyResults> resultsBuf = ObjectBuf.buffer();
		Range thrRange = new Range (0, K - 1);

		// Number of active workers.
		int workerCount = K;

		// Repeat until no more work.
		while (workerCount > 0)
			{
			// Receive results from any process and any worker thread.
			CommStatus status = world.receive (null, thrRange, resultsBuf);
			int workerRank = status.fromRank;
			int workerThr = status.tag;

			// Send next task to that worker thread.
			Range vertexRange = schedule.next (workerThr);
			world.send (workerRank, workerThr, ObjectBuf.buffer (vertexRange));
			if (vertexRange == null)
				{
				-- workerCount;
				}

			// Record results.
			globalResults.addAll (resultsBuf.item);
			}
		}

	/**
	 * Perform the worker section.
	 *
	 * @param  thr  Thread index in the range 0 .. K - 1.
	 *
	 * @param  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private static void worker
		(int thr)
		throws IOException
		{
		// Set up object for sending results to master.
		MaximumParsimonyResults results = new MaximumParsimonyResults (T);

		// Set up buffer for receiving task from master.
		ObjectItemBuf<Range> vertexBuf = ObjectBuf.buffer();
		Range vertexRange;

		// Set up maximum parsimony search algorithm object.
		MaximumParsimonyBnbHyb searcher =
			new MaximumParsimonyBnbHyb (excisedList, bound, results);

		// Repeat until no more work.
		for (;;)
			{
			// Send results of previous partition to the master. Message tag =
			// thread index. A copy is sent so that the master can be reading
			// the copy while the worker is changing the original.
			world.send
				(0,
				 thr,
				 ObjectBuf.buffer (new MaximumParsimonyResults (results)));

			// Get next task from the master. Message tag = thread index. If
			// null, no more work.
			world.receive (0, thr, vertexBuf);
			vertexRange = vertexBuf.item;
			if (vertexRange == null) break;

			// Search the partition.
			results.clear();
			searcher.findTrees (startLevel, vertexRange.lb(), vertexRange.ub());
			}
		}

//	/**
//	 * Helper class for sending maximum parsimony results plus vertex index from
//	 * worker to master.
//	 */
//	private static class WorkerResults
//		extends MaximumParsimonyResults
//		{
//		public int vertex;
//
//		public WorkerResults
//			(int treeStoreLimit)
//			{
//			super (treeStoreLimit);
//			}
//		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java [-Dpj.np=<Kp>] [-Dpj.nt=<Kt>] [-Dpj.schedule=<schedule>] edu.rit.compbio.phyl.PhylogenyParsBnbHyb <infile> <outdir> [<N> [<T>]]");
		System.err.println ("<Kp> = Number of parallel processes (default: 1)");
		System.err.println ("<Kt> = Number of parallel threads per process (default: number of CPUs)");
		System.err.println ("<schedule> = Load balancing schedule (default: dynamic(1))");
		System.err.println ("<infile> = Input DNA sequence list file name");
		System.err.println ("<outdir> = Output directory name");
		System.err.println ("<N> = Number of DNA sequences to use (default: all)");
		System.err.println ("<T> = Number of trees to report (default: 100)");
		System.exit (1);
		}

	}
