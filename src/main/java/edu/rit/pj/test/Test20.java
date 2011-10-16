//******************************************************************************
//
// File:    Test20.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test20
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

package edu.rit.pj.test;

import edu.rit.pj.Comm;
import edu.rit.pj.LongForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.pj.reduction.IntegerOp;

import edu.rit.pj.replica.ReplicatedInteger;

import edu.rit.util.LongRange;
import edu.rit.util.Random;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Map;

/**
 * Class Test20 is a hybrid parallel unit test program for class {@linkplain
 * edu.rit.pj.replica.ReplicatedInteger}. The program runs on a number of
 * processors in the cluster with a number of threads in each process. The
 * program generates a sequence of random integers. Each random integer updates
 * a ReplicatedInteger variable using minimum as the reduction operator; thus,
 * the variable always holds the smallest random integer. Once all processes
 * have finished, each process prints the value of its variable.
 * <P>
 * Usage: java -Dpj.np=<I>np</I> -Dpj.nt=<I>nt</I> edu.rit.pj.test.Test20
 * <I>seed</I> <I>N</I>
 * <BR><I>np</I> = Number of processes
 * <BR><I>nt</I> = Number of CPUs per process
 * <BR><I>seed</I> = Random seed
 * <BR><I>N</I> = Number of random integers
 *
 * @author  Alan Kaminsky
 * @version 12-Sep-2008
 */
public class Test20
	{

// Prevent construction.

	private Test20()
		{
		}

// Global variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static long seed;
	static long N;

	// Replicated, shared reduction variable.
	static ReplicatedInteger littlest;

	// Range of iterations for this process.
	static LongRange range;

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Initialize middleware.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (args.length != 2) usage();
		seed = Long.parseLong (args[0]);
		N = Long.parseLong (args[1]);

		// Set up replicated, shared reduction variable.
		littlest = new ReplicatedInteger
			(/*op          */ IntegerOp.MINIMUM,
			 /*initialValue*/ Integer.MAX_VALUE,
			 /*tag         */ 1);

		// Partition computation among processes.
		range = new LongRange (0, N-1) .subrange (size, rank);

		// Perform computation in parallel threads.
		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				execute (range.lb(), range.ub(), new LongForLoop()
					{
					Random prng = Random.getInstance (seed);
					public void run (long first, long last) throws Exception
						{
						prng.setSeed (seed);
						prng.skip (first);
						for (long i = first; i <= last; ++ i)
							{
							littlest.reduce (prng.nextInt (Integer.MAX_VALUE));
							}
						}
					});
				}
			});

		// Wait until all processes have finished.
		world.barrier();

		// Wait two more seconds to let any in-progress updates finish.
		Thread.sleep (2000L);

		// Print replicated, shared reduction variable.
		System.out.println ("littlest = "+littlest+" ("+rank+")");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<np> -Dpj.nt=<nt> edu.rit.pj.test.Test20 <seed> <N>");
		System.err.println ("<np> = Number of processes");
		System.err.println ("<nt> = Number of CPUs per process");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<N> = Number of random integers");
		System.exit (1);
		}

	}
