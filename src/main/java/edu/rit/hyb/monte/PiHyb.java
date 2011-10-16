//******************************************************************************
//
// File:    PiHyb.java
// Package: edu.rit.hyb.monte
// Unit:    Class edu.rit.hyb.monte.PiHyb
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

package edu.rit.hyb.monte;

import edu.rit.mp.LongBuf;

import edu.rit.mp.buf.LongItemBuf;

import edu.rit.pj.Comm;
import edu.rit.pj.LongForLoop;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.ParallelRegion;

import edu.rit.pj.reduction.LongOp;
import edu.rit.pj.reduction.SharedLong;

import edu.rit.util.LongRange;
import edu.rit.util.Random;

/**
 * Class PiHyb is a hybrid SMP cluster parallel program that calculates an
 * approximate value for &pi; using a Monte Carlo technique. The program
 * generates a number of random points in the unit square (0,0) to (1,1) and
 * counts how many of them lie within a circle of radius 1 centered at the
 * origin. The fraction of the points within the circle is approximately &pi;/4.
 * <P>
 * Usage: java -Dpj.np=<I>Kp</I> -Dpj.nt=<I>Kt</I> edu.rit.hyb.monte.PiHyb
 * <I>seed</I> <I>N</I>
 * <BR><I>Kp</I> = Number of parallel processes
 * <BR><I>Kt</I> = Number of parallel threads per process
 * <BR><I>seed</I> = Random seed
 * <BR><I>N</I> = Number of random points
 * <P>
 * The computation is performed in parallel in multiple processes with multiple
 * threads per process. The program uses class edu.rit.util.Random for its
 * pseudorandom number generator. To improve performance, each thread has its
 * own pseudorandom number generator, and the program uses the reduction pattern
 * to determine the count. The program uses the "sequence splitting" technique
 * with the pseudorandom number generators to yield results identical to the
 * sequential version. The program measures the computation's running time.
 *
 * @author  Alan Kaminsky
 * @version 29-Feb-2008
 */
public class PiHyb
	{

// Prevent construction.

	private PiHyb()
		{
		}

// Program shared variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static long seed;
	static long N;

	// Subrange of points this process will do.
	static long lb;
	static long ub;

	// Number of points within the unit circle.
	static SharedLong count;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long time = -System.currentTimeMillis();

		// Initialize PJ middleware.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Validate command line arguments.
		if (args.length != 2) usage();
		seed = Long.parseLong (args[0]);
		N = Long.parseLong (args[1]);

		// Determine subrange of points this process will do.
		LongRange range = new LongRange (0, N-1) .subrange (size, rank);
		lb = range.lb();
		ub = range.ub();

		// Generate n random points in the unit square, count how many are in
		// the unit circle.
		count = new SharedLong (0);
		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				execute (lb, ub, new LongForLoop()
					{
					// Set up per-thread PRNG and counter.
					Random prng_thread = Random.getInstance (seed);
					long count_thread = 0;

					// Extra padding to avert cache interference.
					long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
					long pad8, pad9, pada, padb, padc, padd, pade, padf;

					// Parallel loop body.
					public void run (long first, long last)
						{
						// Skip PRNG ahead to index <first>.
						prng_thread.setSeed (seed);
						prng_thread.skip (2 * first);

						// Generate random points.
						for (long i = first; i <= last; ++ i)
							{
							double x = prng_thread.nextDouble();
							double y = prng_thread.nextDouble();
							if (x*x + y*y <= 1.0) ++ count_thread;
							}
						}

					public void finish()
						{
						// Reduce per-thread counts into shared count.
						count.addAndGet (count_thread);
						}
					});
				}
			});

		// Reduce all processes' shared counts into process 0.
		LongItemBuf buf = LongBuf.buffer (count.longValue());
		world.reduce (0, buf, LongOp.SUM);

		// Stop timing.
		time += System.currentTimeMillis();

		// Print results.
		if (rank == 0)
			{
			System.out.println
				("pi = 4 * " + buf.item + " / " + N + " = " +
				 (4.0 * buf.item / N));
			}
		System.out.println (time + " msec " + rank);
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<Kp> -Dpj.nt=<Kt> edu.rit.hyb.monte.PiHyb <seed> <N>");
		System.err.println ("<Kp> = Number of parallel processes");
		System.err.println ("<Kt> = Number of parallel threads per process");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<N> = Number of random points");
		System.exit (1);
		}

	}
