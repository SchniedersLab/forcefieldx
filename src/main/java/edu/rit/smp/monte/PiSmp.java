//******************************************************************************
//
// File:    PiSmp.java
// Package: edu.rit.smp.monte
// Unit:    Class edu.rit.smp.monte.PiSmp
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

package edu.rit.smp.monte;

//import edu.rit.pj.Comm;
import edu.rit.pj.LongForLoop;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.ParallelRegion;

import edu.rit.pj.reduction.SharedLong;

import java.util.Random;

/**
 * Class PiSmp is an SMP parallel program that calculates an approximate value
 * for &pi; using a Monte Carlo technique. The program generates a number of
 * random points in the unit square (0,0) to (1,1) and counts how many of them
 * lie within a circle of radius 1 centered at the origin. The fraction of the
 * points within the circle is approximately &pi;/4.
 * <P>
 * Usage: java -Dpj.nt=<I>K</I> edu.rit.smp.monte.PiSmp <I>seed</I> <I>N</I>
 * <BR><I>K</I> = Number of parallel threads
 * <BR><I>seed</I> = Random seed
 * <BR><I>N</I> = Number of random points
 * <P>
 * The computation is performed in parallel in multiple threads. All threads
 * share the same pseudorandom number generator and counter. The program
 * measures the computation's running time.
 *
 * @author  Alan Kaminsky
 * @version 29-Feb-2008
 */
public class PiSmp
	{

// Prevent construction.

	private PiSmp()
		{
		}

// Program shared variables.

	// Command line arguments.
	static long seed;
	static long N;

	// Pseudorandom number generator.
	static Random prng;

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
		//Comm.init (args);

		// Start timing.
		long time = -System.currentTimeMillis();

		// Validate command line arguments.
		if (args.length != 2) usage();
		seed = Long.parseLong (args[0]);
		N = Long.parseLong (args[1]);

		// Set up PRNG.
		prng = new Random (seed);

		// Generate n random points in the unit square, count how many are in
		// the unit circle.
		count = new SharedLong (0);
		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				execute (0, N-1, new LongForLoop()
					{
					public void run (long first, long last)
						{
						for (long i = first; i <= last; ++ i)
							{
							double x = prng.nextDouble();
							double y = prng.nextDouble();
							if (x*x + y*y <= 1.0) count.incrementAndGet();
							}
						}
					});
				}
			});

		// Stop timing.
		time += System.currentTimeMillis();

		// Print results.
		System.out.println
			("pi = 4 * " + count + " / " + N + " = " +
			 (4.0 * count.doubleValue() / N));
		System.out.println (time + " msec");
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.nt=<K> edu.rit.smp.monte.PiSmp <seed> <N>");
		System.err.println ("<K> = Number of parallel threads");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<N> = Number of random points");
		System.exit (1);
		}

	}
