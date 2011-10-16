//******************************************************************************
//
// File:    Test13.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test13
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.util.Random;
import edu.rit.util.Range;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;

/**
 * Class Test13 is a unit test main program for the SMP barrier functionality in
 * package {@linkplain edu.rit.pj}. The program does the O(<I>N</I><SUP>3</SUP>)
 * Floyd's Algorithm on an <I>N</I>x<I>N</I> matrix filled with random data. The
 * final matrix is stored in a binary file using java.io.DataOutput. Running
 * time measurements are for the computational core only, not including
 * initializing the matrix or writing the output file.
 * <P>
 * Usage: java -Dpj.nt=<I>K</I> edu.rit.pj.test.Test13 <I>seed</I> <I>N</I>
 * <I>file</I>
 * <BR><I>K</I> = Number of parallel threads
 * <BR><I>seed</I> = Random seed
 * <BR><I>N</I> = Matrix size (<I>N</I>x<I>N</I>)
 * <BR><I>file</I> = Output file name
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
public class Test13
	{

// Prevent construction.

	private Test13()
		{
		}

// Global variables.

	// Command line arguments.
	static long seed;
	static int N;
	static File file;

	// Distance matrix.
	static double[][] D;

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 3) usage();
		seed = Long.parseLong (args[0]);
		N = Integer.parseInt (args[1]);
		file = new File (args[2]);

		// Set up distance matrix.
		D = new double [N] [N+16]; // 128 bytes extra padding in each row
		Random prng = Random.getInstance (seed);
		for (int i = 0; i < N; ++ i)
			{
			double[] D_i = D[i];
			for (int j = 0; j < N; ++ j)
				{
				D_i[j] = prng.nextDouble();
				}
			}

		// Start timing.
		long time = -System.currentTimeMillis();

		// Parallel region.
		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				// Determine this thread's middle loop range.
				int size = getThreadCount();
				int rank = getThreadIndex();
				Range range = new Range (0, N-1) .subrange (size, rank);
				int lb = range.lb();
				int ub = range.ub();

				// Outer loop, N iterations.
				for (int i = 0; i < N; ++ i)
					{
					double[] D_i = D[i];

					// Middle loop, N iterations in parallel.
					for (int j = lb; j <= ub; ++ j)
						{
						double[] D_j = D[j];

						// Inner loop, N iterations.
						for (int k = 0; k < N; ++ k)
							{
							D_j[k] = Math.min (D_j[k], D_j[i] + D_i[k]);
							}
						}

					// Team threads wait for each other each outer loop
					// iteration.
					barrier();
					}
				}
			});

		// Stop timing.
		time += System.currentTimeMillis();

		// Write output file.
		DataOutputStream out =
			new DataOutputStream
				(new BufferedOutputStream
					(new FileOutputStream (file)));
		for (int i = 0; i < N; ++ i)
			{
			double[] D_i = D[i];
			for (int j = 0; j < N; ++ j)
				{
				out.writeDouble (D_i[j]);
				}
			}
		out.close();

		System.out.println (time);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.nt=<K> edu.rit.pj.test.Test13 <seed> <N> <file>");
		System.err.println ("<K> = Number of parallel threads");
		System.err.println ("<seed> = Random seed");
		System.err.println ("<N> = Matrix size (NxN)");
		System.err.println ("<file> = Output file name");
		System.exit (1);
		}

	}
