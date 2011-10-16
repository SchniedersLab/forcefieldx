//******************************************************************************
//
// File:    Test12.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test12
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

import edu.rit.crypto.blockcipher.AES256CipherSmp;

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import edu.rit.util.Range;

/**
 * Class Test12 is a unit test main program for the SMP barrier functionality in
 * package {@linkplain edu.rit.pj}. The program obtains three numbers,
 * <I>N1</I>, <I>N2</I>, and <I>N3</I>, from the command line. Then the program
 * does the following:
 * <PRE>
 * block = Array of N2 16-byte blocks
 * for i = 0 to N1 - 1
 *     parallel for j = 0 to N2 - 1
 *         for k = 0 to N3 - 1
 *             Encrypt block[j] using AES with an all-zero key
 *     barrier
 * </PRE>
 * <P>
 * Usage: java -Dpj.nt=<I>K</I> edu.rit.pj.test.Test12 <I>N1</I> <I>N2</I>
 * <I>N3</I>
 * <BR><I>K</I> = Number of parallel threads
 * <BR><I>N1</I> = Number of outer loop iterations
 * <BR><I>N2</I> = Number of middle loop iterations
 * <BR><I>N3</I> = Number of inner loop iterations
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
public class Test12
	{

// Prevent construction.

	private Test12()
		{
		}

// Global variables.

	// Command line arguments.
	static int N1;
	static int N2;
	static int N3;

	// Array of blocks.
	static byte[][] block;

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long time = -System.currentTimeMillis();

		// Parse command line arguments.
		if (args.length != 3) usage();
		N1 = Integer.parseInt (args[0]);
		N2 = Integer.parseInt (args[1]);
		N3 = Integer.parseInt (args[2]);

		// Set up blocks. Each block has 16 bytes data plus 128 bytes padding to
		// avert cache interference.
		block = new byte [N2] [144];

		// Parallel region.
		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				// Determine this thread's middle loop range.
				Range range =
					new Range (0, N2-1)
						.subrange (getThreadCount(), getThreadIndex());
				int lb = range.lb();
				int ub = range.ub();

				// Set up cipher object.
				AES256CipherSmp cipher = new AES256CipherSmp (new byte [32]);

				// Outer loop, N1 iterations.
				for (int i = 0; i < N1; ++ i)
					{
					// Middle loop, N2 iterations in parallel.
					for (int j = lb; j <= ub; ++ j)
						{
						byte[] block_j = block[j];

						// Inner loop, N iterations.
						for (int k = 0; k < N3; ++ k)
							{
							cipher.encrypt (block_j, block_j);
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
		System.out.println (time);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.nt=<K> edu.rit.pj.test.Test12 <N1> <N2> <N3>");
		System.err.println ("<K> = Number of parallel threads");
		System.err.println ("<N1> = Number of outer loop iterations");
		System.err.println ("<N2> = Number of middle loop iterations");
		System.err.println ("<N3> = Number of inner loop iterations");
		System.exit (1);
		}

	}
