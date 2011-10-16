//******************************************************************************
//
// File:    Test18.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test18
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

import edu.rit.mp.IntegerBuf;

import edu.rit.pj.Comm;

import edu.rit.util.Range;

/**
 * Class Test18 is a unit test main program for the all-to-all collective
 * communication operation in class {@linkplain edu.rit.pj.Comm}. The program
 * runs on a number of processors in the cluster. Each process sets up arrays
 * of data, fills them in, and prints them; the processes do an all-to-all;
 * and each process prints its final arrays.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.pj.test.Test18 <I>N</I>
 * <BR><I>K</I> = Number of parallel processes
 * <BR><I>N</I> = Array length
 *
 * @author  Alan Kaminsky
 * @version 04-May-2008
 */
public class Test18
	{

// Prevent construction.

	private Test18()
		{
		}

// Global variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static int N;

	// Arrays.
	static int[] senddata;
	static int[] recvdata;

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
		if (args.length != 1) usage();
		N = Integer.parseInt (args[0]);

		// Set up arrays.
		senddata = new int [size*N];
		recvdata = new int [size*N];
		for (int i = 0; i < size*N; ++ i)
			{
			senddata[i] = 100 * rank + i + 1;
			}

		// Print array before all-to-all.
		System.out.print (rank);
		System.out.print (" senddata:");
		for (int i = 0; i < size*N; ++ i)
			{
			System.out.print (' ');
			System.out.print (senddata[i]);
			}
		System.out.println();

		// All-to-all.
		world.allToAll
			(IntegerBuf.sliceBuffers
				(senddata,
				 new Range (0, size*N-1) .subranges (size)),
			 IntegerBuf.sliceBuffers
				(recvdata,
				 new Range (0, size*N-1) .subranges (size)));

		// Print array after all-to-all.
		System.out.print (rank);
		System.out.print (" recvdata:");
		for (int i = 0; i < size*N; ++ i)
			{
			System.out.print (' ');
			System.out.print (recvdata[i]);
			}
		System.out.println();
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.pj.test.Test18 <N>");
		System.err.println ("<K> = Number of parallel processes");
		System.err.println ("<N> = Array length");
		System.exit (1);
		}

	}
