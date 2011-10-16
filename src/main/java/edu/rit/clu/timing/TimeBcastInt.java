//******************************************************************************
//
// File:    TimeBcastInt.java
// Package: edu.rit.clu.timing
// Unit:    Class edu.rit.clu.timing.TimeBcastInt
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

package edu.rit.clu.timing;

import edu.rit.mp.IntegerBuf;

import edu.rit.pj.Comm;

import java.io.IOException;

import java.text.DecimalFormat;

import java.util.Date;

/**
 * Class TimeBcastInt measures the time to broadcast a message on a workstation
 * cluster computer using Parallel Java (PJ). The program runs on <I>K</I> &gt;=
 * 2 processors. The program creates an <I>n</I>-integer message, broadcasts the
 * message repeatedly for a certain number of repetitions, and measures the time
 * to broadcast one message. The program repeats this for certain values of the
 * message size <I>n</I>.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.timing.TimeBcastInt <I>reps</I>
 * <I>n1</I> [ <I>n2</I> . . . ]
 * <BR><I>K</I> = Number of processors, <I>K</I> &gt;= 2
 * <BR><I>reps</I> = Number of repetitions for each value of <I>n</I>
 * <BR><I>n1</I> = First value of <I>n</I>
 * <BR><I>n2</I> = Second value of <I>n</I> . . .
 *
 * @author  Alan Kaminsky
 * @version 13-Dec-2007
 */
public class TimeBcastInt
	{

// Prevent construction.

	private TimeBcastInt()
		{
		}

// Hidden constants.

	private static final DecimalFormat FMT3 = new DecimalFormat ("0.00E0");

// Global variables.

	static Comm world;
	static int size;
	static int rank;

	static int reps;
	static int numn;
	static int[] n;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Initialize PJ.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (args.length < 2) usage();
		reps = Integer.parseInt (args[0]);
		numn = args.length - 1;
		n = new int [numn];
		for (int i = 0; i < numn; ++ i)
			{
			n[i] = Integer.parseInt (args[i+1]);
			}

		// Process 0 prints labels.
		if (rank == 0)
			{
			System.out.println
				("java -Dpj.np=" + size +
				 " edu.rit.clu.timing.TimeBcastInt " + reps);
			System.out.println (new Date());
			System.out.println ("n\ttime1\ttime2\tBcast time (sec)");
			}

		// Test each value of n.
		for (int i = 0; i < numn; ++ i)
			{
			// Create message buffer.
			int n_i = n[i];
			int[] bufarray = new int [n_i];
			IntegerBuf buf = IntegerBuf.buffer (bufarray);

			// Time repetitions without broadcasting messages.
			long time1 = -System.currentTimeMillis();
			for (int j = 0; j < reps; ++ j)
				{
				fill (n_i, bufarray, buf);
				}
			time1 += System.currentTimeMillis();

			// Time repetitions with broadcasting messages.
			long time2 = -System.currentTimeMillis();
			for (int j = 0; j < reps; ++ j)
				{
				fillBroadcast (n_i, bufarray, buf);
				}
			time2 += System.currentTimeMillis();

			// Process 0 prints results.
			if (rank == 0)
				{
				double bcasttime =
					((double)(time2 - time1)) / reps / size / 1000.0;
				System.out.print (n_i);
				System.out.print ('\t');
				System.out.print (time1);
				System.out.print ('\t');
				System.out.print (time2);
				System.out.print ('\t');
				System.out.print (FMT3.format (bcasttime));
				System.out.println();
				}
			}
		}

// Hidden operations.

	/**
	 * Fill the buffer.
	 */
	private static void fill
		(int n_i,
		 int[] bufarray,
		 IntegerBuf buf)
		{
		// Fill buffer.
		for (int k = 0; k < n_i; ++ k)
			{
			bufarray[k] = k;
			}
		}

	/**
	 * Fill the buffer, broadcast messages.
	 */
	private static void fillBroadcast
		(int n_i,
		 int[] bufarray,
		 IntegerBuf buf)
		throws IOException
		{
		// Fill buffer.
		for (int k = 0; k < n_i; ++ k)
			{
			bufarray[k] = k;
			}

		// Broadcast from root 0, size-1, size-2, ..., 1.
		for (int i = size; i > 0; -- i)
			{
			world.broadcast (i % size, buf);
			}
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.timing.TimeBcastInt <reps> <n1> [<n2> ...]");
		System.err.println ("<K> = Number of processors, <K> >= 2");
		System.err.println ("<reps> = Number of repetitions for each value of <n>");
		System.err.println ("<n1> = First value of <n>");
		System.err.println ("<n2> = Second value of <n> . . .");
		System.exit (1);
		}

	}
