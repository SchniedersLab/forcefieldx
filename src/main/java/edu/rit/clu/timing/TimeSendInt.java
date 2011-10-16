//******************************************************************************
//
// File:    TimeSendInt.java
// Package: edu.rit.clu.timing
// Unit:    Class edu.rit.clu.timing.TimeSendInt
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
 * Class TimeSendInt measures the time to send a message on a workstation
 * cluster computer using Parallel Java (PJ). The program runs on two
 * processors. The program creates an <I>n</I>-integer message, sends the
 * message repeatedly for a certain number of repetitions, and measures the time
 * to send one message. The program repeats this for certain values of the
 * message size <I>n</I>.
 * <P>
 * Usage: java -Dpj.np=2 edu.rit.clu.timing.TimeSendInt <I>reps</I> <I>n1</I>
 * [ <I>n2</I> . . . ]
 * <BR><I>reps</I> = Number of repetitions for each value of <I>n</I>
 * <BR><I>n1</I> = First value of <I>n</I>
 * <BR><I>n2</I> = Second value of <I>n</I> . . .
 *
 * @author  Alan Kaminsky
 * @version 29-Nov-2007
 */
public class TimeSendInt
	{

// Prevent construction.

	private TimeSendInt()
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

		// Process 0.
		if (rank == 0)
			{
			System.out.println
				("java -Dpj.np=2 edu.rit.clu.timing.TimeSendInt " + reps);
			System.out.println (new Date());
			System.out.println ("n\ttime1\ttime2\tSend time (sec)");

			// Test each value of n.
			for (int i = 0; i < numn; ++ i)
				{
				// Create message buffer.
				int n_i = n[i];
				int[] bufarray = new int [n_i];
				IntegerBuf buf = IntegerBuf.buffer (bufarray);

				// Time repetitions without sending messages.
				long time1 = -System.currentTimeMillis();
				for (int j = 0; j < reps; ++ j)
					{
					fill (n_i, bufarray, buf);
					}
				time1 += System.currentTimeMillis();

				// Time repetitions with sending messages.
				long time2 = -System.currentTimeMillis();
				for (int j = 0; j < reps; ++ j)
					{
					fillSendReceive (n_i, bufarray, buf);
					}
				time2 += System.currentTimeMillis();

				// Print results.
				double sendtime =
					((double)(time2 - time1)) / ((double) reps) / 2000.0;
				System.out.print (n_i);
				System.out.print ('\t');
				System.out.print (time1);
				System.out.print ('\t');
				System.out.print (time2);
				System.out.print ('\t');
				System.out.print (FMT3.format (sendtime));
				System.out.println();
				}
			}

		// Process 1.
		else if (rank == 1)
			{
			// Test each value of n.
			for (int i = 0; i < numn; ++ i)
				{
				// Create message buffer.
				int n_i = n[i];
				int[] bufarray = new int [n_i];
				IntegerBuf buf = IntegerBuf.buffer (bufarray);

				// Do repetitions without receiving messages.
				for (int j = 0; j < reps; ++ j)
					{
					fill (n_i, bufarray, buf);
					}

				// Do repetitions with receiving messages.
				for (int j = 0; j < reps; ++ j)
					{
					fillReceiveSend (n_i, bufarray, buf);
					}
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
	 * Fill the buffer, send a message, receive a message.
	 */
	private static void fillSendReceive
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

		// Send message.
		world.send (1, buf);

		// Receive message back.
		world.receive (1, buf);
		}

	/**
	 * Fill the buffer, receive a message, send a message.
	 */
	private static void fillReceiveSend
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

		// Receive message.
		world.receive (0, buf);

		// Send message back.
		world.send (0, buf);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=2 edu.rit.clu.timing.TimeSendInt <reps> <n1> [<n2> ...]");
		System.err.println ("<reps> = Number of repetitions for each value of <n>");
		System.err.println ("<n1> = First value of <n>");
		System.err.println ("<n2> = Second value of <n> . . .");
		System.exit (1);
		}

	}
