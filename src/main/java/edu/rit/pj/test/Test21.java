//******************************************************************************
//
// File:    Test21.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test21
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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
import edu.rit.pj.WorkerIntegerForLoop;
import edu.rit.pj.WorkerIntegerStrideForLoop;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;

/**
 * Class Test21 is a unit test main program for classes {@linkplain
 * edu.rit.pj.WorkerTeam WorkerTeam}, {@linkplain edu.rit.pj.WorkerRegion
 * WorkerRegion}, and {@linkplain edu.rit.pj.WorkerIntegerForLoop
 * WorkerIntegerForLoop}. A worker for loop iterates over a given range of
 * indexes. Each iteration prints the loop index and the worker index executing
 * that loop index.
 * <P>
 * Usage: java [ -Dpj.np=<I>K</I> ] [ -Dpj.schedule=<I>schedule</I> ]
 * edu.rit.pj.test.Test21 <I>lb</I> <I>ub</I> <I>stride</I>
 * <BR><I>K</I> = Number of parallel processes
 * <BR><I>schedule</I> = Worker for loop schedule
 * <BR><I>lb</I> = Loop index lower bound, inclusive
 * <BR><I>ub</I> = Loop index upper bound, inclusive
 * <BR><I>stride</I> = Loop stride
 *
 * @author  Alan Kaminsky
 * @version 19-Jan-2010
 */
public class Test21
	{

// Prevent construction.

	private Test21()
		{
		}

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		Comm.init (args);

		if (args.length != 3) usage();
		final int lb = Integer.parseInt (args[0]);
		final int ub = Integer.parseInt (args[1]);
		final int stride = Integer.parseInt (args[2]);

		new WorkerTeam().execute (new WorkerRegion()
			{
			public void run() throws Exception
				{
				final int w = getThreadIndex();
				System.out.printf ("Begin thread %d%n", w);
				if (stride == 1)
					{
					execute (lb, ub, new WorkerIntegerForLoop()
						{
						public void run (int first, int last)
							{
							for (int i = first; i <= last; ++ i)
								{
								System.out.printf
									("i = %d, thread = %d%n", i, w);
								}
							}
						});
					}
				else
					{
					execute (lb, ub, stride, new WorkerIntegerStrideForLoop()
						{
						public void run (int first, int last, int stride)
							{
							for (int i = first; i <= last; i += stride)
								{
								System.out.printf
									("i = %d, thread = %d%n", i, w);
								}
							}
						});
					}
				System.out.printf ("End thread %d%n", w);
				}
			});
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java [-Dpj.np=<K>] [-Dpj.schedule=<schedule>] edu.rit.pj.test.Test21 <lb> <ub> <stride>");
		System.err.println ("<K> = Number of parallel processes");
		System.err.println ("<schedule> = Worker for loop schedule");
		System.err.println ("<lb> = Loop index lower bound, inclusive");
		System.err.println ("<ub> = Loop index upper bound, inclusive");
		System.err.println ("<stride> = Loop stride");
		System.exit (1);
		}

	}
