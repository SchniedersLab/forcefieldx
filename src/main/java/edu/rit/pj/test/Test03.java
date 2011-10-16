//******************************************************************************
//
// File:    Test03.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test03
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

import edu.rit.pj.BarrierAction;
import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerStrideForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import java.util.concurrent.BrokenBarrierException;

/**
 * Class Test03 is a unit test main program for classes {@linkplain
 * edu.rit.pj.ParallelTeam ParallelTeam}, {@linkplain edu.rit.pj.ParallelRegion
 * ParallelRegion}, and {@linkplain edu.rit.pj.IntegerForLoop IntegerForLoop}.
 * A parallel for loop iterates over a given range of indexes. Each iteration
 * prints the loop index and the thread index executing that loop index.
 * <P>
 * Usage: java [ -Dpj.nt=<I>K</I> ] [ -Dpj.schedule=<I>schedule</I> ]
 * edu.rit.pj.test.Test03 <I>lb</I> <I>ub</I> <I>stride</I> <I>wait</I>
 * <BR><I>K</I> = Number of parallel threads
 * <BR><I>schedule</I> = Parallel for loop schedule
 * <BR><I>lb</I> = Loop index lower bound, inclusive
 * <BR><I>ub</I> = Loop index upper bound, inclusive
 * <BR><I>stride</I> = Loop stride
 * <BR><I>wait</I> = Wait option, WAIT or NO_WAIT
 *
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public class Test03
	{

// Prevent construction.

	private Test03()
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
		if (args.length != 4) usage();
		final int lb = Integer.parseInt (args[0]);
		final int ub = Integer.parseInt (args[1]);
		final int stride = Integer.parseInt (args[2]);
		final BarrierAction wait =
			args[3].equals ("WAIT") ? BarrierAction.WAIT :
			args[3].equals ("NO_WAIT") ? BarrierAction.NO_WAIT :
			null;

		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				System.out.println ("Begin thread " + getThreadIndex());
				if (stride == 1)
					{
					execute (lb, ub, new IntegerForLoop()
						{
						public void run (int first, int last)
							{
							for (int i = first; i <= last; ++ i)
								{
								System.out.println
									("i = " + i +
									 ", thread = " + getThreadIndex());
								}
							}
						},
						wait);
					}
				else
					{
					execute (lb, ub, stride, new IntegerStrideForLoop()
						{
						public void run (int first, int last, int stride)
							{
							for (int i = first; i <= last; i += stride)
								{
								System.out.println
									("i = " + i +
									 ", thread = " + getThreadIndex());
								}
							}
						},
						wait);
					}
				System.out.println
					("End thread " + getThreadIndex());
				}
			});
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java [-Dpj.nt=<K>] [-Dpj.schedule=<schedule>] edu.rit.pj.test.Test03 <lb> <ub> <stride> <wait>");
		System.err.println ("<K> = Number of parallel threads");
		System.err.println ("<schedule> = Parallel for loop schedule");
		System.err.println ("<lb> = Loop index lower bound, inclusive");
		System.err.println ("<ub> = Loop index upper bound, inclusive");
		System.err.println ("<stride> = Loop stride");
		System.err.println ("<wait> = Wait option, WAIT or NO_WAIT");
		System.exit (1);
		}

	}
