//******************************************************************************
//
// File:    Test08.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test08
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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerStrideForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.Schedule;

/**
 * Class Test08 is a unit test main program for classes {@linkplain
 * edu.rit.pj.ParallelTeam ParallelTeam}, {@linkplain edu.rit.pj.ParallelRegion
 * ParallelRegion}, and {@linkplain edu.rit.pj.IntegerForLoop IntegerForLoop}.
 * A parallel for loop iterates over a given range of indexes. A portion of the
 * loop body executes in parallel. Another portion of the loop body executes
 * sequentially inside an <TT>ordered()</TT> method call.
 * <P>
 * Usage: java [ -Dpj.nt=<I>K</I> ] [ -Dpj.schedule=<I>schedule</I> ]
 * edu.rit.pj.test.Test08 <I>lb</I> <I>ub</I> <I>stride</I>
 * <BR><I>K</I> = Number of parallel threads
 * <BR><I>schedule</I> = Parallel for loop schedule
 * <BR><I>lb</I> = Loop index lower bound, inclusive
 * <BR><I>ub</I> = Loop index upper bound, inclusive
 * <BR><I>stride</I> = Loop stride
 *
 * @author  Alan Kaminsky
 * @version 31-May-2007
 */
public class Test08
	{

// Prevent construction.

	private Test08()
		{
		}

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		if (args.length != 3) usage();
		final int lb = Integer.parseInt (args[0]);
		final int ub = Integer.parseInt (args[1]);
		final int stride = Integer.parseInt (args[2]);

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
							throws Exception
							{
							for (int i = first; i <= last; ++ i)
								{
								final int ii = i;
								System.out.println
									("i=" + ii + ", thread=" +
									 getThreadIndex() + ", unordered");
								ordered (new ParallelSection()
									{
									public void run()
										{
										System.out.println
											("\t\t\t\ti=" + ii + ", thread=" +
											 getThreadIndex() + ", ordered");
										}
									});
								}
							}
						});
					}
				else
					{
					execute (lb, ub, stride, new IntegerStrideForLoop()
						{
						public void run (int first, int last, int stride)
							throws Exception
							{
							for (int i = first; i <= last; i += stride)
								{
								final int ii = i;
								System.out.println
									("i=" + ii + ", thread=" +
									 getThreadIndex() + ", unordered");
								ordered (new ParallelSection()
									{
									public void run()
										{
										System.out.println
											("\t\t\t\ti=" + ii + ", thread=" +
											 getThreadIndex() + ", ordered");
										}
									});
								}
							}
						});
					}
				System.out.println ("End thread " + getThreadIndex());
				}
			});
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java [-Dpj.nt=<K>] [-Dpj.schedule=<schedule>] edu.rit.pj.test.Test08 <lb> <ub> <stride>");
		System.err.println ("<K> = Number of parallel threads");
		System.err.println ("<schedule> = Parallel for loop schedule");
		System.err.println ("<lb> = Loop index lower bound, inclusive");
		System.err.println ("<ub> = Loop index upper bound, inclusive");
		System.err.println ("<stride> = Loop stride");
		System.exit (1);
		}

	}
