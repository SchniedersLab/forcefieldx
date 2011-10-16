//******************************************************************************
//
// File:    Test06.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test06
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
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;

/**
 * Class Test06 is a unit test main program for classes {@linkplain
 * edu.rit.pj.ParallelTeam ParallelTeam}, {@linkplain edu.rit.pj.ParallelRegion
 * ParallelRegion}, and {@linkplain edu.rit.pj.ParallelSection ParallelSection}.
 * A number of parallel sections are created and executed. Each parallel
 * section's code is slightly different.
 * <P>
 * Usage: java -Dpj.nt=<I>K</I> edu.rit.pj.test.Test06 <I>N</I>
 * <BR><I>K</I> = Number of parallel threads
 * <BR><I>N</I> = Number of parallel sections
 *
 * @author  Alan Kaminsky
 * @version 04-Jun-2007
 */
public class Test06
	{

// Prevent construction.

	private Test06()
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
		if (args.length != 1) usage();
		final int N = Integer.parseInt (args[0]);

		final ParallelSection[] sections = new ParallelSection [N];
		for (int i = 0; i < N; ++ i)
			{
			final int ii = i + 1;
			sections[i] = new ParallelSection()
				{
				public void run()
					{
					System.out.println
						("Parallel section " + ii +
						 ", thread " + getThreadIndex());
					}
				};
			}

		new ParallelTeam().execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				System.out.println ("Begin thread " + getThreadIndex());
				execute (sections);
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
		System.err.println ("Usage: java -Dpj.nt=<K> edu.rit.pj.test.Test06 <N>");
		System.err.println ("<K> = Number of parallel threads");
		System.err.println ("<N> = Number of parallel sections");
		System.exit (1);
		}

	}
