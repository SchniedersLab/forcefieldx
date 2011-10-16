//******************************************************************************
//
// File:    Test24.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test24
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
import edu.rit.pj.WorkerIteration;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;

/**
 * Class Test24 is a unit test main program for classes {@linkplain
 * edu.rit.pj.WorkerTeam WorkerTeam}, {@linkplain edu.rit.pj.WorkerRegion
 * WorkerRegion}, and {@linkplain edu.rit.pj.WorkerIteration WorkerIteration}. A
 * worker iteration iterates over the strings given as command line arguments.
 * Each iteration prints the item (string) and the worker index processing that
 * item.
 * <P>
 * Usage: java [ -Dpj.np=<I>Kp</I> ] edu.rit.pj.test.Test24 <I>strings</I>
 * <BR><I>Kp</I> = Number of parallel processes
 * <BR><I>strings</I> = Strings to be iterated over
 *
 * @author  Alan Kaminsky
 * @version 07-Oct-2010
 */
public class Test24
	{

// Prevent construction.

	private Test24()
		{
		}

// Main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(final String[] args)
		throws Throwable
		{
		Comm.init (args);

		new WorkerTeam().execute (new WorkerRegion()
			{
			public void run() throws Exception
				{
				final int w = getThreadIndex();
				System.out.printf
					("Begin thread %d, rank %d%n",
					 w, Comm.world().rank());
				execute (args, new WorkerIteration<String>()
					{
					public void run (String s) throws Exception
						{
						Thread.sleep (1000L);
						System.out.printf ("s = \"%s\", thread = %d%n", s, w);
						}
					});
				System.out.printf ("End thread %d%n", w);
				}
			});
		}

	}
