//******************************************************************************
//
// File:    Test02.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.Test02
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

package edu.rit.pj.cluster.test;

import edu.rit.pj.Comm;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

/**
 * Class Test02 is a unit test main program for class {@linkplain
 * edu.rit.pj.Comm}. The program runs on the number of nodes and with the number
 * of processes specified. In each process, the program runs the number of
 * threads specified. For further information about specifying the number of
 * nodes, processes, and threads, see class {@linkplain edu.rit.pj.PJProperties
 * PJProperties}. Each thread in each process prints out a "Hello, world"
 * message and echoes the command line arguments.
 * <P>
 * Usage: java -Dpj.nn=<I>Kn</I> -Dpj.np=<I>Kp</I> -Dpj.nt=<I>Kt</I>
 * edu.rit.pj.cluster.test.Test02 [ <I>args</I> ]
 * <BR><I>Kn</I> = Number of nodes
 * <BR><I>Kp</I> = Number of processes
 * <BR><I>Kt</I> = Number of threads per process
 *
 * @author  Alan Kaminsky
 * @version 11-Mar-2009
 */
public class Test02
	{

// Prevent construction.

	private Test02()
		{
		}

// Global variables.

	static Comm world;
	static int size;
	static int rank;
	static int threadSize;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(final String[] args)
		throws Exception
		{
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();
		threadSize = ParallelTeam.getDefaultThreadCount();
		new ParallelTeam(threadSize).execute (new ParallelRegion()
			{
			public void run()
				{
				StringBuilder buf = new StringBuilder();
				buf.append ("Hello, world from thread ");
				buf.append (getThreadIndex());
				buf.append (" of ");
				buf.append (threadSize);
				buf.append (", process ");
				buf.append (rank);
				buf.append (" of ");
				buf.append (size);
				buf.append ("!");
				for (String arg : args)
					{
					buf.append (" ");
					buf.append (arg);
					}
				System.out.println (buf);
				}
			});
		}

	}
