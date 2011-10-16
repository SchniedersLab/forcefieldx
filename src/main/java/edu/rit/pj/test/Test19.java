//******************************************************************************
//
// File:    Test19.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test19
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

import edu.rit.pj.Comm;
import edu.rit.pj.PJProperties;

import edu.rit.mp.IntegerBuf;

/**
 * Class Test19 is a unit test main program for class {@linkplain
 * edu.rit.pj.Comm} and class {@linkplain edu.rit.pj.PJProperties}. The program
 * runs on a number of processors in the cluster. Each process prints out its
 * own PJ property settings.
 * <P>
 * Usage: java -Dpj.nn=<I>nn</I> -Dpj.np=<I>np</I> -Dpj.nt=<I>nt</I>
 * edu.rit.pj.test.Test19
 * <BR><I>nn</I> = Number of backend nodes
 * <BR><I>np</I> = Number of processes
 * <BR><I>nt</I> = Number of CPUs per process
 *
 * @author  Alan Kaminsky
 * @version 21-May-2008
 */
public class Test19
	{

// Prevent construction.

	private Test19()
		{
		}

// Global variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

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

		// Wait for previous process to print.
		if (rank > 0) world.receive (rank-1, IntegerBuf.emptyBuffer());

		// Print.
		System.out.println ("Process " + rank);
		System.out.println ("\tworld.size() = " + size);
		System.out.println ("\tworld.rank() = " + rank);
		System.out.println ("\tworld.host() = " + world.host());
		System.out.println ("\tpj.nn = " + PJProperties.getPjNn());
		System.out.println ("\tpj.np = " + PJProperties.getPjNp());
		System.out.println ("\tpj.nt = " + PJProperties.getPjNt());
		System.out.println ("\tpj.schedule = " + PJProperties.getPjSchedule());
		System.out.println ("\tpj.host = " + PJProperties.getPjHost());
		System.out.println ("\tpj.port = " + PJProperties.getPjPort());
		System.out.println ("\tpj.jobtime = " + PJProperties.getPjJobTime());
		System.out.println ("\tpj.jvmflags = " + PJProperties.getPjJvmFlags());
		System.out.println ("\tpj.prng = " + PJProperties.getPjPrng());

		// Tell next process to print.
		if (rank < size-1) world.send (rank+1, IntegerBuf.emptyBuffer());
		}

	}
