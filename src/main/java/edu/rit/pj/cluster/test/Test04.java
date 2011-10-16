//******************************************************************************
//
// File:    Test04.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.Test04
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

package edu.rit.pj.cluster.test;

import edu.rit.pj.Comm;

import edu.rit.util.Random;

/**
 * Class Test04 is a unit test main program for class {@linkplain
 * edu.rit.pj.Comm}. Each process waits a random number of seconds, then prints
 * a message.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.pj.cluster.test.Test04 [ <I>lb</I> [
 * <I>ub</I> ] ]
 * <BR><I>lb</I> = Lower bound wait time (sec, default 1)
 * <BR><I>ub</I> = Upper bound wait time (sec, default <I>lb</I>+9)
 *
 * @author  Alan Kaminsky
 * @version 23-Apr-2008
 */
public class Test04
	{

// Prevent construction.

	private Test04()
		{
		}

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		Comm.init (args);
		Comm world = Comm.world();
		int size = world.size();
		int rank = world.rank();
		int lb = args.length >= 1 ? Integer.parseInt (args[0]) : 1;
		int ub = args.length >= 2 ? Integer.parseInt (args[1]) : lb+9;
		Random prng = Random.getInstance (System.currentTimeMillis());
		int sec = prng.nextInt (ub-lb+1) + lb;
		Thread.sleep (sec*1000L);
		System.out.println (rank + ": " + sec + " sec");
		}

	}
