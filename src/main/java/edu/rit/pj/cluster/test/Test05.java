//******************************************************************************
//
// File:    Test05.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.Test05
//
// This Java source file is copyright (C) 2012 by Alan Kaminsky. All rights
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
import edu.rit.pj.cluster.JobBackend;

/**
 * Class Test05 is a unit test main program for the <TT>reportComment()</TT>
 * method in class {@linkplain edu.rit.pj.cluster.JobBackend}. Each process
 * reports a comment once a second. The program runs until killed externally.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.pj.cluster.test.Test05
 *
 * @author  Alan Kaminsky
 * @version 24-Jan-2012
 */
public class Test05
	{

// Prevent construction.

	private Test05()
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
		int tick = 0;
		for (;;)
			{
			JobBackend.getJobBackend().setComment
				("Tick "+tick+" from process "+rank+" of "+size);
			++ tick;
			Thread.sleep (1000L);
			}
		}

	}
