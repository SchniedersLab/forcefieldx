//******************************************************************************
//
// File:    Test01.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.Test01
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

package edu.rit.pj.cluster.test;

import edu.rit.pj.PJProperties;
import edu.rit.pj.cluster.JobBackend;

/**
 * Class Test01 is a unit test main program for class {@linkplain
 * edu.rit.pj.cluster.JobBackend}. It prints out certain information about the
 * static job backend object.
 *
 * @author  Alan Kaminsky
 * @version 14-Mar-2007
 */
public class Test01
	{

// Prevent construction.

	private Test01()
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
		JobBackend backend = JobBackend.getJobBackend();
		if (backend == null)
			{
			System.out.println ("No job backend object");
			}
		else
			{
			int K = backend.getK();
			int rank = backend.getRank();
			Thread.sleep (2000L * rank);
			System.out.println ("********************************");
			System.out.println ("K = " + K);
			System.out.println ("rank = " + rank);
			System.out.println ("pj.nt = " + PJProperties.getPjNt());
			System.out.println ("pj.schedule = " + PJProperties.getPjSchedule());
			System.out.println ("pj.np = " + PJProperties.getPjNp());
			System.out.println ("pj.host = " + PJProperties.getPjHost());
			System.out.println ("pj.port = " + PJProperties.getPjPort());
			System.out.println ("pj.jvmflags = " + PJProperties.getPjJvmFlags());
			Thread.sleep (2000L * (K - rank));
			System.out.println ("Process " + rank + " finished");
			}
		}

	}
