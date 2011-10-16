//******************************************************************************
//
// File:    Test02.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.Test02
//
// This Java source file is copyright (C) 2005 by Alan Kaminsky. All rights
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
import edu.rit.pj.ParallelTeam;

/**
 * Class Test02 is a unit test main program for classes {@linkplain
 * edu.rit.pj.ParallelTeam ParallelTeam} and {@linkplain
 * edu.rit.pj.ParallelRegion ParallelRegion}. Each thread in the team prints its
 * thread index on the standard output.
 * <P>
 * Usage: java -Dpj.nt=<I>K</I> edu.rit.pj.test.Test02
 * <BR><I>K</I> = Number of parallel threads
 *
 * @author  Alan Kaminsky
 * @version 29-Aug-2005
 */
public class Test02
	{

// Prevent construction.

	private Test02()
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
		ParallelTeam team = new ParallelTeam();
		team.execute (new ParallelRegion()
			{
			public void start()
				{
				System.out.println ("Start first region");
				}

			public void run() throws Exception
				{
				int i = getThreadIndex();
				Thread.sleep (1000L * i);
				System.out.println
					("Run thread " + i + " first region");
				}

			public void finish()
				{
				System.out.println ("Finish first region");
				}
			});
		team.execute (new ParallelRegion()
			{
			public void start()
				{
				System.out.println ("Start second region");
				}

			public void run() throws Exception
				{
				int i = getThreadIndex();
				Thread.sleep (1000L * i);
				System.out.println
					("Run thread " + i + " second region");
				}

			public void finish()
				{
				System.out.println ("Finish second region");
				}
			});
		}

	}
