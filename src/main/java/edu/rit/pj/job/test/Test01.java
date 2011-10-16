//******************************************************************************
//
// File:    Test01.java
// Package: edu.rit.pj.job
// Unit:    Class edu.rit.pj.job.Test01
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

package edu.rit.pj.job.test;

import edu.rit.io.Stdio;

import edu.rit.pj.job.Job;
import edu.rit.pj.job.JobGenerator;

import java.io.File;

/**
 * Class Test01 is a {@linkplain edu.rit.pj.job.JobGenerator JobGenerator} for
 * unit testing. Each job, when run, sleeps for one second, then prints its job
 * number and a message on the standard output. The standard output is
 * redirected to a file named <TT>"out<I>i</I>"</TT>, where <TT><I>i</I></TT> is
 * the job number.
 *
 * @author  Alan Kaminsky
 * @version 08-Oct-2010
 */
public class Test01
	extends JobGenerator
	{

// Hidden data members.

	private int N;

// Exported constructors.

	/**
	 * Construct a new test job generator.
	 *
	 * @param  N  Number of jobs.
	 */
	public Test01
		(int N)
		{
		this.N = N;
		}

// Hidden operations.

	/**
	 * Returns the number of jobs in the job group, <I>N</I>.
	 *
	 * @return  Number of jobs.
	 */
	protected int jobCount()
		{
		return N;
		}

	/**
	 * Create the job with the given job number. This method must create and
	 * return an instance of class {@linkplain Job} whose job number is
	 * <TT>theJobNumber</TT>.
	 *
	 * @param  theJobNumber  Job number (0 .. <I>N</I>&minus;1).
	 */
	protected Job createJob
		(int theJobNumber)
		{
		Job job = new Job
			(theJobNumber, "Test01", "edu.rit.pj.job.test.Test01");
		job.addArgument ("Job");
		job.addArgument (""+theJobNumber);
		job.addArgument ("Hello");
		job.addArgument ("world");
		job.stdoutToFile (new File ("out"+theJobNumber));
		job.stderrToStdout();
		return job;
		}

// Job main program.

	/**
	 * Test01 job main program. It simply prints its arguments on the standard
	 * output.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		Thread.sleep (1000L);
		for (String arg : args) Stdio.out().printf ("%s ", arg);
		Stdio.out().println();
		}

	}
