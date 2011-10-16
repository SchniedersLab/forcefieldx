//******************************************************************************
//
// File:    JobSchedulerTest.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.JobSchedulerTest
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

import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;

import edu.rit.pj.cluster.JobSchedulerProxy;
import edu.rit.pj.cluster.JobSchedulerRef;

import java.net.InetSocketAddress;

import java.util.Scanner;

/**
 * Class JobSchedulerTest is a unit test main program for class {@linkplain
 * edu.rit.pj.cluster.JobScheduler}. Class JobSchedulerTest mimics a job
 * frontend process. It lets you send messages to a Job Scheduler process by
 * typing commands on the console. It displays messages from the Job Scheduler
 * process on the console.
 * <P>
 * Usage: java edu.rit.pj.cluster.test.JobSchedulerTest <I>host</I> <I>port</I>
 * <BR><I>host</I> = Job Scheduler process host name
 * <BR><I>port</I> = Job Scheduler process port number
 *
 * @author  Alan Kaminsky
 * @version 21-May-2008
 */
public class JobSchedulerTest
	{

// Prevent construction.

	private JobSchedulerTest()
		{
		}

// Unit test main program.

	/**
	 * Unit test main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 2) usage();
		String host = args[0];
		int port = Integer.parseInt (args[1]);

		// Set up proxy for Job Scheduler.
		ChannelGroup channelgroup = new ChannelGroup();
		Channel channel =
			channelgroup.connect (new InetSocketAddress (host, port));
		JobSchedulerRef jobscheduler =
			new JobSchedulerProxy (channelgroup, channel);

		// Set up stub for Job Frontend.
		JobFrontendStub jobfrontend = new JobFrontendStub (channelgroup);
		jobfrontend.start();

		// Read Job Scheduler commands from the standard input.
		Scanner scanner = new Scanner (System.in);
		for (;;)
			{
			String line = scanner.nextLine();
			Scanner linescanner = new Scanner (line);
			String command = linescanner.next();
			if (command.equals ("backendFailed"))
				{
				String name = linescanner.next();
				jobscheduler.backendFailed (jobfrontend, name);
				}
			else if (command.equals ("cancelJob"))
				{
				String errmsg = linescanner.next();
				jobscheduler.cancelJob (jobfrontend, errmsg);
				}
			else if (command.equals ("jobFinished"))
				{
				jobscheduler.jobFinished (jobfrontend);
				}
			else if (command.equals ("renewLease"))
				{
				jobscheduler.renewLease (jobfrontend);
				}
			else if (command.equals ("requestJob"))
				{
				String username = linescanner.next();
				int Nn = linescanner.nextInt();
				int Np = linescanner.nextInt();
				int Nt = linescanner.nextInt();
				jobscheduler.requestJob (jobfrontend, username, Nn, Np, Nt);
				}
			else if (command.equals ("close"))
				{
				jobscheduler.close();
				channelgroup.close();
				System.exit (0);
				}
			else
				{
				commands();
				}
			}
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.pj.cluster.test.JobSchedulerTest <host> <port>");
		System.err.println ("<host> = Job Scheduler process host name");
		System.err.println ("<port> = Job Scheduler process port number");
		System.exit (1);
		}

	/**
	 * Print list of commands.
	 */
	private static void commands()
		{
		System.out.println ("requestJob <username> <Nn> <Np> <Nt>");
		System.out.println ("renewLease");
		System.out.println ("jobFinished");
		System.out.println ("cancelJob <errmsg>");
		System.out.println ("backendFailed <name>");
		System.out.println ("close");
		}

	}
