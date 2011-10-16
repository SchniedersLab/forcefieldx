//******************************************************************************
//
// File:    Test02.java
// Package: edu.rit.mp.test
// Unit:    Class edu.rit.mp.test.Test02
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

package edu.rit.mp.test;

import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;

import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;

import java.text.DecimalFormat;

/**
 * Class Test02 is a main program that sends and receives MP messages using the
 * loopback channel. Each message consists of <I>length</I> integers, starting
 * from 0. The number of messages sent is <I>count</I>. The program prints
 * nothing as messages are sent and received, unless an error occurs. The
 * program prints the time to send all the messages, minus the loop overhead
 * time.
 * <P>
 * Usage: java edu.rit.mp.test.Test02 <I>length</I> <I>count</I>
 * <BR><I>length</I> = Length of each message
 * <BR><I>count</I> = Number of messages
 *
 * @author  Alan Kaminsky
 * @version 18-Jun-2007
 */
public class Test02
	{

	/**
	 * Prevent construction.
	 */
	private Test02()
		{
		}

	/**
	 * Main routine.
	 */
	private void run
		(String[] args)
		throws Throwable
		{
		// Parse command line arguments.
		if (args.length != 2) usage();
		final int length = Integer.parseInt (args[0]);
		final int count = Integer.parseInt (args[1]);

		// Set up channel group.
		final ChannelGroup channelgroup = new ChannelGroup();

		// Get loopback channel.
		final Channel channel = channelgroup.loopbackChannel();

		// Run send and receive sections in parallel threads.
		new ParallelTeam(2).execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				execute
					(new ParallelSection()
						{
						// Send section.
						public void run() throws Exception
							{
							// Measure loop overhead time.
							long t1 = System.currentTimeMillis();
							for (int msgnum = 0; msgnum < count; ++ msgnum)
								{
								}
							long t2 = System.currentTimeMillis();
							System.out.println
								((t2-t1) + " msec loop overhead time");

							// Set up send buffer.
							Buf01 src = new Buf01 (length);

							// Send messages and measure time.
							long t3 = System.currentTimeMillis();
							for (int msgnum = 0; msgnum < count; ++ msgnum)
								{
								channelgroup.send (channel, src);
								}
							long t4 = System.currentTimeMillis();
							long t = t4 - t3 - t2 + t1;
							System.out.println
								(t + " msec message send time");
							double tmsg =
								((double)t) / ((double)count) / 1000.0;
							System.out.println
								(new DecimalFormat("0.00E0").format (tmsg) +
								 " sec per message");
							}
						},

					 new ParallelSection()
						{
						// Receive section.
						public void run() throws Exception
							{
							// Set up receive buffer.
							Buf01 dst = new Buf01 (length);

							// Receive messages.
							for (int msgnum = 0; msgnum < count; ++ msgnum)
								{
								channelgroup.receive (channel, dst);
								}
							}
						});
				}
			});
		}

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		new Test02().run (args);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mp.test.Test02 <tohost> <toport> <length> <count>");
		System.err.println ("<tohost> = Host to which to send messages");
		System.err.println ("<toport> = Port to which to send messages");
		System.err.println ("<length> = Length of each message");
		System.err.println ("<count> = Number of messages");
		System.exit (1);
		}

	}
