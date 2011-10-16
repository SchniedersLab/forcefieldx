//******************************************************************************
//
// File:    Send01.java
// Package: edu.rit.mp.test
// Unit:    Class edu.rit.mp.test.Send01
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

import java.net.InetSocketAddress;

import java.text.DecimalFormat;

/**
 * Class Send01 is a main program that sends MP messages to the {@linkplain
 * Receive01} program. Each message consists of <I>length</I> integers, starting
 * from 0. The number of messages sent is <I>count</I>. The program prints
 * nothing as messages are sent, unless an error occurs. The program prints the
 * time to send all the messages, minus the loop overhead time.
 * <P>
 * Usage: java edu.rit.mp.test.Send01 <I>tohost</I> <I>toport</I> <I>length</I>
 * <I>count</I>
 * <BR><I>tohost</I> = Host to which to send messages
 * <BR><I>toport</I> = Port to which to send messages
 * <BR><I>length</I> = Length of each message
 * <BR><I>count</I> = Number of messages
 *
 * @author  Alan Kaminsky
 * @version 27-Jan-2006
 */
public class Send01
	{

	/**
	 * Prevent construction.
	 */
	private Send01()
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
		if (args.length != 4) usage();
		String tohost = args[0];
		int toport = Integer.parseInt (args[1]);
		int length = Integer.parseInt (args[2]);
		int count = Integer.parseInt (args[3]);

		// Set up item source.
		Buf01 src = new Buf01 (length);

		// Set up channel group.
		ChannelGroup channelgroup = new ChannelGroup();

		// Set up a connection to the far end.
		Channel channel =
			channelgroup.connect (new InetSocketAddress (tohost, toport));

		// Measure loop overhead time.
		long t1 = System.currentTimeMillis();
		for (int msgnum = 0; msgnum < count; ++ msgnum)
			{
			}
		long t2 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec loop overhead time");

		// Send messages and measure time.
		long t3 = System.currentTimeMillis();
		for (int msgnum = 0; msgnum < count; ++ msgnum)
			{
			channelgroup.send (channel, src);
			}
		long t4 = System.currentTimeMillis();
		long t = t4 - t3 - t2 + t1;
		System.out.println (t + " msec message send time");
		double tmsg = ((double)t) / ((double)count) / 1000.0;
		System.out.println
			(new DecimalFormat ("0.00E0") .format (tmsg) + " sec per message");
		}

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		new Send01().run (args);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mp.test.Send01 <tohost> <toport> <length> <count>");
		System.err.println ("<tohost> = Host to which to send messages");
		System.err.println ("<toport> = Port to which to send messages");
		System.err.println ("<length> = Length of each message");
		System.err.println ("<count> = Number of messages");
		System.exit (1);
		}

	}
