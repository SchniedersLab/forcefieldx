//******************************************************************************
//
// File:    Receive03.java
// Package: edu.rit.mp.test
// Unit:    Class edu.rit.mp.test.Receive03
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

import edu.rit.mp.ChannelGroup;
import edu.rit.mp.Status;

import edu.rit.mp.FloatBuf;

import java.net.InetSocketAddress;

/**
 * Class Receive03 is a main program that receives MP messages from the
 * {@linkplain Send03} program. Each message consists of <I>length</I> floats,
 * starting from 0. The number of messages received is <I>count</I>. The program
 * prints nothing as messages are received, unless an error occurs.
 * <P>
 * Usage: java edu.rit.mp.test.Receive03 <I>tohost</I> <I>toport</I>
 * <I>length</I> <I>count</I>
 * <BR><I>tohost</I> = Host to which to send messages
 * <BR><I>toport</I> = Port to which to send messages
 * <BR><I>length</I> = Length of each message
 * <BR><I>count</I> = Number of messages
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
public class Receive03
	{

	/**
	 * Prevent construction.
	 */
	private Receive03()
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

		// Set up item destination.
		float[] items = new float [length];
		FloatBuf dst = FloatBuf.buffer (items);

		// Set up channel group.
		ChannelGroup channelgroup =
			new ChannelGroup (new InetSocketAddress (tohost, toport));
		channelgroup.startListening();

		// Receive messages.
		Status status = null;
		for (int msgnum = 0; msgnum < count; ++ msgnum)
			{
			status = channelgroup.receive (null, dst);
			}

		// Check final messge.
		if (status.length != length)
			{
			System.out.print ("Invalid length = ");
			System.out.print (status.length);
			System.out.println();
			}
		for (int i = 0; i < length; ++ i)
			{
			if (items[i] != i)
				{
				System.out.print ("Invalid items[");
				System.out.print (i);
				System.out.print ("] = ");
				System.out.print (items[i]);
				System.out.println();
				}
			}
		}

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		new Receive03().run (args);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mp.test.Receive03 <tohost> <toport> <length> <count>");
		System.err.println ("<tohost> = Host to which to send messages");
		System.err.println ("<toport> = Port to which to send messages");
		System.err.println ("<length> = Length of each message");
		System.err.println ("<count> = Number of messages");
		System.exit (1);
		}

	}
