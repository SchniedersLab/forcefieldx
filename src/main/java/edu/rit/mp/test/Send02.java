//******************************************************************************
//
// File:    Send02.java
// Package: edu.rit.mp.test
// Unit:    Class edu.rit.mp.test.Send02
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
import edu.rit.mp.ObjectBuf;

import edu.rit.util.Range;

import java.net.InetSocketAddress;

/**
 * Class Send02 is a main program that sends an MP message to the {@linkplain
 * Receive02} program. The message consists of zero or more strings from the
 * command line, sent as serialized objects.
 * <P>
 * Usage: java edu.rit.mp.test.Send02 <I>tohost</I> <I>toport</I>
 * [ <I>string</I> . . . ]
 * <BR><I>tohost</I> = Host to which to send messages
 * <BR><I>toport</I> = Port to which to send messages
 * <BR><I>string</I> = String to send (zero or more)
 *
 * @author  Alan Kaminsky
 * @version 09-Mar-2006
 */
public class Send02
	{

	/**
	 * Prevent construction.
	 */
	private Send02()
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
		if (args.length < 2) usage();
		String tohost = args[0];
		int toport = Integer.parseInt (args[1]);

		// Set up item source.
		ObjectBuf src =
			ObjectBuf.sliceBuffer (args, new Range (2, args.length-1));

		// Set up channel group.
		ChannelGroup channelgroup = new ChannelGroup();

		// Set up a connection to the far end.
		Channel channel =
			channelgroup.connect (new InetSocketAddress (tohost, toport));

		// Send message.
		channelgroup.send (channel, src);
		}

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		new Send02().run (args);
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mp.test.Send02 <tohost> <toport> [ <string> ... ]");
		System.err.println ("<tohost> = Host to which to send messages");
		System.err.println ("<toport> = Port to which to send messages");
		System.err.println ("<string> = String to send (zero or more)");
		System.exit (1);
		}

	}
