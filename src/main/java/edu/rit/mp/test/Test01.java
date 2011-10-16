//******************************************************************************
//
// File:    Test01.java
// Package: edu.rit.mp.test
// Unit:    Class edu.rit.mp.test.Test01
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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
import edu.rit.mp.ConnectListener;

import edu.rit.util.PrintStreamLogger;

import java.net.InetSocketAddress;

/**
 * Class Test01 is a unit test main program that creates and accepts channel
 * connections in the Message Protocol (MP). First run this command in one
 * process to accept connection requests:
 * <P>
 * java edu.rit.mp.test.Test01 <I>host1</I> <I>port1</I>
 * <P>
 * Then run this command in another process to request a connection to the first
 * process:
 * <P>
 * java edu.rit.mp.test.Test01 <I>host2</I> <I>port2</I> <I>host1</I>
 * <I>port1</I>
 *
 * @author  Alan Kaminsky
 * @version 11-Mar-2009
 */
public class Test01
	{

// Prevent construction.

	private Test01()
		{
		}

// Unit test main program.

	/**
	 * Unit test main program.
	 * <P>
	 * Usage: java edu.rit.mp.test.Test01 <I>nearhost</I> <I>nearport</I> [
	 * <I>farhost</I> <I>farport</I> ]
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		if (args.length != 2 && args.length != 4) usage();
		String nearhost = args[0];
		int nearport = Integer.parseInt (args[1]);
		ChannelGroup channelgroup =
			new ChannelGroup
				(new InetSocketAddress (nearhost, nearport));
//		ChannelGroup channelgroup =
//			new ChannelGroup
//				(new InetSocketAddress (nearhost, nearport),
//				 new PrintStreamLogger());
		channelgroup.setConnectListener (new ConnectListener()
			{
			public void nearEndConnected
				(ChannelGroup theChannelGroup,
				 Channel theChannel)
				{
				System.out.println ("nearEndConnected (" + theChannel + ")");
				}
			public void farEndConnected
				(ChannelGroup theChannelGroup,
				 Channel theChannel)
				{
				System.out.println ("farEndConnected (" + theChannel + ")");
				}
			});
		channelgroup.startListening();
		if (args.length == 4)
			{
			String farhost = args[2];
			int farport = Integer.parseInt (args[3]);
			channelgroup.connect (new InetSocketAddress (farhost, farport));
			channelgroup.close();
			}
		else
			{
			Thread.currentThread().join();
			}
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mp.test.Test01 <nearhost> <nearport> [<farhost> <farport>]");
		System.exit (1);
		}

	}
