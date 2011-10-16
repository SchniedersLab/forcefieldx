//******************************************************************************
//
// File:    LoopbackChannel.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.LoopbackChannel
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

package edu.rit.mp;

import java.io.IOException;
import java.io.InterruptedIOException;

import java.net.InetSocketAddress;

/**
 * Class LoopbackChannel provides a channel for sending and receiving messages
 * within the same process in the Message Protocol (MP).
 *
 * @author  Alan Kaminsky
 * @version 13-May-2008
 */
class LoopbackChannel
	extends Channel
	{

// Hidden constructors.

	/**
	 * Construct a new loopback channel.
	 *
	 * @param  theChannelGroup  Enclosing channel group.
	 */
	LoopbackChannel
		(ChannelGroup theChannelGroup)
		{
		super (theChannelGroup);
		}

// Exported operations.

	/**
	 * Obtain the channel group ID of this channel's near end channel group.
	 *
	 * @return  Near end channel group ID.
	 */
	public int nearEndChannelGroupId()
		{
		return myChannelGroup.myChannelGroupId;
		}

	/**
	 * Obtain the channel group ID of this channel's far end channel group.
	 *
	 * @return  Far end channel group ID.
	 */
	public int farEndChannelGroupId()
		{
		return myChannelGroup.myChannelGroupId;
		}

	/**
	 * Obtain this channel's near end address. This is the host and port of the
	 * near end of this channel's connection.
	 *
	 * @return  Near end address.
	 */
	public InetSocketAddress nearEndAddress()
		{
		return new InetSocketAddress (0);
		}

	/**
	 * Obtain this channel's far end address. This is the host and port of the
	 * far end of this channel's connection.
	 *
	 * @return  Far end address.
	 */
	public InetSocketAddress farEndAddress()
		{
		return new InetSocketAddress (0);
		}

// Hidden operations.

	/**
	 * Send a message via this channel. The I/O request object must be newly
	 * constructed with the message tag and source buffer fields filled in. This
	 * method is allowed to return immediately and let the message be sent in a
	 * separate thread. The calling thread should use the I/O request object to
	 * wait for the message send to complete.
	 *
	 * @param  theIORequest  I/O request object.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void send
		(IORequest theIORequest)
		throws IOException
		{
		synchronized (this)
			{
			// Check whether channel is closed.
			if (myWriteState == WRITE_CLOSED)
				{
				throw new IOException
					("edu.rit.mp.LoopbackChannel: Channel closed");
				}
			}

		// Wait until there is a matching receive request.
		IORequest recvRequest;
		try
			{
			recvRequest = myIORequestList.waitForMatch (theIORequest);
			}
		catch (InterruptedException exc)
			{
			IOException exc2 =
				new InterruptedIOException
					("edu.rit.mp.LoopbackChannel: Send interrupted");
			exc2.initCause (exc);
			throw exc2;
			}

		// Copy source buffer to destination buffer.
		recvRequest.myBuf.copy (theIORequest.myBuf);

		// Set up status object.
		recvRequest.myStatus =
			new Status
				(this,
				 theIORequest.myTagLb,
				 theIORequest.myBuf.myLength);

		// Report success.
		theIORequest.reportSuccess();
		recvRequest.reportSuccess();
		}

	}
