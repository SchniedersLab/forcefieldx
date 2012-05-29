//******************************************************************************
//
// File:    NetworkChannel.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.NetworkChannel
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
import java.net.InetSocketAddress;
import java.nio.channels.SocketChannel;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Class NetworkChannel provides a channel for sending and receiving messages
 * over the network in the Message Protocol (MP).
 *
 * @author  Alan Kaminsky
 * @version 23-Apr-2008
 */
class NetworkChannel
	extends Channel
	{

// Hidden data members.

	// Underlying socket channel.
	private SocketChannel mySocketChannel;

	// Far end channel group ID.
	private int myFarChannelGroupId;

	// Queue of outgoing I/O requests.
	private LinkedBlockingQueue<IORequest> myOutgoingQueue;

	// Network channel send and receive threads.
	private NetworkChannelSendThread myNetworkChannelSendThread;
	private NetworkChannelReceiveThread myNetworkChannelReceiveThread;

// Hidden constructors.

	/**
	 * Construct a new network channel.
	 *
	 * @param  theChannelGroup       Enclosing channel group.
	 * @param  theSocketChannel      Underlying socket channel.
	 * @param  theFarChannelGroupId  Far end channel group ID.
	 */
	NetworkChannel
		(ChannelGroup theChannelGroup,
		 SocketChannel theSocketChannel,
		 int theFarChannelGroupId)
		{
		super (theChannelGroup);
		mySocketChannel = theSocketChannel;
		myFarChannelGroupId = theFarChannelGroupId;
		myOutgoingQueue = new LinkedBlockingQueue<IORequest>();
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
		return myFarChannelGroupId;
		}

	/**
	 * Obtain this channel's near end address. This is the host and port of the
	 * near end of this channel's connection.
	 *
	 * @return  Near end address.
	 */
	public InetSocketAddress nearEndAddress()
		{
		return (InetSocketAddress)
			mySocketChannel.socket().getLocalSocketAddress();
		}

	/**
	 * Obtain this channel's far end address. This is the host and port of the
	 * far end of this channel's connection.
	 *
	 * @return  Far end address.
	 */
	public InetSocketAddress farEndAddress()
		{
		return (InetSocketAddress)
			mySocketChannel.socket().getRemoteSocketAddress();
		}

// Hidden operations.

	/**
	 * Start sending and receiving messages via this channel.
	 */
	void start()
		{
		myNetworkChannelSendThread =
			new NetworkChannelSendThread
				(this,
				 mySocketChannel,
				 myOutgoingQueue);
		myNetworkChannelReceiveThread =
			new NetworkChannelReceiveThread
				(this,
				 mySocketChannel);
		}

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
	synchronized void send
		(IORequest theIORequest)
		throws IOException
		{
		// Check whether channel is closed.
		if (myWriteState == WRITE_CLOSED)
			{
			throw new IOException
				("NetworkChannel.send(): Channel closed");
			}

		// Put I/O request in outgoing queue. Send thread will actually send the
		// message.
		myOutgoingQueue.add (theIORequest);
		}

	/**
	 * Perform additional close actions in a subclass.
	 */
	void subclassClose()
		{
		// Close socket channel.
		if (mySocketChannel != null)
			{
			try { mySocketChannel.close(); } catch (IOException exc) {}
			}

		// Interrupt the send and receive threads, which will cause them to
		// terminate.
		if (myNetworkChannelSendThread != null)
			{
			myNetworkChannelSendThread.interrupt();
			}
		if (myNetworkChannelReceiveThread != null)
			{
			myNetworkChannelReceiveThread.interrupt();
			}

		// Enable garbage collection of data members.
		mySocketChannel = null;
		myOutgoingQueue = null;
		myNetworkChannelSendThread = null;
		myNetworkChannelReceiveThread = null;
		}

	/**
	 * Shut down the input side of this network channel.
	 */
	synchronized void shutdownInput()
		{
		myReadState = READ_CLOSED;
		myNetworkChannelReceiveThread = null;
		}

	}
