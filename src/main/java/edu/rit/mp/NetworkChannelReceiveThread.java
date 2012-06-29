//******************************************************************************
//
// File:    NetworkChannelReceiveThread.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.NetworkChannelReceiveThread
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

import java.io.EOFException;
import java.io.IOException;
import java.io.InterruptedIOException;

import java.nio.ByteBuffer;

import java.nio.channels.ClosedChannelException;
import java.nio.channels.SocketChannel;

/**
 * Class NetworkChannelReceiveThread provides a thread for receiving incoming
 * messages for a {@linkplain NetworkChannel}.
 *
 * @author  Alan Kaminsky
 * @version 23-Apr-2008
 */
class NetworkChannelReceiveThread
	extends Thread
	{

// Hidden data members.

	// Enclosing network channel and channel group.
	private NetworkChannel myNetworkChannel;
	private ChannelGroup myChannelGroup;

	// Underlying socket channel.
	private SocketChannel mySocketChannel;

	// Queue of incoming I/O requests.
	private IORequestList myIORequestList;

	// Byte buffer.
	private ByteBuffer myByteBuffer;

// Hidden constructors.

	/**
	 * Construct a new network channel receive thread.
	 *
	 * @param  theNetworkChannel  Enclosing network channel.
	 * @param  theSocketChannel   Underlying socket channel.
	 */
	NetworkChannelReceiveThread
		(NetworkChannel theNetworkChannel,
		 SocketChannel theSocketChannel)
		{
		myNetworkChannel = theNetworkChannel;
		myChannelGroup = theNetworkChannel.myChannelGroup;
		mySocketChannel = theSocketChannel;
		myIORequestList = theNetworkChannel.myIORequestList;
		myByteBuffer = ByteBuffer.allocateDirect (Constants.BUFFER_SIZE);
		setDaemon (true);
		start();
		}

// Exported operations.

	/**
	 * Run this network channel receive thread.
	 */
	public void run()
		{
		IORequest iorequest = null;
		Buf buf = null;
		myByteBuffer.position (0);
		myByteBuffer.limit (0);

		// Processing loop.
		receiveloop: for (;;)
			{
			try
				{
				// Read the next 13-byte message header.
				while (myByteBuffer.remaining() < 13)
					{
					myByteBuffer.compact();
					if (mySocketChannel.read (myByteBuffer) == -1)
						{
						break receiveloop;
						}
					myByteBuffer.flip();
					}

				// Extract message header fields.
				int magic = myByteBuffer.getInt();
				int messagetag = myByteBuffer.getInt();
				byte messagetype = (byte) myByteBuffer.get();
				int messagelength = myByteBuffer.getInt();

				// If the magic number is incorrect, bad error. Close the
				// channel and terminate this thread.
				if (magic != Constants.MAGIC_NUMBER)
					{
					myChannelGroup.myLogger.log
						("edu.rit.mp.NetworkChannelReceiveThread: Invalid magic number received");
					myNetworkChannel.close();
					break receiveloop;
					}
				
				// Wait for a matching I/O request to show up in the incoming
				// queue.
				iorequest =
					myIORequestList.waitForMatch
						(myNetworkChannel, messagetag, messagetype);

				// Message preprocessing.
				buf = iorequest.myBuf;
				buf.preReceive (messagelength);
				int buflength = buf.myMessageLength;
				int n;
				int i = 0;
				int num = Math.min (messagelength, buflength);

				// Repeatedly transfer items from socket channel to byte buffer,
				// then from byte buffer to destination buffer.
				n = buf.receiveItems (i, num, myByteBuffer);
				i += n;
				num -= n;
				while (num > 0)
					{
					myByteBuffer.compact();
					if (mySocketChannel.read (myByteBuffer) == -1)
						{
						throw new EOFException
							("Unexpected end-of-stream while receiving message");
						}
					myByteBuffer.flip();
					n = buf.receiveItems (i, num, myByteBuffer);
					i += n;
					num -= n;
					}

				// If there are more items in the message than in the
				// destination buffer, suck out the extra message items.
				num = messagelength - buflength;
				if (num > 0)
					{
					num -= buf.skipItems (num, myByteBuffer);
					while (num > 0)
						{
						myByteBuffer.compact();
						if (mySocketChannel.read (myByteBuffer) == -1)
							{
							throw new EOFException
								("Unexpected end-of-stream while receiving message");
							}
						myByteBuffer.flip();
						num -= buf.skipItems (num, myByteBuffer);
						}
					}

				// Message postprocessing.
				Status status =
					new Status (myNetworkChannel, messagetag, messagelength);
				buf.postReceive (status, myChannelGroup.myClassLoader);
				iorequest.myStatus = status;

				// Report success to receiving thread.
				iorequest.reportSuccess();

				iorequest = null;
				buf = null;
				}

			catch (IOException exc)
				{
				// Report failure to receiving thread. Terminate this thread.
				if (iorequest != null)
					{
					iorequest.reportFailure (exc);
					}
				break receiveloop;
				}

			catch (InterruptedException exc)
				{
				// Report failure to sending thread. Terminate this thread.
				if (iorequest != null)
					{
					InterruptedIOException exc2 = new InterruptedIOException();
					exc2.initCause (exc);
					iorequest.reportFailure (exc2);
					}
				break receiveloop;
				}

			catch (RuntimeException exc)
				{
				// Report failure to sending thread. Terminate this thread.
				if (iorequest != null)
					{
					iorequest.reportFailure (exc);
					}
				break receiveloop;
				}

			catch (Error exc)
				{
				// Report failure to sending thread. Terminate this thread.
				if (iorequest != null)
					{
					iorequest.reportFailure (exc);
					}
				break receiveloop;
				}
			}

		// This thread is terminating. Enable garbage collection of data
		// members.
		myNetworkChannel.shutdownInput();
		myNetworkChannel = null;
		myChannelGroup = null;
		mySocketChannel = null;
		myIORequestList = null;
		myByteBuffer = null;
		}

	}
