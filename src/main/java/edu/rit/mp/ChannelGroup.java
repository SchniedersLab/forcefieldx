//******************************************************************************
//
// File:    ChannelGroup.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.ChannelGroup
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

package edu.rit.mp;

import edu.rit.util.Logger;
import edu.rit.util.PrintStreamLogger;
import edu.rit.util.Range;
import edu.rit.util.Timer;
import edu.rit.util.TimerTask;
import edu.rit.util.TimerThread;

import java.io.IOException;
import java.io.PrintStream;

import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.Socket;

import java.nio.ByteBuffer;

import java.nio.channels.ClosedChannelException;
import java.nio.channels.ServerSocketChannel;
import java.nio.channels.SocketChannel;

import java.util.LinkedList;
import java.util.List;

/**
 * Class ChannelGroup provides a group of {@linkplain Channel}s for sending and
 * receiving messages in the Message Protocol (MP).
 * <P>
 * <B>Creating Channels</B>
 * <P>
 * A channel group can be used to create channels in two ways: by accepting a
 * connection request from another process, and by requesting a connection to
 * another process.
 * <P>
 * The channel group can be configured to listen to a certain host and port for
 * connection requests. Configure the host and port by by calling the
 * <TT>listen()</TT> method. To start accepting connection requests, call the
 * <TT>startListening()</TT> method.
 * <P>
 * If desired, an application can receive notification of newly created channels
 * by providing a {@linkplain ConnectListener} object to the channel group's
 * <TT>setConnectListener()</TT> method. Specify a connect listener object
 * before calling the <TT>startListening()</TT> method, otherwise the
 * application may not receive some notifications.
 * <P>
 * When a connection request arrives, the channel group sets up a new channel
 * object for communicating over the connection. If a connect listener has been
 * registered, the channel group then passes the channel to the connect
 * listener's <TT>farEndConnected()</TT> method, which does whatever the
 * application needs to record the new channel's presence. The channel does not
 * start sending and receiving messages until after the connect listener's
 * <TT>farEndConnected()</TT> method (if any) has returned.
 * <P>
 * The application can also call the channel group's <TT>connect()</TT> method
 * to request a connection to another host and port. The channel group sets up a
 * new channel object for communicating over the connection. If a connect
 * listener has been registered, the channel group then passes the channel to
 * the connect listener's <TT>nearEndConnected()</TT> method, which does
 * whatever the application needs to record the new channel's presence. The
 * channel does not start sending and receiving messages until after the connect
 * listener's <TT>nearEndConnected()</TT> method (if any) has returned. The
 * <TT>connect()</TT> method also returns the new channel.
 * <P>
 * Once a connection has been set up and a channel object has been created on
 * each side, the applications can use their respective channels to send and
 * receive messages.
 * <P>
 * If a channel group does not need to accept incoming connection requests, the
 * channel group need not listen to any host and port. The channel group can
 * still be used to make outgoing connection requests.
 * <P>
 * <B>Channel Group IDs</B>
 * <P>
 * Each channel group has a channel group ID. The channel group ID is an
 * integer, initially 0, that can be changed by the <TT>setChannelGroupId()</TT>
 * method. The channel group attaches no significance to the channel group ID;
 * it is provided for the use of the application using the channel group.
 * <P>
 * You can query a channel group object to determine its channel group ID. You
 * can also query a channel object to determine the ID of the channel group at
 * the near end of the channel and the ID of the channel group at the far end of
 * the channel.
 * <P>
 * <B>Sending Messages</B>
 * <P>
 * To send a message, the application creates a message buffer (class
 * {@linkplain Buf}) specifying where to get the items to be sent. The
 * application calls the channel group's <TT>send()</TT> method, passing the
 * channel on which to send the message, the message buffer, and the message
 * tag. (If the message tag is not specified, it defaults to 0.) The channel
 * group extracts the items from the message buffer and sends a message over the
 * channel's connection. When the <TT>send()</TT> method returns, the message
 * has been fully sent, but the message may not have been fully received yet.
 * <P>
 * The far end application must receive the message from the channel at the
 * other end of the connection. If no application is receiving the message, the
 * <TT>send()</TT> method may block (because of flow control). This in turn may
 * lead to a deadlock.
 * <P>
 * At most one outgoing message at a time may be in progress on a channel. If a
 * second thread tries to send a message on a channel while a first thread is
 * still sending a message on that channel, the second thread will block until
 * the first thread has finished sending the message.
 * <P>
 * <B>Receiving Messages</B>
 * <P>
 * To receive a message, the application creates a message buffer (class
 * {@linkplain Buf}) specifying where to put the items to be received. The
 * application calls the channel group's <TT>receive()</TT> method, passing the
 * channel from which to receive the message, the desired message tag, and the
 * message buffer. The application can specify "any channel" instead of a
 * specific channel. The application can specify a range of tags or "any tag"
 * instead of a specific tag. Any number of threads can have receive requests
 * pending at the same time.
 * <P>
 * When the channel group receives a message from a channel, the channel group
 * tries to match the message with the pending receive requests. A message
 * matches a receive request if (a) the message's channel is the same as the
 * receive request's channel, or the receive request specified "any channel;"
 * and (b) the message's tag is the same as the receive request's tag, or the
 * message's tag falls within the receive request's range of tags, or the
 * receive request specified "any tag;" and (c) the message's item type is the
 * same as the receive request's item type (as given by the message buffer). The
 * pending receive requests are maintained in FIFO order. If no receive request
 * matches the message, the channel group does not read the message until such
 * time as a matching receive request occurs. If more than one receive request
 * matches the message, the channel group chooses the first matching receive
 * request.
 * <P>
 * Once the channel group has matched the incoming message with the receive
 * request, the channel group reads the items from the message and stores them
 * into the message buffer. If there are fewer items in the message than the
 * length of the message buffer, the extra items at the end of the message
 * buffer are not set to anything. If there are more items in the message than
 * the length of the message buffer, the extra items are read from the message
 * and discarded. Once the message has been read, the <TT>receive()</TT> method
 * returns a {@linkplain Status} object reporting the channel on which the
 * message arrived, the message tag, and the actual number of items in the
 * message (which may or may not be the same as the number of items in the
 * message buffer).
 * <P>
 * If the receive requests do not match properly with the incoming messages, a
 * deadlock may occur.
 * <P>
 * <B>Sending and Receiving Within the Same Process</B>
 * <P>
 * Each channel group has a "loopback" channel that is used to send messages
 * within the same process. To obtain the loopback channel, call the
 * <TT>loopbackChannel()</TT> method. Then one thread can send messages using
 * the loopback channel while a different thread receives messages using the
 * loopback channel. If the same thread both sends and receives using the
 * loopback channel, a deadlock may occur.
 * <P>
 * The loopback channel uses the <TT>copy()</TT> method of class {@linkplain
 * Buf} to transfer data items directly from the source buffer to the
 * destination buffer. The loopback channel does not do any network
 * communication.
 * <P>
 * <B>Non-Blocking Send and Receive Operations</B>
 * <P>
 * The <TT>send()</TT> method described so far does a <B>blocking send</B>
 * operation; the <TT>send()</TT> method does not return until the message has
 * been fully sent. There is also a <B>non-blocking send</B> operation,
 * <TT>sendNoWait()</TT>, which includes an {@linkplain IORequest} argument. The
 * <TT>sendNoWait()</TT> method initiates the send operation and returns
 * immediately. This allows the caller to continue processing while the channel
 * group sends the message in a separate thread. To wait for the message to be
 * fully sent, the caller must call the IORequest object's
 * <TT>waitForFinish()</TT> method.
 * <P>
 * Likewise, the <TT>receive()</TT> method described so far does a <B>blocking
 * receive</B> operation; the <TT>receive()</TT> method does not return until
 * the message has been fully received. There is also a <B>non-blocking
 * receive</B> operation, <TT>receiveNoWait()</TT>, which includes an
 * {@linkplain IORequest} argument. The <TT>receiveNoWait()</TT> method
 * initiates the receive operation and returns immediately. This allows the
 * caller to continue processing while the channel group receives the message in
 * a separate thread. To wait for the message to be fully received, the caller
 * must call the IORequest object's <TT>waitForFinish()</TT> method, which
 * returns a {@linkplain Status} object giving the results of the receive
 * operation.
 *
 * @author  Alan Kaminsky
 * @version 11-Mar-2009
 */
public class ChannelGroup
	{

// Hidden data members.

	// Channel group ID.
	int myChannelGroupId;

	// Server socket channel for accepting incoming connections, or null if not
	// accepting incoming connections.
	ServerSocketChannel myServerSocketChannel;

	// I/O request list for matching incoming messages to receive I/O requests.
	IORequestList myIORequestList;

	// Alternate class loader for use when receiving objects.
	ClassLoader myClassLoader;

	// Loopback channel.
	LoopbackChannel myLoopbackChannel;

	// List of open channels.
	List<Channel> myChannelList;

	// Accepting thread, or null if not accepting.
	AcceptThread myAcceptThread;

	// Registered connect listener, or null if none.
	ConnectListener myConnectListener;

	// For logging error messages.
	Logger myLogger;

	// For timeouts during channel setup.
	TimerThread myTimerThread;

// Hidden helper classes.

	/**
	 * Class AcceptThread provides a thread that accepts incoming connections.
	 *
	 * @author  Alan Kaminsky
	 * @version 11-Mar-2009
	 */
	private class AcceptThread
		extends Thread
		{
		public AcceptThread()
			{
			setDaemon (true);
			start();
			}

		public void run()
			{
			acceptloop : for (;;)
				{
				// Wait for an incoming connection.
				SocketChannel connection = null;
				try
					{
					connection = myServerSocketChannel.accept();
					}
				catch (ClosedChannelException exc)
					{
					myLogger.log
						("ChannelGroup: Channel closed",
						 exc);
					break acceptloop;
					}
				catch (IOException exc)
					{
					myLogger.log
						("ChannelGroup: I/O error while accepting connection",
						 exc);
					break acceptloop;
					}

				// Set up channel over connection.
				if (connection != null)
					{
					try
						{
						farEndConnect (connection);
						}
					catch (IOException exc)
						{
						// Clear thread's interrupted status, otherwise accept()
						// above will throw an exception.
						Thread.interrupted();
						myLogger.log
							("ChannelGroup: I/O error while setting up channel",
							 exc);
						try { connection.close(); } catch (IOException exc2) {}
						}
					}
				}
			}
		}

// Exported constructors.

	/**
	 * Construct a new channel group. The channel group ID is initially 0. The
	 * channel group will not listen for connection requests. To listen for
	 * connection requests at a later time, call the <TT>listen()</TT> method
	 * followed by the <TT>startListening</TT> method.
	 * <P>
	 * The channel group will log error messages on the standard error.
	 */
	public ChannelGroup()
		{
		this (new PrintStreamLogger());
		}

	/**
	 * Construct a new channel group. The channel group ID is initially 0. The
	 * channel group will listen for connection requests on the given host and
	 * port. To start actively listening, call the <TT>startListening()</TT>
	 * method.
	 * <P>
	 * The channel group will log error messages on the standard error.
	 *
	 * @param  theListenAddress  Host and port at which to listen.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theListenAddress</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public ChannelGroup
		(InetSocketAddress theListenAddress)
		throws IOException
		{
		this (theListenAddress, new PrintStreamLogger());
		}

	/**
	 * Construct a new channel group. The channel group ID is initially 0. The
	 * channel group will listen for connection requests using the given server
	 * socket channel. The server socket channel must be bound to a host and
	 * port. To start actively listening, call the <TT>startListening()</TT>
	 * method.
	 * <P>
	 * The channel group will log error messages on the standard error.
	 *
	 * @param  theServerSocketChannel  Server socket channel.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theServerSocketChannel</TT> is
	 *     null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred. Thrown if
	 *     <TT>theServerSocketChannel</TT> is not bound.
	 */
	public ChannelGroup
		(ServerSocketChannel theServerSocketChannel)
		throws IOException
		{
		this (theServerSocketChannel, new PrintStreamLogger());
		}

	/**
	 * Construct a new channel group. The channel group ID is initially 0. The
	 * channel group will not listen for connection requests. To listen for
	 * connection requests at a later time, call the <TT>listen()</TT> method
	 * followed by the <TT>startListening</TT> method.
	 * <P>
	 * The channel group will log error messages using the given logger.
	 *
	 * @param  theLogger  Logger for error messages.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theLogger</TT> is null.
	 */
	public ChannelGroup
		(Logger theLogger)
		{
		if (theLogger == null)
			{
			throw new NullPointerException
				("ChannelGroup(): theLogger is null");
			}
		myIORequestList = new IORequestList();
		myLoopbackChannel = new LoopbackChannel (this);
		myChannelList = new LinkedList<Channel>();
		myChannelList.add (myLoopbackChannel);
		myLogger = theLogger;
		myTimerThread = new TimerThread();
		myTimerThread.setDaemon (true);
		myTimerThread.start();
		}

	/**
	 * Construct a new channel group. The channel group ID is initially 0. The
	 * channel group will listen for connection requests on the given host and
	 * port. To start actively listening, call the <TT>startListening()</TT>
	 * method.
	 * <P>
	 * The channel group will log error messages using the given logger.
	 *
	 * @param  theListenAddress  Host and port at which to listen.
	 * @param  theLogger         Logger for error messages.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theListenAddress</TT> is null.
	 *     Thrown if <TT>theLogger</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public ChannelGroup
		(InetSocketAddress theListenAddress,
		 Logger theLogger)
		throws IOException
		{
		this (theLogger);
		listen (theListenAddress);
		}

	/**
	 * Construct a new channel group. The channel group ID is initially 0. The
	 * channel group will listen for connection requests using the given server
	 * socket channel. The server socket channel must be bound to a host and
	 * port. To start actively listening, call the <TT>startListening()</TT>
	 * method.
	 * <P>
	 * The channel group will log error messages using the given logger.
	 *
	 * @param  theServerSocketChannel  Server socket channel.
	 * @param  theLogger               Logger for error messages.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theServerSocketChannel</TT> is
	 *     null. Thrown if <TT>theLogger</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred. Thrown if
	 *     <TT>theServerSocketChannel</TT> is not bound.
	 */
	public ChannelGroup
		(ServerSocketChannel theServerSocketChannel,
		 Logger theLogger)
		throws IOException
		{
		this (theLogger);
		listen (theServerSocketChannel);
		}

// Exported operations.

	/**
	 * Set this channel group's channel group ID.
	 *
	 * @param  theChannelGroupId  Channel group ID.
	 */
	public void setChannelGroupId
		(int theChannelGroupId)
		{
		myChannelGroupId = theChannelGroupId;
		}

	/**
	 * Obtain this channel group's channel group ID.
	 *
	 * @return  Channel group ID.
	 */
	public int getChannelGroupId()
		{
		return myChannelGroupId;
		}
		
	/**
	 * Obtain this channel group's listen address. This is the near end host and
	 * port to which this channel group is listening for connection requests. If
	 * this channel group is not listening for connection requests, null is
	 * returned.
	 *
	 * @return  Near end address, or null.
	 */
	public synchronized InetSocketAddress listenAddress()
		{
		return
			myServerSocketChannel == null ?
				null :
				(InetSocketAddress)
					myServerSocketChannel.socket().getLocalSocketAddress();
		}

	/**
	 * Listen for connection requests on the given host and port. To start
	 * actively listening, call the <TT>startListening()</TT> method.
	 *
	 * @param  theListenAddress  Host and port at which to listen.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theListenAddress</TT> is null.
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if listening has already started.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public synchronized void listen
		(InetSocketAddress theListenAddress)
		throws IOException
		{
		if (theListenAddress == null)
			{
			throw new NullPointerException
				("ChannelGroup.listen(): theListenAddress is null");
			}
		ServerSocketChannel channel = ServerSocketChannel.open();
		channel.socket().bind (theListenAddress);
		listen (channel);
		}

	/**
	 * Listen for connection requests using the given server socket channel. The
	 * server socket channel must be bound to a host and port. To start actively
	 * listening, call the <TT>startListening()</TT> method.
	 *
	 * @param  theServerSocketChannel  Server socket channel.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theServerSocketChannel</TT> is
	 *     null.
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if listening has already started.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred. Thrown if
	 *     <TT>theServerSocketChannel</TT> is not bound.
	 */
	public synchronized void listen
		(ServerSocketChannel theServerSocketChannel)
		throws IOException
		{
		if (theServerSocketChannel == null)
			{
			throw new NullPointerException
				("ChannelGroup.listen(): theServerSocketChannel is null");
			}
		if (! theServerSocketChannel.socket().isBound())
			{
			throw new IOException
				("ChannelGroup.listen(): theServerSocketChannel is not bound");
			}
		if (myAcceptThread != null)
			{
			throw new IllegalStateException
				("ChannelGroup.listen(): Listening has already started");
			}
		if (myIORequestList == null)
			{
			throw new IOException
				("ChannelGroup.listen(): Channel group closed");
			}

		myServerSocketChannel = theServerSocketChannel;
		}

	/**
	 * Register the given connect listener with this channel group. Thereafter,
	 * this channel group will report each connected channel by calling
	 * <TT>theConnectListener</TT>'s <TT>nearEndConnected()</TT> method (if the
	 * connection request originated in this process) or
	 * <TT>farEndConnected()</TT> method (if the connection request originated
	 * in another process). It is assumed that these methods will not do any
	 * lengthy processing and will not block the calling thread.
	 * <P>
	 * At most one connect listener may be registered. If a connect listener is
	 * already registered, it is replaced with the given connect listener. If
	 * <TT>theConnectListener</TT> is null, any registered connect listener is
	 * discarded, and this channel group will not report connected channels.
	 * <P>
	 * Call the <TT>setConnectListener()</TT> method before calling the
	 * <TT>startListening()</TT> method, otherwise the application may not
	 * receive some connection notifications.
	 *
	 * @param  theConnectListener  Connect listener, or null.
	 */
	public synchronized void setConnectListener
		(ConnectListener theConnectListener)
		{
		myConnectListener = theConnectListener;
		}

	/**
	 * Start actively listening for connection requests.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if a host and port or a server socket
	 *     channel upon which to listen has not been specified. Thrown if
	 *     listening has already started.
	 */
	public synchronized void startListening()
		{
		if (myServerSocketChannel == null)
			{
			throw new IllegalStateException
				("ChannelGroup.startListening(): No server socket channel");
			}
		if (myAcceptThread != null)
			{
			throw new IllegalStateException
				("ChannelGroup.listen(): Listening has already started");
			}

		myAcceptThread = new AcceptThread();
		}

	/**
	 * Create a new channel connected to the given far end host and port. In the
	 * far end computer, there must be a channel group listening to the given
	 * host and port. Once the connection is set up, if a connect listener has
	 * been registered, the channel group calls the connect listener's
	 * <TT>nearEndConnected()</TT> method to report the new channel.
	 *
	 * @param  theFarEndAddress  Host and port of far end channel group.
	 *
	 * @return  New channel.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Channel connect
		(InetSocketAddress theFarEndAddress)
		throws IOException
		{
		synchronized (this)
			{
			if (myIORequestList == null)
				{
				throw new IOException
					("ChannelGroup.connect(): Channel group closed");
				}
			}

		SocketChannel connection = null;
		try
			{
			connection = SocketChannel.open (theFarEndAddress);
			return nearEndConnect (connection);
			}
		catch (IOException exc)
			{
			// Clear thread's interrupted status.
			Thread.interrupted();
			if (connection != null)
				{
				try { connection.close(); } catch (IOException exc2) {}
				}
			throw exc;
			}
		}

	/**
	 * Obtain this channel group's loopback channel. If this channel group is
	 * closed, null is returned.
	 *
	 * @return  Loopback channel, or null.
	 */
	public synchronized Channel loopbackChannel()
		{
		return myLoopbackChannel;
		}

	/**
	 * Send a message to the given channel. The message uses a tag of 0. The
	 * message items come from the given item source buffer.
	 * <P>
	 * The <TT>send()</TT> method does not return until the message has been
	 * fully sent. (The message may not have been fully received yet.)
	 * <P>
	 * The <TT>send()</TT> method assumes that <TT>theChannel</TT> was created
	 * by this channel group. If not, the <TT>send()</TT> method's behavior is
	 * unspecified.
	 *
	 * @param  theChannel  Channel.
	 * @param  theSrc      Item source buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theChannel</TT> is null or
	 *     <TT>theSrc</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void send
		(Channel theChannel,
		 Buf theSrc)
		throws IOException
		{
		IORequest req = new IORequest();
		sendNoWait (theChannel, 0, theSrc, req);
		req.waitForFinish();
		}

	/**
	 * Send a message to the given channel with the given tag. The message items
	 * come from the given item source buffer.
	 * <P>
	 * The <TT>send()</TT> method does not return until the message has been
	 * fully sent. (The message may not have been fully received yet.)
	 * <P>
	 * The <TT>send()</TT> method assumes that <TT>theChannel</TT> was created
	 * by this channel group. If not, the <TT>send()</TT> method's behavior is
	 * unspecified.
	 *
	 * @param  theChannel  Channel.
	 * @param  theTag      Message tag.
	 * @param  theSrc      Item source buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theChannel</TT> is null or
	 *     <TT>theSrc</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void send
		(Channel theChannel,
		 int theTag,
		 Buf theSrc)
		throws IOException
		{
		IORequest req = new IORequest();
		sendNoWait (theChannel, theTag, theSrc, req);
		req.waitForFinish();
		}

	/**
	 * Send (non-blocking) a message to the given channel. The message uses a
	 * tag of 0. The message items come from the given item source buffer.
	 * <TT>theIORequest</TT> is the IORequest object to be associated with the
	 * send operation.
	 * <P>
	 * The <TT>sendNoWait()</TT> method returns immediately. To wait for the
	 * message to be fully sent, call <TT>theIORequest.waitForFinish()</TT>.
	 * <P>
	 * The <TT>sendNoWait()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>sendNoWait()</TT> method's
	 * behavior is unspecified.
	 *
	 * @param  theChannel    Channel.
	 * @param  theSrc        Item source buffer.
	 * @param  theIORequest  IORequest object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theChannel</TT> is null,
	 *     <TT>theSrc</TT> is null, or <TT>theIORequest</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void sendNoWait
		(Channel theChannel,
		 Buf theSrc,
		 IORequest theIORequest)
		throws IOException
		{
		sendNoWait (theChannel, 0, theSrc, theIORequest);
		}

	/**
	 * Send (non-blocking) a message to the given channel with the given tag.
	 * The message items come from the given item source buffer.
	 * <TT>theIORequest</TT> is the IORequest object to be associated with the
	 * send operation.
	 * <P>
	 * The <TT>sendNoWait()</TT> method returns immediately. To wait for the
	 * message to be fully sent, call <TT>theIORequest.waitForFinish()</TT>.
	 * <P>
	 * The <TT>sendNoWait()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>sendNoWait()</TT> method's
	 * behavior is unspecified.
	 *
	 * @param  theChannel    Channel.
	 * @param  theTag        Message tag.
	 * @param  theSrc        Item source buffer.
	 * @param  theIORequest  IORequest object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theChannel</TT> is null,
	 *     <TT>theSrc</TT> is null, or <TT>theIORequest</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void sendNoWait
		(Channel theChannel,
		 int theTag,
		 Buf theSrc,
		 IORequest theIORequest)
		throws IOException
		{
		// Note: This method is not synchronized. Synchronization happens inside
		// theChannel.send().

		// Verify preconditions.
		if (myIORequestList == null)
			{
			throw new IOException
				("ChannelGroup.sendNoWait(): Channel group closed");
			}
		if (theSrc == null)
			{
			throw new NullPointerException
				("ChannelGroup.sendNoWait(): Source buffer is null");
			}

		theIORequest.initialize (theChannel, theTag, theTag, theSrc);
		theChannel.send (theIORequest);
		}

	/**
	 * Receive a message from the given channel. If <TT>theChannel</TT> is null,
	 * a message will be received from any channel in this channel group. The
	 * message must have a tag of 0. The message items are stored in the given
	 * item destination buffer.
	 * <P>
	 * The <TT>receive()</TT> method does not return until the message has been
	 * fully received.
	 * <P>
	 * The <TT>receive()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receive()</TT> method's
	 * behavior is unspecified.
	 *
	 * @param  theChannel  Channel, or null to receive from any channel.
	 * @param  theDst      Item destination buffer.
	 *
	 * @return  Status object giving the outcome of the message reception.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Status receive
		(Channel theChannel,
		 Buf theDst)
		throws IOException
		{
		IORequest req = new IORequest();
		receiveNoWait (theChannel, 0, 0, theDst, req);
		return req.waitForFinish();
		}

	/**
	 * Receive a message from the given channel with the given tag. If
	 * <TT>theChannel</TT> is null, a message will be received from any channel
	 * in this channel group. The message items are stored in the given item
	 * destination buffer.
	 * <P>
	 * The <TT>receive()</TT> method does not return until the message has been
	 * fully received.
	 * <P>
	 * The <TT>receive()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receive()</TT> method's
	 * behavior is unspecified.
	 *
	 * @param  theChannel  Channel, or null to receive from any channel.
	 * @param  theTag      Message tag.
	 * @param  theDst      Item destination buffer.
	 *
	 * @return  Status object giving the outcome of the message reception.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Status receive
		(Channel theChannel,
		 int theTag,
		 Buf theDst)
		throws IOException
		{
		IORequest req = new IORequest();
		receiveNoWait (theChannel, theTag, theTag, theDst, req);
		return req.waitForFinish();
		}

	/**
	 * Receive a message from the given channel with the given range of tags. If
	 * <TT>theChannel</TT> is null, a message will be received from any channel
	 * in this channel group. If <TT>theTagRange</TT> is null, a message will be
	 * received with any tag. The message items are stored in the given item
	 * destination buffer.
	 * <P>
	 * The <TT>receive()</TT> method does not return until the message has been
	 * fully received.
	 * <P>
	 * The <TT>receive()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receive()</TT> method's
	 * behavior is unspecified.
	 *
	 * @param  theChannel   Channel, or null to receive from any channel.
	 * @param  theTagRange  Message tag range, or null to receive any tag.
	 * @param  theDst       Item destination buffer.
	 *
	 * @return  Status object giving the outcome of the message reception.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Status receive
		(Channel theChannel,
		 Range theTagRange,
		 Buf theDst)
		throws IOException
		{
		IORequest req = new IORequest();
		if (theTagRange == null)
			{
			receiveNoWait
				(theChannel, Integer.MIN_VALUE, Integer.MAX_VALUE, theDst, req);
			}
		else
			{
			receiveNoWait
				(theChannel, theTagRange.lb(), theTagRange.ub(), theDst, req);
			}
		return req.waitForFinish();
		}

	/**
	 * Receive (non-blocking) a message from the given channel. If
	 * <TT>theChannel</TT> is null, a message will be received from any channel
	 * in this channel group. The message must have a tag of 0. The message
	 * items are stored in the given item destination buffer.
	 * <TT>theIORequest</TT> is the IORequest object to be associated with the
	 * receive operation.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method returns immediately. To wait for the
	 * message to be fully received, call <TT>theIORequest.waitForFinish()</TT>.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receiveNoWait()</TT>
	 * method's behavior is unspecified.
	 *
	 * @param  theChannel    Channel, or null to receive from any channel.
	 * @param  theDst        Item destination buffer.
	 * @param  theIORequest  IORequest object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null or
	 *     <TT>theIORequest</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void receiveNoWait
		(Channel theChannel,
		 Buf theDst,
		 IORequest theIORequest)
		throws IOException
		{
		receiveNoWait (theChannel, 0, 0, theDst, theIORequest);
		}

	/**
	 * Receive (non-blocking) a message from the given channel with the given
	 * tag. If <TT>theChannel</TT> is null, a message will be received from any
	 * channel in this channel group. The message items are stored in the given
	 * item destination buffer. <TT>theIORequest</TT> is the IORequest object to
	 * be associated with the receive operation.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method returns immediately. To wait for the
	 * message to be fully received, call <TT>theIORequest.waitForFinish()</TT>.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receiveNoWait()</TT>
	 * method's behavior is unspecified.
	 *
	 * @param  theChannel    Channel, or null to receive from any channel.
	 * @param  theTag        Message tag.
	 * @param  theDst        Item destination buffer.
	 * @param  theIORequest  IORequest object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null or
	 *     <TT>theIORequest</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void receiveNoWait
		(Channel theChannel,
		 int theTag,
		 Buf theDst,
		 IORequest theIORequest)
		throws IOException
		{
		receiveNoWait (theChannel, theTag, theTag, theDst, theIORequest);
		}

	/**
	 * Receive (non-blocking) a message from the given channel with the given
	 * range of tags. If <TT>theChannel</TT> is null, a message will be received
	 * from any channel in this channel group. If <TT>theTagRange</TT> is null,
	 * a message will be received with any tag. The message items are stored in
	 * the given item destination buffer. <TT>theIORequest</TT> is the IORequest
	 * object to be associated with the receive operation.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method returns immediately. To wait for the
	 * message to be fully received, call <TT>theIORequest.waitForFinish()</TT>.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receiveNoWait()</TT>
	 * method's behavior is unspecified.
	 *
	 * @param  theChannel    Channel, or null to receive from any channel.
	 * @param  theTagRange   Message tag range, or null to receive any tag.
	 * @param  theDst        Item destination buffer.
	 * @param  theIORequest  IORequest object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null or
	 *     <TT>theIORequest</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void receiveNoWait
		(Channel theChannel,
		 Range theTagRange,
		 Buf theDst,
		 IORequest theIORequest)
		throws IOException
		{
		if (theTagRange == null)
			{
			receiveNoWait
				(theChannel,
				 Integer.MIN_VALUE,
				 Integer.MAX_VALUE,
				 theDst,
				 theIORequest);
			}
		else
			{
			receiveNoWait
				(theChannel,
				 theTagRange.lb(),
				 theTagRange.ub(),
				 theDst,
				 theIORequest);
			}
		}

	/**
	 * Receive (non-blocking) a message from the given channel with the given
	 * tag range. If <TT>theChannel</TT> is null, a message will be received
	 * from any channel in this channel group. The message items are stored in
	 * the given item destination buffer. <TT>theIORequest</TT> is the IORequest
	 * object to be associated with the receive operation.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method returns immediately. To wait for the
	 * message to be fully received, call <TT>theIORequest.waitForFinish()</TT>.
	 * <P>
	 * The <TT>receiveNoWait()</TT> method assumes that <TT>theChannel</TT> was
	 * created by this channel group. If not, the <TT>receiveNoWait()</TT>
	 * method's behavior is unspecified.
	 *
	 * @param  theChannel    Channel, or null to receive from any channel.
	 * @param  theTagLb      Message tag range lower bound.
	 * @param  theTagUb      Message tag range upper bound.
	 * @param  theDst        Item destination buffer.
	 * @param  theIORequest  IORequest object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDst</TT> is null or
	 *     <TT>theIORequest</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void receiveNoWait
		(Channel theChannel,
		 int theTagLb,
		 int theTagUb,
		 Buf theDst,
		 IORequest theIORequest)
		throws IOException
		{
		// Note: This method is not synchronized. Synchronization happens inside
		// myIORequestList.add().

		// Verify preconditions.
		if (myIORequestList == null)
			{
			throw new IOException
				("ChannelGroup.receiveNoWait(): Channel group closed");
			}
		if (theDst == null)
			{
			throw new NullPointerException
				("ChannelGroup.receiveNoWait(): Destination buffer is null");
			}
		if (theChannel != null)
			{
			synchronized (theChannel)
				{
				// Check whether channel is closed.
				if (theChannel.myReadState == Channel.READ_CLOSED)
					{
					throw new IOException
						("ChannelGroup.receiveNoWait(): Channel closed");
					}
				}
			}

		theIORequest.initialize (theChannel, theTagLb, theTagUb, theDst);
		myIORequestList.add (theIORequest);
		}

	/**
	 * Specify an alternate class loader for this channel group. When objects
	 * are received in a message via this channel group, the given class loader
	 * will be used to load the objects' classes. If
	 * <TT>setAlternateClassLoader()</TT> is never called, or if
	 * <TT>theClassLoader</TT> is null, an alternate class loader will not be
	 * used.
	 *
	 * @param  theClassLoader  Alternate class loader, or null.
	 */
	public synchronized void setAlternateClassLoader
		(ClassLoader theClassLoader)
		{
		myClassLoader = theClassLoader;
		}

	/**
	 * Close this channel group. Any pending receive requests will fail with a
	 * {@linkplain ChannelGroupClosedException}.
	 */
	public synchronized void close()
		{
		// Stop listening for connections.
		if (myServerSocketChannel != null)
			{
			try
				{
				myServerSocketChannel.close();
				}
			catch (IOException exc)
				{
				}
			}

		// Close all channels.
		if (myChannelList != null)
			{
			while (! myChannelList.isEmpty())
				{
				myChannelList.get(0).close();
				}
			}

		// Report failure to all pending receive requests.
		if (myIORequestList != null)
			{
			myIORequestList.reportFailure
				(new ChannelGroupClosedException ("Channel group closed"));
			}

		// Enable garbage collection of fields.
		myServerSocketChannel = null;
		myIORequestList = null;
		myClassLoader = null;
		myLoopbackChannel = null;
		myChannelList = null;
		myAcceptThread = null;
		}

	/**
	 * Finalize this channel group.
	 */
	protected void finalize()
		{
		close();
		}

	/**
	 * Dump the state of this channel group on the given print stream. For
	 * debugging.
	 *
	 * @param  out     Print stream.
	 * @param  prefix  String to print at the beginning of each line.
	 */
	public void dump
		(PrintStream out,
		 String prefix)
		{
		out.println (prefix+getClass().getName()+"@"+Integer.toHexString(System.identityHashCode(this)));
		out.println (prefix+"myChannelGroupId = "+myChannelGroupId);
		out.println (prefix+"myServerSocketChannel = "+myServerSocketChannel);
		out.println (prefix+"myIORequestList:");
		myIORequestList.dump (out, prefix+"\t");
		out.println (prefix+"myClassLoader = "+myClassLoader);
		out.println (prefix+"myLoopbackChannel = "+myLoopbackChannel);
		out.println (prefix+"myChannelList:");
		out.println (prefix+"\t"+myChannelList.size()+" entries");
		for (Channel c : myChannelList)
			{
			c.dump (out, prefix+"\t");
			}
		out.println (prefix+"myAcceptThread = "+myAcceptThread);
		out.println (prefix+"myConnectListener = "+myConnectListener);
		out.println (prefix+"myLogger = "+myLogger);
		out.println (prefix+"myTimerThread = "+myTimerThread);
		}

// Hidden operations.

	/**
	 * Create a new network channel using the given socket channel. The
	 * connection request originated from the near end. If this channel group is
	 * closed, null is returned.
	 *
	 * @param  theSocketChannel  Socket channel.
	 *
	 * @return  New channel.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	Channel nearEndConnect
		(SocketChannel theSocketChannel)
		throws IOException
		{
		// Note: This method is not synchronized. Synchronization happens inside
		// createNetworkChannel().

		// Turn on socket's TCP no-delay option.
		Socket socket = theSocketChannel.socket();
		socket.setTcpNoDelay (true);

		// Send channel group ID to far end.
		ByteBuffer buf = ByteBuffer.allocate (4);
		buf.putInt (myChannelGroupId);
		buf.flip();
		if (theSocketChannel.write (buf) != 4)
			{
			throw new IOException
				("ChannelGroup.nearEndConnect(): Cannot send channel group ID");
			}

		// Receive channel group ID from far end with a 30-second timeout.
		buf.clear();
		final Thread thread = Thread.currentThread();
		Timer timer = myTimerThread.createTimer (new TimerTask()
			{
			public void action (Timer theTimer)
				{
				thread.interrupt();
				}
			});
		timer.start (30000L);
		if (theSocketChannel.read (buf) != 4)
			{
			throw new IOException
				("ChannelGroup.nearEndConnect(): Cannot receive channel group ID");
			}
		timer.stop();
		buf.flip();
		int farChannelGroupId = buf.getInt();

		// Set up channel.
		Channel channel =
			createNetworkChannel (theSocketChannel, farChannelGroupId);

		// Inform listener if any.
		if (myConnectListener != null)
			{
			myConnectListener.nearEndConnected (this, channel);
			}

		// Start the channel sending and receiving messages.
		channel.start();

		return channel;
		}

	/**
	 * Create a new network channel using the given socket channel. The
	 * connection request originated from the far end. If this channel group is
	 * closed, null is returned.
	 *
	 * @param  theSocketChannel  Socket channel.
	 *
	 * @return  New channel.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	Channel farEndConnect
		(SocketChannel theSocketChannel)
		throws IOException
		{
		// Note: This method is not synchronized. Synchronization happens inside
		// createNetworkChannel().

		// Turn on socket's TCP no-delay option.
		Socket socket = theSocketChannel.socket();
		socket.setTcpNoDelay (true);

		// Start a 30-second timeout for receiving channel group ID.
		final Thread thread = Thread.currentThread();
		Timer timer = myTimerThread.createTimer (new TimerTask()
			{
			public void action (Timer theTimer)
				{
				thread.interrupt();
				}
			});
		timer.start (30000L);

		try
			{
			// Receive channel group ID from far end.
			ByteBuffer buf = ByteBuffer.allocate (4);
			if (theSocketChannel.read (buf) != 4)
				{
				throw new IOException
					("ChannelGroup.farEndConnect(): Cannot receive channel group ID");
				}
			timer.stop();
			buf.flip();
			int farChannelGroupId = buf.getInt();

			// Send channel group ID to far end.
			buf.clear();
			buf.putInt (myChannelGroupId);
			buf.flip();
			if (theSocketChannel.write (buf) != 4)
				{
				throw new IOException
					("ChannelGroup.farEndConnect(): Cannot send channel group ID");
				}

			// Set up channel.
			Channel channel =
				createNetworkChannel (theSocketChannel, farChannelGroupId);

			// Inform listener if any.
			if (myConnectListener != null)
				{
				myConnectListener.farEndConnected (this, channel);
				}

			// Start the channel sending and receiving messages.
			channel.start();

			return channel;
			}

		// Stop timer when an IOException is thrown.
		catch (IOException exc)
			{
			timer.stop();
			throw exc;
			}
		}

	/**
	 * Create a new network channel using the given socket channel. If this
	 * channel group is closed, null is returned.
	 *
	 * @param  theSocketChannel      Socket channel.
	 * @param  theFarChannelGroupId  Far end channel group ID.
	 *
	 * @return  New channel, or null.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	synchronized Channel createNetworkChannel
		(SocketChannel theSocketChannel,
		 int theFarChannelGroupId)
		throws IOException
		{
		Channel channel = null;
		if (myIORequestList != null)
			{
			channel =
				new NetworkChannel
					(this, theSocketChannel, theFarChannelGroupId);
			myChannelList.add (channel);
			}
		return channel;
		}

	/**
	 * Remove the given channel from this channel group.
	 *
	 * @param  Channel.
	 */
	synchronized void removeChannel
		(Channel channel)
		{
		if (myChannelList != null)
			{
			myChannelList.remove (channel);
			}
		}

	}
