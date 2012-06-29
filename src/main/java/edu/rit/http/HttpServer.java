//******************************************************************************
//
// File:    HttpServer.java
// Package: edu.rit.http
// Unit:    Class edu.rit.http.HttpServer
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

package edu.rit.http;

import edu.rit.util.Logger;
import edu.rit.util.PrintStreamLogger;

import java.io.IOException;

import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;

// For unit test main program
// import java.io.PrintWriter;
// import java.util.Map;

/**
 * Class HttpServer provides a lightweight HTTP/1.0 server. The HTTP server is
 * designed to be embedded inside another application.
 * <P>
 * When constructed, the HTTP server starts a thread listening for connections
 * to a given host and port. When a web browser sets up a connection, the HTTP
 * server calls the <TT>process()</TT> method to process the request. This is an
 * abstract method that must be overridden in a subclass. The <TT>process()</TT>
 * method's arguments are an {@linkplain HttpRequest} object from which the
 * method reads the HTTP request message and an {@linkplain HttpResponse} object
 * to which the method writes the HTTP response message.
 * <P>
 * The HTTP server assumes that the <TT>process()</TT> method will not perform
 * any lengthy processing and will not block the calling thread. If the HTTP
 * request requires lengthy processing or requires blocking, the
 * <TT>process()</TT> method should spawn a separate thread to process the
 * request. A thread pool may prove useful; see package java.util.concurrent.
 * <P>
 * When a client opens a socket connection to the HTTP server, the server places
 * a two-second timeout on reading the socket. If the client does not send an
 * HTTP request message before the timeout, the {@linkplain HttpRequest}'s
 * <TT>isValid()</TT> method returns false.
 *
 * @author  Alan Kaminsky
 * @version 29-Jul-2010
 */
public abstract class HttpServer
	{

// Hidden data members.

	private ServerSocket myServerSocket;
	private AcceptorThread myAcceptorThread;
	private Logger myLogger;

// Hidden helper classes.

	private class AcceptorThread
		extends Thread
		{
		public void run()
			{
			try
				{
				for (;;)
					{
					Socket socket = myServerSocket.accept();
					socket.setSoTimeout (2000);
					HttpRequest request = new HttpRequest (socket);
					HttpResponse response = new HttpResponse (socket);
					try
						{
						process (request, response);
						}
					catch (Throwable exc)
						{
						// Any exception while processing a request: Ignore.
						myLogger.log
							("Exception while processing HTTP request",
							 exc);
						}
					finally
						{
						try
							{
							if (! socket.isClosed()) socket.close();
							}
						catch (Throwable exc)
							{
							// Any exception while closing socket: Ignore.
							myLogger.log
								("Exception while closing HTTP socket",
								 exc);
							}
						}
					socket = null;
					request = null;
					response = null;
					}
				}
			catch (Throwable exc)
				{
				// Any exception while accepting a connection: Terminate thread.
				if (! myServerSocket.isClosed())
					{
					myLogger.log
						("Exception while accepting HTTP connection",
						 exc);
					}
				}
			finally
				{
				myLogger.log ("HTTP server terminating");
				}
			}
		}

// Exported constructors.

	/**
	 * Construct a new HTTP server. The HTTP server will print error messages
	 * on the standard error.
	 *
	 * @param  address  Host and port to which the HTTP server will listen for
	 *                  connections.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public HttpServer
		(InetSocketAddress address)
		throws IOException
		{
		this (address, null);
		}

	/**
	 * Construct a new HTTP server. The HTTP server will print error messages
	 * using the given logger.
	 *
	 * @param  address  Host and port to which the HTTP server will listen for
	 *                  connections.
	 * @param  logger   Error message logger. If null, the HTTP server will
	 *                  print error messages on the standard error.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public HttpServer
		(InetSocketAddress address,
		 Logger logger)
		throws IOException
		{
		myLogger = logger == null ? new PrintStreamLogger() : logger;
		myServerSocket = new ServerSocket();
		myServerSocket.bind (address);
		myAcceptorThread = new AcceptorThread();
		myAcceptorThread.setDaemon (true);
		myAcceptorThread.start();
		}

// Exported operations.

	/**
	 * Obtain the host and port to which this HTTP server is listening for
	 * connections.
	 *
	 * @return  Host and port.
	 */
	public InetSocketAddress getAddress()
		{
		return (InetSocketAddress) myServerSocket.getLocalSocketAddress();
		}

	/**
	 * Close this HTTP server.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void close()
		throws IOException
		{
		myServerSocket.close();
		}

// Hidden operations.

	/**
	 * Process the given HTTP request. The <TT>process()</TT> method must be
	 * overridden in a subclass to read the HTTP request from
	 * <TT>theRequest</TT> and write the HTTP response to <TT>theResponse</TT>.
	 *
	 * @param  theRequest   HTTP request.
	 * @param  theResponse  HTTP response.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	protected abstract void process
		(HttpRequest theRequest,
		 HttpResponse theResponse)
		throws IOException;

// Unit test main program.

//	/**
//	 * Unit test main program. The program listens for connections to
//	 * localhost:8080. The program reads each HTTP request from a web browser
//	 * and merely echoes the request data back to the browser.
//	 * <P>
//	 * Usage: java edu.rit.http.HttpServer
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		HttpServer server =
//			new HttpServer (new InetSocketAddress ("localhost", 8080))
//				{
//				protected void process
//					(HttpRequest request,
//					 HttpResponse response)
//					throws IOException
//					{
//					if (request.isValid())
//						{
//						PrintWriter out = response.getPrintWriter();
//						out.println ("<HTML>");
//						out.println ("<HEAD>");
//						out.println ("</HEAD>");
//						out.println ("<BODY>");
//						out.println ("<UL>");
//						out.println ("<LI>");
//						out.print   ("Method = <TT>\"");
//						out.print   (request.getMethod());
//						out.println ("\"</TT>");
//						out.println ("<LI>");
//						out.print   ("URI = <TT>\"");
//						out.print   (request.getUri());
//						out.println ("\"</TT>");
//						out.println ("<LI>");
//						out.print   ("Version = <TT>\"");
//						out.print   (request.getHttpVersion());
//						out.println ("\"</TT>");
//						for (Map.Entry<String,String> entry :
//									request.getHeaders())
//							{
//							out.println ("<LI>");
//							out.print   ("Header name = <TT>\"");
//							out.print   (entry.getKey());
//							out.print   ("\"</TT>, value = <TT>\"");
//							out.print   (entry.getValue());
//							out.println ("\"</TT>");
//							}
//						out.println ("</UL>");
//						out.println ("</BODY>");
//						out.println ("</HTML>");
//						}
//					else
//						{
//						response.setStatusCode
//							(HttpResponse.Status.STATUS_400_BAD_REQUEST);
//						PrintWriter out = response.getPrintWriter();
//						out.println ("<HTML>");
//						out.println ("<HEAD>");
//						out.println ("</HEAD>");
//						out.println ("<BODY>");
//						out.println ("<P>");
//						out.println ("400 Bad Request");
//						out.println ("<P>");
//						out.println ("You idiot.");
//						out.println ("</BODY>");
//						out.println ("</HTML>");
//						}
//					response.close();
//					}
//				};
//		
//		Thread.currentThread().join();
//		}

	}
