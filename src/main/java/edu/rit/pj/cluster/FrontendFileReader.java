//******************************************************************************
//
// File:    FrontendFileReader.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.FrontendFileReader
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

package edu.rit.pj.cluster;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Class FrontendFileReader provides an object that reads sequential files in
 * the job frontend process.
 *
 * @author  Alan Kaminsky
 * @version 05-Nov-2006
 */
public class FrontendFileReader
	{

// Hidden data members.

	private JobFrontend myJobFrontend;

	// Mapping from frontend file descriptor to file handler.
	private Map<Integer,FileHandler> myFileHandlerForFFD =
		new HashMap<Integer,FileHandler>();

	// Next frontend file descriptor.
	private int myNextFFD = 2;

// Hidden helper classes.

	/**
	 * Class FileHandler is an object that performs each file operation in a
	 * separate thread, so as not to block the job frontend's message processing
	 * thread.
	 *
	 * @author  Alan Kaminsky
	 * @version 05-Nov-2006
	 */
	private class FileHandler
		extends Thread
		{
		private LinkedBlockingQueue<Invocation> myQueue =
			new LinkedBlockingQueue<Invocation>();

		private InputStream myInputStream;
		private byte[] myBuffer = new byte [0];

		private abstract class Invocation
			{
			public abstract boolean invoke()
				throws IOException;
			}

		/**
		 * Construct a new file handler.
		 */
		public FileHandler()
			{
			setDaemon (true);
			start();
			}

		/**
		 * Construct a new file handler to read the given input stream.
		 *
		 * @param  theInputStream  Input stream.
		 */
		public FileHandler
			(InputStream theInputStream)
			{
			myInputStream = theInputStream;
			setDaemon (true);
			start();
			}

		/**
		 * Run this file handler.
		 */
		public void run()
			{
			try
				{
				while (myQueue.take().invoke());
				}
			catch (Throwable exc)
				{
				myJobFrontend.terminateCancelJobOther (exc);
				}
			}

		/**
		 * Open the given input file for reading.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  bfd            Backend file descriptor.
		 * @param  file           File.
		 */
		public void inputFileOpen
			(JobBackendRef theJobBackend,
			 int bfd,
			 File file)
			{
			myQueue.offer
				(new InputFileOpenInvocation (theJobBackend, bfd, file));
			}

		private class InputFileOpenInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int bfd;
			private File file;

			public InputFileOpenInvocation
				(JobBackendRef theJobBackend,
				 int bfd,
				 File file)
				{
				this.theJobBackend = theJobBackend;
				this.bfd = bfd;
				this.file = file;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeInputFileOpen (theJobBackend, bfd, file);
				}
			}

		private boolean invokeInputFileOpen
			(JobBackendRef theJobBackend,
			 int bfd,
			 File file)
			throws IOException
			{
			int ffd = 0;
			IOException result = null;
			boolean more = false;
			try
				{
				myInputStream = new FileInputStream (file);
				synchronized (myFileHandlerForFFD)
					{
					ffd = myNextFFD ++;
					myFileHandlerForFFD.put (ffd, this);
					}
				more = true;
				}
			catch (IOException exc)
				{
				result = exc;
				}
			theJobBackend.inputFileOpenResult (myJobFrontend, bfd, ffd, result);
			return more;
			}

		/**
		 * Read bytes from the given input file. <TT>ffd</TT> = 1 refers to the
		 * job's standard input stream; other values refer to a previously
		 * opened file.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  ffd            Frontend file descriptor.
		 * @param  len            Number of bytes to read.
		 */
		public void inputFileRead
			(JobBackendRef theJobBackend,
			 int ffd,
			 int len)
			{
			myQueue.offer
				(new InputFileReadInvocation (theJobBackend, ffd, len));
			}

		private class InputFileReadInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int ffd;
			private int len;

			public InputFileReadInvocation
				(JobBackendRef theJobBackend,
				 int ffd,
				 int len)
				{
				this.theJobBackend = theJobBackend;
				this.ffd = ffd;
				this.len = len;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeInputFileRead (theJobBackend, ffd, len);
				}
			}

		private boolean invokeInputFileRead
			(JobBackendRef theJobBackend,
			 int ffd,
			 int len)
			throws IOException
			{
			int resultlen = 0;
			IOException resultexc = null;
			boolean more = false;
			try
				{
				if (myBuffer.length < len) myBuffer = new byte [len];
				resultlen = myInputStream.read (myBuffer, 0, len);
				more = true;
				}
			catch (IOException exc)
				{
				resultexc = exc;
				try { myInputStream.close(); } catch (IOException exc2) {}
				synchronized (myFileHandlerForFFD)
					{
					myFileHandlerForFFD.remove (ffd);
					}
				}
			theJobBackend.inputFileReadResult
				(myJobFrontend, ffd, myBuffer, resultlen, resultexc);
			return more;
			}

		/**
		 * Skip bytes from the given input file.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  ffd            Frontend file descriptor.
		 * @param  len            Number of bytes to skip.
		 */
		public void inputFileSkip
			(JobBackendRef theJobBackend,
			 int ffd,
			 long len)
			{
			myQueue.offer
				(new InputFileSkipInvocation (theJobBackend, ffd, len));
			}

		private class InputFileSkipInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int ffd;
			private long len;

			public InputFileSkipInvocation
				(JobBackendRef theJobBackend,
				 int ffd,
				 long len)
				{
				this.theJobBackend = theJobBackend;
				this.ffd = ffd;
				this.len = len;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeInputFileSkip (theJobBackend, ffd, len);
				}
			}

		private boolean invokeInputFileSkip
			(JobBackendRef theJobBackend,
			 int ffd,
			 long len)
			throws IOException
			{
			long resultlen = 0L;
			IOException resultexc = null;
			boolean more = false;
			try
				{
				resultlen = myInputStream.skip (len);
				more = true;
				}
			catch (IOException exc)
				{
				resultexc = exc;
				try { myInputStream.close(); } catch (IOException exc2) {}
				synchronized (myFileHandlerForFFD)
					{
					myFileHandlerForFFD.remove (ffd);
					}
				}
			theJobBackend.inputFileSkipResult
				(myJobFrontend, ffd, resultlen, resultexc);
			return more;
			}

		/**
		 * Close the given input file.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  ffd            Frontend file descriptor.
		 */
		public void inputFileClose
			(JobBackendRef theJobBackend,
			 int ffd)
			{
			myQueue.offer
				(new InputFileCloseInvocation (theJobBackend, ffd));
			}

		private class InputFileCloseInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int ffd;

			public InputFileCloseInvocation
				(JobBackendRef theJobBackend,
				 int ffd)
				{
				this.theJobBackend = theJobBackend;
				this.ffd = ffd;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeInputFileClose (theJobBackend, ffd);
				}
			}

		private boolean invokeInputFileClose
			(JobBackendRef theJobBackend,
			 int ffd)
			throws IOException
			{
			IOException result = null;
			try
				{
				myInputStream.close();
				}
			catch (IOException exc)
				{
				result = exc;
				}
			synchronized (myFileHandlerForFFD)
				{
				myFileHandlerForFFD.remove (ffd);
				}
			theJobBackend.inputFileCloseResult (myJobFrontend, ffd, result);
			return false;
			}
		}

// Exported constructors.

	/**
	 * Construct a new frontend file reader.
	 *
	 * @param  theJobFrontend  Job Frontend.
	 */
	public FrontendFileReader
		(JobFrontend theJobFrontend)
		{
		myJobFrontend = theJobFrontend;

		// Set up frontend file descriptor 1 (stdin).
		myFileHandlerForFFD.put (1, new FileHandler (System.in));
		}

// Exported operations.

	/**
	 * Open the given input file for reading.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  bfd            Backend file descriptor.
	 * @param  file           File.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void inputFileOpen
		(JobBackendRef theJobBackend,
		 int bfd,
		 File file)
		throws IOException
		{
		new FileHandler().inputFileOpen (theJobBackend, bfd, file);
		}

	/**
	 * Read bytes from the given input file. <TT>ffd</TT> = 1 refers to the
	 * job's standard input stream; other values refer to a previously opened
	 * file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 * @param  len            Number of bytes to read.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void inputFileRead
		(JobBackendRef theJobBackend,
		 int ffd,
		 int len)
		throws IOException
		{
		FileHandler handler = null;
		synchronized (myFileHandlerForFFD)
			{
			handler = myFileHandlerForFFD.get (ffd);
			}
		if (handler != null)
			{
			handler.inputFileRead (theJobBackend, ffd, len);
			}
		else
			{
			theJobBackend.inputFileReadResult
				(myJobFrontend, ffd, null, -1,
				 new IOException ("File closed, ffd=" + ffd));
			}
		}

	/**
	 * Skip bytes from the given input file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 * @param  len            Number of bytes to skip.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void inputFileSkip
		(JobBackendRef theJobBackend,
		 int ffd,
		 long len)
		throws IOException
		{
		FileHandler handler = null;
		synchronized (myFileHandlerForFFD)
			{
			handler = myFileHandlerForFFD.get (ffd);
			}
		if (handler != null)
			{
			handler.inputFileSkip (theJobBackend, ffd, len);
			}
		else
			{
			theJobBackend.inputFileSkipResult
				(myJobFrontend, ffd, 0L,
				 new IOException ("File closed, ffd=" + ffd));
			}
		}

	/**
	 * Close the given input file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void inputFileClose
		(JobBackendRef theJobBackend,
		 int ffd)
		throws IOException
		{
		FileHandler handler = null;
		synchronized (myFileHandlerForFFD)
			{
			handler = myFileHandlerForFFD.get (ffd);
			}
		if (handler != null)
			{
			handler.inputFileClose (theJobBackend, ffd);
			}
		else
			{
			theJobBackend.inputFileCloseResult
				(myJobFrontend, ffd,
				 new IOException ("File closed, ffd=" + ffd));
			}
		}

	}
