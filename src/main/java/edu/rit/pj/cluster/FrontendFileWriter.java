//******************************************************************************
//
// File:    FrontendFileWriter.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.FrontendFileWriter
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

package edu.rit.pj.cluster;

import edu.rit.mp.ByteBuf;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import java.util.HashMap;
import java.util.Map;

import java.util.concurrent.LinkedBlockingQueue;

/**
 * Class FrontendFileWriter provides an object that writes sequential files in
 * the job frontend process.
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2007
 */
public class FrontendFileWriter
	{

// Hidden data members.

	private JobFrontend myJobFrontend;

	// Mapping from frontend file descriptor to file handler.
	private Map<Integer,FileHandler> myFileHandlerForFFD =
		new HashMap<Integer,FileHandler>();

	// Next frontend file descriptor.
	private int myNextFFD = 3;

// Hidden helper classes.

	/**
	 * Class FileHandler is an object that performs each file operation in a
	 * separate thread, so as not to block the job frontend's message processing
	 * thread.
	 *
	 * @author  Alan Kaminsky
	 * @version 20-Nov-2006
	 */
	private class FileHandler
		extends Thread
		{
		private LinkedBlockingQueue<Invocation> myQueue =
			new LinkedBlockingQueue<Invocation>();

		private OutputStream myOutputStream;

		private byte[] myByteArray = new byte [0];
		private ByteBuf myByteBuf = ByteBuf.buffer (myByteArray);

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
		 * Construct a new file handler to write the given output stream.
		 *
		 * @param  theOutputStream  Output stream.
		 */
		public FileHandler
			(OutputStream theOutputStream)
			{
			myOutputStream = theOutputStream;
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
		 * Open the given output file for writing or appending.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  bfd            Backend file descriptor.
		 * @param  file           File.
		 * @param  append         True to append, false to overwrite.
		 */
		public void outputFileOpen
			(JobBackendRef theJobBackend,
			 int bfd,
			 File file,
			 boolean append)
			{
			myQueue.offer
				(new OutputFileOpenInvocation
					(theJobBackend, bfd, file, append));
			}

		private class OutputFileOpenInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int bfd;
			private File file;
			private boolean append;

			public OutputFileOpenInvocation
				(JobBackendRef theJobBackend,
				 int bfd,
				 File file,
				 boolean append)
				{
				this.theJobBackend = theJobBackend;
				this.bfd = bfd;
				this.file = file;
				this.append = append;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeOutputFileOpen (theJobBackend, bfd, file, append);
				}
			}

		private boolean invokeOutputFileOpen
			(JobBackendRef theJobBackend,
			 int bfd,
			 File file,
			 boolean append)
			throws IOException
			{
			int ffd = 0;
			IOException result = null;
			boolean more = false;
			try
				{
				myOutputStream = new FileOutputStream (file, append);
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
			theJobBackend.outputFileOpenResult
				(myJobFrontend, bfd, ffd, result);
			return more;
			}

		/**
		 * Write the given bytes to the given output file. <TT>ffd</TT> = 1
		 * refers to the job's standard output stream; <TT>ffd</TT> = 2 refers
		 * to the job's standard error stream; other values refer to a
		 * previously opened file.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  ffd            Frontend file descriptor.
		 * @param  len            Number of bytes to write.
		 */
		public void outputFileWrite
			(JobBackendRef theJobBackend,
			 int ffd,
			 int len)
			{
			myQueue.offer
				(new OutputFileWriteInvocation
					(theJobBackend, ffd, len));
			}

		private class OutputFileWriteInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int ffd;
			private int len;

			public OutputFileWriteInvocation
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
				return invokeOutputFileWrite
					(theJobBackend, ffd, len);
				}
			}

		private boolean invokeOutputFileWrite
			(JobBackendRef theJobBackend,
			 int ffd,
			 int len)
			throws IOException
			{
			IOException result = null;
			boolean more = false;
			try
				{
				if (myByteArray.length < len)
					{
					myByteArray = new byte [len];
					myByteBuf = ByteBuf.buffer (myByteArray);
					}
				((JobBackendProxy) theJobBackend).receive (ffd, myByteBuf);
				myOutputStream.write (myByteArray, 0, len);
				more = true;
				}
			catch (IOException exc)
				{
				result = exc;
				try { myOutputStream.close(); } catch (IOException exc2) {}
				synchronized (myFileHandlerForFFD)
					{
					myFileHandlerForFFD.remove (ffd);
					}
				}
			theJobBackend.outputFileWriteResult (myJobFrontend, ffd, result);
			return more;
			}

		/**
		 * Flush accumulated bytes to the given output file.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  ffd            Frontend file descriptor.
		 */
		public void outputFileFlush
			(JobBackendRef theJobBackend,
			 int ffd)
			{
			myQueue.offer
				(new OutputFileFlushInvocation (theJobBackend, ffd));
			}

		private class OutputFileFlushInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int ffd;

			public OutputFileFlushInvocation
				(JobBackendRef theJobBackend,
				 int ffd)
				{
				this.theJobBackend = theJobBackend;
				this.ffd = ffd;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeOutputFileFlush (theJobBackend, ffd);
				}
			}

		private boolean invokeOutputFileFlush
			(JobBackendRef theJobBackend,
			 int ffd)
			throws IOException
			{
			IOException result = null;
			boolean more = false;
			try
				{
				myOutputStream.flush();
				more = true;
				}
			catch (IOException exc)
				{
				result = exc;
				try { myOutputStream.close(); } catch (IOException exc2) {}
				synchronized (myFileHandlerForFFD)
					{
					myFileHandlerForFFD.remove (ffd);
					}
				}
			theJobBackend.outputFileFlushResult (myJobFrontend, ffd, result);
			return more;
			}

		/**
		 * Close the given output file.
		 *
		 * @param  theJobBackend  Job Backend that is calling this method.
		 * @param  ffd            Frontend file descriptor.
		 */
		public void outputFileClose
			(JobBackendRef theJobBackend,
			 int ffd)
			{
			myQueue.offer
				(new OutputFileCloseInvocation (theJobBackend, ffd));
			}

		private class OutputFileCloseInvocation
			extends Invocation
			{
			private JobBackendRef theJobBackend;
			private int ffd;

			public OutputFileCloseInvocation
				(JobBackendRef theJobBackend,
				 int ffd)
				{
				this.theJobBackend = theJobBackend;
				this.ffd = ffd;
				}

			public boolean invoke()
				throws IOException
				{
				return invokeOutputFileClose (theJobBackend, ffd);
				}
			}

		private boolean invokeOutputFileClose
			(JobBackendRef theJobBackend,
			 int ffd)
			throws IOException
			{
			IOException result = null;
			try
				{
				myOutputStream.close();
				}
			catch (IOException exc)
				{
				result = exc;
				}
			synchronized (myFileHandlerForFFD)
				{
				myFileHandlerForFFD.remove (ffd);
				}
			theJobBackend.outputFileCloseResult (myJobFrontend, ffd, result);
			return false;
			}
		}

// Exported constructors.

	/**
	 * Construct a new frontend file writer.
	 *
	 * @param  theJobFrontend  Job Frontend.
	 */
	public FrontendFileWriter
		(JobFrontend theJobFrontend)
		{
		myJobFrontend = theJobFrontend;

		// Set up frontend file descriptor 1 (stdout) and 2 (stderr).
		myFileHandlerForFFD.put (1, new FileHandler (System.out));
		myFileHandlerForFFD.put (2, new FileHandler (System.err));
		}

// Exported operations.

	/**
	 * Open the given output file for writing or appending.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  bfd            Backend file descriptor.
	 * @param  file           File.
	 * @param  append         True to append, false to overwrite.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileOpen
		(JobBackendRef theJobBackend,
		 int bfd,
		 File file,
		 boolean append)
		throws IOException
		{
		new FileHandler().outputFileOpen (theJobBackend, bfd, file, append);
		}

	/**
	 * Write the given bytes to the given output file. <TT>ffd</TT> = 1 refers
	 * to the job's standard output stream; <TT>ffd</TT> = 2 refers to the job's
	 * standard error stream; other values refer to a previously opened file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 * @param  len            Number of bytes to write.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileWrite
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
			handler.outputFileWrite (theJobBackend, ffd, len);
			}
		else
			{
			theJobBackend.outputFileWriteResult
				(myJobFrontend, ffd,
				 new IOException ("File closed, ffd=" + ffd));
			}
		}

	/**
	 * Flush accumulated bytes to the given output file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileFlush
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
			handler.outputFileFlush (theJobBackend, ffd);
			}
		else
			{
			theJobBackend.outputFileFlushResult
				(myJobFrontend, ffd,
				 new IOException ("File closed, ffd=" + ffd));
			}
		}

	/**
	 * Close the given output file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileClose
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
			handler.outputFileClose (theJobBackend, ffd);
			}
		else
			{
			theJobBackend.outputFileCloseResult
				(myJobFrontend, ffd,
				 new IOException ("File closed, ffd=" + ffd));
			}
		}

	}
