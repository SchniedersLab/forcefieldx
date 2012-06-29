//******************************************************************************
//
// File:    BackendFileWriter.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.BackendFileWriter
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

import edu.rit.io.LineBufferedOutputStream;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import java.util.HashMap;
import java.util.Map;

/**
 * Class BackendFileWriter provides an object that writes sequential files in
 * the job backend process.
 * <P>
 * <I>Note:</I> Class BackendFileWriter is not multiple thread safe; it assumes
 * it is being called by a synchronized method in the job backend.
 *
 * @author  Alan Kaminsky
 * @version 05-Nov-2006
 */
public class BackendFileWriter
	{

// Exported data members.

	/**
	 * Print stream for printing on the job frontend's standard output.
	 */
	public final PrintStream out;

	/**
	 * Print stream for printing on the job frontend's standard error.
	 */
	public final PrintStream err;

// Hidden data members.

	private JobFrontendRef myJobFrontend;
	private JobBackendRef myJobBackend;

	// Mapping from backend file descriptor to backend file output stream.
	private Map<Integer,BackendFileOutputStream> myOutputStreamForBFD =
		new HashMap<Integer,BackendFileOutputStream>();

	// Mapping from frontend file descriptor to backend file output stream.
	private Map<Integer,BackendFileOutputStream> myOutputStreamForFFD =
		new HashMap<Integer,BackendFileOutputStream>();

	// Next backend file descriptor.
	private int myNextBFD = 1;

// Exported constructors.

	/**
	 * Construct a new backend file writer.
	 *
	 * @param  theJobFrontend  Job Frontend.
	 * @param  theJobBackend   Job Backend.
	 */
	public BackendFileWriter
		(JobFrontendRef theJobFrontend,
		 JobBackendRef theJobBackend)
		{
		myJobFrontend = theJobFrontend;
		myJobBackend = theJobBackend;

		// Set up output streams for stdout and stderr.
		BackendFileOutputStream outstream =
			new BackendFileOutputStream
				(myJobFrontend, myJobBackend, 1);
		BackendFileOutputStream errstream =
			new BackendFileOutputStream
				(myJobFrontend, myJobBackend, 2);
		myOutputStreamForFFD.put (1, outstream);
		myOutputStreamForFFD.put (2, errstream);

		// Set up print streams for stdout and stderr.
		out =
			new PrintStream
				(new LineBufferedOutputStream (outstream),
				 true);
		err =
			new PrintStream
				(new LineBufferedOutputStream (errstream),
				 true);
		}

// Exported operations.

	/**
	 * Open a backend file output stream on the given file.
	 *
	 * @param  file    File.
	 * @param  append  True to append, false to overwrite.
	 *
	 * @return  Backend file output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public BackendFileOutputStream open
		(File file,
		 boolean append)
		throws IOException
		{
		BackendFileOutputStream stream = null;
		int bfd = 0;
		int ffd = 0;

		synchronized (this)
			{
			stream = new BackendFileOutputStream (myJobFrontend, myJobBackend);
			bfd = myNextBFD ++;
			myOutputStreamForBFD.put (bfd, stream);
			}

		ffd = stream.open (bfd, file, append);

		synchronized (this)
			{
			myOutputStreamForFFD.put (ffd, stream);
			}
		return stream;
		}

	/**
	 * Report the result of opening the given output file.
	 *
	 * @param  theJobFrontend  Job Frontend that is calling this method.
	 * @param  bfd             Backend file descriptor.
	 * @param  ffd             Frontend file descriptor if success.
	 * @param  exc             Null if success, exception if failure.
	 */
	public void outputFileOpenResult
		(JobFrontendRef theJobFrontend,
		 int bfd,
		 int ffd,
		 IOException exc)
		{
		BackendFileOutputStream stream = null;
		synchronized (this)
			{
			stream = myOutputStreamForBFD.remove (bfd);
			}
		if (stream != null)
			{
			stream.putResult (ffd, exc);
			}
		}

	/**
	 * Report the result of writing the given output file.
	 *
	 * @param  theJobFrontend  Job Frontend that is calling this method.
	 * @param  ffd             Frontend file descriptor.
	 * @param  exc             Null if success, exception if failure.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileWriteResult
		(JobFrontendRef theJobFrontend,
		 int ffd,
		 IOException exc)
		{
		BackendFileOutputStream stream = null;
		synchronized (this)
			{
			stream = myOutputStreamForFFD.get (ffd);
			}
		if (stream != null)
			{
			stream.putResult (ffd, exc);
			}
		}

	/**
	 * Report the result of flushing the given output file.
	 *
	 * @param  theJobFrontend  Job Frontend that is calling this method.
	 * @param  ffd             Frontend file descriptor.
	 * @param  exc             Null if success, exception if failure.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileFlushResult
		(JobFrontendRef theJobFrontend,
		 int ffd,
		 IOException exc)
		{
		BackendFileOutputStream stream = null;
		synchronized (this)
			{
			stream = myOutputStreamForFFD.get (ffd);
			}
		if (stream != null)
			{
			stream.putResult (ffd, exc);
			}
		}

	/**
	 * Report the result of closing the given output file.
	 *
	 * @param  theJobFrontend  Job Frontend that is calling this method.
	 * @param  ffd             Frontend file descriptor.
	 * @param  exc             Null if success, exception if failure.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileCloseResult
		(JobFrontendRef theJobFrontend,
		 int ffd,
		 IOException exc)
		{
		BackendFileOutputStream stream = null;
		synchronized (this)
			{
			stream = myOutputStreamForFFD.remove (ffd);
			}
		if (stream != null)
			{
			stream.putResult (ffd, exc);
			}
		}

	}
