//******************************************************************************
//
// File:    JobFrontendStub.java
// Package: edu.rit.pj.cluster.test
// Unit:    Class edu.rit.pj.cluster.test.JobFrontendStub
//
// This Java source file is copyright (C) 2012 by Alan Kaminsky. All rights
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

package edu.rit.pj.cluster.test;

import java.io.File;
import java.net.InetSocketAddress;

import edu.rit.mp.ChannelGroup;
import edu.rit.mp.Status;
import edu.rit.mp.buf.ObjectItemBuf;
import edu.rit.pj.cluster.JobBackendRef;
import edu.rit.pj.cluster.JobFrontendMessage;
import edu.rit.pj.cluster.JobFrontendRef;
import edu.rit.pj.cluster.JobSchedulerRef;
import edu.rit.pj.cluster.Message;

/**
 * Class JobFrontendStub provides a stub for the PJ job frontend process. It
 * receives job frontend messages from a channel and prints each message on the
 * console.
 *
 * @author  Alan Kaminsky
 * @version 24-Jan-2012
 */
public class JobFrontendStub
	extends Thread
	implements JobFrontendRef
	{

// Hidden data members.

	private ChannelGroup myChannelGroup;

// Exported constructors.

	/**
	 * Construct a new job frontend stub. The stub will receive messages from
	 * the given channel group.
	 *
	 * @param  theChannelGroup  Channel group.
	 */
	public JobFrontendStub
		(ChannelGroup theChannelGroup)
		{
		myChannelGroup = theChannelGroup;
		}

// Exported operations.

	/**
	 * Run this job frontend stub.
	 */
	public void run()
		{
		ObjectItemBuf<JobFrontendMessage> buf =
			new ObjectItemBuf<JobFrontendMessage>();
		Status status = null;
		JobFrontendMessage message = null;

		for (;;)
			{
			try
				{
				// Receive a message from any channel.
				status = myChannelGroup.receive (null, null, buf);
				message = buf.item;

				// Process message.
				if (status.tag == Message.FROM_JOB_BACKEND)
					{
					message.invoke (this, (JobBackendRef) null);
					}
				else if (status.tag == Message.FROM_JOB_SCHEDULER)
					{
					message.invoke (this, (JobSchedulerRef) null);
					}
				else
					{
					System.out.println ("Bad tag = " + status.tag);
					}

				// Enable garbage collection of no-longer-needed objects while
				// waiting to receive next message.
				buf.item = null;
				status = null;
				message = null;
				}

			catch (Throwable exc)
				{
				// Log any exception but otherwise ignore it.
				exc.printStackTrace (System.out);
				}
			}
		}

	/**
	 * Assign a backend process to the job.
	 *
	 * @param  theJobScheduler  Job Scheduler that is calling this method.
	 * @param  name             Backend node name.
	 * @param  host             Host name for SSH remote login.
	 * @param  jvm              Full pathname of Java Virtual Machine.
	 * @param  classpath        Java class path for PJ Library.
	 * @param  jvmflags         Array of JVM command line flags.
	 * @param  Nt               Number of CPUs assigned to the process.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void assignBackend
		(JobSchedulerRef theJobScheduler,
		 String name,
		 String host,
		 String jvm,
		 String classpath,
		 String[] jvmflags,
		 int Nt)
		{
		StringBuilder b = new StringBuilder();
		b.append ("assignBackend (theJobScheduler, \"");
		b.append (name);
		b.append ("\", \"");
		b.append (host);
		b.append ("\", \"");
		b.append (jvm);
		b.append ("\", \"");
		b.append (classpath);
		b.append ("\", {");
		for (int i = 0; i < jvmflags.length; ++ i)
			{
			if (i > 0) b.append (',');
			b.append ('"');
			b.append (jvmflags[i]);
			b.append ('"');
			}
		b.append ("}, ");
		b.append (Nt);
		b.append (")");
		System.out.println (b);
		}

	/**
	 * Assign a job number to the job. The host name for the job frontend's
	 * middleware channel group is also specified.
	 *
	 * @param  theJobScheduler  Job Scheduler that is calling this method.
	 * @param  jobnum           Job number.
	 * @param  pjhost           Host name for middleware channel group.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void assignJobNumber
		(JobSchedulerRef theJobScheduler,
		 int jobnum,
		 String pjhost)
		{
		System.out.println
			("assignJobNumber (theJobScheduler, " +
			 jobnum + ", \"" +
			 pjhost + "\")");
		}

	/**
	 * Cancel the job.
	 *
	 * @param  theJobScheduler  Job Scheduler that is calling this method.
	 * @param  errmsg           Error message string.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void cancelJob
		(JobSchedulerRef theJobScheduler,
		 String errmsg)
		{
		System.out.println
			("cancelJob (theJobScheduler, \"" +
			 errmsg + "\")");
		}

	/**
	 * Renew the lease on the job.
	 *
	 * @param  theJobScheduler  Job Scheduler that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void renewLease
		(JobSchedulerRef theJobScheduler)
		{
		System.out.println
			("renewLease (theJobScheduler)");
		}

	/**
	 * Report that a backend process has finished executing the job.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void backendFinished
		(JobBackendRef theJobBackend)
		{
		System.out.println
			("backendFinished (theJobBackend)");
		}

	/**
	 * Report that a backend process is ready to commence executing the job.
	 *
	 * @param  theJobBackend
	 *     Job Backend that is calling this method.
	 * @param  rank
	 *     Rank of the job backend process.
	 * @param  middlewareAddress
	 *     Host/port to which the job backend process is listening for
	 *     middleware messages.
	 * @param  worldAddress
	 *     Host/port to which the job backend process is listening for the world
	 *     communicator.
	 * @param  frontendAddress
	 *     Host/port to which the job backend process is listening for the
	 *     frontend communicator, or null if the frontend communicator does not
	 *     exist.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void backendReady
		(JobBackendRef theJobBackend,
		 int rank,
		 InetSocketAddress middlewareAddress,
		 InetSocketAddress worldAddress,
		 InetSocketAddress frontendAddress)
		{
		System.out.println
			("backendReady (theJobBackend, " +
			 rank + ", " +
			 middlewareAddress + ", " +
			 worldAddress + ", " +
			 frontendAddress + ")");
		}

	/**
	 * Cancel the job.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  errmsg         Error message string.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void cancelJob
		(JobBackendRef theJobBackend,
		 String errmsg)
		{
		System.out.println
			("cancelJob (theJobBackend, \"" +
			 errmsg + "\")");
		}

	/**
	 * Renew the lease on the job.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void renewLease
		(JobBackendRef theJobBackend)
		{
		System.out.println
			("renewLease (theJobBackend)");
		}

	/**
	 * Request the given resource from this job frontend's class loader.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  resourceName   Resource name.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void requestResource
		(JobBackendRef theJobBackend,
		 String resourceName)
		{
		System.out.println
			("requestResource (theJobBackend, \"" +
			 resourceName + "\")");
		}

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
		{
		System.out.println
			("outputFileOpen (theJobBackend, " +
			 bfd + ", " +
			 file + ", " +
			 append + ")");
		}

	/**
	 * Write the given bytes to the given output file. <TT>ffd</TT> = 1 refers
	 * to the job's standard output stream; <TT>ffd</TT> = 2 refers to the job's
	 * standard error stream; other values refer to a previously opened file.
	 *
	 * @param  theJobBackend  Job Backend that is calling this method.
	 * @param  ffd            Frontend file descriptor.
	 * @param  buf            Array of bytes to write.
	 * @param  off            Index of first byte to write.
	 * @param  len            Number of bytes to write.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void outputFileWrite
		(JobBackendRef theJobBackend,
		 int ffd,
		 byte[] buf,
		 int off,
		 int len)
		{
		System.out.println
			("outputFileWrite (theJobBackend, " +
			 ffd + ", byte[]{\"" +
			 new String (buf) + "\"}, " +
			 off + ", " +
			 len + ")");
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
		{
		System.out.println
			("outputFileFlush (theJobBackend, " +
			 ffd + ")");
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
		{
		System.out.println
			("outputFileClose (theJobBackend, " +
			 ffd + ")");
		}

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
		{
		System.out.println
			("inputFileOpen (theJobBackend, " +
			 bfd + ", " +
			 file + ")");
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
		{
		System.out.println
			("inputFileRead (theJobBackend, " +
			 ffd + ", " +
			 len + ")");
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
		{
		System.out.println
			("inputFileSkip (theJobBackend, " +
			 ffd + ", " +
			 len + ")");
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
		{
		System.out.println
			("inputFileClose (theJobBackend, " +
			 ffd + ")");
		}

	/**
	 * Report a comment for a process.
	 *
	 * @param  theJobBackend  Job backend that is calling this method.
	 * @param  rank           Process rank.
	 * @param  comment        Comment string.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void reportComment
		(JobBackendRef theJobBackend,
		 int rank,
		 String comment)
		{
		System.out.println
			("reportComment (theJobBackend, " +
			 rank + ", \"" +
			 comment + "\")");
		}

	/**
	 * Close communication with this Job Frontend.
	 */
	public void close()
		{
		System.out.println
			("close()");
		}

	}
