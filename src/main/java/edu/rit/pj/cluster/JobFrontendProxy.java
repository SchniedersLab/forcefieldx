//******************************************************************************
//
// File:    JobFrontendProxy.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobFrontendProxy
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

package edu.rit.pj.cluster;

import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;

import edu.rit.mp.ByteBuf;

import edu.rit.util.Range;

import java.io.File;
import java.io.IOException;

import java.net.InetSocketAddress;

/**
 * Class JobFrontendProxy provides a proxy object for sending messages to a PJ
 * job frontend process.
 *
 * @author  Alan Kaminsky
 * @version 20-Jun-2012
 */
public class JobFrontendProxy
	extends Proxy
	implements JobFrontendRef
	{

// Exported constructors.

	/**
	 * Construct a new job frontend proxy. The proxy will use the given channel
	 * in the given channel group to send messages to the job frontend process.
	 *
	 * @param  theChannelGroup  Channel group.
	 * @param  theChannel       Channel.
	 */
	public JobFrontendProxy
		(ChannelGroup theChannelGroup,
		 Channel theChannel)
		{
		super (theChannelGroup, theChannel);
		}

// Exported operations.

	/**
	 * Assign a backend process to the job.
	 *
	 * @param  theJobScheduler  Job Scheduler that is calling this method.
	 * @param  name             Backend node name.
	 * @param  host             Host name for SSH remote login.
	 * @param  jvm              Full pathname of Java Virtual Machine.
	 * @param  classpath        Java class path for PJ Library.
	 * @param  jvmflags         Array of JVM command line flags.
	 * @param  shellCommand     Shell command string.
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
		 String shellCommand,
		 int Nt)
		throws IOException
		{
		send
			(JobFrontendMessage.assignBackend
				(theJobScheduler, name, host, jvm, classpath, jvmflags,
				 shellCommand, Nt));
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
		throws IOException
		{
		send
			(JobFrontendMessage.assignJobNumber
				(theJobScheduler, jobnum, pjhost));
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
		throws IOException
		{
		send (JobFrontendMessage.cancelJob (theJobScheduler, errmsg));
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
		throws IOException
		{
		send (JobFrontendMessage.renewLease (theJobScheduler));
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
		throws IOException
		{
		send (JobFrontendMessage.backendFinished (theJobBackend));
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
		throws IOException
		{
		send
			(JobFrontendMessage.backendReady
				(theJobBackend, rank, middlewareAddress,
				 worldAddress, frontendAddress));
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
		throws IOException
		{
		send (JobFrontendMessage.cancelJob (theJobBackend, errmsg));
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
		throws IOException
		{
		send (JobFrontendMessage.renewLease (theJobBackend));
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
		throws IOException
		{
		send (JobFrontendMessage.requestResource (theJobBackend, resourceName));
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
		throws IOException
		{
		send
			(JobFrontendMessage.outputFileOpen
				(theJobBackend, bfd, file, append));
		}

	/**
	 * Write the given bytes to the given output file. <TT>fd</TT> = 1 refers to
	 * the job's standard output stream; <TT>fd</TT> = 2 refers to the job's
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
		throws IOException
		{
		send
			(JobFrontendMessage.outputFileWrite (theJobBackend, ffd, len));
		send
			(ffd,
			 ByteBuf.sliceBuffer (buf, new Range (off, off+len-1)));
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
		send (JobFrontendMessage.outputFileFlush (theJobBackend, ffd));
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
		send (JobFrontendMessage.outputFileClose (theJobBackend, ffd));
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
		throws IOException
		{
		send (JobFrontendMessage.inputFileOpen (theJobBackend, bfd, file));
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
		send (JobFrontendMessage.inputFileRead (theJobBackend, ffd, len));
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
		send (JobFrontendMessage.inputFileSkip (theJobBackend, ffd, len));
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
		send (JobFrontendMessage.inputFileClose (theJobBackend, ffd));
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
		throws IOException
		{
		send (JobFrontendMessage.reportComment (theJobBackend, rank, comment));
		}

	}
