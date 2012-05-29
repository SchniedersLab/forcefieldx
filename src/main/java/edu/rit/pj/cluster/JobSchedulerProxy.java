//******************************************************************************
//
// File:    JobSchedulerProxy.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobSchedulerProxy
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

import java.io.IOException;

import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;

/**
 * Class JobSchedulerProxy provides a proxy object for sending messages to a PJ
 * job scheduler process.
 *
 * @author  Alan Kaminsky
 * @version 24-Jan-2012
 */
public class JobSchedulerProxy
	extends Proxy
	implements JobSchedulerRef
	{

// Exported constructors.

	/**
	 * Construct a new job scheduler proxy. The proxy will use the given channel
	 * in the given channel group to send messages to the job scheduler process.
	 *
	 * @param  theChannelGroup  Channel group.
	 * @param  theChannel       Channel.
	 */
	public JobSchedulerProxy
		(ChannelGroup theChannelGroup,
		 Channel theChannel)
		{
		super (theChannelGroup, theChannel);
		}

// Exported operations.

	/**
	 * Report that a backend node failed.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  name            Backend node name.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void backendFailed
		(JobFrontendRef theJobFrontend,
		 String name)
		throws IOException
		{
		send (JobSchedulerMessage.backendFailed (theJobFrontend, name));
		}

	/**
	 * Cancel a job.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  errmsg          Error message string.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void cancelJob
		(JobFrontendRef theJobFrontend,
		 String errmsg)
		throws IOException
		{
		send (JobSchedulerMessage.cancelJob (theJobFrontend, errmsg));
		}

	/**
	 * Report that a job finished.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void jobFinished
		(JobFrontendRef theJobFrontend)
		throws IOException
		{
		send (JobSchedulerMessage.jobFinished (theJobFrontend));
		}

	/**
	 * Renew the lease on a job.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void renewLease
		(JobFrontendRef theJobFrontend)
		throws IOException
		{
		send (JobSchedulerMessage.renewLease (theJobFrontend));
		}

	/**
	 * Report a comment for a process.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  rank            Process rank.
	 * @param  comment         Comment string.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void reportComment
		(JobFrontendRef theJobFrontend,
		 int rank,
		 String comment)
		throws IOException
		{
		send
			(JobSchedulerMessage.reportComment
				(theJobFrontend, rank, comment));
		}

	/**
	 * Request that a job be scheduled.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  username        User name.
	 * @param  Nn              Number of backend nodes.
	 * @param  Np              Number of processes.
	 * @param  Nt              Number of CPUs per process. 0 means "all CPUs."
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void requestJob
		(JobFrontendRef theJobFrontend,
		 String username,
		 int Nn,
		 int Np,
		 int Nt)
		throws IOException
		{
		send
			(JobSchedulerMessage.requestJob
				(theJobFrontend, username, Nn, Np, Nt));
		}

	}
