//******************************************************************************
//
// File:    JobSchedulerRef.java
// Package: edu.rit.pj.cluster
// Unit:    Interface edu.rit.pj.cluster.JobSchedulerRef
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

/**
 * Interface JobSchedulerRef specifies the interface for the PJ Job Scheduler
 * Daemon process.
 *
 * @author  Alan Kaminsky
 * @version 24-Jan-2012
 */
public interface JobSchedulerRef
	{

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
		throws IOException;

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
		throws IOException;

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
		throws IOException;

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
		throws IOException;

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
		throws IOException;

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
		throws IOException;

	/**
	 * Close communication with this Job Scheduler.
	 */
	public void close();

	}
