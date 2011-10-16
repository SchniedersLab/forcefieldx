//******************************************************************************
//
// File:    JobInfo.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobInfo
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

package edu.rit.pj.cluster;

import edu.rit.util.Timer;

/**
 * Class JobInfo provides a record of information about one job in a parallel
 * computer in the PJ cluster middleware.
 *
 * @author  Alan Kaminsky
 * @version 21-May-2008
 */
public class JobInfo
	{

// Exported enumerations.

	/**
	 * The state of a job.
	 */
	public static enum State
		{
		/**
		 * The job is waiting to run.
		 */
		WAITING ("Waiting"),

		/**
		 * The job is running.
		 */
		RUNNING ("Running");

		private final String stringForm;

		/**
		 * Construct a new State value.
		 *
		 * @param  stringForm  String form.
		 */
		State
			(String stringForm)
			{
			this.stringForm = stringForm;
			}

		/**
		 * Returns a string version of this State value.
		 *
		 * @return  String version.
		 */
		public String toString()
			{
			return stringForm;
			}
		}

// Exported data members.

	/**
	 * The job number.
	 */
	public int jobnum;

	/**
	 * The job's state.
	 */
	public State state;

	/**
	 * The time when the job entered its current state (milliseconds since
	 * midnight 01-Jan-1970 GMT).
	 */
	public long stateTime;

	/**
	 * The job's user name.
	 */
	public String username;

	/**
	 * The number of backend nodes in the job.
	 */
	public int Nn;

	/**
	 * The number of processes in the job.
	 */
	public int Np;

	/**
	 * The number of CPUs per process in the job.
	 */
	public int Nt;

	/**
	 * The number of processes that have been assigned to the job so far.
	 */
	public int count;

	/**
	 * Array of backend nodes for each process assigned to the job in rank
	 * order. The array has <TT>Np</TT> total elements. The first <TT>count</TT>
	 * elements have been assigned.
	 */
	public BackendInfo[] backend;

	/**
	 * Number of CPUs assigned to each process in the job in rank order. The
	 * array has <TT>Np</TT> total elements. The first <TT>count</TT> elements
	 * have been assigned.
	 */
	public int[] cpus;

	/**
	 * The number of nodes that have been assigned to the job so far.
	 */
	public int nodeCount;

	/**
	 * Reference to the job frontend process.
	 */
	public JobFrontendRef frontend;

	/**
	 * Lease renewal timer.
	 */
	public Timer renewTimer;

	/**
	 * Lease expiration timer.
	 */
	public Timer expireTimer;

	/**
	 * Maximum job time timer.
	 */
	public Timer jobTimer;

// Exported constructors.

	/**
	 * Construct a new job information record.
	 *
	 * @param  jobnum       The job number.
	 * @param  state        The job's state.
	 * @param  stateTime    The time when the job entered its current state.
	 * @param  username     The job's user name.
	 * @param  Nn           The number of backend nodes in the job.
	 * @param  Np           The number of processes in the job.
	 * @param  Nt           The number of CPUs per process in the job.
	 * @param  count        The number of processes that have been assigned to
	 *                      the job so far.
	 * @param  backend      Array of backends assigned to the job in rank order.
	 * @param  cpus         Array of CPUs for each process in rank order.
	 * @param  nodeCount    The number of nodes that have been assigned to the
	 *                      job so far.
	 * @param  frontend     Reference to the job frontend process.
	 * @param  renewTimer   Lease renewal timer.
	 * @param  expireTimer  Lease expiration timer.
	 * @param  jobTimer     Maximum job time timer.
	 */
	public JobInfo
		(int jobnum,
		 State state,
		 long stateTime,
		 String username,
		 int Nn,
		 int Np,
		 int Nt,
		 int count,
		 BackendInfo[] backend,
		 int[] cpus,
		 int nodeCount,
		 JobFrontendRef frontend,
		 Timer renewTimer,
		 Timer expireTimer,
		 Timer jobTimer)
		{
		this.jobnum = jobnum;
		this.state = state;
		this.stateTime = stateTime;
		this.username = username;
		this.Nn = Nn;
		this.Np = Np;
		this.Nt = Nt;
		this.count = count;
		this.backend = backend;
		this.cpus = cpus;
		this.nodeCount = nodeCount;
		this.frontend = frontend;
		this.renewTimer = renewTimer;
		this.expireTimer = expireTimer;
		this.jobTimer = jobTimer;
		}

	}
