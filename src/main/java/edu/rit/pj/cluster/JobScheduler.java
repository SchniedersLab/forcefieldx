//******************************************************************************
//
// File:    JobScheduler.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobScheduler
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

import edu.rit.http.HttpRequest;
import edu.rit.http.HttpResponse;
import edu.rit.http.HttpServer;

import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;
import edu.rit.mp.ChannelGroupClosedException;
import edu.rit.mp.ConnectListener;
import edu.rit.mp.Status;

import edu.rit.mp.ObjectBuf;

import edu.rit.mp.buf.ObjectItemBuf;

import edu.rit.pj.Version;

import edu.rit.util.Logger;
import edu.rit.util.PrintStreamLogger;
import edu.rit.util.Timer;
import edu.rit.util.TimerTask;
import edu.rit.util.TimerThread;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;

import java.net.InetSocketAddress;

import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Scanner;

/**
 * Class JobScheduler is the main program for the PJ Job Scheduler Daemon
 * process for a parallel computer.
 * <P>
 * Run the Job Scheduler Daemon on the cluster's frontend processor by typing
 * this command:
 * <P>
 * java edu.rit.pj.cluster.JobScheduler <I>configfile</I>
 * <BR><I>configfile</I> = Configuration file name
 * <P>
 * For further information about the configuration file, see class {@linkplain
 * Configuration}.
 *
 * @author  Alan Kaminsky
 * @version 20-Jun-2012
 */
public class JobScheduler
	implements JobSchedulerRef
	{

// Hidden data members.

	// Cluster name.
	private String myClusterName;

	// Log file.
	private Logger myLog;

	// Web interface host and port.
	private String myWebHost;
	private int myWebPort;

	// Job Scheduler host and port.
	private String mySchedulerHost;
	private int mySchedulerPort;

	// Job frontend host.
	private String myFrontendHost;

	// Maximum job time, or 0 if no maximum.
	private int myJobTime;

	// Mapping from backend processor name to backend info.
	private Map<String,BackendInfo> myNameToBackendMap =
		new HashMap<String,BackendInfo>();

	// Array of backend info records.
	private BackendInfo[] myBackendInfo;
	private int myBackendCount;

	// Next backend number to assign to a job.
	private int myNextBackendNumber = 0;

	// Next job number.
	private int myNextJobNumber = 1;

	// Mapping from job frontend to job info.
	private Map<JobFrontendRef,JobInfo> myFrontendToJobMap =
		new HashMap<JobFrontendRef,JobInfo>();

	// Queue of running jobs.
	private List<JobInfo> myRunningJobList =
		new LinkedList<JobInfo>();

	// Queue of waiting jobs.
	private List<JobInfo> myWaitingJobList =
		new LinkedList<JobInfo>();

	// Timer thread for lease renewals and expirations.
	private TimerThread myLeaseTimerThread;

	// Channel group for communicating with job frontend processes.
	private ChannelGroup myChannelGroup;

	// Server for web interface.
	private HttpServer myHttpServer;

	// Total compute time (msec) of all jobs.
	private long myTotalComputeTime;

	// Date and time when Job Scheduler started.
	private long myStartDateTime;

// Hidden constructors.

	/**
	 * Construct a new Job Scheduler Daemon.
	 *
	 * @param  configfile  Configuration file name.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private JobScheduler
		(String configfile)
		throws IOException
		{
		long now = System.currentTimeMillis();
		myStartDateTime = now;

		// Parse configuration file.
		Configuration config = new Configuration (configfile);
		myClusterName = config.getClusterName();
		myLog =
			new PrintStreamLogger
				(new PrintStream
					(new FileOutputStream (config.getLogFile(), true),
					 true));
		myWebHost = config.getWebHost();
		myWebPort = config.getWebPort();
		mySchedulerHost = config.getSchedulerHost();
		mySchedulerPort = config.getSchedulerPort();
		myFrontendHost = config.getFrontendHost();
		myJobTime = config.getJobTime();
		myBackendCount = config.getBackendCount();
		myBackendInfo = new BackendInfo [myBackendCount];
		for (int i = 0; i < myBackendCount; ++ i)
			{
			BackendInfo backendinfo = config.getBackendInfo (i);
			myNameToBackendMap.put (backendinfo.name, backendinfo);
			myBackendInfo[i] = backendinfo;
			}

		// Log startup.
		myLog.log (now, "Started " + Version.PJ_VERSION);

		// Set up shutdown hook.
		Runtime.getRuntime().addShutdownHook (new Thread()
			{
			public void run()
				{
				shutdown();
				}
			});

		// Set up lease timer thread.
		myLeaseTimerThread = new TimerThread();
		myLeaseTimerThread.setDaemon (true);
		myLeaseTimerThread.start();

		// Set up channel group.
		myChannelGroup =
			new ChannelGroup
				(new InetSocketAddress (mySchedulerHost, mySchedulerPort),
				 myLog);
		myLog.log (now, "Job Scheduler at " + myChannelGroup.listenAddress());
		myChannelGroup.setConnectListener (new ConnectListener()
			{
			public void nearEndConnected
				(ChannelGroup theChannelGroup,
				 Channel theChannel)
				{
				}
			public void farEndConnected
				(ChannelGroup theChannelGroup,
				 Channel theChannel)
				{
				createJob (theChannel);
				}
			});

		// Set up server for web interface.
		myHttpServer =
			new HttpServer (new InetSocketAddress (myWebHost, myWebPort), myLog)
				{
				protected void process
					(HttpRequest request,
					 HttpResponse response)
					throws IOException
					{
					processHttpRequest (request, response);
					}
				};
		myLog.log (now, "Web interface at " + myHttpServer.getAddress());

		// Log backend nodes.
		for (BackendInfo backend : myBackendInfo)
			{
			myLog.log
				(now,
				 "Backend " + backend.name + " at " + backend.host +
					", " + backend.totalCpus +
					" CPU" + (backend.totalCpus==1?"":"s"));
			}

		// Start accepting jobs.
		myChannelGroup.startListening();
		}

// Hidden operations.

	/**
	 * Create a job associated with the given channel.
	 *
	 * @param  theChannel  Channel for talking to Job Frontend process.
	 */
	private synchronized void createJob
		(Channel theChannel)
		{
		// Create Job Frontend proxy object for the channel.
		JobFrontendRef frontend =
			new JobFrontendProxy (myChannelGroup, theChannel);
		theChannel.info (frontend);

		// Create job information record.
		JobInfo jobinfo = getJobInfo (frontend);

		// Start lease timers.
		jobinfo.renewTimer.start
			(Constants.LEASE_RENEW_INTERVAL,
			 Constants.LEASE_RENEW_INTERVAL);
		jobinfo.expireTimer.start
			(Constants.LEASE_EXPIRE_INTERVAL);
		}

	/**
	 * Run this Job Scheduler.
	 */
	private void run()
		{
		ObjectItemBuf<JobSchedulerMessage> buf =
			ObjectBuf.buffer ((JobSchedulerMessage) null);
		Status status = null;
		JobSchedulerMessage message = null;
		JobFrontendRef frontend = null;

		receiveloop : for (;;)
			{
			// Receive a message from any channel.
			try
				{
				status = myChannelGroup.receive (null, null, buf);
				}
			catch (ChannelGroupClosedException exc)
				{
				// Normal termination.
				break receiveloop;
				}
			catch (Throwable exc)
				{
				myLog.log ("Exception while receiving message", exc);
				break receiveloop;
				}
			message = buf.item;

			// Get job frontend proxy associated with channel.
			frontend = (JobFrontendRef) status.channel.info();

			// Process message.
			try
				{
				message.invoke (this, frontend);
				}
			catch (Throwable exc)
				{
				myLog.log ("Exception while processing message", exc);
				}

			// Enable garbage collection of no-longer-needed objects while
			// waiting to receive next message.
			buf.item = null;
			status = null;
			message = null;
			frontend = null;
			}
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
	public synchronized void backendFailed
		(JobFrontendRef theJobFrontend,
		 String name)
		throws IOException
		{
		BackendInfo backendinfo = myNameToBackendMap.get (name);
		if (backendinfo != null)
			{
			long now = System.currentTimeMillis();
			myLog.log (now, "Backend " + name + " failed");
//			if (backendinfo.state != BackendInfo.State.FAILED)
//				{
//				/*TBD*/ Cancel any reserved or running job
//				backendinfo.state = BackendInfo.State.FAILED;
//				backendinfo.stateTime = now;
//				backendinfo.job = null;
//				assignResourcesToJobs (now);
//				}
			}
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
	public synchronized void cancelJob
		(JobFrontendRef theJobFrontend,
		 String errmsg)
		throws IOException
		{
		JobInfo jobinfo = getJobInfo (theJobFrontend);
		doCancelJob (System.currentTimeMillis(), jobinfo, errmsg);
		}

	/**
	 * Report that a job finished.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public synchronized void jobFinished
		(JobFrontendRef theJobFrontend)
		throws IOException
		{
		JobInfo jobinfo = getJobInfo (theJobFrontend);
		doFinishJob (System.currentTimeMillis(), jobinfo);
		}

	/**
	 * Renew the lease on a job.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public synchronized void renewLease
		(JobFrontendRef theJobFrontend)
		throws IOException
		{
		JobInfo jobinfo = getJobInfo (theJobFrontend);
		jobinfo.expireTimer.start (Constants.LEASE_EXPIRE_INTERVAL);
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
	public synchronized void reportComment
		(JobFrontendRef theJobFrontend,
		 int rank,
		 String comment)
		{
		JobInfo jobinfo = getJobInfo (theJobFrontend);
		jobinfo.comment[rank] = comment;
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
	public synchronized void requestJob
		(JobFrontendRef theJobFrontend,
		 String username,
		 int Nn,
		 int Np,
		 int Nt)
		throws IOException
		{
		JobInfo jobinfo = getJobInfo (theJobFrontend);
		long now = System.currentTimeMillis();
		myLog.log
			(now,
			 "Job " + jobinfo.jobnum + " queued, username=" + username +
				", nn=" + Nn + ", np=" + Np + ", nt=" + Nt);

		// Record job parameters.
		jobinfo.username = username;
		jobinfo.Nn = Math.min (Nn, Np);
		jobinfo.Np = Np;
		jobinfo.Nt = Nt;
		jobinfo.backend = new BackendInfo [Np];
		jobinfo.cpus = new int [Np];
		jobinfo.comment = new String [Np];
		for (int i = 0; i < Np; ++ i) jobinfo.comment[i] = "";

		// If the cluster doesn't have enough resources, cancel the job.
		if (! enoughResourcesForJob (jobinfo.Nn, jobinfo.Np, jobinfo.Nt))
			{
			doCancelJobTooFewResources (now, jobinfo);
			return;
			}

		// Add job to queue of waiting jobs.
		myWaitingJobList.add (jobinfo);

		// Inform job frontend of job number.
		theJobFrontend.assignJobNumber (this, jobinfo.jobnum, myFrontendHost);

		// Assign idle nodes to waiting jobs.
		assignResourcesToJobs (now);
		}

	/**
	 * Close communication with this Job Scheduler.
	 */
	public void close()
		{
		}

// More hidden operations.

	/**
	 * Take action when a job's lease renewal timer times out.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private synchronized void renewTimeout
		(Timer theTimer,
		 JobFrontendRef theJobFrontend)
		throws IOException
		{
		if (theTimer.isTriggered())
			{
			theJobFrontend.renewLease (this);
			}
		}

	/**
	 * Take action when a job's lease expiration timer times out.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private synchronized void expireTimeout
		(Timer theTimer,
		 JobFrontendRef theJobFrontend)
		throws IOException
		{
		if (theTimer.isTriggered())
			{
			JobInfo jobinfo = getJobInfo (theJobFrontend);
			doCancelJob
				(System.currentTimeMillis(),
				 jobinfo,
				 "Job frontend lease expired");
			}
		}

	/**
	 * Take action when a job's maximum job time timer times out.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private synchronized void jobTimeout
		(Timer theTimer,
		 JobFrontendRef theJobFrontend)
		throws IOException
		{
		if (theTimer.isTriggered())
			{
			JobInfo jobinfo = getJobInfo (theJobFrontend);
			String errmsg =
				"Maximum job time (" + myJobTime + " seconds) exceeded";
			jobinfo.frontend.cancelJob (this, errmsg);
			doCancelJob (System.currentTimeMillis(), jobinfo, errmsg);
			}
		}

	/**
	 * Get the job info record associated with the given job frontend.
	 *
	 * @param  frontend  Job frontend.
	 *
	 * @return  Job info record.
	 */
	private JobInfo getJobInfo
		(JobFrontendRef frontend)
		{
		final JobFrontendRef fe = frontend;
		JobInfo jobinfo = myFrontendToJobMap.get (frontend);
		if (jobinfo == null)
			{
			jobinfo = new JobInfo
				(/*jobnum   */ myNextJobNumber ++,
				 /*state    */ JobInfo.State.WAITING,
				 /*stateTime*/ System.currentTimeMillis(),
				 /*username */ null,
				 /*Nn       */ 0,
				 /*Np       */ 0,
				 /*Nt       */ 0,
				 /*count    */ 0,
				 /*backend  */ null,
				 /*cpus     */ null,
				 /*nodeCount*/ 0,
				 /*frontend */ fe,
				 /*renewTimer*/
					myLeaseTimerThread.createTimer (new TimerTask()
						{
						public void action (Timer theTimer)
							{
							try
								{
								renewTimeout (theTimer, fe);
								}
							catch (Throwable exc)
								{
								myLog.log (exc);
								}
							}
						}),
				 /*expireTimer*/
					myLeaseTimerThread.createTimer (new TimerTask()
						{
						public void action (Timer theTimer)
							{
							try
								{
								expireTimeout (theTimer, fe);
								}
							catch (Throwable exc)
								{
								myLog.log (exc);
								}
							}
						}),
				 /*jobTimer*/
					myLeaseTimerThread.createTimer (new TimerTask()
						{
						public void action (Timer theTimer)
							{
							try
								{
								jobTimeout (theTimer, fe);
								}
							catch (Throwable exc)
								{
								myLog.log (exc);
								}
							}
						}));
			myFrontendToJobMap.put (frontend, jobinfo);
			}
		return jobinfo;
		}

	/**
	 * Finish the given job.
	 *
	 * @param  now      Current time.
	 * @param  jobinfo  Job info record.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void doFinishJob
		(long now,
		 JobInfo jobinfo)
		throws IOException
		{
		myLog.log (now, "Job " + jobinfo.jobnum + " finished");
		doCleanupJob (now, jobinfo);
		}

	/**
	 * Cancel the given job.
	 *
	 * @param  now      Current time.
	 * @param  jobinfo  Job info record.
	 * @param  errmsg   Error message.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void doCancelJob
		(long now,
		 JobInfo jobinfo,
		 String errmsg)
		throws IOException
		{
		myLog.log (now, "Job " + jobinfo.jobnum + " canceled: " + errmsg);
		doCleanupJob (now, jobinfo);
		}

	/**
	 * Cancel the given job because of too few resources.
	 *
	 * @param  now      Current time.
	 * @param  jobinfo  Job info record.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void doCancelJobTooFewResources
		(long now,
		 JobInfo jobinfo)
		throws IOException
		{
		String errmsg;
		if (jobinfo.Nt == 0)
			{
			errmsg =
				"Too few resources available to assign " +
				jobinfo.Nn + " node" + (jobinfo.Nn==1?"":"s") + " and " +
				jobinfo.Np + " process" + (jobinfo.Np==1?"":"es");
			}
		else
			{
			errmsg =
				"Too few resources available to assign " +
				jobinfo.Nn + " node" + (jobinfo.Nn==1?"":"s") + ", " +
				jobinfo.Np + " process" + (jobinfo.Np==1?"":"es") + ", and " +
				jobinfo.Nt + " CPU" + (jobinfo.Nt==1?"":"s") + " per process";
			}
		jobinfo.frontend.cancelJob (this, errmsg);
		doCancelJob (now, jobinfo, errmsg);
		}

	/**
	 * Clean up the given job.
	 *
	 * @param  now      Current time.
	 * @param  jobinfo  Job info record.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void doCleanupJob
		(long now,
		 JobInfo jobinfo)
		throws IOException
		{
		// Stop lease timers.
		jobinfo.renewTimer.stop();
		jobinfo.expireTimer.stop();
		jobinfo.jobTimer.stop();

		// Stop communication with job frontend.
		jobinfo.frontend.close();

		// Remove job from queues.
		myFrontendToJobMap.remove (jobinfo.frontend);
		myRunningJobList.remove (jobinfo);
		myWaitingJobList.remove (jobinfo);

		// Make each of the job's nodes idle (but not failed nodes).
		for (int i = 0; i < jobinfo.count; ++ i)
			{
			BackendInfo backendinfo = jobinfo.backend[i];
			if (backendinfo.state != BackendInfo.State.FAILED)
				{
				backendinfo.state = BackendInfo.State.IDLE;
				backendinfo.stateTime = now;
				backendinfo.job = null;
				}
			}

		// Update total compute time.
		myTotalComputeTime += (now - jobinfo.stateTime);

		// Assign idle nodes to waiting jobs.
		assignResourcesToJobs (now);
		}

	/**
	 * Assign idle nodes to waiting jobs.
	 *
	 * @param  now  Current time.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void assignResourcesToJobs
		(long now)
		throws IOException
		{
		// List of jobs to be canceled.
		List<JobInfo> cancelList = new LinkedList<JobInfo>();

		// Decide what to do with each waiting job.
		Iterator<JobInfo> iter = myWaitingJobList.iterator();
		jobLoop : while (iter.hasNext())
			{
			JobInfo jobinfo = iter.next();

			// If the cluster doesn't have enough resources, don't try to
			// reserve any.
			if (! enoughResourcesForJob (jobinfo.Nn, jobinfo.Np, jobinfo.Nt))
				{
				iter.remove();
				cancelList.add (jobinfo);
				continue jobLoop;
				}

			// Used to decide how many processes for each node.
			int Np_div_Nn = jobinfo.Np / jobinfo.Nn;
			int Np_rem_Nn = jobinfo.Np % jobinfo.Nn;

			// Reserve idle nodes for this job until there are no more idle
			// nodes or this job has all the nodes it needs.
			int be = myNextBackendNumber;
			do
				{
				// Decide how many processes for this node.
				int Nproc = Np_div_Nn;
				if (jobinfo.nodeCount < Np_rem_Nn) ++ Nproc;

				// Reserve this node only if it is idle and it has enough CPUs.
				BackendInfo backendinfo = myBackendInfo[be];
				if (backendinfo.state == BackendInfo.State.IDLE &&
						backendinfo.totalCpus >= Nproc)
					{
					// Reserve node.
					backendinfo.state = BackendInfo.State.RESERVED;
					backendinfo.stateTime = now;
					backendinfo.job = jobinfo;

					// Used to decide how many CPUs for each process.
					int Nt_div_Nproc = backendinfo.totalCpus / Nproc;
					int Nt_rem_Nproc = backendinfo.totalCpus % Nproc;

					// Assign Np processes.
					for (int i = 0; i < Nproc; ++ i)
						{
						// Decide how many CPUs for this process.
						int Ncpus = jobinfo.Nt;
						if (Ncpus == 0)
							{
							Ncpus = Nt_div_Nproc;
							if (i < Nt_rem_Nproc) ++ Ncpus;
							}

						// Log information.
						myLog.log
							(now,
							 "Job " + jobinfo.jobnum + " assigned " +
							 backendinfo.name + ", rank=" + jobinfo.count +
							 ", CPUs=" + Ncpus);

						// Record information about process.
						jobinfo.backend[jobinfo.count] = backendinfo;
						jobinfo.cpus[jobinfo.count] = Ncpus;
						++ jobinfo.count;

						// Inform Job Frontend.
						jobinfo.frontend.assignBackend
							(/*theJobScheduler*/ this,
							 /*name           */ backendinfo.name,
							 /*host           */ backendinfo.host,
							 /*jvm            */ backendinfo.jvm,
							 /*classpath      */ backendinfo.classpath,
							 /*jvmflags       */ backendinfo.jvmflags,
							 /*shellCommand   */ backendinfo.shellCommand,
							 /*Nt             */ Ncpus);
						}

					// Assign one node.
					++ jobinfo.nodeCount;
					}

				// Consider next node.
				be = (be + 1) % myBackendCount;
				}
			while (be != myNextBackendNumber && jobinfo.count < jobinfo.Np);
			myNextBackendNumber = be;

			// If this job now has Np processes, start running this job.
			if (jobinfo.count == jobinfo.Np)
				{
				// Log information.
				myLog.log (now, "Job " + jobinfo.jobnum + " started");

				// Mark job as running.
				iter.remove();
				myRunningJobList.add (jobinfo);
				jobinfo.state = JobInfo.State.RUNNING;
				jobinfo.stateTime = now;

				// Mark all the job's nodes as running.
				for (BackendInfo backendinfo : jobinfo.backend)
					{
					backendinfo.state = BackendInfo.State.RUNNING;
					backendinfo.stateTime = now;
					}

				// If the Job Scheduler is imposing a maximum job time, start
				// job timer.
				if (myJobTime > 0)
					{
					jobinfo.jobTimer.start (myJobTime * 1000L);
					}
				}

			// If this job does not yet have Np processes, don't schedule any
			// further jobs.
			else
				{
				break jobLoop;
				}
			}

		// Cancel jobs for which there are insufficient resources.
		for (JobInfo jobinfo : cancelList)
			{
			doCancelJobTooFewResources (now, jobinfo);
			}
		}

	/**
	 * Determine if there are enough resources to run a job.
	 *
	 * @param  Nn  Number of backend nodes required.
	 * @param  Np  Number of processes required.
	 * @param  Nt  Number of CPUs per process required. 0 means "all CPUs."
	 *
	 * @return  True if there are enough resources, false if not.
	 */
	private boolean enoughResourcesForJob
		(int Nn,
		 int Np,
		 int Nt)
		{
		// Determine worst-case processes per node.
		int Ppn = (Np + Nn - 1) / Nn;

		// If number of CPUs per process is "all CPUs," assume one CPU per
		// process.
		if (Nt == 0) Nt = 1;

		// Count how many nodes meet the requirements.
		int nodeCount = 0;
		for (BackendInfo backendinfo : myBackendInfo)
			{
			// The node must not have failed.
			if (backendinfo.state != BackendInfo.State.FAILED &&

			// The node must have at least Ppn*Nt CPUs.
					backendinfo.totalCpus >= Ppn*Nt)
				{
				// The node meets the requirements.
				++ nodeCount;
				}
			}

		// Return outcome.
		return nodeCount >= Nn;
		}

	/**
	 * Process the given HTTP request.
	 *
	 * @param  request   HTTP request.
	 * @param  response  HTTP response.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void processHttpRequest
		(HttpRequest request,
		 HttpResponse response)
		throws IOException
		{
		long now = System.currentTimeMillis();

		// Reject an invalid HTTP request.
		if (! request.isValid())
			{
			response.setStatusCode
				(HttpResponse.Status.STATUS_400_BAD_REQUEST);
			PrintWriter out = response.getPrintWriter();
			printStatusHtmlStart (out, now);
			out.println ("<P>");
			out.println ("400 Bad Request");
			printStatusHtmlEnd (out);
			}

		// Reject all methods except GET.
		else if (! request.getMethod().equals (HttpRequest.GET_METHOD))
			{
			response.setStatusCode
				(HttpResponse.Status.STATUS_501_NOT_IMPLEMENTED);
			PrintWriter out = response.getPrintWriter();
			printStatusHtmlStart (out, now);
			out.println ("<P>");
			out.println ("501 Not Implemented");
			printStatusHtmlEnd (out);
			}

		// Print the status document.
		else if (request.getUri().equals ("/") ||
					request.getUri().equals ("/?"))
			{
			PrintWriter out = response.getPrintWriter();
			printStatusHtmlStart (out, now);
			printStatusHtmlBody (out, now);
			printStatusHtmlEnd (out);
			}

		// Print the debug document.
		else if (request.getUri().equals ("/debug"))
			{
			PrintWriter out = response.getPrintWriter();
			printDebugHtmlStart (out, now);
			printDebugHtmlBody (out);
			printStatusHtmlEnd (out);
			}

		// Print the detailed job status document.
		else if (request.getUri().startsWith ("/job/"))
			{
			String jobString = request.getUri().substring (5);
			try
				{
				int jobNum = Integer.parseInt (jobString);
				PrintWriter out = response.getPrintWriter();
				printJobDetailHtmlStart (out, now, jobNum);
				printJobDetailHtmlBody (out, now, jobNum);
				printStatusHtmlEnd (out);
				}
			catch (NumberFormatException exc)
				{
				PrintWriter out = response.getPrintWriter();
				printErrorHtmlStart (out);
				out.printf ("<P>Invalid job number \"%s\"</P>\n", jobString);
				printErrorHtmlEnd (out);
				}
			}

		// Reject all other URIs.
		else
			{
			response.setStatusCode
				(HttpResponse.Status.STATUS_404_NOT_FOUND);
			PrintWriter out = response.getPrintWriter();
			printErrorHtmlStart (out);
			out.println ("<P>404 Not Found</P>");
			printErrorHtmlEnd (out);
			}

		// Send the response.
		response.close();
		}

	/**
	 * Print the start of the status HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 * @param  now  Current time.
	 */
	private void printStatusHtmlStart
		(PrintWriter out,
		 long now)
		{
		out.println ("<HTML>");
		out.println ("<HEAD>");
		out.print   ("<TITLE>");
		out.print   (myClusterName);
		out.println ("</TITLE>");
		out.print   ("<META HTTP-EQUIV=\"refresh\" CONTENT=\"20;url=");
		printWebInterfaceURL (out);
		out.println ("\">");
		out.println ("<STYLE TYPE=\"text/css\">");
		out.println ("<!--");
		out.println ("* {font-family: Arial, Helvetica, Sans-Serif;}");
		out.println ("body {font-size: small;}");
		out.println ("h1 {font-size: 140%; font-weight: bold;}");
		out.println ("table {font-size: 100%;}");
		out.println ("-->");
		out.println ("</STYLE>");
		out.println ("</HEAD>");
		out.println ("<BODY>");
		out.print   ("<H1>");
		out.print   (myClusterName);
		out.println ("</H1>");
		out.println ("<P>");
		out.print   ("<FORM ACTION=\"");
		printWebInterfaceURL (out);
		out.println ("\" METHOD=\"get\">");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"center\">");
		out.print   ("<INPUT TYPE=\"submit\" VALUE=\"Refresh\">");
		out.println ("</TD>");
		out.println ("<TD WIDTH=20> </TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"center\">");
		out.print   (new Date (now));
		out.print   (" -- ");
		out.print   (Version.PJ_VERSION);
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</FORM>");
		}

	/**
	 * Print the body of the status HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 * @param  now  Current time.
	 */
	private synchronized void printStatusHtmlBody
		(PrintWriter out,
		 long now)
		{
		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"center\" VALIGN=\"top\">");

		out.println ("Nodes");
		out.println ("<TABLE BORDER=1 CELLPADDING=3 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");

		out.println ("<TABLE BORDER=0 CELLPADDING=3 CELLSPACING=0>");
		printBackendLabels (out);
		int i = 0;
		for (BackendInfo backend : myBackendInfo)
			{
			printBackendInfo (out, now, backend, i);
			++ i;
			}
		out.println ("</TABLE>");

		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");

		out.println ("</TD>");
		out.println ("<TD WIDTH=40> </TD>");
		out.println ("<TD ALIGN=\"center\" VALIGN=\"top\">");

		out.println ("Jobs");
		out.println ("<TABLE BORDER=1 CELLPADDING=3 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");

		out.println ("<TABLE BORDER=0 CELLPADDING=3 CELLSPACING=0>");
		printJobLabels (out);
		i = 0;
		for (JobInfo job : myRunningJobList)
			{
			printJobInfo (out, now, job, i);
			++ i;
			}
		for (JobInfo job : myWaitingJobList)
			{
			printJobInfo (out, now, job, i);
			++ i;
			}
		out.println ("</TABLE>");

		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");

		printTotalComputeTime (out);
		out.print ("<BR>");
		printJobCount (out);
		out.println ("<BR>Since " + new Date (myStartDateTime));

		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		}

	/**
	 * Print the job count.
	 *
	 * @param  out  Print writer.
	 */
	private void printJobCount
		(PrintWriter out)
		{
		if (myNextJobNumber == 2)
			{
			out.print ("1 job");
			}
		else
			{
			out.print (myNextJobNumber-1);
			out.print (" jobs");
			}
		out.println (" served");
		}

	/**
	 * Print the total CPU time.
	 *
	 * @param  out  Print writer.
	 */
	private void printTotalComputeTime
		(PrintWriter out)
		{
		if (myTotalComputeTime < 1000000L)
			{
			out.print (myTotalComputeTime / 1000L);
			}
		else if (myTotalComputeTime < 1000000000L)
			{
			out.print ("Over ");
			out.print (myTotalComputeTime / 1000000L);
			out.print (" thousand");
			}
		else if (myTotalComputeTime < 1000000000000L)
			{
			out.print ("Over ");
			out.print (myTotalComputeTime / 1000000000L);
			out.print (" million");
			}
		else if (myTotalComputeTime < 1000000000000000L)
			{
			out.print ("Over ");
			out.print (myTotalComputeTime / 1000000000000L);
			out.print (" billion");
			}
		else
			{
			out.print ("Over ");
			out.print (myTotalComputeTime / 1000000000000000L);
			out.print (" trillion");
			}
		out.println (" CPU seconds served");
		}

	/**
	 * Print the end of the status HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printStatusHtmlEnd
		(PrintWriter out)
		{
		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("Job queue web interface:&nbsp;&nbsp;");
		out.println ("</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<A HREF=\"");
		printWebInterfaceURL (out);
		out.print   ("\">");
		printWebInterfaceURL (out);
		out.println ("</A>");
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("Powered by Parallel Java:&nbsp;&nbsp;");
		out.println ("</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("<A HREF=\"http://www.cs.rit.edu/~ark/pj.shtml\">http://www.cs.rit.edu/~ark/pj.shtml</A>");
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("Developed by Alan Kaminsky:&nbsp;&nbsp;");
		out.println ("</TD>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.println ("<A HREF=\"http://www.cs.rit.edu/~ark/\">http://www.cs.rit.edu/~ark/</A>");
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</BODY>");
		out.println ("</HTML>");
		}

	/**
	 * Print the web interface URL on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printWebInterfaceURL
		(PrintWriter out)
		{
		out.printf ("http://%s:%d/", myWebHost, myWebPort);
		}

	/**
	 * Print the URL for the given job number on the given print writer.
	 *
	 * @param  out     Print writer.
	 * @param  jobNum  Job number.
	 */
	private void printJobNumberURL
		(PrintWriter out,
		 int jobNum)
		{
		out.printf ("http://%s:%d/job/%d", myWebHost, myWebPort, jobNum);
		}

	/**
	 * Print a link for the given job number on the given print writer.
	 *
	 * @param  out     Print writer.
	 * @param  jobNum  Job number.
	 */
	private void printJobNumberLink
		(PrintWriter out,
		 int jobNum)
		{
		out.printf ("<A HREF=\"http://%s:%d/job/%d\">&nbsp;%d&nbsp;</A>",
			myWebHost, myWebPort, jobNum, jobNum);
		}

	/**
	 * Print the difference between the given times on the given print writer.
	 *
	 * @param  out   Print writer.
	 * @param  now   Time now.
	 * @param  then  Time then.
	 */
	private void printDeltaTime
		(PrintWriter out,
		 long now,
		 long then)
		{
		out.print ((now - then + 500L) / 1000L);
		out.print (" sec");
		}

	/**
	 * Print the backend labels on the given print writer.
	 *
	 * @param  out      Print writer.
	 */
	private void printBackendLabels
		(PrintWriter out)
		{
		out.println ("<TR BGCOLOR=\"#E8E8E8\">");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Node</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>CPUs</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Status</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Job</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Time</I>");
		out.println ("</TD>");
		out.println ("</TR>");
		}

	/**
	 * Print the given backend info on the given print writer.
	 *
	 * @param  out      Print writer.
	 * @param  now      Current time.
	 * @param  backend  Backend info.
	 * @param  i        Even = white background, odd = gray background.
	 */
	private void printBackendInfo
		(PrintWriter out,
		 long now,
		 BackendInfo backend,
		 int i)
		{
		out.print   ("<TR BGCOLOR=\"#");
		out.print   (i%2==0 ? "FFFFFF" : "E8E8E8");
		out.println ("\">");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (backend.name);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (backend.totalCpus);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (backend.state == BackendInfo.State.FAILED)
			{
			out.print ("<FONT COLOR=\"#FF0000\"><B>");
			out.print (backend.state);
			out.print ("</B></FONT>");
			}
		else
			{
			out.print (backend.state);
			}
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (backend.job != null)
			{
			printJobNumberLink (out, backend.job.jobnum);
			}
		else
			{
			out.print ("&nbsp;");
			}
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (backend.job != null)
			{
			printDeltaTime (out, now, backend.job.stateTime);
			}
		else
			{
			out.print ("&nbsp;");
			}
		out.println ("</TD>");
		out.println ("</TR>");
		}

	/**
	 * Print the job labels on the given print writer.
	 *
	 * @param  out      Print writer.
	 */
	private void printJobLabels
		(PrintWriter out)
		{
		out.println ("<TR BGCOLOR=\"#E8E8E8\">");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Job</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>User</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>nn</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>np</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>nt</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Rank</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Node</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>CPUs</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Status</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Time</I>");
		out.println ("</TD>");
		out.println ("</TR>");
		}

	/**
	 * Print the given job info on the given print writer.
	 *
	 * @param  out  Print writer.
	 * @param  now  Current time.
	 * @param  job  Job info.
	 * @param  i    Even = white background, odd = gray background.
	 */
	private void printJobInfo
		(PrintWriter out,
		 long now,
		 JobInfo job,
		 int i)
		{
		boolean first;
		out.print   ("<TR BGCOLOR=\"#");
		out.print   (i%2==0 ? "FFFFFF" : "E8E8E8");
		out.println ("\">");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		printJobNumberLink (out, job.jobnum);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (job.username);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (job.Nn);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (job.Np);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (job.Nt == 0 ? "all" : ""+job.Nt);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (job.count == 0)
			{
			out.print ("&nbsp;");
			}
		else
			{
			for (int j = 0; j < job.count; ++ j)
				{
				if (j > 0) out.print ("<BR>");
				out.print (j);
				}
			}
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (job.count == 0)
			{
			out.print ("&nbsp;");
			}
		else
			{
			for (int j = 0; j < job.count; ++ j)
				{
				if (j > 0) out.print ("<BR>");
				out.print (job.backend[j].name);
				}
			}
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (job.count == 0)
			{
			out.print ("&nbsp;");
			}
		else
			{
			for (int j = 0; j < job.count; ++ j)
				{
				if (j > 0) out.print ("<BR>");
				out.print (job.cpus[j]);
				}
			}
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   (job.state);
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		printDeltaTime (out, now, job.stateTime);
		out.println ("</TD>");
		out.println ("</TR>");
		}

	/**
	 * Print the start of the debug HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printDebugHtmlStart
		(PrintWriter out,
		 long now)
		{
		out.println ("<HTML>");
		out.println ("<HEAD>");
		out.print   ("<TITLE>");
		out.print   (myClusterName);
		out.println ("</TITLE>");
		out.println ("<STYLE TYPE=\"text/css\">");
		out.println ("<!--");
		out.println ("* {font-family: Arial, Helvetica, Sans-Serif;}");
		out.println ("body {font-size: small;}");
		out.println ("h1 {font-size: 140%; font-weight: bold;}");
		out.println ("table {font-size: 100%;}");
		out.println ("-->");
		out.println ("</STYLE>");
		out.println ("</HEAD>");
		out.println ("<BODY>");
		out.print   ("<H1>");
		out.print   (myClusterName);
		out.println ("</H1>");
		out.println ("<P>");
		out.print   (new Date (now));
		out.print   (" -- ");
		out.print   (Version.PJ_VERSION);
		out.println ("</P>");
		}

	/**
	 * Print the body of the debug HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printDebugHtmlBody
		(PrintWriter out)
		{
		out.println ("<P>");
		out.println ("<HR/>");
		out.println ("<H3>Thread Dump</H3>");
		out.println ("</P>");
		Map<Thread,StackTraceElement[]> traces = Thread.getAllStackTraces();
		for (Map.Entry<Thread,StackTraceElement[]> entry : traces.entrySet())
			{
			Thread thread = entry.getKey();
			out.println ("<P>");
			out.print   ("Name: ");
			out.print   (thread.getName());
			out.println ("&nbsp;&nbsp;&nbsp;&nbsp;");
			out.print   (" ID: ");
			out.print   (thread.getId());
			out.println ("&nbsp;&nbsp;&nbsp;&nbsp;");
			out.print   (" Daemon: ");
			out.print   (thread.isDaemon() ? "yes" : "no");
			out.println ("&nbsp;&nbsp;&nbsp;&nbsp;");
			out.print   (" State: ");
			out.print   (thread.getState());
			out.println ("&nbsp;&nbsp;&nbsp;&nbsp;");
			out.print   (" Priority: ");
			out.print   (thread.getPriority());
			out.println ("&nbsp;&nbsp;&nbsp;&nbsp;");
			out.print   (" Thread Group: ");
			out.print   (thread.getThreadGroup().getName());
			out.println ();
			for (StackTraceElement element : entry.getValue())
				{
				out.print   ("<BR/>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;");
				out.println (element);
				}
			out.println ("</P>");
			}
		out.println ("<P>");
		out.println ("<HR/>");
		out.println ("</P>");
		}

	/**
	 * Print the start of the detailed job status HTML document on the given
	 * print writer.
	 *
	 * @param  out     Print writer.
	 * @param  now     Current time.
	 * @param  jobNum  Job number.
	 */
	private void printJobDetailHtmlStart
		(PrintWriter out,
		 long now,
		 int jobNum)
		{
		out.println ("<HTML>");
		out.println ("<HEAD>");
		out.print   ("<TITLE>");
		out.print   (myClusterName);
		out.println ("</TITLE>");
		out.print   ("<META HTTP-EQUIV=\"refresh\" CONTENT=\"20;url=");
		printJobNumberURL (out, jobNum);
		out.println ("\">");
		out.println ("<STYLE TYPE=\"text/css\">");
		out.println ("<!--");
		out.println ("* {font-family: Arial, Helvetica, Sans-Serif;}");
		out.println ("body {font-size: small;}");
		out.println ("h1 {font-size: 140%; font-weight: bold;}");
		out.println ("table {font-size: 100%;}");
		out.println ("-->");
		out.println ("</STYLE>");
		out.println ("</HEAD>");
		out.println ("<BODY>");
		out.print   ("<H1>");
		out.print   (myClusterName);
		out.println ("</H1>");
		out.println ("<P>");
		out.print   ("<FORM ACTION=\"");
		printJobNumberURL (out, jobNum);
		out.println ("\" METHOD=\"get\">");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"center\">");
		out.print   ("<INPUT TYPE=\"submit\" VALUE=\"Refresh\">");
		out.println ("</TD>");
		out.println ("<TD WIDTH=20> </TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"center\">");
		out.print   (new Date (now));
		out.print   (" -- ");
		out.print   (Version.PJ_VERSION);
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</FORM>");
		}

	/**
	 * Print the body of the detailed job status HTML document on the given
	 * print writer.
	 *
	 * @param  out     Print writer.
	 * @param  now     Current time.
	 * @param  jobNum  Job number.
	 */
	private synchronized void printJobDetailHtmlBody
		(PrintWriter out,
		 long now,
		 int jobNum)
		{
		JobInfo jobInfo = null;

		// Find job info.
		for (JobInfo job : myRunningJobList)
			{
			if (job.jobnum == jobNum)
				{
				jobInfo = job;
				break;
				}
			}
		if (jobInfo == null)
			{
			for (JobInfo job : myWaitingJobList)
				{
				if (job.jobnum == jobNum)
					{
					jobInfo = job;
					break;
					}
				}
			}

		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\"><B>Job:</B></TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.printf  ("<TD ALIGN=\"left\" VALIGN=\"top\"><B>%d</B></TD>",
			jobNum);
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">User:</TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.printf  ("<TD ALIGN=\"left\" VALIGN=\"top\">%s</TD>",
			jobInfo == null ? " " : jobInfo.username);
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Nodes (nn):</TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.printf  ("<TD ALIGN=\"left\" VALIGN=\"top\">%s</TD>",
			jobInfo == null ? " " : ""+jobInfo.Nn);
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Processes (np):</TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.printf  ("<TD ALIGN=\"left\" VALIGN=\"top\">%s</TD>",
			jobInfo == null ? " " : ""+jobInfo.Np);
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Threads (nt):</TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.printf  ("<TD ALIGN=\"left\" VALIGN=\"top\">%s</TD>",
			jobInfo == null ? " " : jobInfo.Nt == 0 ? "All" : ""+jobInfo.Nt);
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Status:</TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.printf  ("<TD ALIGN=\"left\" VALIGN=\"top\">%s</TD>",
			jobInfo == null ? "Not in queue" : jobInfo.state);
		out.println ("</TR>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">Time:</TD>");
		out.println ("<TD WIDTH=10> </TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		if (jobInfo == null)
			{
			out.print (" ");
			}
		else
			{
			printDeltaTime (out, now, jobInfo.stateTime);
			}
		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</P>");

		if (jobInfo == null || jobInfo.count == 0) return;

		out.println ("<P>");
		out.println ("<TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"center\" VALIGN=\"top\">");

		out.println ("Processes");
		out.println ("<TABLE BORDER=1 CELLPADDING=3 CELLSPACING=0>");
		out.println ("<TR>");
		out.println ("<TD ALIGN=\"left\" VALIGN=\"top\">");

		out.println ("<TABLE BORDER=0 CELLPADDING=3 CELLSPACING=0>");
		printJobDetailProcessLabels (out);
		for (int i = 0; i < jobInfo.count; ++ i)
			{
			printJobDetailProcessInfo (out, jobInfo, i);
			}
		out.println ("</TABLE>");

		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");

		out.println ("</TD>");
		out.println ("</TR>");
		out.println ("</TABLE>");
		out.println ("</P>");
		}

	/**
	 * Print the detailed job status process labels on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printJobDetailProcessLabels
		(PrintWriter out)
		{
		out.println ("<TR BGCOLOR=\"#E8E8E8\">");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Rank</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Node</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>CPUs</I>");
		out.println ("</TD>");
		out.print   ("<TD ALIGN=\"left\" VALIGN=\"top\">");
		out.print   ("<I>Comment</I>");
		out.println ("</TD>");
		out.println ("</TR>");
		}

	/**
	 * Print the detailed job status process information on the given print
	 * writer.
	 *
	 * @param  out      Print writer.
	 * @param  jobInfo  Job info.
	 * @param  rank     Process rank.
	 */
	private void printJobDetailProcessInfo
		(PrintWriter out,
		 JobInfo jobInfo,
		 int rank)
		{
		out.printf ("<TR BGCOLOR=\"#%s\">\n",
			rank%2 == 0 ? "FFFFFF" : "E8E8E8");
		out.printf ("<TD ALIGN=\"left\" VALIGN=\"top\">%d&nbsp;&nbsp;</TD>\n",
			rank);
		out.printf ("<TD ALIGN=\"left\" VALIGN=\"top\">%s&nbsp;&nbsp;</TD>\n",
			jobInfo.backend[rank].name);
		out.printf ("<TD ALIGN=\"left\" VALIGN=\"top\">%d&nbsp;&nbsp;</TD>\n",
			jobInfo.cpus[rank]);
		out.printf ("<TD ALIGN=\"left\" VALIGN=\"top\">%s</TD>\n",
			jobInfo.comment[rank]);
		out.println ("</TR>");
		}

	/**
	 * Print the start of the error HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printErrorHtmlStart
		(PrintWriter out)
		{
		out.println ("<HTML>");
		out.println ("<HEAD>");
		out.print   ("<TITLE>");
		out.print   (myClusterName);
		out.println ("</TITLE>");
		out.println ("<STYLE TYPE=\"text/css\">");
		out.println ("<!--");
		out.println ("* {font-family: Arial, Helvetica, Sans-Serif;}");
		out.println ("body {font-size: small;}");
		out.println ("h1 {font-size: 140%; font-weight: bold;}");
		out.println ("table {font-size: 100%;}");
		out.println ("-->");
		out.println ("</STYLE>");
		out.println ("</HEAD>");
		out.println ("<BODY>");
		}

	/**
	 * Print the end of the error HTML document on the given print writer.
	 *
	 * @param  out  Print writer.
	 */
	private void printErrorHtmlEnd
		(PrintWriter out)
		{
		out.println ("</BODY>");
		out.println ("</HTML>");
		}

	/**
	 * Shut down this Job Scheduler.
	 */
	private void shutdown()
		{
		if (myChannelGroup != null)
			{
			myChannelGroup.close();
			}
		if (myHttpServer != null)
			{
			try { myHttpServer.close(); } catch (IOException exc) {}
			}
		myLog.log ("Stopped");
		}

// Main program.

	/**
	 * Job Scheduler main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		if (args.length != 1)
			{
			System.err.println
				("Usage: java edu.rit.pj.cluster.JobScheduler <configfile>");
			System.exit (1);
			}

		JobScheduler scheduler = new JobScheduler (args[0]);
		scheduler.run();
		}

	}
