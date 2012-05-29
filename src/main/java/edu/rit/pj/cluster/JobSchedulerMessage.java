//******************************************************************************
//
// File:    JobSchedulerMessage.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobSchedulerMessage
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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class JobSchedulerMessage provides a message sent to a Job Scheduler process
 * (interface {@linkplain JobSchedulerRef}) in the PJ cluster middleware.
 *
 * @author  Alan Kaminsky
 * @version 24-Jan-2012
 */
public abstract class JobSchedulerMessage
	extends Message
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = -7379945472003527741L;

// Exported constructors.

	/**
	 * Construct a new job scheduler message.
	 */
	public JobSchedulerMessage()
		{
		}

	/**
	 * Construct a new job scheduler message with the given message tag.
	 *
	 * @param  theTag  Message tag to use when sending this message.
	 */
	public JobSchedulerMessage
		(int theTag)
		{
		super (theTag);
		}

// Exported operations.

	/**
	 * Construct a new "backend failed" message.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  name            Backend node name.
	 *
	 * @return  "Backend failed" message.
	 */
	public static JobSchedulerMessage backendFailed
		(JobFrontendRef theJobFrontend,
		 String name)
		{
		return new BackendFailedMessage (theJobFrontend, name);
		}

	/**
	 * Construct a new "cancel job" message.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  errmsg          Error message string.
	 *
	 * @return  "Cancel job" message.
	 */
	public static JobSchedulerMessage cancelJob
		(JobFrontendRef theJobFrontend,
		 String errmsg)
		{
		return new CancelJobMessage (theJobFrontend, errmsg);
		}

	/**
	 * Construct a new "job finished" message.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @return  "Job finished" message.
	 */
	public static JobSchedulerMessage jobFinished
		(JobFrontendRef theJobFrontend)
		{
		return new JobFinishedMessage (theJobFrontend);
		}

	/**
	 * Construct a new "renew lease" message.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 *
	 * @return  "Renew lease" message.
	 */
	public static JobSchedulerMessage renewLease
		(JobFrontendRef theJobFrontend)
		{
		return new RenewLeaseMessage (theJobFrontend);
		}

	/**
	 * Construct a new "report comment" message.
	 *
	 * @param  theJobFrontend  Job frontend that is calling this method.
	 * @param  rank            Process rank.
	 * @param  comment         Comment string.
	 *
	 * @return  "Report comment" message.
	 */
	public static JobSchedulerMessage reportComment
		(JobFrontendRef theJobFrontend,
		 int rank,
		 String comment)
		{
		return new ReportCommentMessage (theJobFrontend, rank, comment);
		}

	/**
	 * Construct a new "request job" message.
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
	public static JobSchedulerMessage requestJob
		(JobFrontendRef theJobFrontend,
		 String username,
		 int Nn,
		 int Np,
		 int Nt)
		{
		return new RequestJobMessage (theJobFrontend, username, Nn, Np, Nt);
		}

	/**
	 * Invoke the method corresponding to this job scheduler message on the
	 * given Job Scheduler object. The method arguments come from the fields of
	 * this job scheduler message object.
	 *
	 * @param  theJobScheduler  Job Scheduler on which to invoke the method.
	 * @param  theJobFrontend   Job Frontend that is calling the method.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void invoke
		(JobSchedulerRef theJobScheduler,
		 JobFrontendRef theJobFrontend)
		throws IOException
		{
		throw new UnsupportedOperationException();
		}

	/**
	 * Write this job scheduler message to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		}

	/**
	 * Read this job scheduler message from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		}

// Hidden subclasses.

	/**
	 * Class BackendFailedMessage provides the Job Scheduler "backend failed"
	 * message in the PJ cluster middleware.
	 *
	 * @author  Alan Kaminsky
	 * @version 13-Oct-2006
	 */
	private static class BackendFailedMessage
		extends JobSchedulerMessage
		{
		private static final long serialVersionUID = 6495614788809259018L;

		private String name;

		public BackendFailedMessage()
			{
			}

		public BackendFailedMessage
			(JobFrontendRef theJobFrontend,
			 String name)
			{
			super (Message.FROM_JOB_FRONTEND);
			this.name = name;
			}

		public void invoke
			(JobSchedulerRef theJobScheduler,
			 JobFrontendRef theJobFrontend)
			throws IOException
			{
			theJobScheduler.backendFailed (theJobFrontend, name);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeUTF (name);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException
			{
			name = in.readUTF();
			}
		}

	/**
	 * Class CancelJobMessage provides the Job Scheduler "cancel job" message in
	 * the PJ cluster middleware.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Oct-2006
	 */
	private static class CancelJobMessage
		extends JobSchedulerMessage
		{
		private static final long serialVersionUID = 2902818757044365344L;

		private String errmsg;

		public CancelJobMessage()
			{
			}

		public CancelJobMessage
			(JobFrontendRef theJobFrontend,
			 String errmsg)
			{
			super (Message.FROM_JOB_FRONTEND);
			this.errmsg = errmsg;
			}

		public void invoke
			(JobSchedulerRef theJobScheduler,
			 JobFrontendRef theJobFrontend)
			throws IOException
			{
			theJobScheduler.cancelJob (theJobFrontend, errmsg);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeUTF (errmsg);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException
			{
			errmsg = in.readUTF();
			}
		}

	/**
	 * Class JobFinishedMessage provides the Job Scheduler "job finished"
	 * message in the PJ cluster middleware.
	 *
	 * @author  Alan Kaminsky
	 * @version 13-Oct-2006
	 */
	private static class JobFinishedMessage
		extends JobSchedulerMessage
		{
		private static final long serialVersionUID = -1179228962545666153L;

		public JobFinishedMessage()
			{
			}

		public JobFinishedMessage
			(JobFrontendRef theJobFrontend)
			{
			super (Message.FROM_JOB_FRONTEND);
			}

		public void invoke
			(JobSchedulerRef theJobScheduler,
			 JobFrontendRef theJobFrontend)
			throws IOException
			{
			theJobScheduler.jobFinished (theJobFrontend);
			}
		}

	/**
	 * Class RenewLeaseMessage provides the Job Scheduler "renew lease" message
	 * in the PJ cluster middleware.
	 *
	 * @author  Alan Kaminsky
	 * @version 12-Oct-2006
	 */
	private static class RenewLeaseMessage
		extends JobSchedulerMessage
		{
		private static final long serialVersionUID = 8547605668292095227L;

		public RenewLeaseMessage()
			{
			}

		public RenewLeaseMessage
			(JobFrontendRef theJobFrontend)
			{
			super (Message.FROM_JOB_FRONTEND);
			}

		public void invoke
			(JobSchedulerRef theJobScheduler,
			 JobFrontendRef theJobFrontend)
			throws IOException
			{
			theJobScheduler.renewLease (theJobFrontend);
			}
		}

	/**
	 * Class ReportCommentMessage provides the Job Scheduler "report comment"
	 * message in the PJ cluster middleware.
	 *
	 * @author  Alan Kaminsky
	 * @version 24-Jan-2012
	 */
	private static class ReportCommentMessage
		extends JobSchedulerMessage
		{
		private static final long serialVersionUID = -7431990305653172900L;

		private int rank;
		private String comment;

		public ReportCommentMessage()
			{
			}

		public ReportCommentMessage
			(JobFrontendRef theJobFrontend,
			 int rank,
			 String comment)
			{
			super (Message.FROM_JOB_FRONTEND);
			this.rank = rank;
			this.comment = comment == null ? "" : comment;
			}

		public void invoke
			(JobSchedulerRef theJobScheduler,
			 JobFrontendRef theJobFrontend)
			throws IOException
			{
			theJobScheduler.reportComment (theJobFrontend, rank, comment);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeInt (rank);
			out.writeUTF (comment);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException
			{
			rank = in.readInt();
			comment = in.readUTF();
			}
		}

	/**
	 * Class RequestJobMessage provides the Job Scheduler "request job" message
	 * in the PJ cluster middleware.
	 *
	 * @author  Alan Kaminsky
	 * @version 21-May-2008
	 */
	private static class RequestJobMessage
		extends JobSchedulerMessage
		{
		private static final long serialVersionUID = -6712799261136980645L;

		private String username;
		private int Nn;
		private int Np;
		private int Nt;

		public RequestJobMessage()
			{
			}

		public RequestJobMessage
			(JobFrontendRef theJobFrontend,
			 String username,
			 int Nn,
			 int Np,
			 int Nt)
			{
			super (Message.FROM_JOB_FRONTEND);
			this.username = username;
			this.Nn = Nn;
			this.Np = Np;
			this.Nt = Nt;
			}

		public void invoke
			(JobSchedulerRef theJobScheduler,
			 JobFrontendRef theJobFrontend)
			throws IOException
			{
			theJobScheduler.requestJob (theJobFrontend, username, Nn, Np, Nt);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeUTF (username);
			out.writeInt (Nn);
			out.writeInt (Np);
			out.writeInt (Nt);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException
			{
			username = in.readUTF();
			Nn = in.readInt();
			Np = in.readInt();
			Nt = in.readInt();
			}
		}

	}
