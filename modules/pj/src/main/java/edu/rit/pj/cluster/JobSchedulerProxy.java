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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
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
 * @author Alan Kaminsky
 * @version 24-Jan-2012
 */
public class JobSchedulerProxy
        extends Proxy
        implements JobSchedulerRef {

// Exported constructors.
    /**
     * Construct a new job scheduler proxy. The proxy will use the given channel
     * in the given channel group to send messages to the job scheduler process.
     *
     * @param theChannelGroup Channel group.
     * @param theChannel Channel.
     */
    public JobSchedulerProxy(ChannelGroup theChannelGroup,
            Channel theChannel) {
        super(theChannelGroup, theChannel);
    }

// Exported operations.
    /**
     * {@inheritDoc}
     *
     * Report that a backend node failed.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void backendFailed(JobFrontendRef theJobFrontend,
            String name)
            throws IOException {
        send(JobSchedulerMessage.backendFailed(theJobFrontend, name));
    }

    /**
     * {@inheritDoc}
     *
     * Cancel a job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void cancelJob(JobFrontendRef theJobFrontend,
            String errmsg)
            throws IOException {
        send(JobSchedulerMessage.cancelJob(theJobFrontend, errmsg));
    }

    /**
     * {@inheritDoc}
     *
     * Report that a job finished.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void jobFinished(JobFrontendRef theJobFrontend)
            throws IOException {
        send(JobSchedulerMessage.jobFinished(theJobFrontend));
    }

    /**
     * {@inheritDoc}
     *
     * Renew the lease on a job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void renewLease(JobFrontendRef theJobFrontend)
            throws IOException {
        send(JobSchedulerMessage.renewLease(theJobFrontend));
    }

    /**
     * {@inheritDoc}
     *
     * Report a comment for a process.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void reportComment(JobFrontendRef theJobFrontend,
            int rank,
            String comment)
            throws IOException {
        send(JobSchedulerMessage.reportComment(theJobFrontend, rank, comment));
    }

    /**
     * {@inheritDoc}
     *
     * Request that a job be scheduled.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void requestJob(JobFrontendRef theJobFrontend,
            String username,
            int Nn,
            int Np,
            int Nt)
            throws IOException {
        send(JobSchedulerMessage.requestJob(theJobFrontend, username, Nn, Np, Nt));
    }

}
