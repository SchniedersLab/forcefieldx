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

import java.io.File;
import java.io.IOException;
import java.net.InetSocketAddress;

import edu.rit.mp.ByteBuf;
import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;
import edu.rit.util.Range;

/**
 * Class JobFrontendProxy provides a proxy object for sending messages to a PJ
 * job frontend process.
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2012
 */
public class JobFrontendProxy
        extends Proxy
        implements JobFrontendRef {

// Exported constructors.
    /**
     * Construct a new job frontend proxy. The proxy will use the given channel
     * in the given channel group to send messages to the job frontend process.
     *
     * @param theChannelGroup Channel group.
     * @param theChannelGroup Channel group.
     * @param theChannel Channel.
     */
    public JobFrontendProxy(ChannelGroup theChannelGroup,
            Channel theChannel) {
        super(theChannelGroup, theChannel);
    }

// Exported operations.
    /**
     * {@inheritDoc}
     *
     * Assign a backend process to the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void assignBackend(JobSchedulerRef theJobScheduler,
            String name,
            String host,
            String jvm,
            String classpath,
            String[] jvmflags,
            String shellCommand,
            int Nt)
            throws IOException {
        send(JobFrontendMessage.assignBackend(theJobScheduler, name, host, jvm, classpath, jvmflags,
                shellCommand, Nt));
    }

    /**
     * {@inheritDoc}
     *
     * Assign a job number to the job. The host name for the job frontend's
     * middleware channel group is also specified.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void assignJobNumber(JobSchedulerRef theJobScheduler,
            int jobnum,
            String pjhost)
            throws IOException {
        send(JobFrontendMessage.assignJobNumber(theJobScheduler, jobnum, pjhost));
    }

    /**
     * {@inheritDoc}
     *
     * Cancel the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void cancelJob(JobSchedulerRef theJobScheduler,
            String errmsg)
            throws IOException {
        send(JobFrontendMessage.cancelJob(theJobScheduler, errmsg));
    }

    /**
     * Renew the lease on the job.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void renewLease(JobSchedulerRef theJobScheduler)
            throws IOException {
        send(JobFrontendMessage.renewLease(theJobScheduler));
    }

    /**
     * {@inheritDoc}
     *
     * Report that a backend process has finished executing the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void backendFinished(JobBackendRef theJobBackend)
            throws IOException {
        send(JobFrontendMessage.backendFinished(theJobBackend));
    }

    /**
     * {@inheritDoc}
     *
     * Report that a backend process is ready to commence executing the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void backendReady(JobBackendRef theJobBackend,
            int rank,
            InetSocketAddress middlewareAddress,
            InetSocketAddress worldAddress,
            InetSocketAddress frontendAddress)
            throws IOException {
        send(JobFrontendMessage.backendReady(theJobBackend, rank, middlewareAddress,
                worldAddress, frontendAddress));
    }

    /**
     * {@inheritDoc}
     *
     * Cancel the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void cancelJob(JobBackendRef theJobBackend,
            String errmsg)
            throws IOException {
        send(JobFrontendMessage.cancelJob(theJobBackend, errmsg));
    }

    /**
     * {@inheritDoc}
     *
     * Renew the lease on the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void renewLease(JobBackendRef theJobBackend)
            throws IOException {
        send(JobFrontendMessage.renewLease(theJobBackend));
    }

    /**
     * {@inheritDoc}
     *
     * Request the given resource from this job frontend's class loader.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void requestResource(JobBackendRef theJobBackend,
            String resourceName)
            throws IOException {
        send(JobFrontendMessage.requestResource(theJobBackend, resourceName));
    }

    /**
     * {@inheritDoc}
     *
     * Open the given output file for writing or appending.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file,
            boolean append)
            throws IOException {
        send(JobFrontendMessage.outputFileOpen(theJobBackend, bfd, file, append));
    }

    /**
     * {@inheritDoc}
     *
     * Write the given bytes to the given output file. <TT>fd</TT> = 1 refers to
     * the job's standard output stream; <TT>fd</TT> = 2 refers to the job's
     * standard error stream; other values refer to a previously opened file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileWrite(JobBackendRef theJobBackend,
            int ffd,
            byte[] buf,
            int off,
            int len)
            throws IOException {
        send(JobFrontendMessage.outputFileWrite(theJobBackend, ffd, len));
        send(ffd,
                ByteBuf.sliceBuffer(buf, new Range(off, off + len - 1)));
    }

    /**
     * {@inheritDoc}
     *
     * Flush accumulated bytes to the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileFlush(JobBackendRef theJobBackend,
            int ffd)
            throws IOException {
        send(JobFrontendMessage.outputFileFlush(theJobBackend, ffd));
    }

    /**
     * {@inheritDoc}
     *
     * Close the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileClose(JobBackendRef theJobBackend,
            int ffd)
            throws IOException {
        send(JobFrontendMessage.outputFileClose(theJobBackend, ffd));
    }

    /**
     * {@inheritDoc}
     *
     * Open the given input file for reading.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file)
            throws IOException {
        send(JobFrontendMessage.inputFileOpen(theJobBackend, bfd, file));
    }

    /**
     * {@inheritDoc}
     *
     * Read bytes from the given input file. <TT>ffd</TT> = 1 refers to the
     * job's standard input stream; other values refer to a previously opened
     * file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileRead(JobBackendRef theJobBackend,
            int ffd,
            int len)
            throws IOException {
        send(JobFrontendMessage.inputFileRead(theJobBackend, ffd, len));
    }

    /**
     * {@inheritDoc}
     *
     * Skip bytes from the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileSkip(JobBackendRef theJobBackend,
            int ffd,
            long len)
            throws IOException {
        send(JobFrontendMessage.inputFileSkip(theJobBackend, ffd, len));
    }

    /**
     * {@inheritDoc}
     *
     * Close the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileClose(JobBackendRef theJobBackend,
            int ffd)
            throws IOException {
        send(JobFrontendMessage.inputFileClose(theJobBackend, ffd));
    }

    /**
     * {@inheritDoc}
     *
     * Report a comment for a process.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void reportComment(JobBackendRef theJobBackend,
            int rank,
            String comment)
            throws IOException {
        send(JobFrontendMessage.reportComment(theJobBackend, rank, comment));
    }

}
