//******************************************************************************
//
// File:    JobBackendProxy.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobBackendProxy
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
import java.net.InetSocketAddress;
import java.util.Properties;

import edu.rit.mp.ByteBuf;
import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;
import edu.rit.util.ByteSequence;
import edu.rit.util.Range;

/**
 * Class JobBackendProxy provides a proxy object for sending messages to a PJ
 * job backend process.
 *
 * @author Alan Kaminsky
 * @version 20-Nov-2006
 */
public class JobBackendProxy
        extends Proxy
        implements JobBackendRef {

// Exported constructors.
    /**
     * Construct a new job backend proxy. The proxy will use the given channel
     * in the given channel group to send messages to the job backend process.
     *
     * @param theChannelGroup Channel group.
     * @param theChannel Channel.
     */
    public JobBackendProxy(ChannelGroup theChannelGroup,
            Channel theChannel) {
        super(theChannelGroup, theChannel);
    }

// Exported operations.
    /**
     * {@inheritDoc}
     *
     * Cancel the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void cancelJob(JobFrontendRef theJobFrontend,
            String errmsg)
            throws IOException {
        send(JobBackendMessage.cancelJob(theJobFrontend, errmsg));
    }

    /**
     * Commence the job.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param middlewareAddress Array of hosts/ports for middleware messages.
     * The first <I>K</I>
     * elements are for the job backend processes in rank order, the
     * <I>K</I>+1st element is for the job frontend process. If the
     * @param worldAddress Array of hosts/ports for the world communicator. The
     * <I>K</I>
     * elements are for the job backend processes in rank order.
     * @param frontendAddress Array of hosts/ports for the frontend
     * communicator. The first
     * <I>K</I> elements are for the job backend processes in rank order, the
     * <I>K</I>+1st element is for the job frontend process. If the frontend
     * communicator does not exist, <code>frontendAddress</code> is null.
     * @param properties Java system properties.
     * @param mainClassName Fully qualified class name of the Java main program
     * class to execute.
     * @param args Array of 0 or more Java command line arguments.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void commenceJob(JobFrontendRef theJobFrontend,
            InetSocketAddress[] middlewareAddress,
            InetSocketAddress[] worldAddress,
            InetSocketAddress[] frontendAddress,
            Properties properties,
            String mainClassName,
            String[] args)
            throws IOException {
        send(JobBackendMessage.commenceJob(theJobFrontend, middlewareAddress, worldAddress,
                frontendAddress, properties, mainClassName, args));
    }

    /**
     * {@inheritDoc}
     *
     * Report that the job finished.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void jobFinished(JobFrontendRef theJobFrontend)
            throws IOException {
        send(JobBackendMessage.jobFinished(theJobFrontend));
    }

    /**
     * {@inheritDoc}
     *
     * Renew the lease on the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void renewLease(JobFrontendRef theJobFrontend)
            throws IOException {
        send(JobBackendMessage.renewLease(theJobFrontend));
    }

    /**
     * Report the content for a previously-requested resource.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param resourceName Resource name.
     * @param content Resource content, or null if resource not found.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void reportResource(JobFrontendRef theJobFrontend,
            String resourceName,
            byte[] content)
            throws IOException {
        send(JobBackendMessage.reportResource(theJobFrontend, resourceName, content));
    }

    /**
     * {@inheritDoc}
     *
     * Report the content for a previously-requested resource.
     *
     * @exception IOException Thrown if an I/O error occurred.
     * @param theJobFrontend a {@link edu.rit.pj.cluster.JobFrontendRef} object
     * @param resourceName a {@link java.lang.String} object
     * @param content a {@link edu.rit.util.ByteSequence} object
     * @throws java.io.IOException if any.
     */
    public void reportResource(JobFrontendRef theJobFrontend,
            String resourceName,
            ByteSequence content)
            throws IOException {
        send(JobBackendMessage.reportResource(theJobFrontend, resourceName, content));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of opening the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileOpenResult(JobFrontendRef theJobFrontend,
            int bfd,
            int ffd,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.outputFileOpenResult(theJobFrontend, bfd, ffd, exc));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of writing the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileWriteResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.outputFileWriteResult(theJobFrontend, ffd, exc));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of flushing the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileFlushResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.outputFileFlushResult(theJobFrontend, ffd, exc));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of closing the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void outputFileCloseResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.outputFileCloseResult(theJobFrontend, ffd, exc));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of opening the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileOpenResult(JobFrontendRef theJobFrontend,
            int bfd,
            int ffd,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.inputFileOpenResult(theJobFrontend, bfd, ffd, exc));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of reading the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileReadResult(JobFrontendRef theJobFrontend,
            int ffd,
            byte[] buf,
            int len,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.inputFileReadResult(theJobFrontend, ffd, len, exc));
        if (len > 0) {
            send(ffd, ByteBuf.sliceBuffer(buf, new Range(0, len - 1)));
        }
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of skipping the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileSkipResult(JobFrontendRef theJobFrontend,
            int ffd,
            long len,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.inputFileSkipResult(theJobFrontend, ffd, len, exc));
    }

    /**
     * {@inheritDoc}
     *
     * Report the result of closing the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileCloseResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc)
            throws IOException {
        send(JobBackendMessage.inputFileCloseResult(theJobFrontend, ffd, exc));
    }

}
