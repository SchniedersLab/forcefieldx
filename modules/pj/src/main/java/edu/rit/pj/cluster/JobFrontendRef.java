//******************************************************************************
//
// File:    JobFrontendRef.java
// Package: edu.rit.pj.cluster
// Unit:    Interface edu.rit.pj.cluster.JobFrontendRef
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

/**
 * Interface JobFrontendRef specifies the interface for the PJ job frontend
 * process.
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2012
 */
public interface JobFrontendRef {

// Exported operations.
    /**
     * Assign a backend process to the job.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @param name Backend node name.
     * @param host Host name for SSH remote login.
     * @param jvm Full pathname of Java Virtual Machine.
     * @param jvmflags Array of JVM command line flags.
     * @param classpath Java class path for PJ Library.
     * @param shellCommand Shell command string.
     * @param Nt Number of CPUs assigned to the process.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void assignBackend(JobSchedulerRef theJobScheduler,
            String name,
            String host,
            String jvm,
            String classpath,
            String[] jvmflags,
            String shellCommand,
            int Nt)
            throws IOException;

    /**
     * Assign a job number to the job. The host name for the job frontend's
     * middleware channel group is also specified.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @param jobnum Job number.
     * @param pjhost Host name for middleware channel group.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void assignJobNumber(JobSchedulerRef theJobScheduler,
            int jobnum,
            String pjhost)
            throws IOException;

    /**
     * Cancel the job.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @param errmsg Error message string.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void cancelJob(JobSchedulerRef theJobScheduler,
            String errmsg)
            throws IOException;

    /**
     * Renew the lease on the job.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void renewLease(JobSchedulerRef theJobScheduler)
            throws IOException;

    /**
     * Report that a backend process has finished executing the job.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void backendFinished(JobBackendRef theJobBackend)
            throws IOException;

    /**
     * Report that a backend process is ready to commence executing the job.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param rank Rank of the job backend process.
     * @param middlewareAddress Host/port to which the job backend process is
     * listening for middleware messages.
     * @param worldAddress Host/port to which the job backend process is
     * listening for the world communicator.
     * @param frontendAddress Host/port to which the job backend process is
     * listening for the frontend communicator, or null if the frontend
     * communicator does not exist.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void backendReady(JobBackendRef theJobBackend,
            int rank,
            InetSocketAddress middlewareAddress,
            InetSocketAddress worldAddress,
            InetSocketAddress frontendAddress)
            throws IOException;

    /**
     * Cancel the job.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param errmsg Error message string.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void cancelJob(JobBackendRef theJobBackend,
            String errmsg)
            throws IOException;

    /**
     * Renew the lease on the job.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void renewLease(JobBackendRef theJobBackend)
            throws IOException;

    /**
     * Request the given resource from this job frontend's class loader.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param resourceName Resource name.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void requestResource(JobBackendRef theJobBackend,
            String resourceName)
            throws IOException;

    /**
     * Open the given output file for writing or appending.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param file File.
     * @param append True to append, false to overwrite.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void outputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file,
            boolean append)
            throws IOException;

    /**
     * Write the given bytes to the given output file. <code>ffd</code> = 1 refers
     * to the job's standard output stream; <code>ffd</code> = 2 refers to the job's
     * standard error stream; other values refer to a previously opened file.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param buf Array of bytes to write.
     * @param off Index of first byte to write.
     * @param len Number of bytes to write.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void outputFileWrite(JobBackendRef theJobBackend,
            int ffd,
            byte[] buf,
            int off,
            int len)
            throws IOException;

    /**
     * Flush accumulated bytes to the given output file.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void outputFileFlush(JobBackendRef theJobBackend,
            int ffd)
            throws IOException;

    /**
     * Close the given output file.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void outputFileClose(JobBackendRef theJobBackend,
            int ffd)
            throws IOException;

    /**
     * Open the given input file for reading.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param file File.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void inputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file)
            throws IOException;

    /**
     * Read bytes from the given input file. <code>ffd</code> = 1 refers to the
     * job's standard input stream; other values refer to a previously opened
     * file.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes to read.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void inputFileRead(JobBackendRef theJobBackend,
            int ffd,
            int len)
            throws IOException;

    /**
     * Skip bytes from the given input file.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes to skip.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void inputFileSkip(JobBackendRef theJobBackend,
            int ffd,
            long len)
            throws IOException;

    /**
     * Close the given input file.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void inputFileClose(JobBackendRef theJobBackend,
            int ffd)
            throws IOException;

    /**
     * Report a comment for a process.
     *
     * @param theJobBackend Job backend that is calling this method.
     * @param rank Process rank.
     * @param comment Comment string.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void reportComment(JobBackendRef theJobBackend,
            int rank,
            String comment)
            throws IOException;

    /**
     * Close communication with this Job Frontend.
     */
    public void close();

}
