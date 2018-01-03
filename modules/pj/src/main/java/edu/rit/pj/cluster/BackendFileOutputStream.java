//******************************************************************************
//
// File:    BackendFileOutputStream.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.BackendFileOutputStream
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj.cluster;

import java.io.File;
import java.io.IOException;
import java.io.InterruptedIOException;
import java.io.OutputStream;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Class BackendFileOutputStream provides an object in a job backend process
 * that writes a sequential file in the job frontend process. A backend file
 * output stream is not constructed directly, rather it is created by a factory
 * method in class {@linkplain BackendFileWriter}.
 * <P>
 * <I>Note:</I> Class BackendFileOutputStream does not do any buffering. Each
 * method call sends a message to and receives a message from the job frontend.
 * Consider layering a BufferedOutputStream on top of the
 * BackendFileOutputStream.
 *
 * @author Alan Kaminsky
 * @version 05-Nov-2006
 */
public class BackendFileOutputStream
        extends OutputStream {

// Hidden data members.
    private JobFrontendRef myJobFrontend;
    private JobBackendRef myJobBackend;

    // Queue of results from job frontend.
    private LinkedBlockingQueue<Result> myResultQueue
            = new LinkedBlockingQueue<Result>();

    private static class Result {

        public int ffd;
        public IOException exc;

        public Result(int ffd,
                IOException exc) {
            this.ffd = ffd;
            this.exc = exc;
        }
    }

    // Frontend file descriptor.
    private int ffd;

// Hidden constructors.
    /**
     * Construct a new backend file output stream. Call the <TT>open()</TT>
     * method to open the file and obtain the frontend file descriptor.
     *
     * @param theJobFrontend Job Frontend.
     * @param theJobBackend Job Backend.
     */
    BackendFileOutputStream(JobFrontendRef theJobFrontend,
            JobBackendRef theJobBackend) {
        this.myJobFrontend = theJobFrontend;
        this.myJobBackend = theJobBackend;
    }

    /**
     * Construct a new backend file output stream. Use the given frontend file
     * descriptor.
     *
     * @param theJobFrontend Job Frontend.
     * @param theJobBackend Job Backend.
     * @param ffd Frontend file descriptor.
     */
    BackendFileOutputStream(JobFrontendRef theJobFrontend,
            JobBackendRef theJobBackend,
            int ffd) {
        this.myJobFrontend = theJobFrontend;
        this.myJobBackend = theJobBackend;
        this.ffd = ffd;
    }

// Exported operations.
    /**
     * Write the given byte to this output stream. Only the least significant
     * eight bits are written.
     *
     * @param b Byte.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void write(int b)
            throws IOException {
        write(new byte[]{(byte) b});
    }

    /**
     * Write the given byte array to this output stream.
     *
     * @param buf Byte array.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void write(byte[] buf)
            throws IOException {
        write(buf, 0, buf.length);
    }

    /**
     * {@inheritDoc}
     *
     * Write a portion of the given byte array to this output stream.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>off</TT> &lt; 0, <TT>len</TT>
     * &lt; 0, or <TT>off+len</TT> &gt; <TT>buf.length</TT>.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void write(byte[] buf,
            int off,
            int len)
            throws IOException {
        if (off < 0 || len < 0 || off + len > buf.length) {
            throw new IndexOutOfBoundsException();
        }
        verifyOpen();
        myJobFrontend.outputFileWrite(myJobBackend, ffd, buf, off, len);
        getResult();
    }

    /**
     * Flush this output stream.
     *
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void flush()
            throws IOException {
        verifyOpen();
        myJobFrontend.outputFileFlush(myJobBackend, ffd);
        getResult();
    }

    /**
     * Close this output stream.
     *
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void close()
            throws IOException {
        verifyOpen();
        try {
            myJobFrontend.outputFileClose(myJobBackend, ffd);
            getResult();
        } finally {
            ffd = 0;
        }
    }

// Hidden operations.
    /**
     * Request the Job Frontend to open the file.
     *
     * @param bfd Backend file descriptor.
     * @param file File.
     * @param append True to append, false to overwrite.
     *
     * @return Frontend file descriptor.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    int open(int bfd,
            File file,
            boolean append)
            throws IOException {
        myJobFrontend.outputFileOpen(myJobBackend, bfd, file, append);
        this.ffd = getResult().ffd;
        return this.ffd;
    }

    /**
     * Get the next result from the result queue. Throw an IOException if
     * necessary.
     *
     * @return Result object.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private Result getResult()
            throws IOException {
        try {
            Result result = myResultQueue.take();
            if (result.exc != null) {
                throw result.exc;
            }
            return result;
        } catch (InterruptedException exc) {
            IOException exc2 = new InterruptedIOException("I/O interrupted");
            exc2.initCause(exc);
            throw exc2;
        }
    }

    /**
     * Put the given result into the result queue.
     *
     * @param ffd Frontend file descriptor.
     * @param exc Null if success, exception if failure.
     */
    void putResult(int ffd,
            IOException exc) {
        myResultQueue.offer(new Result(ffd, exc));
    }

    /**
     * Verify that this file is open.
     *
     * @exception IOException Thrown if this file is not open.
     */
    private void verifyOpen()
            throws IOException {
        if (ffd == 0) {
            throw new IOException("File closed");
        }
    }

}
