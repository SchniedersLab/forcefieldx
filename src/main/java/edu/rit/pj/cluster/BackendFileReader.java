//******************************************************************************
//
// File:    BackendFileReader.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.BackendFileReader
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

import java.util.HashMap;
import java.util.Map;

/**
 * Class BackendFileReader provides an object that reads sequential files in the
 * job backend process.
 * <P>
 * <I>Note:</I> Class BackendFileReader is not multiple thread safe; it assumes
 * it is being called by a synchronized method in the job backend.
 *
 * @author Alan Kaminsky
 * @version 20-Nov-2006
 */
public class BackendFileReader {

// Exported data members.
    /**
     * Input stream for reading from the job frontend's standard input.
     */
    public final BackendFileInputStream in;

// Hidden data members.
    private JobFrontendRef myJobFrontend;
    private JobBackendRef myJobBackend;

    // Mapping from backend file descriptor to backend file input stream.
    private Map<Integer, BackendFileInputStream> myInputStreamForBFD
            = new HashMap<Integer, BackendFileInputStream>();

    // Mapping from frontend file descriptor to backend file input stream.
    private Map<Integer, BackendFileInputStream> myInputStreamForFFD
            = new HashMap<Integer, BackendFileInputStream>();

    // Next backend file descriptor.
    private int myNextBFD = 1;

// Exported constructors.
    /**
     * Construct a new backend file reader.
     *
     * @param theJobFrontend Job Frontend.
     * @param theJobBackend Job Backend.
     */
    public BackendFileReader(JobFrontendRef theJobFrontend,
            JobBackendRef theJobBackend) {
        myJobFrontend = theJobFrontend;
        myJobBackend = theJobBackend;

        // Set up input stream for stdin.
        in = new BackendFileInputStream(myJobFrontend, myJobBackend, 1);
        myInputStreamForFFD.put(1, in);
    }

// Exported operations.
    /**
     * Open a backend file input stream on the given file.
     *
     * @param file File.
     *
     * @return Backend file input stream.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public BackendFileInputStream open(File file)
            throws IOException {
        BackendFileInputStream stream = null;
        int bfd = 0;
        int ffd = 0;

        synchronized (this) {
            stream = new BackendFileInputStream(myJobFrontend, myJobBackend);
            bfd = myNextBFD++;
            myInputStreamForBFD.put(bfd, stream);
        }

        ffd = stream.open(bfd, file);

        synchronized (this) {
            myInputStreamForFFD.put(ffd, stream);
        }
        return stream;
    }

    /**
     * Report the result of opening the given input file.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param ffd Frontend file descriptor if success.
     * @param exc Null if success, exception if failure.
     */
    public void inputFileOpenResult(JobFrontendRef theJobFrontend,
            int bfd,
            int ffd,
            IOException exc) {
        BackendFileInputStream stream = null;
        synchronized (this) {
            stream = myInputStreamForBFD.remove(bfd);
        }
        if (stream != null) {
            stream.putResult(ffd, 0, 0L, exc);
        }
    }

    /**
     * Report the result of reading the given input file.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes read, or -1 if EOF.
     * @param exc Null if success, exception if failure.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileReadResult(JobFrontendRef theJobFrontend,
            int ffd,
            int len,
            IOException exc) {
        BackendFileInputStream stream = null;
        synchronized (this) {
            stream = myInputStreamForFFD.get(ffd);
        }
        if (stream != null) {
            stream.putResult(ffd, len, 0L, exc);
        }
    }

    /**
     * Report the result of skipping the given input file.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes skipped.
     * @param exc Null if success, exception if failure.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileSkipResult(JobFrontendRef theJobFrontend,
            int ffd,
            long len,
            IOException exc) {
        BackendFileInputStream stream = null;
        synchronized (this) {
            stream = myInputStreamForFFD.get(ffd);
        }
        if (stream != null) {
            stream.putResult(ffd, 0, len, exc);
        }
    }

    /**
     * Report the result of closing the given input file.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param exc Null if success, exception if failure.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void inputFileCloseResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc) {
        BackendFileInputStream stream = null;
        synchronized (this) {
            stream = myInputStreamForFFD.remove(ffd);
        }
        if (stream != null) {
            stream.putResult(ffd, 0, 0L, exc);
        }
    }

}
