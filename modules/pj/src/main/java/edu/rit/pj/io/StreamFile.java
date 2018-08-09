//******************************************************************************
//
// File:    StreamFile.java
// Package: edu.rit.pj.io
// Unit:    Class edu.rit.pj.io.StreamFile
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
package edu.rit.pj.io;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import edu.rit.pj.cluster.JobBackend;

/**
 * Class StreamFile represents a file that resides in the user's account in the
 * job frontend process of a PJ cluster parallel program. Operations are
 * provided to open an input stream or an output stream to read or write the
 * file in the frontend processor.
 *
 * @author Alan Kaminsky
 * @version 06-Nov-2006
 */
public class StreamFile {

// Hidden data members.
    private File myFile;

// Exported constructors.
    /**
     * Construct a new stream file that refers to the given file in the frontend
     * processor.
     *
     * @param theFile File.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theFile</TT> is null.
     */
    public StreamFile(File theFile) {
        if (theFile == null) {
            throw new NullPointerException("StreamFile(): theFile is null");
        }
        myFile = theFile;
    }

// Exported operations.
    /**
     * Obtain the file in the frontend processor to which this stream file
     * refers.
     *
     * @return File.
     */
    public File getFile() {
        return myFile;
    }

    /**
     * Open an output stream for writing this stream file. If the file does not
     * exist, it is created; if the file exists, it is overwritten.
     * <P>
     * When called from a job backend process in a cluster parallel program, the
     * returned output stream communicates with the job frontend process to
     * write the file in the frontend processor. Otherwise, the returned output
     * stream is a normal file output stream to write the file directly.
     * <P>
     * <I>Note:</I> The returned output stream does not do any buffering. Each
     * method call sends a message to and receives a message from the job
     * frontend. Consider layering a BufferedOutputStream on top of the returned
     * output stream.
     *
     * @return Output stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public OutputStream getOutputStream()
            throws IOException {
        return getOutputStream(false);
    }

    /**
     * Open an output stream for writing this stream file. If the file does not
     * exist, it is created; if the file exists and the <TT>append</TT> flag is
     * false, the file is overwritten; if the file exists and the
     * <TT>append</TT> flag is true, data is written after the end of the file.
     * <P>
     * When called from a job backend process in a cluster parallel program, the
     * returned output stream communicates with the job frontend process to
     * write the file in the frontend processor. Otherwise, the returned output
     * stream is a normal file output stream to write the file directly.
     * <P>
     * <I>Note:</I> The returned output stream does not do any buffering. Each
     * method call sends a message to and receives a message from the job
     * frontend. Consider layering a BufferedOutputStream on top of the returned
     * output stream.
     *
     * @param append True to append, false to overwrite.
     * @return Output stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public OutputStream getOutputStream(boolean append)
            throws IOException {
        JobBackend backend = JobBackend.getJobBackend();
        if (backend == null) {
            return new FileOutputStream(myFile, append);
        } else {
            return backend.getFileWriter().open(myFile, append);
        }
    }

    /**
     * Open an input stream for reading this stream file.
     * <P>
     * When called from a job backend process in a cluster parallel program, the
     * returned input stream communicates with the job frontend process to read
     * the file in the frontend processor. Otherwise, the returned input stream
     * is a normal file input stream to read the file directly.
     * <P>
     * <I>Note:</I> The returned input stream does not do any buffering. Each
     * method call sends a message to and receives a message from the job
     * frontend. Consider layering a BufferedInputStream on top of the returned
     * input stream.
     *
     * @return Input stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public InputStream getInputStream()
            throws IOException {
        JobBackend backend = JobBackend.getJobBackend();
        if (backend == null) {
            return new FileInputStream(myFile);
        } else {
            return backend.getFileReader().open(myFile);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Determine if this stream file is equal to the given object.
     */
    public boolean equals(Object obj) {
        return (obj instanceof StreamFile)
                && (this.myFile.equals(((StreamFile) obj).myFile));
    }

    /**
     * Returns a hash code for this stream file.
     *
     * @return Hash code.
     */
    public int hashCode() {
        return myFile.hashCode();
    }

    /**
     * Returns a string version of this stream file.
     *
     * @return String version.
     */
    public String toString() {
        return myFile.toString();
    }

}
