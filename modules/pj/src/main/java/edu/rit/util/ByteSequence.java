//******************************************************************************
//
// File:    ByteSequence.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.ByteSequence
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
package edu.rit.util;

import java.io.DataOutput;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Class ByteSequence provides an abstraction for a sequence of bytes. The
 * contents of the byte sequence are specified at construction time; the
 * contents may come from a byte array, an input stream, or another byte
 * sequence. You can obtain the byte sequence's contents as a byte array or
 * write the byte sequence's contents to an output stream.
 *
 * @author Alan Kaminsky
 * @version 02-Nov-2006
 */
public class ByteSequence {

// Hidden data members.
    private LinkedList<byte[]> myChunkList = new LinkedList<byte[]>();
    private LinkedList<Integer> myLengthList = new LinkedList<Integer>();
    private int myTotalLength;

// Exported constructors.
    /**
     * Construct a new byte sequence whose contents are a copy of the given byte
     * array.
     *
     * @param buf Byte array to copy.
     * @exception NullPointerException Thrown if <TT>buf</TT> is null.
     */
    public ByteSequence(byte[] buf) {
        this(buf, 0, buf.length);
    }

    /**
     * Construct a new byte sequence whose contents are a copy of a portion of
     * the given byte array.
     *
     * @param buf Byte array to copy.
     * @param off Index of first byte to copy.
     * @param len Number of bytes to copy.
     * @exception NullPointerException Thrown if <TT>buf</TT> is null.
     * @exception IndexOutOfBoundsException Thrown if <TT>off</TT> &lt; 0,
     * <TT>len</TT> &lt; 0, or
     * <TT>off+len</TT> &gt; <TT>buf.length</TT>.
     */
    public ByteSequence(byte[] buf,
            int off,
            int len) {
        if (off < 0 || len < 0 || off + len > buf.length) {
            throw new IndexOutOfBoundsException();
        }
        byte[] chunk = new byte[len];
        System.arraycopy(buf, off, chunk, 0, len);
        myChunkList.add(chunk);
        myLengthList.add(len);
        myTotalLength = len;
    }

    /**
     * Construct a new byte sequence whose contents come from the given input
     * stream. Bytes are read from <TT>theInputStream</TT> into the byte
     * sequence until the end-of-stream is encountered, then
     * <TT>theInputStream</TT> is closed. If <TT>theInputStream</TT> is null,
     * the byte sequence's length is 0.
     *
     * @param theInputStream Input stream, or null.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public ByteSequence(InputStream theInputStream)
            throws IOException {
        if (theInputStream != null) {
            for (;;) {
                byte[] chunk = new byte[4096];
                int length = theInputStream.read(chunk);
                if (length == -1) {
                    break;
                }
                myChunkList.add(chunk);
                myLengthList.add(length);
                myTotalLength += length;
                if (myTotalLength < 0) {
                    throw new IOException("ByteSequence(): Input stream too long");
                }
            }
            theInputStream.close();
        }
    }

    /**
     * Construct a new byte sequence whose contents are a copy of the given byte
     * sequence.
     *
     * @param theByteSequence Byte sequence to copy.
     * @exception NullPointerException Thrown if <TT>theByteSequence</TT> is
     * null.
     */
    public ByteSequence(ByteSequence theByteSequence) {
        byte[] chunk = theByteSequence.toByteArray();
        myChunkList.add(chunk);
        myLengthList.add(chunk.length);
        myTotalLength = chunk.length;
    }

// Exported operations.
    /**
     * Obtain the length of this byte sequence.
     *
     * @return Number of bytes.
     */
    public int length() {
        return myTotalLength;
    }

    /**
     * Obtain a byte array with a copy of this byte sequence's contents. A new
     * byte array of the proper size is created and returned.
     *
     * @return Contents.
     */
    public byte[] toByteArray() {
        byte[] result = new byte[myTotalLength];
        copy(result);
        return result;
    }

    /**
     * Copy this byte sequence's contents into the given byte array. Bytes are
     * copied into <TT>buf</TT> starting at index 0. The number of bytes copied
     * is <TT>buf.length</TT> or this byte sequence's length, whichever is
     * smaller.
     *
     * @param buf Buffer to hold the copy.
     * @return Actual number of bytes copied.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null.
     */
    public int copy(byte[] buf) {
        return copy(buf, 0, buf.length);
    }

    /**
     * Copy this byte sequence's contents into a portion of the given byte
     * array. Bytes are copied into <TT>buf</TT> starting at index <TT>off</TT>.
     * The number of bytes copied is <TT>len</TT> or this byte sequence's
     * length, whichever is smaller.
     *
     * @param buf Buffer to hold the copy.
     * @param off Index in <TT>buf</TT> at which to start copying.
     * @param len Maximum number of bytes to copy.
     * @return Actual number of bytes copied.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null.
     * @exception IndexOutOfBoundsException Thrown if <TT>off</TT> &lt; 0,
     * <TT>len</TT> &lt; 0, or
     * <TT>off+len</TT> &gt; <TT>buf.length</TT>.
     */
    public int copy(byte[] buf,
            int off,
            int len) {
        if (off < 0 || len < 0 || off + len > buf.length) {
            throw new IndexOutOfBoundsException();
        }
        int total = 0;
        Iterator<byte[]> chunkiter = myChunkList.iterator();
        Iterator<Integer> lengthiter = myLengthList.iterator();
        while (len > 0 && chunkiter.hasNext()) {
            byte[] chunk = chunkiter.next();
            int length = lengthiter.next();
            int n = Math.min(length, len);
            System.arraycopy(chunk, 0, buf, off, n);
            off += n;
            len -= n;
            total += n;
        }
        return total;
    }

    /**
     * Write this byte sequence's contents to the given output stream.
     *
     * @param theOutputStream Output stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void write(OutputStream theOutputStream)
            throws IOException {
        Iterator<byte[]> chunkiter = myChunkList.iterator();
        Iterator<Integer> lengthiter = myLengthList.iterator();
        while (chunkiter.hasNext()) {
            byte[] chunk = chunkiter.next();
            int length = lengthiter.next();
            theOutputStream.write(chunk, 0, length);
        }
    }

    /**
     * Write this byte sequence's contents to the given data output stream.
     *
     * @param theOutputStream Data output stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void write(DataOutput theOutputStream)
            throws IOException {
        Iterator<byte[]> chunkiter = myChunkList.iterator();
        Iterator<Integer> lengthiter = myLengthList.iterator();
        while (chunkiter.hasNext()) {
            byte[] chunk = chunkiter.next();
            int length = lengthiter.next();
            theOutputStream.write(chunk, 0, length);
        }
    }

}
