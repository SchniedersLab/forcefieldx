//******************************************************************************
//
// File:    SharedFloatArrayBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.SharedFloatArrayBuf
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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
package edu.rit.mp.buf;

import java.nio.ByteBuffer;

import edu.rit.mp.Buf;
import edu.rit.mp.FloatBuf;
import edu.rit.pj.reduction.FloatOp;
import edu.rit.pj.reduction.Op;
import edu.rit.pj.reduction.SharedFloatArray;
import edu.rit.util.Range;

/**
 * Class SharedFloatArrayBuf provides a buffer for a multiple thread safe array
 * of float items sent or received using the Message Protocol (MP). The array
 * element stride may be 1 or greater than 1. While an instance of class
 * SharedFloatArrayBuf may be constructed directly, normally you will use a
 * factory method in class {@linkplain edu.rit.mp.FloatBuf
 * FloatBuf}. See that class for further information.
 *
 * @author Alan Kaminsky
 * @version 26-Oct-2007
 */
public class SharedFloatArrayBuf
        extends FloatBuf {

// Hidden data members.
    SharedFloatArray myArray;
    Range myRange;
    int myArrayOffset;
    int myStride;

// Exported constructors.
    /**
     * Construct a new shared float array buffer.
     *
     * @param theArray Shared array.
     * @param theRange Range of array elements to include in the buffer.
     */
    public SharedFloatArrayBuf(SharedFloatArray theArray,
            Range theRange) {
        super(theRange.length());
        myArray = theArray;
        myRange = theRange;
        myArrayOffset = theRange.lb();
        myStride = theRange.stride();
    }

// Exported operations.
    /**
     * {@inheritDoc}
     *
     * Obtain the given item from this buffer.
     * <P>
     * The <TT>get()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     */
    public float get(int i) {
        return myArray.get(myArrayOffset + i * myStride);
    }

    /**
     * {@inheritDoc}
     *
     * Store the given item in this buffer.
     * <P>
     * The <TT>put()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     */
    public void put(int i,
            float item) {
        myArray.set(myArrayOffset + i * myStride, item);
    }

    /**
     * {@inheritDoc}
     *
     * Create a buffer for performing parallel reduction using the given binary
     * operation. The results of the reduction are placed into this buffer.
     * @exception ClassCastException (unchecked exception) Thrown if this
     * buffer's element data type and the given binary operation's argument data
     * type are not the same.
     */
    public Buf getReductionBuf(Op op) {
        return new SharedFloatArrayReductionBuf(myArray, myRange, (FloatOp) op);
    }

// Hidden operations.
    /**
     * {@inheritDoc}
     *
     * Send as many items as possible from this buffer to the given byte buffer.
     * <P>
     * The <TT>sendItems()</TT> method must not block the calling thread; if it
     * does, all message I/O in MP will be blocked.
     */
    protected int sendItems(int i,
            ByteBuffer buffer) {
        int index = i;
        int off = myArrayOffset + i * myStride;
        while (index < myLength && buffer.remaining() >= 4) {
            buffer.putFloat(myArray.get(off));
            ++index;
            off += myStride;
        }
        return index - i;
    }

    /**
     * {@inheritDoc}
     *
     * Receive as many items as possible from the given byte buffer to this
     * buffer.
     * <P>
     * The <TT>receiveItems()</TT> method must not block the calling thread; if
     * it does, all message I/O in MP will be blocked.
     */
    protected int receiveItems(int i,
            int num,
            ByteBuffer buffer) {
        int index = i;
        int off = myArrayOffset + i * myStride;
        int max = Math.min(i + num, myLength);
        while (index < max && buffer.remaining() >= 4) {
            myArray.set(off, buffer.getFloat());
            ++index;
            off += myStride;
        }
        return index - i;
    }

}
