//******************************************************************************
//
// File:    SharedUnsigned16BitIntegerArrayBuf_1.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.SharedUnsigned16BitIntegerArrayBuf_1
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

import edu.rit.mp.Buf;
import edu.rit.mp.Unsigned16BitIntegerBuf;

import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.reduction.Op;
import edu.rit.pj.reduction.SharedIntegerArray;

import edu.rit.util.Range;

import java.nio.ByteBuffer;

/**
 * Class SharedUnsigned16BitIntegerArrayBuf_1 provides a buffer for a multiple
 * thread safe array of unsigned 16-bit integer items sent or received using the
 * Message Protocol (MP). The array element stride must be 1. While an instance
 * of class SharedUnsigned16BitIntegerArrayBuf_1 may be constructed directly,
 * normally you will use a factory method in class {@linkplain
 * edu.rit.mp.Unsigned16BitIntegerBuf Unsigned16BitIntegerBuf}. See that class
 * for further information.
 *
 * @author Alan Kaminsky
 * @version 26-Oct-2007
 */
public class SharedUnsigned16BitIntegerArrayBuf_1
        extends SharedUnsigned16BitIntegerArrayBuf {

// Exported constructors.
    /**
     * Construct a new shared unsigned 16-bit integer array buffer.
     *
     * @param theArray Shared array.
     * @param theRange Range of array elements to include in the buffer. The
     * stride is assumed to be 1.
     */
    public SharedUnsigned16BitIntegerArrayBuf_1(SharedIntegerArray theArray,
            Range theRange) {
        super(theArray, theRange);
    }

// Exported operations.
    /**
     * Obtain the given item from this buffer.
     * <P>
     * The <TT>get()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <TT>length()</TT>-1.
     *
     * @return Item at index <TT>i</TT>.
     */
    public int get(int i) {
        return myArray.get(myArrayOffset + i);
    }

    /**
     * Store the given item in this buffer.
     * <P>
     * The <TT>put()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <TT>length()</TT>-1.
     * @param item Item to be stored at index <TT>i</TT>.
     */
    public void put(int i,
            int item) {
        myArray.set(myArrayOffset + i, item);
    }

    /**
     * Create a buffer for performing parallel reduction using the given binary
     * operation. The results of the reduction are placed into this buffer.
     *
     * @param op Binary operation.
     *
     * @exception ClassCastException (unchecked exception) Thrown if this
     * buffer's element data type and the given binary operation's argument data
     * type are not the same.
     */
    public Buf getReductionBuf(Op op) {
        return new SharedUnsigned16BitIntegerArrayReductionBuf_1(myArray, myRange, (IntegerOp) op);
    }

// Hidden operations.
    /**
     * Send as many items as possible from this buffer to the given byte buffer.
     * <P>
     * The <TT>sendItems()</TT> method must not block the calling thread; if it
     * does, all message I/O in MP will be blocked.
     *
     * @param i Index of first item to send, in the range 0 ..
     * <TT>length</TT>-1.
     * @param buffer Byte buffer.
     *
     * @return Number of items sent.
     */
    protected int sendItems(int i,
            ByteBuffer buffer) {
        int index = i;
        int off = myArrayOffset + i;
        while (index < myLength && buffer.remaining() >= 2) {
            buffer.putShort((short) myArray.get(off));
            ++index;
            ++off;
        }
        return index - i;
    }

    /**
     * Receive as many items as possible from the given byte buffer to this
     * buffer.
     * <P>
     * The <TT>receiveItems()</TT> method must not block the calling thread; if
     * it does, all message I/O in MP will be blocked.
     *
     * @param i Index of first item to receive, in the range 0 ..
     * <TT>length</TT>-1.
     * @param num Maximum number of items to receive.
     * @param buffer Byte buffer.
     *
     * @return Number of items received.
     */
    protected int receiveItems(int i,
            int num,
            ByteBuffer buffer) {
        int index = i;
        int off = myArrayOffset + i;
        int max = Math.min(i + num, myLength);
        while (index < max && buffer.remaining() >= 2) {
            myArray.set(off, buffer.getShort() & 0xFFFF);
            ++index;
            ++off;
        }
        return index - i;
    }

}
