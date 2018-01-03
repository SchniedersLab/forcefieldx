//******************************************************************************
//
// File:    LongArrayBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.LongArrayBuf
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
import edu.rit.mp.LongBuf;
import edu.rit.pj.reduction.LongOp;
import edu.rit.pj.reduction.Op;
import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class LongArrayBuf provides a buffer for an array of long items sent or
 * received using the Message Protocol (MP). The array element stride may be 1
 * or greater than 1. While an instance of class LongArrayBuf may be constructed
 * directly, normally you will use a factory method in class
 * {@linkplain edu.rit.mp.LongBuf LongBuf}. See that class for further
 * information.
 *
 * @author Alan Kaminsky
 * @version 25-Oct-2007
 */
public class LongArrayBuf
        extends LongBuf {

// Hidden data members.
    long[] myArray;
    Range myRange;
    int myArrayOffset;
    int myStride;

// Exported constructors.
    /**
     * Construct a new long array buffer.
     *
     * @param theArray Array.
     * @param theRange Range of array elements to include in the buffer.
     */
    public LongArrayBuf(long[] theArray,
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
    public long get(int i) {
        return myArray[myArrayOffset + i * myStride];
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
            long item) {
        myArray[myArrayOffset + i * myStride] = item;
    }

    /**
     * {@inheritDoc}
     *
     * Copy items from the given buffer to this buffer. The number of items
     * copied is this buffer's length or <TT>theSrc</TT>'s length, whichever is
     * smaller. If <TT>theSrc</TT> is this buffer, the <TT>copy()</TT> method
     * does nothing.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>theSrc</TT>'s item data type is not the same as this buffer's item
     * data type.
     */
    public void copy(Buf theSrc) {
        if (theSrc == this) {
        } else if (theSrc instanceof LongArrayBuf) {
            LongArrayBuf src = (LongArrayBuf) theSrc;
            Arrays.copy(src.myArray, src.myRange, this.myArray, this.myRange);
        } else {
            LongBuf.defaultCopy((LongBuf) theSrc, this);
        }
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
        return new LongArrayReductionBuf(myArray, myRange, (LongOp) op);
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
        while (index < myLength && buffer.remaining() >= 8) {
            buffer.putLong(myArray[off]);
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
        while (index < max && buffer.remaining() >= 8) {
            myArray[off] = buffer.getLong();
            ++index;
            off += myStride;
        }
        return index - i;
    }

}
