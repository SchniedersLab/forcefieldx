//******************************************************************************
//
// File:    ObjectArrayBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.ObjectArrayBuf
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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
import edu.rit.mp.ObjectBuf;

import edu.rit.pj.reduction.ObjectOp;
import edu.rit.pj.reduction.Op;

import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class ObjectArrayBuf provides a buffer for an array of object items sent or
 * received using the Message Protocol (MP). The array element stride may be 1
 * or greater than 1. While an instance of class ObjectArrayBuf may be
 * constructed directly, normally you will use a factory method in class
 * {@linkplain edu.rit.mp.ObjectBuf ObjectBuf}. See that class for further
 * information.
 *
 * @param <T> Data type of the objects in the buffer.
 *
 * @author Alan Kaminsky
 * @version 23-Mar-2009
 */
public class ObjectArrayBuf<T>
        extends ObjectBuf<T> {

// Hidden data members.
    T[] myArray;
    Range myRange;
    int myArrayOffset;
    int myStride;

// Exported constructors.
    /**
     * Construct a new object array buffer.
     *
     * @param theArray Array.
     * @param theRange Range of array elements to include in the buffer.
     */
    public ObjectArrayBuf(T[] theArray,
            Range theRange) {
        super(theRange.length());
        myArray = theArray;
        myRange = theRange;
        myArrayOffset = theRange.lb();
        myStride = theRange.stride();
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
    public T get(int i) {
        return myArray[myArrayOffset + i * myStride];
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
            T item) {
        myArray[myArrayOffset + i * myStride] = item;
        reset();
    }

    /**
     * Copy items from the given buffer to this buffer. The number of items
     * copied is this buffer's length or <TT>theSrc</TT>'s length, whichever is
     * smaller. If <TT>theSrc</TT> is this buffer, the <TT>copy()</TT> method
     * does nothing.
     *
     * @param theSrc Source of items to copy into this buffer.
     *
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>theSrc</TT>'s item data type is not the same as this buffer's item
     * data type.
     */
    public void copy(Buf theSrc) {
        if (theSrc == this) {
        } else if (theSrc instanceof ObjectArrayBuf) {
            ObjectArrayBuf<T> src = (ObjectArrayBuf<T>) theSrc;
            Arrays.copy(src.myArray, src.myRange, this.myArray, this.myRange);
            reset();
        } else {
            ObjectBuf.defaultCopy((ObjectBuf<T>) theSrc, this);
            reset();
        }
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
        return new ObjectArrayReductionBuf<T>(myArray, myRange, (ObjectOp<T>) op, this);
    }

}
