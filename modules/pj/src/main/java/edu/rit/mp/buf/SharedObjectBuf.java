//******************************************************************************
//
// File:    SharedObjectBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.SharedObjectBuf
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
import edu.rit.pj.reduction.SharedObject;

/**
 * Class SharedObjectBuf provides a buffer for a shared object item sent or
 * received using the Message Protocol (MP). While an instance of class
 * SharedObjectBuf may be constructed directly, normally you will use a factory
 * method in class {@linkplain edu.rit.mp.ObjectBuf ObjectBuf}. See that class
 * for further information.
 *
 * @param <T> Data type of the objects in the buffer.
 * @author Alan Kaminsky
 * @version 01-Apr-2012
 */
public class SharedObjectBuf<T>
        extends ObjectBuf<T> {

// Hidden data members.
    SharedObject<T> myItem;

// Exported constructors.
    /**
     * Construct a new shared object buffer.
     *
     * @param item SharedObject object that wraps the item.
     */
    public SharedObjectBuf(SharedObject<T> item) {
        super(1);
        myItem = item;
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
    public T get(int i) {
        return myItem.get();
    }

    /**
     * Store the given item in this buffer.
     * <P>
     * The <TT>put()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <TT>length()</TT>-1.
     * @param item Item to be stored at index <TT>i</TT>.
     * @param item Item to be stored at index <TT>i</TT>.
     */
    public void put(int i,
            T item) {
        myItem.set(item);
        mySerializedItems = null;
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
        return new SharedObjectReductionBuf<T>(myItem, (ObjectOp<T>) op, this);
    }

}
