//******************************************************************************
//
// File:    ObjectItemReductionBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.ObjectItemReductionBuf
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

import edu.rit.mp.ObjectBuf;
import edu.rit.mp.Buf;

import edu.rit.pj.reduction.ObjectOp;
import edu.rit.pj.reduction.Op;

import java.nio.ByteBuffer;

/**
 * Class ObjectItemReductionBuf provides a reduction buffer for class
 * {@linkplain ObjectItemBuf}.
 *
 * @param <T> Data type of the objects in the buffer.
 *
 * @author Alan Kaminsky
 * @version 23-Mar-2009
 */
class ObjectItemReductionBuf<T>
        extends ObjectBuf<T> {

// Hidden data members.
    ObjectItemBuf<T> myBuf;
    ObjectOp<T> myOp;

// Exported constructors.
    /**
     * Construct a new object item reduction buffer.
     *
     * @param buf Buffer containing the item.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>op</TT> is null.
     */
    public ObjectItemReductionBuf(ObjectItemBuf<T> buf,
            ObjectOp<T> op) {
        super(1);
        if (op == null) {
            throw new NullPointerException("ObjectItemReductionBuf(): op is null");
        }
        myBuf = buf;
        myOp = op;
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
        return myBuf.item;
    }

    /**
     * Store the given item in this buffer. The item at index <TT>i</TT> in this
     * buffer is combined with the given <TT>item</TT> using the binary
     * operation.
     * <P>
     * The <TT>put()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <TT>length()</TT>-1.
     * @param item Item to be stored at index <TT>i</TT>.
     */
    public void put(int i,
            T item) {
        myBuf.item = myOp.op(myBuf.item, item);
        reset();
        myBuf.reset();
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
        throw new UnsupportedOperationException();
    }

}
