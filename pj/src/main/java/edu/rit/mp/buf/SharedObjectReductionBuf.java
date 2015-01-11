//******************************************************************************
//
// File:    SharedObjectReductionBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.SharedObjectReductionBuf
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
 * Class SharedObjectReductionBuf provides a reduction buffer for class
 * {@linkplain SharedObjectBuf}.
 *
 * @param <T> Data type of the objects in the buffer.
 *
 * @author Alan Kaminsky
 * @version 01-Apr-2012
 */
class SharedObjectReductionBuf<T>
        extends SharedObjectBuf<T> {

// Hidden data members.
    ObjectOp<T> myOp;
    SharedObjectBuf<T> myBuf;

// Exported constructors.
    /**
     * Construct a new shared object reduction buffer.
     *
     * @param item SharedObject object that wraps the item.
     * @param op Binary operation.
     * @param theBuf Underlying object buffer.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>op</TT> is null.
     */
    public SharedObjectReductionBuf(SharedObject<T> item,
            ObjectOp<T> op,
            SharedObjectBuf<T> theBuf) {
        super(item);
        if (op == null) {
            throw new NullPointerException("SharedObjectReductionBuf(): op is null");
        }
        myOp = op;
        myBuf = theBuf;
    }

// Exported operations.
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
        myItem.reduce(item, myOp);
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
