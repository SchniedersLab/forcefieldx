//******************************************************************************
//
// File:    DoubleItemReductionBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.DoubleItemReductionBuf
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
import edu.rit.mp.DoubleBuf;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.Op;

/**
 * Class DoubleItemReductionBuf provides a reduction buffer for class
 * {@linkplain DoubleItemBuf}.
 *
 * @author Alan Kaminsky
 * @version 25-Oct-2007
 */
class DoubleItemReductionBuf
        extends DoubleBuf {

// Hidden data members.
    DoubleItemBuf myBuf;
    DoubleOp myOp;

// Exported constructors.
    /**
     * Construct a new double item reduction buffer.
     *
     * @param buf Buffer containing the item.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>op</TT> is null.
     */
    public DoubleItemReductionBuf(DoubleItemBuf buf,
            DoubleOp op) {
        super(1);
        if (op == null) {
            throw new NullPointerException("DoubleItemReductionBuf(): op is null");
        }
        myBuf = buf;
        myOp = op;
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
    public double get(int i) {
        return myBuf.item;
    }

    /**
     * {@inheritDoc}
     *
     * Store the given item in this buffer. The item at index <TT>i</TT> in this
     * buffer is combined with the given <TT>item</TT> using the binary
     * operation.
     * <P>
     * The <TT>put()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     */
    public void put(int i,
            double item) {
        myBuf.item = myOp.op(myBuf.item, item);
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
        throw new UnsupportedOperationException();
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
        if (buffer.remaining() >= 8) {
            buffer.putDouble(myBuf.item);
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Receive as many items as possible from the given byte buffer to this
     * buffer. As the items are received, they are combined with the items in
     * this buffer using the binary operation.
     * <P>
     * The <TT>receiveItems()</TT> method must not block the calling thread; if
     * it does, all message I/O in MP will be blocked.
     */
    protected int receiveItems(int i,
            int num,
            ByteBuffer buffer) {
        if (num >= 1 && buffer.remaining() >= 8) {
            myBuf.item = myOp.op(myBuf.item, buffer.getDouble());
            return 1;
        } else {
            return 0;
        }
    }

}
