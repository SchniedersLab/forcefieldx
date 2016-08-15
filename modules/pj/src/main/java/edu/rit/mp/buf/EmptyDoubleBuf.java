//******************************************************************************
//
// File:    EmptyDoubleBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.EmptyDoubleBuf
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
import edu.rit.mp.DoubleBuf;

import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.Op;

import java.nio.ByteBuffer;

/**
 * Class EmptyDoubleBuf provides a double buffer that contains no items for
 * messages using the Message Protocol (MP). When a message is sent from an
 * EmptyDoubleBuf, the message item type is <TT>double</TT> and the message
 * length is 0. When a message is received into an EmptyDoubleBuf, the message
 * item type must be <TT>double</TT>, but all items in the message are
 * discarded.
 *
 * @author Alan Kaminsky
 * @version 19-Nov-2007
 */
public class EmptyDoubleBuf
        extends DoubleBuf {

// Exported constructors.
    /**
     * Construct a new empty double buffer.
     */
    public EmptyDoubleBuf() {
        super(0);
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
        throw new IndexOutOfBoundsException();
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
            double item) {
        throw new IndexOutOfBoundsException();
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
    }

    /**
     * {@inheritDoc}
     *
     * Create a buffer for performing parallel reduction using the given binary
     * operation. The results of the reduction are placed into this buffer.
     * <P>
     * Operations performed on the returned reduction buffer have the same
     * effect as operations performed on this buffer, except whenever a source
     * item <I>S</I> is put into a destination item <I>D</I> in this buffer,
     * <I>D</I> is set to <I>D op S</I>, that is, the reduction of <I>D</I> and
     * <I>S</I> using the given binary operation (rather than just setting
     * <I>D</I> to <I>S</I>).
     * @exception ClassCastException (unchecked exception) Thrown if this
     * buffer's element data type and the given binary operation's argument data
     * type are not the same.
     */
    public Buf getReductionBuf(Op op) {
        DoubleOp doubleop = (DoubleOp) op;
        return this;
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
        return 0;
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
        return 0;
    }

}
