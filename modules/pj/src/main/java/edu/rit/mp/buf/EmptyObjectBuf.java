//******************************************************************************
//
// File:    EmptyObjectBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.EmptyObjectBuf
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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
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

/**
 * Class EmptyObjectBuf provides an object buffer that contains no items for
 * messages using the Message Protocol (MP). When a message is sent from an
 * EmptyObjectBuf, the message item type is Object and the message length is 0.
 * When a message is received into an EmptyObjectBuf, the message item type must
 * be Object, but all items in the message are discarded.
 *
 * @author Alan Kaminsky
 * @version 03-Jul-2008
 */
public class EmptyObjectBuf
        extends ObjectBuf<Object> {

// Exported constructors.
    /**
     * Construct a new empty object buffer.
     */
    public EmptyObjectBuf() {
        super(0);
    }

// Exported operations.
    /**
     * {@inheritDoc}
     *
     * Obtain the given item from this buffer.
     * <P>
     * The <code>get()</code> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     */
    public Object get(int i) {
        throw new IndexOutOfBoundsException();
    }

    /**
     * {@inheritDoc}
     *
     * Store the given item in this buffer.
     * <P>
     * The <code>put()</code> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     */
    public void put(int i,
            Object item) {
        throw new IndexOutOfBoundsException();
    }

    /**
     * {@inheritDoc}
     *
     * Copy items from the given buffer to this buffer. The number of items
     * copied is this buffer's length or <code>theSrc</code>'s length, whichever is
     * smaller. If <code>theSrc</code> is this buffer, the <code>copy()</code> method
     * does nothing.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <code>theSrc</code>'s item data type is not the same as this buffer's item
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
        ObjectOp objectop = (ObjectOp) op;
        return this;
    }

}
