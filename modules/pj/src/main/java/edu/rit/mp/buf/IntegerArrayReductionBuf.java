//******************************************************************************
//
// File:    IntegerArrayReductionBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.IntegerArrayReductionBuf
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

import java.nio.ByteBuffer;

import edu.rit.mp.Buf;
import edu.rit.mp.IntegerBuf;
import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.reduction.Op;
import edu.rit.pj.reduction.ReduceArrays;
import edu.rit.util.Range;

/**
 * Class IntegerArrayReductionBuf provides a reduction buffer for class
 * {@linkplain IntegerArrayBuf}.
 *
 * @author Alan Kaminsky
 * @version 25-Oct-2007
 */
class IntegerArrayReductionBuf
        extends IntegerArrayBuf {

// Hidden data members.
    IntegerOp myOp;

// Exported constructors.
    /**
     * Construct a new integer array reduction buffer.
     *
     * @param theArray Array.
     * @param theRange Range of array elements to include in the buffer.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>op</code> is null.
     */
    public IntegerArrayReductionBuf(int[] theArray,
            Range theRange,
            IntegerOp op) {
        super(theArray, theRange);
        if (op == null) {
            throw new NullPointerException("IntegerArrayReductionBuf(): op is null");
        }
        myOp = op;
    }

// Exported operations.
    /**
     * {@inheritDoc}
     *
     * Store the given item in this buffer.
     * <P>
     * The <code>put()</code> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     */
    public void put(int i,
            int item) {
        int off = myArrayOffset + i * myStride;
        myArray[off] = myOp.op(myArray[off], item);
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
        if (theSrc == this) {
        } else if (theSrc instanceof IntegerArrayBuf) {
            IntegerArrayBuf src = (IntegerArrayBuf) theSrc;
            ReduceArrays.reduce(src.myArray, src.myRange, this.myArray, this.myRange, myOp);
        } else {
            IntegerBuf.defaultCopy((IntegerBuf) theSrc, this);
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
        throw new UnsupportedOperationException();
    }

// Hidden operations.
    /**
     * {@inheritDoc}
     *
     * Receive as many items as possible from the given byte buffer to this
     * buffer.
     * <P>
     * The <code>receiveItems()</code> method must not block the calling thread; if
     * it does, all message I/O in MP will be blocked.
     */
    protected int receiveItems(int i,
            int num,
            ByteBuffer buffer) {
        int index = i;
        int off = myArrayOffset + i * myStride;
        int max = Math.min(i + num, myLength);
        while (index < max && buffer.remaining() >= 4) {
            myArray[off] = myOp.op(myArray[off], buffer.getInt());
            ++index;
            off += myStride;
        }
        return index - i;
    }

}
