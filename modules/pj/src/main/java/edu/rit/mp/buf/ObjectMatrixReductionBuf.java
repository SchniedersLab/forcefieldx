//******************************************************************************
//
// File:    ObjectMatrixReductionBuf.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.ObjectMatrixReductionBuf
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
import edu.rit.pj.reduction.ReduceArrays;
import edu.rit.util.Range;

/**
 * Class ObjectMatrixReductionBuf provides a reduction buffer for class
 * {@linkplain ObjectMatrixBuf}.
 *
 * @param <T> Data type of the objects in the buffer.
 *
 * @author Alan Kaminsky
 * @version 05-Apr-2009
 */
class ObjectMatrixReductionBuf<T>
        extends ObjectMatrixBuf<T> {

// Hidden data members.
    ObjectOp<T> myOp;
    ObjectMatrixBuf<T> myBuf;

// Exported constructors.
    /**
     * Construct a new object matrix reduction buffer. It is assumed that the
     * rows and columns of <TT>theMatrix</TT> are allocated and that each row of
     * <TT>theMatrix</TT> has the same number of columns.
     *
     * @param theMatrix Matrix.
     * @param theRowRange Range of rows to include.
     * @param theColRange Range of columns to include.
     * @param op Binary operation.
     * @param theBuf Underlying object matrix buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>op</TT> is null.
     */
    public ObjectMatrixReductionBuf(T[][] theMatrix,
            Range theRowRange,
            Range theColRange,
            ObjectOp<T> op,
            ObjectMatrixBuf<T> theBuf) {
        super(theMatrix, theRowRange, theColRange);
        if (op == null) {
            throw new NullPointerException("ObjectMatrixReductionBuf(): op is null");
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
     * @param item Item to be stored at index <TT>i</TT>.
     */
    public void put(int i,
            T item) {
        int row = i2r(i) * myRowStride + myLowerRow;
        int col = i2c(i) * myColStride + myLowerCol;
        myMatrix[row][col] = myOp.op(myMatrix[row][col], item);
        reset();
        myBuf.reset();
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
        } else if (theSrc instanceof ObjectMatrixBuf) {
            ObjectMatrixBuf<T> src = (ObjectMatrixBuf<T>) theSrc;
            ReduceArrays.reduce(src.myMatrix, src.myRowRange, src.myColRange,
                    this.myMatrix, this.myRowRange, this.myColRange,
                    myOp);
            reset();
            myBuf.reset();
        } else {
            ObjectBuf.defaultCopy((ObjectBuf<T>) theSrc, this);
            reset();
            myBuf.reset();
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

}
