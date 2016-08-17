//******************************************************************************
//
// File:    Signed16BitIntegerMatrixBuf_1.java
// Package: edu.rit.mp.buf
// Unit:    Class edu.rit.mp.buf.Signed16BitIntegerMatrixBuf_1
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
import edu.rit.mp.Signed16BitIntegerBuf;

import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.reduction.Op;

import edu.rit.util.Arrays;
import edu.rit.util.Range;

import java.nio.ByteBuffer;
import java.nio.ShortBuffer;

/**
 * Class Signed16BitIntegerMatrixBuf_1 provides a buffer for a matrix of signed
 * 16-bit integer items sent or received using the Message Protocol (MP). The
 * matrix row and column strides must both be 1. While an instance of class
 * Signed16BitIntegerMatrixBuf_1 may be constructed directly, normally you will
 * use a factory method in class {@linkplain edu.rit.mp.Signed16BitIntegerBuf
 * Signed16BitIntegerBuf}. See that class for further information.
 *
 * @author Alan Kaminsky
 * @version 05-Apr-2009
 */
public class Signed16BitIntegerMatrixBuf_1
        extends Signed16BitIntegerMatrixBuf {

// Exported constructors.
    /**
     * Construct a new signed 16-bit integer matrix buffer. It is assumed that
     * the rows and columns of <TT>theMatrix</TT> are allocated and that each
     * row of <TT>theMatrix</TT> has the same number of columns.
     *
     * @param theMatrix Matrix.
     * @param theRowRange Range of rows to include. The stride is assumed to be
     * 1.
     * @param theColRange Range of columns to include. The stride is assumed to
     * be 1.
     */
    public Signed16BitIntegerMatrixBuf_1(int[][] theMatrix,
            Range theRowRange,
            Range theColRange) {
        super(theMatrix, theRowRange, theColRange);
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
    public int get(int i) {
        return myMatrix[i2r(i) + myLowerRow][i2c(i) + myLowerCol];
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
            int item) {
        myMatrix[i2r(i) + myLowerRow][i2c(i) + myLowerCol] = item;
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
        return new Signed16BitIntegerMatrixReductionBuf_1(myMatrix, myRowRange, myColRange, (IntegerOp) op);
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
        ShortBuffer shortbuffer = buffer.asShortBuffer();
        int n = 0;
        int r = i2r(i);
        int row = r + myLowerRow;
        int c = i2c(i);
        int col = c + myLowerCol;
        int ncols = Math.min(myColCount - c, shortbuffer.remaining());
        while (r < myRowCount && ncols > 0) {
            int[] myMatrix_row = myMatrix[row];
            while (c < ncols) {
                shortbuffer.put((short) myMatrix_row[col]);
                ++c;
                ++col;
            }
            n += ncols;
            ++r;
            ++row;
            c = 0;
            col = myLowerCol;
            ncols = Math.min(myColCount, shortbuffer.remaining());
        }
        buffer.position(buffer.position() + 2 * n);
        return n;
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
        ShortBuffer shortbuffer = buffer.asShortBuffer();
        num = Math.min(num, shortbuffer.remaining());
        int n = 0;
        int r = i2r(i);
        int row = r + myLowerRow;
        int c = i2c(i);
        int col = c + myLowerCol;
        int ncols = Math.min(myColCount - c, num);
        while (r < myRowCount && ncols > 0) {
            int[] myMatrix_row = myMatrix[row];
            for (c = 0; c < ncols; ++c) {
                myMatrix_row[col] = shortbuffer.get();
                ++col;
            }
            num -= ncols;
            n += ncols;
            ++r;
            ++row;
            col = myLowerCol;
            ncols = Math.min(myColCount, num);
        }
        buffer.position(buffer.position() + 2 * n);
        return n;
    }

}
