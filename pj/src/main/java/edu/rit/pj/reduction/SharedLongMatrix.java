//******************************************************************************
//
// File:    SharedLongMatrix.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedLongMatrix
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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
package edu.rit.pj.reduction;

import java.util.concurrent.atomic.AtomicLongArray;

/**
 * Class SharedLongMatrix provides a matrix reduction variable with elements of
 * type <TT>long</TT>.
 * <P>
 * Class SharedLongMatrix is multiple thread safe. The methods use lock-free
 * atomic compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedLongMatrix is implemented using class
 * java.util.concurrent.atomic.AtomicLongArray.
 *
 * @author Alan Kaminsky
 * @version 12-Feb-2010
 */
public class SharedLongMatrix {

// Hidden data members.
    private AtomicLongArray[] myMatrix;

// Exported constructors.
    /**
     * Construct a new long matrix reduction variable with the given number of
     * rows and columns. Each matrix element is initially 0.
     *
     * @param rows Number of rows.
     * @param cols Number of columns.
     *
     * @exception NegativeArraySizeException (unchecked exception) Thrown if
     * <TT>rows</TT> &lt; 0 or <TT>cols</TT>
     * &lt; 0.
     */
    public SharedLongMatrix(int rows,
            int cols) {
        myMatrix = new AtomicLongArray[rows];
        for (int r = 0; r < rows; ++r) {
            myMatrix[r] = new AtomicLongArray(cols);
        }
    }

    /**
     * Construct a new long matrix reduction variable whose elements are copied
     * from the given matrix. It is assumed that all rows of the
     * <TT>matrix</TT> have the same number of columns.
     *
     * @param matrix Matrix to copy.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>matrix</TT> is null or any row of
     * <TT>matrix</TT> is null.
     */
    public SharedLongMatrix(long[][] matrix) {
        myMatrix = new AtomicLongArray[matrix.length];
        for (int r = 0; r < matrix.length; ++r) {
            myMatrix[r] = new AtomicLongArray(matrix[r]);
        }
    }

// Exported operations.
    /**
     * Returns the number of rows in this matrix reduction variable.
     *
     * @return Rows.
     */
    public int rows() {
        return myMatrix.length;
    }

    /**
     * Returns the number of columns in this matrix reduction variable.
     *
     * @return Columns.
     */
    public int cols() {
        return myMatrix[0].length();
    }

    /**
     * Returns this matrix reduction variable's current value at the given row
     * and column.
     *
     * @param r Row index.
     * @param c Column index.
     *
     * @return Current value.
     */
    public long get(int r,
            int c) {
        return myMatrix[r].get(c);
    }

    /**
     * Set this matrix reduction variable at the given row and column to the
     * given value.
     *
     * @param r Row index.
     * @param c Column index.
     * @param value New value.
     */
    public void set(int r,
            int c,
            long value) {
        myMatrix[r].set(c, value);
    }

    /**
     * Set this matrix reduction variable at the given row and column to the
     * given value and return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     * @param value New value.
     *
     * @return Previous value.
     */
    public long getAndSet(int r,
            int c,
            long value) {
        return myMatrix[r].getAndSet(c, value);
    }

    /**
     * Atomically set this matrix reduction variable at the given row and column
     * to the given updated value if the current value equals the expected
     * value.
     *
     * @param r Row index.
     * @param c Column index.
     * @param expect Expected value.
     * @param update Updated value.
     *
     * @return True if the update happened, false otherwise.
     */
    public boolean compareAndSet(int r,
            int c,
            long expect,
            long update) {
        return myMatrix[r].compareAndSet(c, expect, update);
    }

    /**
     * Atomically set this matrix reduction variable at the given row and column
     * to the given updated value if the current value equals the expected
     * value. May fail spuriously.
     *
     * @param r Row index.
     * @param c Column index.
     * @param expect Expected value.
     * @param update Updated value.
     *
     * @return True if the update happened, false otherwise.
     */
    public boolean weakCompareAndSet(int r,
            int c,
            long expect,
            long update) {
        return myMatrix[r].weakCompareAndSet(c, expect, update);
    }

    /**
     * Add one to this matrix reduction variable at the given row and column and
     * return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     *
     * @return Previous value.
     */
    public long getAndIncrement(int r,
            int c) {
        return myMatrix[r].getAndIncrement(c);
    }

    /**
     * Subtract one from this matrix reduction variable at the given row and
     * column and return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     *
     * @return Previous value.
     */
    public long getAndDecrement(int r,
            int c) {
        return myMatrix[r].getAndDecrement(c);
    }

    /**
     * Add the given value to this matrix reduction variable at the given row
     * and column and return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     * @param value Value to add.
     *
     * @return Previous value.
     */
    public long getAndAdd(int r,
            int c,
            long value) {
        return myMatrix[r].getAndAdd(c, value);
    }

    /**
     * Add one to this matrix reduction variable at the given row and column and
     * return the new value.
     *
     * @param r Row index.
     * @param c Column index.
     *
     * @return New value.
     */
    public long incrementAndGet(int r,
            int c) {
        return myMatrix[r].incrementAndGet(c);
    }

    /**
     * Subtract one from this matrix reduction variable at the given row and
     * column and return the new value.
     *
     * @param r Row index.
     * @param c Column index.
     *
     * @return New value.
     */
    public long decrementAndGet(int r,
            int c) {
        return myMatrix[r].decrementAndGet(c);
    }

    /**
     * Add the given value to this matrix reduction variable at the given row
     * and column and return the new value.
     *
     * @param r Row index.
     * @param c Column index.
     * @param value Value to add.
     *
     * @return New value.
     */
    public long addAndGet(int r,
            int c,
            long value) {
        return myMatrix[r].addAndGet(c, value);
    }

    /**
     * Combine this matrix reduction variable at the given row and column with
     * the given value using the given operation. (This matrix <TT>[r,c]</TT>)
     * is set to (this matrix <TT>[r,c]</TT>) <I>op</I> (<TT>value</TT>), then
     * (this matrix <TT>[r,c]</TT>) is returned.
     *
     * @param r Row index.
     * @param c Column index.
     * @param value Value.
     * @param op Binary operation.
     *
     * @return (This matrix <TT>[r,c]</TT>) <I>op</I> (<TT>value</TT>).
     */
    public long reduce(int r,
            int c,
            long value,
            LongOp op) {
        AtomicLongArray myMatrix_r = myMatrix[r];
        for (;;) {
            long oldvalue = myMatrix_r.get(c);
            long newvalue = op.op(oldvalue, value);
            if (myMatrix_r.compareAndSet(c, oldvalue, newvalue)) {
                return newvalue;
            }
        }
    }

    /**
     * Combine this matrix reduction variable with the given matrix using the
     * given operation. For every row <TT>r</TT> and column <TT>c</TT> in this
     * matrix, (this matrix <TT>[r,c]</TT>) is set to (this matrix
     * <TT>[r,c]</TT>) <I>op</I> (<TT>src[r,c]</TT>).
     * <P>
     * The <TT>reduce()</TT> method is multiple thread safe <I>on a per-element
     * basis.</I> Each individual matrix element is updated atomically, but the
     * matrix as a whole is not updated atomically.
     *
     * @param src Source matrix.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>src</TT> is null. Thrown if
     * <TT>op</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * matrix index would be out of bounds.
     */
    public void reduce(long[][] src,
            LongOp op) {
        reduce(0, 0, src, 0, 0, rows(), cols(), op);
    }

    /**
     * Combine a portion of this matrix reduction variable with a portion of the
     * given matrix using the given operation. For each row index <TT>r</TT>
     * from 0 to <TT>rowlen-1</TT> inclusive, and for each column index
     * <TT>c</TT> from 0 to <TT>collen-1</TT> inclusive, (this matrix
     * <TT>[dstrow+r,dstcol+c]</TT>) is set to (this matrix
     * <TT>[dstrow+r,dstcol+c]</TT>) <I>op</I>
     * (<TT>src[srcrow+r,srccol+c]</TT>).
     * <P>
     * The <TT>reduce()</TT> method is multiple thread safe <I>on a per-element
     * basis.</I> Each individual matrix element is updated atomically, but the
     * matrix as a whole is not updated atomically.
     *
     * @param dstrow Row index of first element to update in this matrix.
     * @param dstcol Column index of first element to update in this matrix.
     * @param src Source matrix.
     * @param srcrow Row index of first element to update from in the source
     * matrix.
     * @param srccol Column index of first element to update from in the source
     * matrix.
     * @param rowlen Number of rows to update.
     * @param collen Number of columns to update.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>src</TT> is null. Thrown if
     * <TT>op</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>rowlen</TT> &lt; 0. Thrown if
     * <TT>collen</TT> &lt; 0. Thrown if any matrix index would be out of
     * bounds.
     */
    public void reduce(int dstrow,
            int dstcol,
            long[][] src,
            int srcrow,
            int srccol,
            int rowlen,
            int collen,
            LongOp op) {
        if (rowlen < 0
                || collen < 0
                || dstrow < 0 || dstrow + rowlen > rows()
                || dstcol < 0 || dstcol + collen > cols()
                || srcrow < 0 || srcrow + rowlen > src.length
                || srccol < 0 || srccol + collen > src[0].length) {
            throw new IndexOutOfBoundsException();
        }

        for (int r = 0; r < rowlen; ++r) {
            AtomicLongArray myMatrix_r = myMatrix[dstrow + r];
            long[] src_r = src[srcrow + r];
            for (int c = 0; c < collen; ++c) {
                int dstcol_c = dstcol + c;
                long src_r_c = src_r[srccol + c];
                updateLoop:
                for (;;) {
                    long oldvalue = myMatrix_r.get(dstcol_c);
                    long newvalue = op.op(oldvalue, src_r_c);
                    if (myMatrix_r.compareAndSet(dstcol_c, oldvalue, newvalue)) {
                        break updateLoop;
                    }
                }
            }
        }
    }

}
