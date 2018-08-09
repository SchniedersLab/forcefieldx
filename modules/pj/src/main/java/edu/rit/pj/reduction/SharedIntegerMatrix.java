//******************************************************************************
//
// File:    SharedIntegerMatrix.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedIntegerMatrix
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
package edu.rit.pj.reduction;

import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * Class SharedIntegerMatrix provides a matrix reduction variable with elements
 * of type <TT>int</TT>.
 * <P>
 * Class SharedIntegerMatrix is multiple thread safe. The methods use lock-free
 * atomic compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedIntegerMatrix is implemented using class
 * java.util.concurrent.atomic.AtomicIntegerArray.
 *
 * @author Alan Kaminsky
 * @version 18-Feb-2010
 */
public class SharedIntegerMatrix {

// Hidden data members.
    private AtomicIntegerArray[] myMatrix;

// Exported constructors.
    /**
     * Construct a new integer matrix reduction variable with the given number
     * of rows and columns. Each matrix element is initially 0.
     *
     * @param rows Number of rows.
     * @param cols Number of columns.
     * @exception NegativeArraySizeException (unchecked exception) Thrown if
     * <TT>rows</TT> &lt; 0 or <TT>cols</TT>
     * &lt; 0.
     */
    public SharedIntegerMatrix(int rows,
            int cols) {
        myMatrix = new AtomicIntegerArray[rows];
        for (int r = 0; r < rows; ++r) {
            myMatrix[r] = new AtomicIntegerArray(cols);
        }
    }

    /**
     * Construct a new integer matrix reduction variable whose elements are
     * copied from the given matrix. It is assumed that all rows of the
     * <TT>matrix</TT> have the same number of columns.
     *
     * @param matrix Matrix to copy.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>matrix</TT> is null or any row of
     * <TT>matrix</TT> is null.
     */
    public SharedIntegerMatrix(int[][] matrix) {
        myMatrix = new AtomicIntegerArray[matrix.length];
        for (int r = 0; r < matrix.length; ++r) {
            myMatrix[r] = new AtomicIntegerArray(matrix[r]);
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
     * @return Current value.
     */
    public int get(int r,
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
            int value) {
        myMatrix[r].set(c, value);
    }

    /**
     * Set this matrix reduction variable at the given row and column to the
     * given value and return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     * @param value New value.
     * @return Previous value.
     */
    public int getAndSet(int r,
            int c,
            int value) {
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
     * @return True if the update happened, false otherwise.
     */
    public boolean compareAndSet(int r,
            int c,
            int expect,
            int update) {
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
     * @return True if the update happened, false otherwise.
     */
    @SuppressWarnings("deprecation")
    public boolean weakCompareAndSet(int r,
            int c,
            int expect,
            int update) {
        return myMatrix[r].weakCompareAndSet(c, expect, update);
    }

    /**
     * Add one to this matrix reduction variable at the given row and column and
     * return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     * @return Previous value.
     */
    public int getAndIncrement(int r,
            int c) {
        return myMatrix[r].getAndIncrement(c);
    }

    /**
     * Subtract one from this matrix reduction variable at the given row and
     * column and return the previous value.
     *
     * @param r Row index.
     * @param c Column index.
     * @return Previous value.
     */
    public int getAndDecrement(int r,
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
     * @return Previous value.
     */
    public int getAndAdd(int r,
            int c,
            int value) {
        return myMatrix[r].getAndAdd(c, value);
    }

    /**
     * Add one to this matrix reduction variable at the given row and column and
     * return the new value.
     *
     * @param r Row index.
     * @param c Column index.
     * @return New value.
     */
    public int incrementAndGet(int r,
            int c) {
        return myMatrix[r].incrementAndGet(c);
    }

    /**
     * Subtract one from this matrix reduction variable at the given row and
     * column and return the new value.
     *
     * @param r Row index.
     * @param c Column index.
     * @return New value.
     */
    public int decrementAndGet(int r,
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
     * @return New value.
     */
    public int addAndGet(int r,
            int c,
            int value) {
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
     * @return (This matrix <TT>[r,c]</TT>) <I>op</I> (<TT>value</TT>).
     */
    public int reduce(int r,
            int c,
            int value,
            IntegerOp op) {
        AtomicIntegerArray myMatrix_r = myMatrix[r];
        for (;;) {
            int oldvalue = myMatrix_r.get(c);
            int newvalue = op.op(oldvalue, value);
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
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>src</TT> is null. Thrown if
     * <TT>op</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * matrix index would be out of bounds.
     */
    public void reduce(int[][] src,
            IntegerOp op) {
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
     * @param srcrow Row index of first element to update from in the source
     * matrix.
     * @param srccol Column index of first element to update from in the source
     * matrix.
     * @param rowlen Number of rows to update.
     * @param collen Number of columns to update.
     * @param op Binary operation.
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
            int[][] src,
            int srcrow,
            int srccol,
            int rowlen,
            int collen,
            IntegerOp op) {
        if (rowlen < 0
                || collen < 0
                || dstrow < 0 || dstrow + rowlen > rows()
                || dstcol < 0 || dstcol + collen > cols()
                || srcrow < 0 || srcrow + rowlen > src.length
                || srccol < 0 || srccol + collen > src[0].length) {
            throw new IndexOutOfBoundsException();
        }

        for (int r = 0; r < rowlen; ++r) {
            AtomicIntegerArray myMatrix_r = myMatrix[dstrow + r];
            int[] src_r = src[srcrow + r];
            for (int c = 0; c < collen; ++c) {
                int dstcol_c = dstcol + c;
                int src_r_c = src_r[srccol + c];
                updateLoop:
                for (;;) {
                    int oldvalue = myMatrix_r.get(dstcol_c);
                    int newvalue = op.op(oldvalue, src_r_c);
                    if (myMatrix_r.compareAndSet(dstcol_c, oldvalue, newvalue)) {
                        break updateLoop;
                    }
                }
            }
        }
    }

}
