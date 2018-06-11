//******************************************************************************
//
// File:    SharedBooleanArray.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedBooleanArray
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
package edu.rit.pj.reduction;

import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * Class SharedBooleanArray provides an array reduction variable with elements
 * of type <TT>boolean</TT>.
 * <P>
 * Class SharedBooleanArray is multiple thread safe. The methods use lock-free
 * atomic compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedBooleanArray is implemented using class
 * java.util.concurrent.atomic.AtomicIntegerArray. Each Boolean array element is
 * stored as an <TT>int</TT> whose values are restricted to the range of type
 * <TT>boolean</TT> (0 = <TT>false</TT>, 1 = <TT>true</TT>).
 *
 * @author Alan Kaminsky
 * @version 07-Jun-2007
 */
public class SharedBooleanArray {

// Hidden data members.
    private AtomicIntegerArray myArray;

// Exported constructors.
    /**
     * Construct a new Boolean array reduction variable with the given length.
     * Each array element is initially <TT>false</TT>.
     *
     * @param len Length.
     * @exception NegativeArraySizeException (unchecked exception) Thrown if
     * <TT>len</TT> &lt; 0.
     */
    public SharedBooleanArray(int len) {
        myArray = new AtomicIntegerArray(len);
    }

    /**
     * Construct a new Boolean array reduction variable whose elements are
     * copied from the given array.
     *
     * @param array Array to copy.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>array</TT> is null.
     */
    public SharedBooleanArray(boolean[] array) {
        int n = array.length;
        int[] intarray = new int[n];
        for (int i = 0; i < n; ++i) {
            intarray[i] = array[i] ? 1 : 0;
        }
        myArray = new AtomicIntegerArray(intarray);
    }

// Exported operations.
    /**
     * Returns this array reduction variable's length.
     *
     * @return Length.
     */
    public int length() {
        return myArray.length();
    }

    /**
     * Returns this array reduction variable's current value at the given index.
     *
     * @param i Index.
     * @return Current value.
     */
    public boolean get(int i) {
        return myArray.get(i) != 0;
    }

    /**
     * Set this array reduction variable at the given index to the given value.
     *
     * @param i Index.
     * @param value New value.
     */
    public void set(int i,
            boolean value) {
        myArray.set(i, value ? 1 : 0);
    }

    /**
     * Set this array reduction variable at the given index to the given value
     * and return the previous value.
     *
     * @param i Index.
     * @param value New value.
     * @return Previous value.
     */
    public boolean getAndSet(int i,
            boolean value) {
        return myArray.getAndSet(i, value ? 1 : 0) != 0;
    }

    /**
     * Atomically set this array reduction variable at the given index to the
     * given updated value if the current value equals the expected value.
     *
     * @param i Index.
     * @param expect Expected value.
     * @param update Updated value.
     * @return True if the update happened, false otherwise.
     */
    public boolean compareAndSet(int i,
            boolean expect,
            boolean update) {
        return myArray.compareAndSet(i, expect ? 1 : 0, update ? 1 : 0);
    }

    /**
     * Atomically set this array reduction variable at the given index to the
     * given updated value if the current value equals the expected value. May
     * fail spuriously.
     *
     * @param i Index.
     * @param expect Expected value.
     * @param update Updated value.
     * @return True if the update happened, false otherwise.
     */
    @SuppressWarnings("deprecation")
    public boolean weakCompareAndSet(int i,
            boolean expect,
            boolean update) {
        return myArray.weakCompareAndSet(i, expect ? 1 : 0, update ? 1 : 0);
    }

    /**
     * Combine this array reduction variable at the given index with the given
     * value using the given operation. (This array <TT>[i]</TT>) is set to
     * (this array <TT>[i]</TT>) <I>op</I> (<TT>value</TT>), then (this array
     * <TT>[i]</TT>) is returned.
     *
     * @param i Index.
     * @param value Value.
     * @param op Binary operation.
     * @return (This array <TT>[i]</TT>) <I>op</I> (<TT>value</TT>).
     */
    public boolean reduce(int i,
            boolean value,
            BooleanOp op) {
        for (;;) {
            int oldvalue = myArray.get(i);
            int newvalue = op.op(oldvalue != 0, value) ? 1 : 0;
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return newvalue != 0;
            }
        }
    }

    /**
     * Combine this array reduction variable with the given array using the
     * given operation. For each index <TT>i</TT> from 0 to this array's
     * length-1, (this array <TT>[i]</TT>) is set to (this array <TT>[i]</TT>)
     * <I>op</I> (<TT>src[i]</TT>).
     * <P>
     * The <TT>reduce()</TT> method is multiple thread safe <I>on a per-element
     * basis.</I> Each individual array element is updated atomically, but the
     * array as a whole is not updated atomically.
     *
     * @param src Source array.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>src</TT> is null. Thrown if
     * <TT>op</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * array index would be out of bounds.
     */
    public void reduce(boolean[] src,
            BooleanOp op) {
        reduce(0, src, 0, myArray.length(), op);
    }

    /**
     * Combine a portion of this array reduction variable with a portion of the
     * given array using the given operation. For each index <TT>i</TT> from 0
     * to <TT>len</TT>-1, (this array <TT>[dstoff+i]</TT>) is set to (this array
     * <TT>[dstoff+i]</TT>) <I>op</I> (<TT>src[srcoff+i]</TT>).
     * <P>
     * The <TT>reduce()</TT> method is multiple thread safe <I>on a per-element
     * basis.</I> Each individual array element is updated atomically, but the
     * array as a whole is not updated atomically.
     *
     * @param dstoff Index of first element to update in this array.
     * @param src Source array.
     * @param srcoff Index of first element to update from in the source array.
     * @param srcoff Index of first element to update from in the source array.
     * @param len Number of array elements to update.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>src</TT> is null. Thrown if
     * <TT>op</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>len</TT> &lt; 0. Thrown if any array index would be out of bounds.
     */
    public void reduce(int dstoff,
            boolean[] src,
            int srcoff,
            int len,
            BooleanOp op) {
        if (len < 0
                || dstoff < 0 || dstoff + len > myArray.length()
                || srcoff < 0 || srcoff + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        while (len > 0) {
            updateLoop:
            for (;;) {
                int oldvalue = myArray.get(dstoff);
                int newvalue = op.op(oldvalue != 0, src[srcoff]) ? 1 : 0;
                if (myArray.compareAndSet(dstoff, oldvalue, newvalue)) {
                    break updateLoop;
                }
            }
            ++dstoff;
            ++srcoff;
            --len;
        }
    }

    /**
     * Returns a string version of this array reduction variable.
     *
     * @return String version.
     */
    public String toString() {
        return myArray.toString();
    }

}
