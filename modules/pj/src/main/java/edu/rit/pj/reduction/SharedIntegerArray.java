//******************************************************************************
//
// File:    SharedIntegerArray.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedIntegerArray
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
package edu.rit.pj.reduction;

import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * Class SharedIntegerArray provides an array reduction variable with elements
 * of type <code>int</code>.
 * <P>
 * Class SharedIntegerArray is multiple thread safe. The methods use lock-free
 * atomic compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedIntegerArray is implemented using class
 * java.util.concurrent.atomic.AtomicIntegerArray.
 *
 * @author Alan Kaminsky
 * @version 24-Aug-2007
 */
public class SharedIntegerArray {

// Hidden data members.
    private AtomicIntegerArray myArray;

// Exported constructors.
    /**
     * Construct a new integer array reduction variable with the given length.
     * Each array element is initially 0.
     *
     * @param len Length.
     * @exception NegativeArraySizeException (unchecked exception) Thrown if
     * <code>len</code> &lt; 0.
     */
    public SharedIntegerArray(int len) {
        myArray = new AtomicIntegerArray(len);
    }

    /**
     * Construct a new integer array reduction variable whose elements are
     * copied from the given array.
     *
     * @param array Array to copy.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>array</code> is null.
     */
    public SharedIntegerArray(int[] array) {
        myArray = new AtomicIntegerArray(array);
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
    public int get(int i) {
        return myArray.get(i);
    }

    /**
     * Set this array reduction variable at the given index to the given value.
     *
     * @param i Index.
     * @param value New value.
     */
    public void set(int i,
            int value) {
        myArray.set(i, value);
    }

    /**
     * Set this array reduction variable at the given index to the given value
     * and return the previous value.
     *
     * @param i Index.
     * @param value New value.
     * @return Previous value.
     */
    public int getAndSet(int i,
            int value) {
        return myArray.getAndSet(i, value);
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
            int expect,
            int update) {
        return myArray.compareAndSet(i, expect, update);
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
            int expect,
            int update) {
        return myArray.weakCompareAndSet(i, expect, update);
    }

    /**
     * Add one to this array reduction variable at the given index and return
     * the previous value.
     *
     * @param i Index.
     * @return Previous value.
     */
    public int getAndIncrement(int i) {
        return myArray.getAndIncrement(i);
    }

    /**
     * Subtract one from this array reduction variable at the given index and
     * return the previous value.
     *
     * @param i Index.
     * @return Previous value.
     */
    public int getAndDecrement(int i) {
        return myArray.getAndDecrement(i);
    }

    /**
     * Add the given value to this array reduction variable at the given index
     * and return the previous value.
     *
     * @param i Index.
     * @param value Value to add.
     * @return Previous value.
     */
    public int getAndAdd(int i,
            int value) {
        return myArray.getAndAdd(i, value);
    }

    /**
     * Add one to this array reduction variable at the given index and return
     * the new value.
     *
     * @param i Index.
     * @return New value.
     */
    public int incrementAndGet(int i) {
        return myArray.incrementAndGet(i);
    }

    /**
     * Subtract one from this array reduction variable at the given index and
     * return the new value.
     *
     * @param i Index.
     * @return New value.
     */
    public int decrementAndGet(int i) {
        return myArray.decrementAndGet(i);
    }

    /**
     * Add the given value to this array reduction variable at the given index
     * and return the new value.
     *
     * @param i Index.
     * @param value Value to add.
     * @return New value.
     */
    public int addAndGet(int i,
            int value) {
        return myArray.addAndGet(i, value);
    }

    /**
     * Combine this array reduction variable at the given index with the given
     * value using the given operation. (This array <code>[i]</code>) is set to
     * (this array <code>[i]</code>) <I>op</I> (<code>value</code>), then (this array
     * <code>[i]</code>) is returned.
     *
     * @param i Index.
     * @param value Value.
     * @param op Binary operation.
     * @return (This array <code>[i]</code>) <I>op</I> (<code>value</code>).
     */
    public int reduce(int i,
            int value,
            IntegerOp op) {
        for (;;) {
            int oldvalue = myArray.get(i);
            int newvalue = op.op(oldvalue, value);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return newvalue;
            }
        }
    }

    /**
     * Combine this array reduction variable with the given array using the
     * given operation. For each index <code>i</code> from 0 to this array's
     * length-1, (this array <code>[i]</code>) is set to (this array <code>[i]</code>)
     * <I>op</I> (<code>src[i]</code>).
     * <P>
     * The <code>reduce()</code> method is multiple thread safe <I>on a per-element
     * basis.</I> Each individual array element is updated atomically, but the
     * array as a whole is not updated atomically.
     *
     * @param src Source array.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>op</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * array index would be out of bounds.
     */
    public void reduce(int[] src,
            IntegerOp op) {
        reduce(0, src, 0, myArray.length(), op);
    }

    /**
     * Combine a portion of this array reduction variable with a portion of the
     * given array using the given operation. For each index <code>i</code> from 0
     * to <code>len</code>-1, (this array <code>[dstoff+i]</code>) is set to (this array
     * <code>[dstoff+i]</code>) <I>op</I> (<code>src[srcoff+i]</code>).
     * <P>
     * The <code>reduce()</code> method is multiple thread safe <I>on a per-element
     * basis.</I> Each individual array element is updated atomically, but the
     * array as a whole is not updated atomically.
     *
     * @param dstoff Index of first element to update in this array.
     * @param src Source array.
     * @param srcoff Index of first element to update from in the source array.
     * @param len Number of array elements to update.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>op</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>len</code> &lt; 0. Thrown if any array index would be out of bounds.
     */
    public void reduce(int dstoff,
            int[] src,
            int srcoff,
            int len,
            IntegerOp op) {
        if (len < 0
                || dstoff < 0 || dstoff + len > myArray.length()
                || srcoff < 0 || srcoff + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        while (len > 0) {
            updateLoop:
            for (;;) {
                int oldvalue = myArray.get(dstoff);
                int newvalue = op.op(oldvalue, src[srcoff]);
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
