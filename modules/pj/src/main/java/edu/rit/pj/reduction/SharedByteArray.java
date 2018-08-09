//******************************************************************************
//
// File:    SharedByteArray.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedByteArray
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
 * Class SharedByteArray provides an array reduction variable with elements of
 * type <TT>byte</TT>.
 * <P>
 * Class SharedByteArray is multiple thread safe. The methods use lock-free
 * atomic compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedByteArray is implemented using class
 * java.util.concurrent.atomic.AtomicIntegerArray. Each byte array element is
 * stored as an <TT>int</TT> whose values are restricted to the range of type
 * <TT>byte</TT>.
 *
 * @author Alan Kaminsky
 * @version 24-Aug-2007
 */
public class SharedByteArray {

// Hidden data members.
    private AtomicIntegerArray myArray;

// Exported constructors.
    /**
     * Construct a new byte array reduction variable with the given length. Each
     * array element is initially 0.
     *
     * @param len Length.
     * @exception NegativeArraySizeException (unchecked exception) Thrown if
     * <TT>len</TT> &lt; 0.
     */
    public SharedByteArray(int len) {
        myArray = new AtomicIntegerArray(len);
    }

    /**
     * Construct a new byte array reduction variable whose elements are copied
     * from the given array.
     *
     * @param array Array to copy.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>array</TT> is null.
     */
    public SharedByteArray(byte[] array) {
        int n = array.length;
        int[] intarray = new int[n];
        for (int i = 0; i < n; ++i) {
            intarray[i] = array[i];
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
    public byte get(int i) {
        return (byte) myArray.get(i);
    }

    /**
     * Set this array reduction variable at the given index to the given value.
     *
     * @param i Index.
     * @param value New value.
     */
    public void set(int i,
            byte value) {
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
    public byte getAndSet(int i,
            byte value) {
        return (byte) myArray.getAndSet(i, value);
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
            byte expect,
            byte update) {
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
            byte expect,
            byte update) {
        return myArray.weakCompareAndSet(i, expect, update);
    }

    /**
     * Add one to this array reduction variable at the given index and return
     * the previous value.
     *
     * @param i Index.
     * @return Previous value.
     */
    public byte getAndIncrement(int i) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = (byte) (oldvalue + 1);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return oldvalue;
            }
        }
    }

    /**
     * Subtract one from this array reduction variable at the given index and
     * return the previous value.
     *
     * @param i Index.
     * @return Previous value.
     */
    public byte getAndDecrement(int i) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = (byte) (oldvalue - 1);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return oldvalue;
            }
        }
    }

    /**
     * Add the given value to this array reduction variable at the given index
     * and return the previous value.
     *
     * @param i Index.
     * @param value Value to add.
     * @return Previous value.
     */
    public byte getAndAdd(int i,
            byte value) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = (byte) (oldvalue + value);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return oldvalue;
            }
        }
    }

    /**
     * Add one to this array reduction variable at the given index and return
     * the new value.
     *
     * @param i Index.
     * @return New value.
     */
    public byte incrementAndGet(int i) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = (byte) (oldvalue + 1);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return newvalue;
            }
        }
    }

    /**
     * Subtract one from this array reduction variable at the given index and
     * return the new value.
     *
     * @param i Index.
     * @return New value.
     */
    public byte decrementAndGet(int i) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = (byte) (oldvalue - 1);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return newvalue;
            }
        }
    }

    /**
     * Add the given value to this array reduction variable at the given index
     * and return the new value.
     *
     * @param i Index.
     * @param value Value to add.
     * @return New value.
     */
    public byte addAndGet(int i,
            byte value) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = (byte) (oldvalue + value);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return newvalue;
            }
        }
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
    public byte reduce(int i,
            byte value,
            ByteOp op) {
        for (;;) {
            byte oldvalue = (byte) myArray.get(i);
            byte newvalue = op.op(oldvalue, value);
            if (myArray.compareAndSet(i, oldvalue, newvalue)) {
                return newvalue;
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
    public void reduce(byte[] src,
            ByteOp op) {
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
            byte[] src,
            int srcoff,
            int len,
            ByteOp op) {
        if (len < 0
                || dstoff < 0 || dstoff + len > myArray.length()
                || srcoff < 0 || srcoff + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        while (len > 0) {
            updateLoop:
            for (;;) {
                byte oldvalue = (byte) myArray.get(dstoff);
                byte newvalue = op.op(oldvalue, src[srcoff]);
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
