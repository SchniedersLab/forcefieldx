//******************************************************************************
//
// File:    SharedBoolean.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedBoolean
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

import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Class SharedBoolean provides a reduction variable for a value of type
 * <TT>boolean</TT>.
 * <P>
 * Class SharedBoolean is multiple thread safe. The methods use lock-free atomic
 * compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedBoolean is implemented using class
 * java.util.concurrent.atomic.AtomicBoolean. The Boolean value is stored as an
 * <TT>int</TT> with 0 = false and 1 = true.
 *
 * @author Alan Kaminsky
 * @version 06-Jun-2007
 */
public class SharedBoolean {

// Hidden data members.
    private AtomicBoolean myValue;

// Exported constructors.
    /**
     * Construct a new Boolean reduction variable with the initial value false.
     */
    public SharedBoolean() {
        myValue = new AtomicBoolean();
    }

    /**
     * Construct a new Boolean reduction variable with the given initial value.
     *
     * @param initialValue Initial value.
     */
    public SharedBoolean(boolean initialValue) {
        myValue = new AtomicBoolean(initialValue);
    }

// Exported operations.
    /**
     * Returns this reduction variable's current value.
     *
     * @return Current value.
     */
    public boolean get() {
        return myValue.get();
    }

    /**
     * Set this reduction variable to the given value.
     *
     * @param value New value.
     */
    public void set(boolean value) {
        myValue.set(value);
    }

    /**
     * Set this reduction variable to the given value and return the previous
     * value.
     *
     * @param value New value.
     *
     * @return Previous value.
     */
    public boolean getAndSet(boolean value) {
        return myValue.getAndSet(value);
    }

    /**
     * Atomically set this reduction variable to the given updated value if the
     * current value equals the expected value.
     *
     * @param expect Expected value.
     * @param update Updated value.
     *
     * @return True if the update happened, false otherwise.
     */
    public boolean compareAndSet(boolean expect,
            boolean update) {
        return myValue.compareAndSet(expect, update);
    }

    /**
     * Atomically set this reduction variable to the given updated value if the
     * current value equals the expected value. May fail spuriously.
     *
     * @param expect Expected value.
     * @param update Updated value.
     *
     * @return True if the update happened, false otherwise.
     */
    public boolean weakCompareAndSet(boolean expect,
            boolean update) {
        return myValue.weakCompareAndSet(expect, update);
    }

    /**
     * Combine this reduction variable with the given value using the given
     * operation. The result is stored back into this reduction variable and is
     * returned.
     *
     * @param value Value.
     * @param op Binary operation.
     *
     * @return (This variable) <I>op</I> (<TT>value</TT>).
     */
    public boolean reduce(boolean value,
            BooleanOp op) {
        for (;;) {
            boolean oldvalue = myValue.get();
            boolean newvalue = op.op(oldvalue, value);
            if (myValue.compareAndSet(oldvalue, newvalue)) {
                return newvalue;
            }
        }
    }

    /**
     * Returns a string version of this reduction variable.
     *
     * @return String version.
     */
    public String toString() {
        return Boolean.toString(get());
    }

}
