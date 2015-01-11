//******************************************************************************
//
// File:    SharedDouble.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedDouble
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

import java.util.concurrent.atomic.AtomicLong;

/**
 * Class SharedDouble provides a reduction variable for a value of type
 * <TT>double</TT>.
 * <P>
 * Class SharedDouble is multiple thread safe. The methods use lock-free atomic
 * compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedDouble is implemented using class
 * java.util.concurrent.atomic.AtomicLong. The double value is stored as a
 * <TT>long</TT> whose bit pattern is the same as the double value.
 *
 * @author Alan Kaminsky
 * @version 07-Jun-2007
 */
public class SharedDouble
        extends Number {

// Hidden data members.
    private AtomicLong myValue;

// Exported constructors.
    /**
     * Construct a new double reduction variable with the initial value 0.
     */
    public SharedDouble() {
        myValue = new AtomicLong(Double.doubleToLongBits(0.0));
    }

    /**
     * Construct a new double reduction variable with the given initial value.
     *
     * @param initialValue Initial value.
     */
    public SharedDouble(double initialValue) {
        myValue = new AtomicLong(Double.doubleToLongBits(initialValue));
    }

// Exported operations.
    /**
     * Returns this reduction variable's current value.
     *
     * @return Current value.
     */
    public double get() {
        return Double.longBitsToDouble(myValue.get());
    }

    /**
     * Set this reduction variable to the given value.
     *
     * @param value New value.
     */
    public void set(double value) {
        myValue.set(Double.doubleToLongBits(value));
    }

    /**
     * Set this reduction variable to the given value and return the previous
     * value.
     *
     * @param value New value.
     *
     * @return Previous value.
     */
    public double getAndSet(double value) {
        return Double.longBitsToDouble(myValue.getAndSet(Double.doubleToLongBits(value)));
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
    public boolean compareAndSet(double expect,
            double update) {
        return myValue.compareAndSet(Double.doubleToLongBits(expect), Double.doubleToLongBits(update));
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
    public boolean weakCompareAndSet(double expect,
            double update) {
        return myValue.weakCompareAndSet(Double.doubleToLongBits(expect), Double.doubleToLongBits(update));
    }

    /**
     * Add one to this reduction variable and return the previous value.
     *
     * @return Previous value.
     */
    public double getAndIncrement() {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = oldvalue + 1.0;
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
                return oldvalue;
            }
        }
    }

    /**
     * Subtract one from this reduction variable and return the previous value.
     *
     * @return Previous value.
     */
    public double getAndDecrement() {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = oldvalue - 1.0;
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
                return oldvalue;
            }
        }
    }

    /**
     * Add the given value to this reduction variable and return the previous
     * value.
     *
     * @param value Value to add.
     *
     * @return Previous value.
     */
    public double getAndAdd(double value) {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = oldvalue + value;
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
                return oldvalue;
            }
        }
    }

    /**
     * Add one to this reduction variable and return the new value.
     *
     * @return New value.
     */
    public double incrementAndGet() {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = oldvalue + 1.0;
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
                return newvalue;
            }
        }
    }

    /**
     * Subtract one from this reduction variable and return the new value.
     *
     * @return New value.
     */
    public double decrementAndGet() {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = oldvalue - 1.0;
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
                return newvalue;
            }
        }
    }

    /**
     * Add the given value to this reduction variable and return the new value.
     *
     * @param value Value to add.
     *
     * @return New value.
     */
    public double addAndGet(double value) {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = oldvalue + value;
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
                return newvalue;
            }
        }
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
    public double reduce(double value,
            DoubleOp op) {
        for (;;) {
            long oldvalueLong = myValue.get();
            double oldvalue = Double.longBitsToDouble(oldvalueLong);
            double newvalue = op.op(oldvalue, value);
            long newvalueLong = Double.doubleToLongBits(newvalue);
            if (myValue.compareAndSet(oldvalueLong, newvalueLong)) {
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
        return Double.toString(get());
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <TT>int</TT>.
     *
     * @return Current value.
     */
    public int intValue() {
        return (int) get();
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <TT>long</TT>.
     *
     * @return Current value.
     */
    public long longValue() {
        return (long) get();
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <TT>float</TT>.
     *
     * @return Current value.
     */
    public float floatValue() {
        return (float) get();
    }

    /**
     * Returns this reduction variable's current value converted to type
     * <TT>double</TT>.
     *
     * @return Current value.
     */
    public double doubleValue() {
        return (double) get();
    }

}
