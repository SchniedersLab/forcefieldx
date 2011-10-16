//******************************************************************************
//
// File:    SharedShort.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.SharedShort
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

import java.util.concurrent.atomic.AtomicInteger;

/**
 * Class SharedShort provides a reduction variable for a value of type
 * <TT>short</TT>.
 * <P>
 * Class SharedShort is multiple thread safe. The methods use lock-free atomic
 * compare-and-set.
 * <P>
 * <I>Note:</I> Class SharedShort is implemented using class
 * java.util.concurrent.atomic.AtomicInteger. The short value is stored as an
 * <TT>int</TT> whose values are restricted to the range of type <TT>short</TT>.
 *
 * @author  Alan Kaminsky
 * @version 07-Jun-2007
 */
public class SharedShort
	extends Number
	{

// Hidden data members.

	private AtomicInteger myValue;

// Exported constructors.

	/**
	 * Construct a new short reduction variable with the initial value 0.
	 */
	public SharedShort()
		{
		myValue = new AtomicInteger();
		}

	/**
	 * Construct a new short reduction variable with the given initial
	 * value.
	 *
	 * @param  initialValue  Initial value.
	 */
	public SharedShort
		(short initialValue)
		{
		myValue = new AtomicInteger (initialValue);
		}

// Exported operations.

	/**
	 * Returns this reduction variable's current value.
	 *
	 * @return  Current value.
	 */
	public short get()
		{
		return (short) myValue.get();
		}

	/**
	 * Set this reduction variable to the given value.
	 *
	 * @param  value  New value.
	 */
	public void set
		(short value)
		{
		myValue.set (value);
		}

	/**
	 * Set this reduction variable to the given value and return the previous
	 * value.
	 *
	 * @param  value  New value.
	 *
	 * @return  Previous value.
	 */
	public short getAndSet
		(short value)
		{
		return (short) myValue.getAndSet (value);
		}

	/**
	 * Atomically set this reduction variable to the given updated value if the
	 * current value equals the expected value.
	 *
	 * @param  expect  Expected value.
	 * @param  update  Updated value.
	 *
	 * @return  True if the update happened, false otherwise.
	 */
	public boolean compareAndSet
		(short expect,
		 short update)
		{
		return myValue.compareAndSet (expect, update);
		}

	/**
	 * Atomically set this reduction variable to the given updated value if the
	 * current value equals the expected value. May fail spuriously.
	 *
	 * @param  expect  Expected value.
	 * @param  update  Updated value.
	 *
	 * @return  True if the update happened, false otherwise.
	 */
	public boolean weakCompareAndSet
		(short expect,
		 short update)
		{
		return myValue.weakCompareAndSet (expect, update);
		}

	/**
	 * Add one to this reduction variable and return the previous value.
	 *
	 * @return  Previous value.
	 */
	public short getAndIncrement()
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = (short) (oldvalue + 1);
			if (myValue.compareAndSet (oldvalue, newvalue)) return oldvalue;
			}
		}

	/**
	 * Subtract one from this reduction variable and return the previous value.
	 *
	 * @return  Previous value.
	 */
	public short getAndDecrement()
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = (short) (oldvalue - 1);
			if (myValue.compareAndSet (oldvalue, newvalue)) return oldvalue;
			}
		}

	/**
	 * Add the given value to this reduction variable and return the previous
	 * value.
	 *
	 * @param  value  Value to add.
	 *
	 * @return  Previous value.
	 */
	public short getAndAdd
		(short value)
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = (short) (oldvalue + value);
			if (myValue.compareAndSet (oldvalue, newvalue)) return oldvalue;
			}
		}

	/**
	 * Add one to this reduction variable and return the new value.
	 *
	 * @return  New value.
	 */
	public short incrementAndGet()
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = (short) (oldvalue + 1);
			if (myValue.compareAndSet (oldvalue, newvalue)) return newvalue;
			}
		}

	/**
	 * Subtract one from this reduction variable and return the new value.
	 *
	 * @return  New value.
	 */
	public short decrementAndGet()
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = (short) (oldvalue - 1);
			if (myValue.compareAndSet (oldvalue, newvalue)) return newvalue;
			}
		}

	/**
	 * Add the given value to this reduction variable and return the new value.
	 *
	 * @param  value  Value to add.
	 *
	 * @return  New value.
	 */
	public short addAndGet
		(short value)
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = (short) (oldvalue + value);
			if (myValue.compareAndSet (oldvalue, newvalue)) return newvalue;
			}
		}

	/**
	 * Combine this reduction variable with the given value using the given
	 * operation. The result is stored back into this reduction variable and is
	 * returned.
	 *
	 * @param  value  Value.
	 * @param  op     Binary operation.
	 *
	 * @return  (This variable) <I>op</I> (<TT>value</TT>).
	 */
	public short reduce
		(short value,
		 ShortOp op)
		{
		for (;;)
			{
			short oldvalue = (short) myValue.get();
			short newvalue = op.op (oldvalue, value);
			if (myValue.compareAndSet (oldvalue, newvalue)) return newvalue;
			}
		}

	/**
	 * Returns a string version of this reduction variable.
	 *
	 * @return  String version.
	 */
	public String toString()
		{
		return Integer.toString (get());
		}

	/**
	 * Returns this reduction variable's current value converted to type
	 * <TT>int</TT>.
	 *
	 * @return  Current value.
	 */
	public int intValue()
		{
		return (int) get();
		}

	/**
	 * Returns this reduction variable's current value converted to type
	 * <TT>long</TT>.
	 *
	 * @return  Current value.
	 */
	public long longValue()
		{
		return (long) get();
		}

	/**
	 * Returns this reduction variable's current value converted to type
	 * <TT>float</TT>.
	 *
	 * @return  Current value.
	 */
	public float floatValue()
		{
		return (float) get();
		}

	/**
	 * Returns this reduction variable's current value converted to type
	 * <TT>double</TT>.
	 *
	 * @return  Current value.
	 */
	public double doubleValue()
		{
		return (double) get();
		}

	}
