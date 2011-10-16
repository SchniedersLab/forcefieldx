//******************************************************************************
//
// File:    UniformPrng.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.UniformPrng
//
// This Java source file is copyright (C) 2011 by Alan Kaminsky. All rights
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

package edu.rit.numeric;

import edu.rit.util.Random;

/**
 * Class UniformPrng provides a pseudorandom number generator (PRNG) that
 * generates random numbers with a uniform(<I>a,b</I>) distribution. The
 * probability density function is
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;<I>f</I>(<I>x</I>) = 1/(<I>b</I> &minus; <I>a</I>), <I>a</I> &le; <I>x</I> &le; <I>b</I>
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;<I>f</I>(<I>x</I>) = 0, otherwise
 * <P>
 * Every call of the <TT>next()</TT> method results in one call of the
 * underlying uniform PRNG's <TT>nextDouble()</TT> method.
 *
 * @author  Alan Kaminsky
 * @version 04-Aug-2011
 */
public class UniformPrng
	extends DoublePrng
	{

// Hidden data members.

	private double a;
	private double b_minus_a;

// Exported constructors.

	/**
	 * Construct a new uniform PRNG.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 * @param  a               Interval lower bound.
	 * @param  b               Interval upper bound.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>a</I> &ge; <I>b</I>.
	 */
	public UniformPrng
		(Random theUniformPrng,
		 double a,
		 double b)
		{
		super (theUniformPrng);
		if (a >= b)
			{
			throw new IllegalArgumentException
				("UniformPrng(): a ("+a+") >= b ("+b+") illegal");
			}
		this.a = a;
		this.b_minus_a = b - a;
		}

// Exported operations.

	/**
	 * Returns the next random number.
	 *
	 * @return  Random number.
	 */
	public double next()
		{
		return myUniformPrng.nextDouble()*b_minus_a + a;
		}

	}
