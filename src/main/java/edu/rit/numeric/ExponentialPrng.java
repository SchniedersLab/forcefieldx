//******************************************************************************
//
// File:    ExponentialPrng.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.ExponentialPrng
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
 * Class ExponentialPrng provides a pseudorandom number generator (PRNG) that
 * generates random numbers with an exponential distribution. The probability
 * density function is
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;<I>f</I>(<I>x</I>) = <I>&lambda;</I>e<SUP>&minus;<I>&lambda;</I><I>x</I></SUP>, <I>x</I> &ge; 0
 * <BR>&nbsp;&nbsp;&nbsp;&nbsp;<I>f</I>(<I>x</I>) = 0, otherwise
 * <BR>The distribution's mean is 1/<I>&lambda;</I> and its standard deviation
 * is 1/<I>&lambda;</I><SUP>2</SUP>.
 * <P>
 * An exponential distribution is often used to model arrivals or departures in
 * a discrete event simulation. The mean arrival or departure rate is
 * <I>&lambda;</I>; the mean interarrival or interdeparture time is
 * 1/<I>&lambda;</I>.
 * <P>
 * Every call of the <TT>next()</TT> method results in one call of the
 * underlying uniform PRNG's <TT>nextDouble()</TT> method.
 *
 * @author  Alan Kaminsky
 * @version 01-Aug-2011
 */
public class ExponentialPrng
	extends DoublePrng
	{

// Hidden data members.

	private double lambda;

// Exported constructors.

	/**
	 * Construct a new exponential PRNG.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 * @param  lambda          Mean rate <I>&lambda;</I> &gt; 0.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>&lambda;</I> &le; 0.
	 */
	public ExponentialPrng
		(Random theUniformPrng,
		 double lambda)
		{
		super (theUniformPrng);
		if (lambda <= 0)
			{
			throw new IllegalArgumentException
				("ExponentialPrng(): lambda = "+lambda+" illegal");
			}
		this.lambda = lambda;
		}

// Exported operations.

	/**
	 * Returns the next random number.
	 *
	 * @return  Random number.
	 */
	public double next()
		{
		return -Math.log(myUniformPrng.nextDouble())/lambda;
		}

	}
