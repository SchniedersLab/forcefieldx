//******************************************************************************
//
// File:    BernoulliPrng.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.BernoulliPrng
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
 * Class BernoulliPrng provides a pseudorandom number generator (PRNG) that
 * generates random values of type <TT>boolean</TT>. The probability of the
 * value true (heads) is <I>p</I>, specified to the constructor. The probability
 * of the value false (tails) is <I>q</I> = 1&nbsp;&minus;&nbsp;<I>p</I>.
 *
 * @author  Alan Kaminsky
 * @version 07-Jun-2011
 */
public class BernoulliPrng
	{

// Hidden data members.

	private final Random myUniformPrng;
	private final double myP;

// Exported constructors.

	/**
	 * Construct a new Bernoulli PRNG with heads probability = 0.5.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 */
	public BernoulliPrng
		(Random theUniformPrng)
		{
		this (theUniformPrng, 0.5);
		}

	/**
	 * Construct a new Bernoulli PRNG with the given heads probability.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 * @param  p               Heads probability, 0 &le; <I>p</I> &le; 1.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>p</I> &lt; 0 or <I>p</I> &gt; 1.
	 */
	public BernoulliPrng
		(Random theUniformPrng,
		 double p)
		{
		if (theUniformPrng == null)
			{
			throw new NullPointerException
				("BernoulliPrng(): theUniformPrng is null");
			}
		if (0.0 > p || p > 1.0)
			{
			throw new IllegalArgumentException
				("BernoulliPrng(): p = "+p+" illegal");
			}
		myUniformPrng = theUniformPrng;
		myP = p;
		}

// Exported operations.

	/**
	 * Returns the next random value. True (heads) is returned with probability
	 * <I>p</I>, false (tails) is returned with probability
	 * 1&nbsp;&minus;&nbsp;<I>p</I>.
	 *
	 * @return  Random value.
	 */
	public boolean next()
		{
		return myUniformPrng.nextDouble() < myP;
		}

	}
