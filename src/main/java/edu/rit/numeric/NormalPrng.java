//******************************************************************************
//
// File:    NormalPrng.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.NormalPrng
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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
 * Class NormalPrng provides a pseudorandom number generator (PRNG) that
 * generates random numbers with a normal distribution.
 * <P>
 * Every two calls of the <TT>next()</TT> method result in two calls of the
 * underlying uniform PRNG's <TT>nextDouble()</TT> method. More precisely, the
 * first of a pair of <TT>next()</TT> calls does two <TT>nextDouble()</TT>
 * calls, the second of a pair of <TT>next()</TT> calls does not call
 * <TT>nextDouble()</TT>.
 * <P>
 * Class NormalPrng uses the Box-Muller method to generate a standard normal
 * distribution. Let <I>x</I><SUB>1</SUB> and <I>x</I><SUB>2</SUB> be drawn from
 * a uniform distribution between 0 and 1. Then <I>y</I><SUB>1</SUB> and
 * <I>y</I><SUB>2</SUB>, defined as follows, are drawn from a standard normal
 * distribution.
 * <CENTER>
 * <I>y</I><SUB>1</SUB> = sqrt(&minus;2 ln <I>x</I><SUB>1</SUB>) cos 2&pi;<I>x</I><SUB>2</SUB>
 * <BR><I>y</I><SUB>2</SUB> = sqrt(&minus;2 ln <I>x</I><SUB>1</SUB>) sin 2&pi;<I>x</I><SUB>2</SUB>
 * </CENTER>
 * And (&mu;&nbsp;+&nbsp;&sigma;<I>y</I><SUB>1</SUB>) and
 * (&mu;&nbsp;+&nbsp;&sigma;<I>y</I><SUB>2</SUB>) are drawn from a normal
 * distribution with mean &mu; and standard deviation &sigma;.
 * <P>
 * <I>Note:</I> While slower than other techniques, this technique does a fixed
 * number of <TT>nextDouble()</TT> calls to generate each normally distributed
 * random number. This behavior is important when parallel programs generate
 * random numbers.
 *
 * @author  Alan Kaminsky
 * @version 10-Jun-2008
 */
public class NormalPrng
	extends DoublePrng
	{

// Hidden data members.

	private double myMean;
	private double myStdev;

	private double u1 = -1.0;
	private double u2;

	private static final double TWO_PI = 2.0*Math.PI;

// Exported constructors.

	/**
	 * Construct a new standard normal PRNG. The mean is 0 and the standard
	 * deviation is 1.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 */
	public NormalPrng
		(Random theUniformPrng)
		{
		this (theUniformPrng, 0.0, 1.0);
		}

	/**
	 * Construct a new normal PRNG with the given mean and standard deviation.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 * @param  theMean         Mean of the normal distribution.
	 * @param  theStdev        Standard deviation of the normal distribution.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 */
	public NormalPrng
		(Random theUniformPrng,
		 double theMean,
		 double theStdev)
		{
		super (theUniformPrng);
		myMean = theMean;
		myStdev = theStdev;
		}

// Exported operations.

	/**
	 * Returns the next random number.
	 *
	 * @return  Random number.
	 */
	public double next()
		{
		double y;
		if (u1 == -1.0)
			{
			double x1 = myUniformPrng.nextDouble();
			double x2 = myUniformPrng.nextDouble();
			u1 = Math.sqrt(-2.0*Math.log(x1));
			u2 = TWO_PI*x2;
			y = u1*Math.cos(u2);
			}
		else
			{
			y = u1*Math.sin(u2);
			u1 = -1.0;
			}
		return myMean + myStdev*y;
		}

	}
