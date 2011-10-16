//******************************************************************************
//
// File:    DoublePrng.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.DoublePrng
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
 * Class DoublePrng is the abstract base class for a pseudorandom number
 * generator (PRNG) that generates random numbers of type <TT>double</TT>. The
 * generated random numbers' probability distribution is determined by the
 * subclass.
 *
 * @author  Alan Kaminsky
 * @version 10-Jun-2008
 */
public abstract class DoublePrng
	{

// Hidden data members.

	/**
	 * The underlying uniform PRNG.
	 */
	protected final Random myUniformPrng;

// Exported constructors.

	/**
	 * Construct a new double PRNG.
	 *
	 * @param  theUniformPrng  The underlying uniform PRNG.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theUniformPrng</TT> is null.
	 */
	public DoublePrng
		(Random theUniformPrng)
		{
		if (theUniformPrng == null)
			{
			throw new NullPointerException
				("DoublePrng(): theUniformPrng is null");
			}
		myUniformPrng = theUniformPrng;
		}

// Exported operations.

	/**
	 * Returns the next random number.
	 *
	 * @return  Random number.
	 */
	public abstract double next();

	}
