//******************************************************************************
//
// File:    MDFunction.java
// Package: edu.rit.numeric
// Unit:    Interface edu.rit.numeric.MDFunction
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

/**
 * Interface MDFunction specifies the interface for a function whose argument is
 * a vector of real values and whose result is a real value; a multidimensional
 * function.
 *
 * @author  Alan Kaminsky
 * @version 05-Oct-2010
 */
public interface MDFunction
	{

// Exported operations.

	/**
	 * Returns the length of the argument vector, <I>N</I>.
	 *
	 * @return  <I>N</I>.
	 */
	public int argumentLength();

	/**
	 * Evaluate this function with the given argument vector.
	 *
	 * @param  x  Argument vector. Must be an <I>N</I>-element array.
	 *
	 * @return  Function value at <TT>x</TT>.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if any argument in <TT>x</TT> is outside
	 *     the allowed set of values for this function.
	 * @exception  RangeException
	 *     (unchecked exception) Thrown if the function value is outside the
	 *     range of type <TT>double</TT>.
	 */
	public double f
		(double[] x);

	}
