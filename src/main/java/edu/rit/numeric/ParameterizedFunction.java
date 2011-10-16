//******************************************************************************
//
// File:    ParameterizedFunction.java
// Package: edu.rit.numeric
// Unit:    Interface edu.rit.numeric.ParameterizedFunction
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

package edu.rit.numeric;

/**
 * Interface ParameterizedFunction specifies the interface for a real-valued
 * function of a real-valued argument, with additional parameters.
 *
 * @author  Alan Kaminsky
 * @version 06-Oct-2010
 */
public interface ParameterizedFunction
	{

// Exported operations.

	/**
	 * Returns the number of parameters, <I>M</I>.
	 *
	 * @return  Number of parameters.
	 */
	public int parameterLength();

	/**
	 * Returns the value of this function evaluated at the given argument, using
	 * the given parameters.
	 *
	 * @param  x  Argument.
	 * @param  p  Array of <I>M</I> parameters.
	 *
	 * @return  <I>f</I> (<I>x</I>, parameters).
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if the argument <TT>x</TT> or any
	 *     parameter in <TT>p</TT> is outside the allowed set of values for this
	 *     function.
	 * @exception  RangeException
	 *     (unchecked exception) Thrown if the result is outside the range of
	 *     type <TT>double</TT>.
	 */
	public double f
		(double x,
		 double[] p);

	}
