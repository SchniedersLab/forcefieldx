//******************************************************************************
//
// File:    VectorFunction.java
// Package: edu.rit.numeric
// Unit:    Interface edu.rit.numeric.VectorFunction
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

/**
 * Interface VectorFunction specifies the interface for a function whose
 * argument is a vector of real values and whose result is a vector of real
 * values.
 *
 * @author  Alan Kaminsky
 * @version 09-Jun-2008
 */
public interface VectorFunction
	{

// Exported operations.

	/**
	 * Returns the length of the result vector, <I>M</I>.
	 *
	 * @return  <I>M</I>.
	 */
	public int resultLength();

	/**
	 * Returns the length of the argument vector, <I>N</I>.
	 *
	 * @return  <I>N</I>.
	 */
	public int argumentLength();

	/**
	 * Evaluate this function with the given argument vector. The result is
	 * stored in the vector <TT>y</TT>. Specifically, for <I>i</I> = 0 to
	 * <I>M</I>&minus;1, <I>y</I><SUB><I>i</I></SUB> =
	 * <I>f</I><SUB><I>i</I></SUB>(<B>x</B>).
	 *
	 * @param  x  Argument vector (input). Must be an <I>N</I>-element array.
	 * @param  y  Result vector (output). Must be an <I>M</I>-element array.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if any argument in <TT>x</TT> is outside
	 *     the allowed set of values for this function.
	 * @exception  RangeException
	 *     (unchecked exception) Thrown if any element of the result vector is
	 *     outside the range of type <TT>double</TT>.
	 */
	public void f
		(double[] x,
		 double[] y);

	/**
	 * Calculate this function's Jacobian matrix with the given argument vector.
	 * The result is stored in the matrix <TT>dydx</TT>. Specifically, for
	 * <I>i</I> = 0 to <I>M</I>&minus;1 and <I>j</I> = 0 to <I>N</I>&minus;1,
	 * <I>dydx</I><SUB><I>i,j</I></SUB> =
	 * &part;<I>f</I><SUB><I>i</I></SUB>(<B>x</B>)&nbsp;&frasl;&nbsp;&part;<I>x</I><SUB><I>j</I></SUB>&nbsp;.
	 *
	 * @param  x     Argument vector (input). Must be an <I>N</I>-element array.
	 * @param  dydx  Jacobian matrix (output). Must be an
	 *               <I>M</I>&times;<I>N</I>-element matrix.
	 *
	 * @exception  DomainException
	 *     (unchecked exception) Thrown if any argument in <TT>x</TT> is outside
	 *     the allowed set of values for this function.
	 * @exception  RangeException
	 *     (unchecked exception) Thrown if any element of the Jacobian matrix is
	 *     outside the range of type <TT>double</TT>.
	 */
	public void df
		(double[] x,
		 double[][] dydx);

	}
