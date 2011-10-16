//******************************************************************************
//
// File:    SampledXYSeries.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.SampledXYSeries
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
 * Class SampledXYSeries provides a series of (<I>x,y</I>) pairs of real values
 * (type <TT>double</TT>) where the Y values are computed by sampling a
 * {@linkplain Function} for a series of X values.
 * <P>
 * Class SampledXYSeries is implemented to minimize space rather than time. Each
 * time <TT>x()</TT> is called, the <I>x</I> value is computed anew. Each time
 * <TT>y()</TT> is called, the <I>x</I> value is computed anew, then the
 * <I>y</I> value is computed anew by evaluating the function at <I>x</I>.
 *
 * @author  Alan Kaminsky
 * @version 06-Jul-2007
 */
public class SampledXYSeries
	extends XYSeries
	{

// Hidden data members.

	private Function myFunction;
	private double myXinit;
	private double myDelta;
	private int myLength;

// Exported constructors.

	/**
	 * Construct a new sampled series sampling values of the given function. The
	 * <I>y</I> values at indexes 0, 1, 2, ... <TT>len-1</TT> are computed by
	 * evaluating <TT>theFunction</TT> at the <I>x</I> values <TT>xinit</TT>,
	 * <TT>xinit+delta</TT>, <TT>xinit+2*delta</TT>, ...
	 * <TT>xinit+(len-1)*delta</TT>.
	 *
	 * @param  theFunction  Function to sample.
	 * @param  xinit        Initial <I>x</I> value.
	 * @param  delta        Difference between successive <I>x</I> values.
	 * @param  len          Length (number of (<I>x,y</I>) pairs).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFunction</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>len</TT> is less than 0.
	 */
	public SampledXYSeries
		(Function theFunction,
		 double xinit,
		 double delta,
		 int len)
		{
		super();
		if (theFunction == null)
			{
			throw new NullPointerException
				("SampledXYSeries(): theFunction is null");
			}
		if (len < 0)
			{
			throw new IllegalArgumentException
				("SampledXYSeries(): len = " + len + " illegal");
			}
		myFunction = theFunction;
		myXinit = xinit;
		myDelta = delta;
		myLength = len;
		}

// Exported operations.

	/**
	 * Returns the number of values in this series.
	 *
	 * @return  Length.
	 */
	public int length()
		{
		return myLength;
		}

	/**
	 * Returns the given X value in this series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The X value in this series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public double x
		(int i)
		{
		return myXinit + i * myDelta;
		}

	/**
	 * Returns the given Y value in this series.
	 *
	 * @param  i  Index.
	 *
	 * @return  The Y value in this series at index <TT>i</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range
	 *     <TT>0</TT> .. <TT>length()-1</TT>.
	 */
	public double y
		(int i)
		{
		return myFunction.f (x (i));
		}

	}
