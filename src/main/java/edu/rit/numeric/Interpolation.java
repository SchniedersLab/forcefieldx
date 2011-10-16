//******************************************************************************
//
// File:    Interpolation.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Interpolation
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
 * Class Interpolation provides an object for interpolating in an {@linkplain
 * XYSeries} of real values (type <TT>double</TT>). Linear interpolation is
 * used.
 * <P>
 * Class Interpolation implements interface {@linkplain Function}. An instance
 * of class Interpolation can be used as a function object.
 *
 * @author  Alan Kaminsky
 * @version 22-Jul-2007
 */
public class Interpolation
	implements Function
	{

// Hidden data members.

	// Series of (x,y) pairs in which to interpolate.
	private ListXYSeries mySeries;

	// Length of mySeries, minus 2.
	private int myLengthMinus2;

	// Index of the lower (x,y) pair of the interval in which the last
	// interpolation occurred.
	private int myIndex;

	// Last interpolation interval was (x1,y1) .. (x2,y2).
	private double x1;
	private double y1;
	private double x2;
	private double y2;

// Exported constructors.

	/**
	 * Construct a new interpolation object that will interpolate between values
	 * in the given X-Y series. The X-Y series must have at least two elements.
	 * <P>
	 * <I>Note:</I> A copy of the given series' elements is made. Changing
	 * <TT>theSeries</TT> will not affect this interpolation object.
	 *
	 * @param  theSeries  X-Y series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theSeries</TT> has fewer than two
	 *     elements.
	 */
	public Interpolation
		(XYSeries theSeries)
		{
		if (theSeries.length() < 2)
			{
			throw new IllegalArgumentException
				("Interpolation(): theSeries length < 2");
			}
		mySeries = new ListXYSeries().add (theSeries);
		myLengthMinus2 = mySeries.length() - 2;
		myIndex = 0;
		x1 = mySeries.x(0);
		y1 = mySeries.y(0);
		x2 = mySeries.x(1);
		y2 = mySeries.y(1);
		}

// Exported operations.

	/**
	 * Using linear interpolation, compute the Y value for the given X value. It
	 * is assumed that the X values in the underlying X-Y series form a strictly
	 * increasing sequence.
	 * <P>
	 * If <TT>x</TT> is less than the smallest X value in the underlying X-Y
	 * series, the Y value is computed by extrapolating the first interval. If
	 * <TT>x</TT> is greater than the largest X value in the underlying X-Y
	 * series, the Y value is computed by extrapolating the last interval.
	 *
	 * @param  x  X value.
	 *
	 * @return  Interpolated or extrapolated Y value.
	 */
	public double f
		(double x)
		{
		// Scan forward if necessary to find the correct interval for x.
		while (myIndex < myLengthMinus2 && x >= x2)
			{
			++ myIndex;
			x1 = x2;
			y1 = y2;
			x2 = mySeries.x (myIndex + 1);
			y2 = mySeries.y (myIndex + 1);
			}

		// Scan backward if necessary to find the correct interval for x.
		while (myIndex > 0 && x < x1)
			{
			-- myIndex;
			x2 = x1;
			y2 = y1;
			x1 = mySeries.x (myIndex);
			y1 = mySeries.y (myIndex);
			}

		// Interpolate on x.
		double dx = (x - x1) / (x2 - x1);
		return (1.0 - dx) * y1 + dx * y2;
		}

	/**
	 * Using linear interpolation, compute the X value for the given Y value. It
	 * is assumed that the Y values in the underlying X-Y series form a strictly
	 * increasing sequence.
	 * <P>
	 * If <TT>y</TT> is less than the smallest Y value in the underlying X-Y
	 * series, the X value is computed by extrapolating the first interval. If
	 * <TT>y</TT> is greater than the largest Y value in the underlying X-Y
	 * series, the X value is computed by extrapolating the last interval.
	 *
	 * @param  y  Y value.
	 *
	 * @return  Interpolated or extrapolated X value.
	 */
	public double fInv
		(double y)
		{
		// Scan forward if necessary to find the correct interval for y.
		while (myIndex < myLengthMinus2 && y >= y2)
			{
			++ myIndex;
			x1 = x2;
			y1 = y2;
			x2 = mySeries.x (myIndex + 1);
			y2 = mySeries.y (myIndex + 1);
			}

		// Scan backward if necessary to find the correct interval for y.
		while (myIndex > 0 && y < y1)
			{
			-- myIndex;
			x2 = x1;
			y2 = y1;
			x1 = mySeries.x (myIndex);
			y1 = mySeries.y (myIndex);
			}

		// Interpolate on y.
		double dy = (y - y1) / (y2 - y1);
		return (1.0 - dy) * x1 + dy * x2;
		}

	}
