//******************************************************************************
//
// File:    Interpolation.java
// Package: edu.rit.numeric
// Unit:    Class edu.rit.numeric.Interpolation
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

import edu.rit.util.Sorting;

/**
 * Class Interpolation provides an object for interpolating in an {@linkplain
 * XYSeries} of real values (type <TT>double</TT>). Linear interpolation is
 * used. The X-Y series must have at least two elements; the X values must be
 * distinct; but the X values need not be in any particular order. When doing
 * interpolations, the (X,Y) pairs are arranged in ascending order of X values.
 * <P>
 * Class Interpolation implements interface {@linkplain Function}. An instance
 * of class Interpolation can be used as a function object.
 *
 * @author  Alan Kaminsky
 * @version 18-Aug-2011
 */
public class Interpolation
	implements Function
	{

// Hidden data members.

	// X and Y values in which to interpolate. X values in ascending order.
	private double[] xData;
	private double[] yData;

	// Number of data values, minus 2.
	private int NM2;

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
	 * in the given X-Y series. The X-Y series must have at least two elements;
	 * the X values must be distinct; but the X values need not be in any
	 * particular order.
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
	 *     elements. Thrown if the X values in <TT>theSeries</TT> are not
	 *     distinct.
	 */
	public Interpolation
		(XYSeries theSeries)
		{
		int N = theSeries.length();
		if (N < 2)
			{
			throw new IllegalArgumentException
				("Interpolation(): theSeries length < 2");
			}
		xData = new double [N];
		yData = new double [N];
		for (int i = 0; i < N; ++ i)
			{
			xData[i] = theSeries.x(i);
			yData[i] = theSeries.y(i);
			}
		Sorting.sort (xData, new Sorting.Double()
			{
			public void swap (double[] x, int a, int b)
				{
				double tmp;
				tmp = xData[a];
				xData[a] = xData[b];
				xData[b] = tmp;
				tmp = yData[a];
				yData[a] = yData[b];
				yData[b] = tmp;
				}
			});
		NM2 = N - 2;
		for (int i = 0; i <= NM2; ++ i)
			{
			if (xData[i] == xData[i+1])
				{
				throw new IllegalArgumentException
					("Interpolation(): Duplicate X value: "+xData[i]);
				}
			}
		myIndex = 0;
		x1 = xData[0];
		y1 = yData[0];
		x2 = xData[1];
		y2 = yData[1];
		}

// Exported operations.

	/**
	 * Using linear interpolation, compute the Y value for the given X value.
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
		while (myIndex < NM2 && x >= x2)
			{
			++ myIndex;
			x1 = x2;
			y1 = y2;
			x2 = xData[myIndex + 1];
			y2 = yData[myIndex + 1];
			}

		// Scan backward if necessary to find the correct interval for x.
		while (myIndex > 0 && x < x1)
			{
			-- myIndex;
			x2 = x1;
			y2 = y1;
			x1 = xData[myIndex];
			y1 = yData[myIndex];
			}

		// Interpolate on x.
		double dx = (x - x1)/(x2 - x1);
		return (1.0 - dx)*y1 + dx*y2;
		}

	}
