//******************************************************************************
//
// File:    LogarithmicAxis.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.LogarithmicAxis
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

package edu.rit.numeric.plot.impl;

import java.awt.Stroke;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.Paint;

import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

/**
 * Class LogarithmicAxis provides a logarithmic axis on a plot. The axis goes
 * from a starting value to an ending value in a logarithmic fashion. The
 * starting and ending values must be integer powers of 10. Each major division
 * corresponds to one decade (factor of 10). Each major division may be
 * subdivided into any number of minor divisions; the minor divisions are
 * distributed linearly between the upper and lower bounds of the major
 * division.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2007
 */
public class LogarithmicAxis
	extends NumericalAxis
	{

// Hidden data members.

	/**
	 * Base-10 logarithm of starting value.
	 */
	private int myStart;

	/**
	 * Base-10 logarithm of ending value.
	 */
	private int myEnd;

	/**
	 * min (myStart, myEnd).
	 */
	private int myMin;

	/**
	 * max (myStart, myEnd).
	 */
	private int myMax;

	/**
	 * Ending value - starting value.
	 */
	private int myRange;

	/**
	 * myLength / myRange.
	 */
	private double myLengthOverMyRange;

	/**
	 * Number of major divisions.
	 */
	private int myMajorDivisionCount;

	/**
	 * Range of a major division.
	 */
	private double myMajorDivisionRange;

	/**
	 * Number of minor divisions per major division.
	 */
	private int myMinorDivisionCount;

	/**
	 * Point where the perpendicular axis crosses this axis.
	 */
	private int myCrossing;

// Exported constructors.

	/**
	 * Construct a new logarithmic axis. The axis is drawn with the default
	 * stroke (solid, width=1) and paint (black). There are no tick marks.
	 *
	 * @param  theStart
	 *     Base-10 logarithm of the starting value for this axis.
	 * @param  theEnd
	 *     Base-10 logarithm of the ending value for this axis. It can be
	 *     greater than or less than <TT>theStart</TT>.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Base-10 logarithm of the point at which the perpendicular axis
	 *     crosses this axis.
	 * @param  theLength
	 *     Length of this axis on the display.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theStart</TT> equals
	 *     <TT>theEnd</TT>. Thrown if <TT>theMinorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LogarithmicAxis
		(int theStart,
		 int theEnd,
		 int theMinorDivisionCount,
		 int theCrossing,
		 double theLength)
		{
		this
			(theStart, theEnd,
			 theMinorDivisionCount,
			 theCrossing, theLength,
			 DEFAULT_STROKE, DEFAULT_PAINT,
			 null, null);
		}

	/**
	 * Construct a new logarithmic axis. The axis is drawn with the given stroke
	 * and paint. There are no tick marks.
	 *
	 * @param  theStart
	 *     Base-10 logarithm of the starting value for this axis.
	 * @param  theEnd
	 *     Base-10 logarithm of the ending value for this axis. It can be
	 *     greater than or less than <TT>theStart</TT>.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Base-10 logarithm of the point at which the perpendicular axis
	 *     crosses this axis.
	 * @param  theLength
	 *     Length of this axis on the display.
	 * @param  theStroke
	 *     Stroke for drawing this axis on the display.
	 * @param  thePaint
	 *     Paint for drawing this axis on the display.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStroke</TT> is null or
	 *     <TT>thePaint</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theStart</TT> equals
	 *     <TT>theEnd</TT>. Thrown if <TT>theMinorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LogarithmicAxis
		(int theStart,
		 int theEnd,
		 int theMinorDivisionCount,
		 int theCrossing,
		 double theLength,
		 Stroke theStroke,
		 Paint thePaint)
		{
		this
			(theStart, theEnd,
			 theMinorDivisionCount,
			 theCrossing, theLength,
			 theStroke, thePaint,
			 null, null);
		}

	/**
	 * Construct a new logarithmic axis. The axis is drawn with the default
	 * stroke (solid, width=1) and paint (black). There are tick marks on either
	 * or both sides of the axis.
	 *
	 * @param  theStart
	 *     Base-10 logarithm of the starting value for this axis.
	 * @param  theEnd
	 *     Base-10 logarithm of the ending value for this axis. It can be
	 *     greater than or less than <TT>theStart</TT>.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Base-10 logarithm of the point at which the perpendicular axis
	 *     crosses this axis.
	 * @param  theLength
	 *     Length of this axis on the display.
	 * @param  theTicksBelowOrLeft
	 *     Ticks to draw below or to the left of the axis (for a horizontal or
	 *     vertical axis, respectively). If null, no ticks are drawn below or to
	 *     the left.
	 * @param  theTicksAboveOrRight
	 *     Ticks to draw above or to the right of the axis (for a horizontal or
	 *     vertical axis, respectively). If null, no ticks are drawn above or to
	 *     the right.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theStart</TT> equals
	 *     <TT>theEnd</TT>. Thrown if <TT>theMinorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LogarithmicAxis
		(int theStart,
		 int theEnd,
		 int theMinorDivisionCount,
		 int theCrossing,
		 double theLength,
		 Ticks theTicksBelowOrLeft,
		 Ticks theTicksAboveOrRight)
		{
		this
			(theStart, theEnd,
			 theMinorDivisionCount,
			 theCrossing, theLength,
			 DEFAULT_STROKE, DEFAULT_PAINT,
			 theTicksBelowOrLeft, theTicksAboveOrRight);
		}

	/**
	 * Construct a new logarithmic axis. The axis is drawn with the given stroke
	 * and paint. There are tick marks on either or both sides of the axis.
	 *
	 * @param  theStart
	 *     Base-10 logarithm of the starting value for this axis.
	 * @param  theEnd
	 *     Base-10 logarithm of the ending value for this axis. It can be
	 *     greater than or less than <TT>theStart</TT>.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Base-10 logarithm of the point at which the perpendicular axis
	 *     crosses this axis.
	 * @param  theLength
	 *     Length of this axis on the display.
	 * @param  theStroke
	 *     Stroke for drawing this axis on the display.
	 * @param  thePaint
	 *     Paint for drawing this axis on the display.
	 * @param  theTicksBelowOrLeft
	 *     Ticks to draw below or to the left of the axis (for a horizontal or
	 *     vertical axis, respectively). If null, no ticks are drawn below or to
	 *     the left.
	 * @param  theTicksAboveOrRight
	 *     Ticks to draw above or to the right of the axis (for a horizontal or
	 *     vertical axis, respectively). If null, no ticks are drawn above or to
	 *     the right.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStroke</TT> is null or
	 *     <TT>thePaint</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theStart</TT> equals
	 *     <TT>theEnd</TT>. Thrown if <TT>theMinorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LogarithmicAxis
		(int theStart,
		 int theEnd,
		 int theMinorDivisionCount,
		 int theCrossing,
		 double theLength,
		 Stroke theStroke,
		 Paint thePaint,
		 Ticks theTicksBelowOrLeft,
		 Ticks theTicksAboveOrRight)
		{
		super
			(theLength, theStroke, thePaint,
			 theTicksBelowOrLeft, theTicksAboveOrRight);
		if
			(theStart == theEnd ||
			 theMinorDivisionCount < 1 ||
			 Math.min (theStart, theEnd) > theCrossing ||
			 theCrossing > Math.max (theStart, theEnd))
			{
			throw new IllegalArgumentException();
			}
		myStart = theStart;
		myEnd = theEnd;
		myMin = Math.min (theStart, theEnd);
		myMax = Math.max (theStart, theEnd);
		myRange = theEnd - theStart;
		myLengthOverMyRange = myLength / myRange;
		myMajorDivisionCount = myMax - myMin;
		myMajorDivisionRange = myRange / myMajorDivisionCount;
		myMinorDivisionCount = theMinorDivisionCount;
		myCrossing = theCrossing;
		}

// Exported operations.

	/**
	 * Returns the starting value of this axis. The starting value is the value
	 * at the left or bottom end of the axis for a horizontal or vertical
	 * orientation, respectively.
	 */
	public double getStart()
		{
		return antilog10 (myStart);
		}

	/**
	 * Returns the ending value of this axis. The ending value is the value at
	 * the right or top end of the axis for a horizontal or vertical
	 * orientation, respectively. Note that the ending value may be greater than
	 * or less than the starting value.
	 */
	public double getEnd()
		{
		return antilog10 (myEnd);
		}

	/**
	 * Returns true if the given value falls within the bounds of this axis,
	 * false otherwise.
	 */
	public boolean includesValue
		(double value)
		{
		double logvalue = log10 (value);
		return myMin <= logvalue && logvalue <= myMax;
		}

	/**
	 * Returns the distance on the display from the start of this axis to the
	 * given value.
	 */
	public double getDisplayDistance
		(double value)
		{
		double logvalue = log10 (value);
		if (logvalue == Double.NEGATIVE_INFINITY)
			{
			return -10.0 * myLengthOverMyRange;
			}
		else
			{
			return (logvalue - myStart) * myLengthOverMyRange;
			}
		}

	/**
	 * Returns the number of major divisions on this axis.
	 */
	public int getMajorDivisionCount()
		{
		return myMajorDivisionCount;
		}

	/**
	 * Returns the value of the given major division.
	 *
	 * @param  i  Major division index in the range 0 ..
	 *            <TT>getMajorDivisionCount()</TT>.
	 */
	public double getMajorDivision
		(int i)
		{
		return antilog10 (myStart + i * myMajorDivisionRange);
		}

	/**
	 * Returns the number of minor divisions on this axis.
	 */
	public int getMinorDivisionCount()
		{
		return myMajorDivisionCount * myMinorDivisionCount;
		}

	/**
	 * Returns the value of the given minor division.
	 *
	 * @param  i  Minor division index in the range 0 ..
	 *            <TT>getMinorDivisionCount()</TT>.
	 */
	public double getMinorDivision
		(int i)
		{
		int decade = i / myMinorDivisionCount;
		int div = i % myMinorDivisionCount;
		double lb = getMajorDivision (decade);
		return
			div == 0 ?
				lb :
				10.0 * lb / myMinorDivisionCount * div;
		}

	/**
	 * Returns the value of the point where the perpendicular axis crosses this
	 * axis.
	 */
	public double getCrossing()
		{
		return antilog10 (myCrossing);
		}

	/**
	 * Determine a scale for a logarithmic axis. The inputs, <TT>min</TT> and
	 * <TT>max</TT>, are the smallest and largest values to be displayed along
	 * the axis; both must be greater than 0. The output is an array of two
	 * integers containing the smallest and the largest exponent of 10 that
	 * encompasses the range <TT>min</TT> through <TT>max</TT>. These values can
	 * then be used to construct an instance of class LogarithmicAxis.
	 *
	 * @param  min  Minimum value to be displayed along the axis.
	 * @param  max  Maximum value to be displayed along the axis.
	 *
	 * @return  Array of minimum and maximum exponents of 10 for the axis.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>min</TT> &lt;= 0 or <TT>max</TT>
	 *     &lt;= 0.
	 */
	public static int[] autoscale
		(double min,
		 double max)
		{
		if (min <= 0.0 || max <= 0.0)
			{
			throw new IllegalArgumentException();
			}
		return new int[]
			{(int) (Math.floor (log10 (min))),
			 (int) (Math.ceil  (log10 (max)))};
		}

// Hidden operations.

	/**
	 * Returns the base-10 logarithm of x.
	 */
	private static double log10
		(double x)
		{
		return Math.log (x) / LOG10;
		}

	/**
	 * Returns the base-10 antilogarithm of x, that is, 10**x.
	 */
	private static double antilog10
		(double x)
		{
		return Math.exp (x * LOG10);
		}

	private static final double LOG10 = Math.log (10.0);

	}
