//******************************************************************************
//
// File:    LinearAxis.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.LinearAxis
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
 * Class LinearAxis provides a linear axis on a plot. The axis goes from a
 * starting value to an ending value in a linear fashion.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2007
 */
public class LinearAxis
	extends NumericalAxis
	{

// Hidden data members.

	/**
	 * Starting value.
	 */
	private double myStart;

	/**
	 * Ending value.
	 */
	private double myEnd;

	/**
	 * min (myStart, myEnd).
	 */
	private double myMin;

	/**
	 * max (myStart, myEnd).
	 */
	private double myMax;

	/**
	 * Ending value - starting value.
	 */
	private double myRange;

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
	 * Number of minor divisions.
	 */
	private int myMinorDivisionCount;

	/**
	 * Range of a minor division.
	 */
	private double myMinorDivisionRange;

	/**
	 * Point where the perpendicular axis crosses this axis.
	 */
	private double myCrossing;

// Exported constructors.

	/**
	 * Construct a new linear axis. The axis is drawn with the default stroke
	 * (solid, width=1) and paint (black). There are no tick marks.
	 *
	 * @param  theStart
	 *     Starting value for this axis.
	 * @param  theEnd
	 *     Ending value for this axis. It can be greater than or less than
	 *     <TT>theStart</TT>.
	 * @param  theMajorDivisionCount
	 *     Number of major divisions for this axis.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Point at which the perpendicular axis crosses this axis.
	 * @param  theLength
	 *     Length of this axis on the display.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theStart</TT> equals
	 *     <TT>theEnd</TT>. Thrown if <TT>theMajorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theMinorDivisionCount</TT> is less than 1.
	 *     Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LinearAxis
		(double theStart,
		 double theEnd,
		 int theMajorDivisionCount,
		 int theMinorDivisionCount,
		 double theCrossing,
		 double theLength)
		{
		this
			(theStart, theEnd,
			 theMajorDivisionCount, theMinorDivisionCount,
			 theCrossing, theLength,
			 DEFAULT_STROKE, DEFAULT_PAINT,
			 null, null);
		}

	/**
	 * Construct a new linear axis. The axis is drawn with the given stroke and
	 * paint. There are no tick marks.
	 *
	 * @param  theStart
	 *     Starting value for this axis.
	 * @param  theEnd
	 *     Ending value for this axis. It can be greater than or less than
	 *     <TT>theStart</TT>.
	 * @param  theMajorDivisionCount
	 *     Number of major divisions for this axis.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Point at which the perpendicular axis crosses this axis.
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
	 *     <TT>theEnd</TT>. Thrown if <TT>theMajorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theMinorDivisionCount</TT> is less than 1.
	 *     Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LinearAxis
		(double theStart,
		 double theEnd,
		 int theMajorDivisionCount,
		 int theMinorDivisionCount,
		 double theCrossing,
		 double theLength,
		 Stroke theStroke,
		 Paint thePaint)
		{
		this
			(theStart, theEnd,
			 theMajorDivisionCount, theMinorDivisionCount,
			 theCrossing, theLength,
			 theStroke, thePaint,
			 null, null);
		}

	/**
	 * Construct a new linear axis. The axis is drawn with the default stroke
	 * (solid, width=1) and paint (black). There are tick marks on either or
	 * both sides of the axis.
	 *
	 * @param  theStart
	 *     Starting value for this axis.
	 * @param  theEnd
	 *     Ending value for this axis. It can be greater than or less than
	 *     <TT>theStart</TT>.
	 * @param  theMajorDivisionCount
	 *     Number of major divisions for this axis.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Point at which the perpendicular axis crosses this axis.
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
	 *     <TT>theEnd</TT>. Thrown if <TT>theMajorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theMinorDivisionCount</TT> is less than 1.
	 *     Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LinearAxis
		(double theStart,
		 double theEnd,
		 int theMajorDivisionCount,
		 int theMinorDivisionCount,
		 double theCrossing,
		 double theLength,
		 Ticks theTicksBelowOrLeft,
		 Ticks theTicksAboveOrRight)
		{
		this
			(theStart, theEnd,
			 theMajorDivisionCount, theMinorDivisionCount,
			 theCrossing, theLength,
			 DEFAULT_STROKE, DEFAULT_PAINT,
			 theTicksBelowOrLeft, theTicksAboveOrRight);
		}

	/**
	 * Construct a new linear axis. The axis is drawn with the given stroke and
	 * paint. There are tick marks on either or both sides of the axis.
	 *
	 * @param  theStart
	 *     Starting value for this axis.
	 * @param  theEnd
	 *     Ending value for this axis. It can be greater than or less than
	 *     <TT>theStart</TT>.
	 * @param  theMajorDivisionCount
	 *     Number of major divisions for this axis.
	 * @param  theMinorDivisionCount
	 *     Number of minor divisions <I>per major division</I> for this axis.
	 * @param  theCrossing
	 *     Point at which the perpendicular axis crosses this axis.
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
	 *     <TT>theEnd</TT>. Thrown if <TT>theMajorDivisionCount</TT> is less
	 *     than 1. Thrown if <TT>theMinorDivisionCount</TT> is less than 1.
	 *     Thrown if <TT>theCrossing</TT> does not fall between
	 *     <TT>theStart</TT> and <TT>theEnd</TT> inclusive. Thrown if
	 *     <TT>theLength</TT> is less than or equal to 0.
	 */
	public LinearAxis
		(double theStart,
		 double theEnd,
		 int theMajorDivisionCount,
		 int theMinorDivisionCount,
		 double theCrossing,
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
			 theMajorDivisionCount < 1 ||
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
		myMajorDivisionCount = theMajorDivisionCount;
		myMajorDivisionRange = myRange / theMajorDivisionCount;
		myMinorDivisionCount = theMajorDivisionCount * theMinorDivisionCount;
		myMinorDivisionRange = myRange / myMinorDivisionCount;
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
		return myStart;
		}

	/**
	 * Returns the ending value of this axis. The ending value is the value at
	 * the right or top end of the axis for a horizontal or vertical
	 * orientation, respectively. Note that the ending value may be greater than
	 * or less than the starting value.
	 */
	public double getEnd()
		{
		return myEnd;
		}

	/**
	 * Returns true if the given value falls within the bounds of this axis,
	 * false otherwise.
	 */
	public boolean includesValue
		(double value)
		{
		return myMin <= value && value <= myMax;
		}

	/**
	 * Returns the distance on the display from the start of this axis to the
	 * given value.
	 */
	public double getDisplayDistance
		(double value)
		{
		return (value - myStart) * myLengthOverMyRange;
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
		return myStart + i * myMajorDivisionRange;
		}

	/**
	 * Returns the number of minor divisions on this axis.
	 */
	public int getMinorDivisionCount()
		{
		return myMinorDivisionCount;
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
		return myStart + i * myMinorDivisionRange;
		}

	/**
	 * Returns the value of the point where the perpendicular axis crosses this
	 * axis.
	 */
	public double getCrossing()
		{
		return myCrossing;
		}

	/**
	 * Determine a pleasing scale for a linear axis. The input, <TT>max</TT>, is
	 * the largest value to be displayed along the axis. The output is the
	 * smallest number of the form 1x10<SUP>k</SUP>, 2x10<SUP>k</SUP>, or
	 * 5x10<SUP>k</SUP>, k an integer, that is greater than or equal to
	 * <TT>max</TT>. This can then be used as the maximum value for the axis.
	 * <P>
	 * <I>Note:</I> The <TT>autoscale()</TT> method works for negative inputs.
	 * If <TT>max</TT> &lt; 0, then <TT>autoscale(max) = -autoscale(-max)</TT>.
	 * <P>
	 * <I>Note:</I> <TT>autoscale(0) = 0</TT>.
	 *
	 * @param  max  Maximum value to be displayed along the axis.
	 *
	 * @return  Maximum value for the axis.
	 */
	public static double autoscale
		(double max)
		{
		if (max == 0.0)
			{
			return 0.0;
			}
		else if (max < 0.0)
			{
			return -autoscale (-max);
			}
		else
			{
			double tentothek = 1.0;
			double twotimestentothek = 2.0;
			double fivetimestentothek = 5.0;
			double tentimestentothek = 10.0;
			while (tentimestentothek < max)
				{
				tentothek *= 10.0;
				twotimestentothek *= 10.0;
				fivetimestentothek *= 10.0;
				tentimestentothek *= 10.0;
				}
			for (;;)
				{
				if (fivetimestentothek < max && max <= tentimestentothek)
					{
					return tentimestentothek;
					}
				else if (twotimestentothek < max && max <= fivetimestentothek)
					{
					return fivetimestentothek;
					}
				else if (tentothek < max && max <= twotimestentothek)
					{
					return twotimestentothek;
					}
				else
					{
					tentothek /= 10.0;
					twotimestentothek /= 10.0;
					fivetimestentothek /= 10.0;
					tentimestentothek /= 10.0;
					}
				}
			}
		}

	}
