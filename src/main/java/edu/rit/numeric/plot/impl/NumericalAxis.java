//******************************************************************************
//
// File:    NumericalAxis.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.NumericalAxis
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

import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;

import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;

/**
 * Class NumericalAxis is the abstract base class for a numerical axis on a
 * plot. A numerical axis represents a numerical quantity; it is divided into
 * major divisions; each major division is divided into minor divisions. A
 * numerical axis may have tick marks at the major and minor divisions. A
 * numerical axis may have tick labels at the major divisions.
 *
 * @author  Alan Kaminsky
 * @version 18-Oct-2007
 */
public abstract class NumericalAxis
	extends Axis
	{

// Hidden data members.

	/**
	 * Ticks below or to the left of the axis.
	 */
	protected Ticks myTicksBelowOrLeft;

	/**
	 * Ticks above or to the right of the axis.
	 */
	protected Ticks myTicksAboveOrRight;

// Hidden constructors.

	/**
	 * Construct a new numerical axis. The axis is drawn with the given length,
	 * stroke, paint, and ticks.
	 *
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
	 *     (unchecked exception) Thrown if <TT>theLength</TT> is less than or
	 *     equal to 0.
	 */
	protected NumericalAxis
		(double theLength,
		 Stroke theStroke,
		 Paint thePaint,
		 Ticks theTicksBelowOrLeft,
		 Ticks theTicksAboveOrRight)
		{
		super (theLength, theStroke, thePaint);
		myTicksBelowOrLeft = theTicksBelowOrLeft;
		myTicksAboveOrRight = theTicksAboveOrRight;
		}

// Exported operations.

	/**
	 * Returns the ticks below or to the left of this axis.
	 */
	public Ticks getTicksBelowOrLeft()
		{
		return myTicksBelowOrLeft;
		}

	/**
	 * Returns the ticks above or to the right of this axis.
	 */
	public Ticks getTicksAboveOrRight()
		{
		return myTicksAboveOrRight;
		}

	/**
	 * Draw this axis in the given graphics context starting at coordinates (0,
	 * 0) and moving horizontally to the right.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void drawHorizontal
		(Graphics2D g2d)
		{
		// Save old settings.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		AffineTransform oldTransform = g2d.getTransform();

		// Draw tick marks if any.
		if (myTicksBelowOrLeft != null || myTicksAboveOrRight != null)
			{
			int major = getMajorDivisionCount();
			int minor = getMinorDivisionCount();
			int minorPerMajor = minor / major;
			for (int i = 0; i <= minor; ++ i)
				{
				double x = getMinorDivision (i);
				g2d.translate (getDisplayDistance (x), 0);
				if (myTicksBelowOrLeft != null)
					{
					if (i % minorPerMajor == 0)
						{
						myTicksBelowOrLeft.drawBelow (g2d, x);
						}
					else
						{
						myTicksBelowOrLeft.drawBelow (g2d);
						}
					}
				if (myTicksAboveOrRight != null)
					{
					if (i % minorPerMajor == 0)
						{
						myTicksAboveOrRight.drawAbove (g2d, x);
						}
					else
						{
						myTicksAboveOrRight.drawAbove (g2d);
						}
					}
				g2d.setTransform (oldTransform);
				}
			}

		// Draw main axis.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, myLength, 0));

		// Restore old settings.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Draw this axis in the given graphics context starting at coordinates (0,
	 * 0) and moving vertically up.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void drawVertical
		(Graphics2D g2d)
		{
		// Save old settings.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		AffineTransform oldTransform = g2d.getTransform();

		// Draw tick marks if any.
		if (myTicksBelowOrLeft != null || myTicksAboveOrRight != null)
			{
			int major = getMajorDivisionCount();
			int minor = getMinorDivisionCount();
			int minorPerMajor = minor / major;
			for (int i = 0; i <= minor; ++ i)
				{
				double y = getMinorDivision (i);
				g2d.translate (0, -getDisplayDistance (y));
				if (myTicksBelowOrLeft != null)
					{
					if (i % minorPerMajor == 0)
						{
						myTicksBelowOrLeft.drawLeft (g2d, y);
						}
					else
						{
						myTicksBelowOrLeft.drawLeft (g2d);
						}
					}
				if (myTicksAboveOrRight != null)
					{
					if (i % minorPerMajor == 0)
						{
						myTicksAboveOrRight.drawRight (g2d, y);
						}
					else
						{
						myTicksAboveOrRight.drawRight (g2d);
						}
					}
				g2d.setTransform (oldTransform);
				}
			}

		// Draw main axis.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, 0, - myLength));

		// Restore old settings.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Returns true if the given value falls within the bounds of this axis,
	 * false otherwise.
	 */
	public abstract boolean includesValue
		(double value);

	/**
	 * Returns the distance on the display from the start of this axis to the
	 * given value.
	 */
	public abstract double getDisplayDistance
		(double value);

	/**
	 * Returns the starting value of this axis. The starting value is the value
	 * at the left or bottom end of the axis for a horizontal or vertical
	 * orientation, respectively.
	 */
	public abstract double getStart();

	/**
	 * Returns the ending value of this axis. The ending value is the value at
	 * the right or top end of the axis for a horizontal or vertical
	 * orientation, respectively. Note that the ending value may be greater than
	 * or less than the starting value.
	 */
	public abstract double getEnd();

	/**
	 * Returns the number of major divisions on this axis.
	 */
	public abstract int getMajorDivisionCount();

	/**
	 * Returns the value of the given major division.
	 *
	 * @param  i  Major division index in the range 0 ..
	 *            <TT>getMajorDivisionCount()</TT>.
	 */
	public abstract double getMajorDivision
		(int i);

	/**
	 * Returns the number of minor divisions on this axis.
	 */
	public abstract int getMinorDivisionCount();

	/**
	 * Returns the value of the given minor division.
	 *
	 * @param  i  Minor division index in the range 0 ..
	 *            <TT>getMinorDivisionCount()</TT>.
	 */
	public abstract double getMinorDivision
		(int i);

	/**
	 * Returns the value of the point where the perpendicular axis crosses this
	 * axis.
	 */
	public abstract double getCrossing();

	}
