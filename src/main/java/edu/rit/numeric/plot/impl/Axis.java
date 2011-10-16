//******************************************************************************
//
// File:    Axis.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.Axis
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

import edu.rit.numeric.plot.Strokes;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;

import java.awt.geom.Rectangle2D;

/**
 * Class Axis is the abstract base class for an axis on a plot. It also provides
 * static methods to help determine the scale for an axis.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2007
 */
public abstract class Axis
	{

// Exported constants.

	/**
	 * The default axis stroke (solid, width=1).
	 */
	public static final Stroke DEFAULT_STROKE = Strokes.solid (1);

	/**
	 * The default axis paint (black).
	 */
	public static final Paint DEFAULT_PAINT = Color.black;

// Hidden data members.

	/**
	 * Length of this axis on the display.
	 */
	protected double myLength;

	/**
	 * Stroke for drawing this axis on the display.
	 */
	protected Stroke myStroke;

	/**
	 * Paint for drawing this axis on the display.
	 */
	protected Paint myPaint;

// Hidden constructors.

	/**
	 * Construct a new axis. The axis is drawn with the given length, stroke,
	 * and paint.
	 *
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
	 *     (unchecked exception) Thrown if <TT>theLength</TT> is less than or
	 *     equal to 0.
	 */
	protected Axis
		(double theLength,
		 Stroke theStroke,
		 Paint thePaint)
		{
		if (theStroke == null || thePaint == null)
			{
			throw new NullPointerException();
			}
		if (theLength <= 0.0)
			{
			throw new IllegalArgumentException();
			}
		myLength = theLength;
		myStroke = theStroke;
		myPaint = thePaint;
		}

// Exported operations.

	/**
	 * Returns the length of this axis on the display.
	 */
	public double getLength()
		{
		return myLength;
		}

	/**
	 * Returns the stroke for drawing this axis on the display.
	 */
	public Stroke getStroke()
		{
		return myStroke;
		}

	/**
	 * Returns the paint for drawing this axis on the display.
	 */
	public Paint getPaint()
		{
		return myPaint;
		}

	/**
	 * Draw this axis in the given graphics context starting at coordinates (0,
	 * 0) and moving horizontally to the right.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public abstract void drawHorizontal
		(Graphics2D g2d);

	/**
	 * Draw this axis in the given graphics context starting at coordinates (0,
	 * 0) and moving vertically up.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public abstract void drawVertical
		(Graphics2D g2d);

	}
