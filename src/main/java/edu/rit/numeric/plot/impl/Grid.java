//******************************************************************************
//
// File:    Grid.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.Grid
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
import java.awt.Stroke;
import java.awt.Paint;

import java.awt.geom.Line2D;

/**
 * Class Grid provides a grid for a plot. The grid specifies the stroke and
 * paint to use for the grid lines, which are drawn at each major and minor
 * division along each axis.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2007
 */
public class Grid
	{

// Exported constants.

	/**
	 * The default gridline stroke (solid, width=0.2).
	 */
	public static final Stroke DEFAULT_STROKE = Strokes.solid (0.2);

	/**
	 * The default gridline paint (black).
	 */
	public static final Paint DEFAULT_PAINT = Color.black;

// Hidden data members.

	/**
	 * Gridline stroke and paint for major divisions.
	 */
	private Stroke myMajorStroke;
	private Paint myMajorPaint;

	/**
	 * Gridline stroke and paint for minor divisions.
	 */
	private Stroke myMinorStroke;
	private Paint myMinorPaint;

// Exported constructors.

	/**
	 * Construct a new grid that uses the default gridline stroke and paint for
	 * the major divisions. Minor divisions do not appear in the grid.
	 */
	public Grid()
		{
		this (DEFAULT_STROKE, DEFAULT_PAINT, null, null);
		}

	/**
	 * Construct a new grid that uses the given gridline stroke and paint for
	 * the major divisions. Minor divisions do not appear in the grid.
	 *
	 * @param  theMajorStroke
	 *     Gridline stroke for major divisions, or null not to draw gridlines
	 *     for major divisions.
	 * @param  theMajorPaint
	 *     Gridline paint for major divisions, or null not to draw gridlines
	 *     for major divisions.
	 */
	public Grid
		(Stroke theMajorStroke,
		 Paint theMajorPaint)
		{
		this (theMajorStroke, theMajorPaint, null, null);
		}

	/**
	 * Construct a new grid that uses the given gridline stroke and paint for
	 * the major and minor divisions.
	 *
	 * @param  theMajorStroke
	 *     Gridline stroke for major divisions, or null not to draw gridlines
	 *     for major divisions.
	 * @param  theMajorPaint
	 *     Gridline paint for major divisions, or null not to draw gridlines
	 *     for major divisions.
	 * @param  theMinorStroke
	 *     Gridline stroke for minor divisions, or null not to draw gridlines
	 *     for minor divisions.
	 * @param  theMinorPaint
	 *     Gridline paint for minor divisions, or null not to draw gridlines
	 *     for minor divisions.
	 */
	public Grid
		(Stroke theMajorStroke,
		 Paint theMajorPaint,
		 Stroke theMinorStroke,
		 Paint theMinorPaint)
		{
		myMajorStroke = theMajorStroke;
		myMajorPaint = theMajorPaint;
		myMinorStroke = theMinorStroke;
		myMinorPaint = theMinorPaint;
		}

// Hidden operations.

	/**
	 * Draw horizontal gridlines at the divisions of the given axis. The axis
	 * starts at <I>y</I> = 0. The gridlines start at <I>x</I> = 0 and are the
	 * given length.
	 *
	 * @param  g2d   Graphics context.
	 * @param  axis  Axis.
	 * @param  len   Gridline length.
	 */
	void drawHorizontalGridlines
		(Graphics2D g2d,
		 NumericalAxis axis,
		 double len)
		{
		boolean drawMajor = myMajorStroke != null && myMajorPaint != null;
		boolean drawMinor = myMinorStroke != null && myMinorPaint != null;

		// Early return if no gridlines.
		if (! drawMajor && ! drawMinor) return;

		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();

		// Do minor divisions.
		if (drawMinor)
			{
			// Set stroke and paint.
			g2d.setStroke (myMinorStroke);
			g2d.setPaint (myMinorPaint);

			// Draw gridlines at the axis's minor divisions.
			int n = axis.getMinorDivisionCount();
			for (int i = 0; i <= n; ++ i)
				{
				double y = axis.getDisplayDistance (axis.getMinorDivision (i));
				g2d.draw (new Line2D.Double (0, -y, len, -y));
				}
			}

		// Do major divisions.
		if (drawMajor)
			{
			// Set stroke and paint.
			g2d.setStroke (myMajorStroke);
			g2d.setPaint (myMajorPaint);

			// Draw gridlines at the axis's minor divisions.
			int n = axis.getMajorDivisionCount();
			for (int i = 0; i <= n; ++ i)
				{
				double y = axis.getDisplayDistance (axis.getMajorDivision (i));
				g2d.draw (new Line2D.Double (0, -y, len, -y));
				}
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Draw vertical gridlines at the divisions of the given axis. The axis
	 * starts at <I>x</I> = 0. The gridlines start at <I>y</I> = 0 and are the
	 * given length.
	 *
	 * @param  g2d   Graphics context.
	 * @param  axis  Axis.
	 * @param  len   Gridline length.
	 */
	void drawVerticalGridlines
		(Graphics2D g2d,
		 NumericalAxis axis,
		 double len)
		{
		boolean drawMajor = myMajorStroke != null && myMajorPaint != null;
		boolean drawMinor = myMinorStroke != null && myMinorPaint != null;

		// Early return if no gridlines.
		if (! drawMajor && ! drawMinor) return;

		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();

		// Do minor divisions.
		if (drawMinor)
			{
			// Set stroke and paint.
			g2d.setStroke (myMinorStroke);
			g2d.setPaint (myMinorPaint);

			// Draw gridlines at the axis's major divisions.
			int n = axis.getMinorDivisionCount();
			for (int i = 0; i <= n; ++ i)
				{
				double x = axis.getDisplayDistance (axis.getMinorDivision (i));
				g2d.draw (new Line2D.Double (x, 0, x, -len));
				}
			}

		// Do major divisions.
		if (drawMajor)
			{
			// Set stroke and paint.
			g2d.setStroke (myMajorStroke);
			g2d.setPaint (myMajorPaint);

			// Draw gridlines at the axis's major divisions.
			int n = axis.getMajorDivisionCount();
			for (int i = 0; i <= n; ++ i)
				{
				double x = axis.getDisplayDistance (axis.getMajorDivision (i));
				g2d.draw (new Line2D.Double (x, 0, x, -len));
				}
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	}
