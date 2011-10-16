//******************************************************************************
//
// File:    XYPlotSeries.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.XYPlotSeries
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

import edu.rit.numeric.XYSeries;

import edu.rit.numeric.plot.Dots;

import java.awt.BasicStroke;
import java.awt.Color;

/**
 * Class XYPlotSeries provides a plot series for an {@linkplain XYPlot}. An
 * {@linkplain edu.rit.numeric.XYSeries XYSeries} specifies the data points to
 * be plotted. A dot is drawn at each data point in the series. A line segment
 * is drawn from each data point to the next data point in the series.
 * <P>
 * Class XYPlotSeries provides operations to tell where to plot data points as
 * well as which {@linkplain edu.rit.swing.Drawable Drawable} object to use to
 * plot the points. Class {@linkplain edu.rit.numeric.plot.Dots Dots} provides
 * several shapes of drawable objects for plotting data points.
 * <P>
 * Class XYPlotSeries provides operations to tell where to plot line segments as
 * well as which stroke and color objects to use to plot the line segments.
 * Class {@linkplain edu.rit.numeric.plot.Strokes Strokes} provides several
 * kinds of strokes for drawing lines.
 *
 * @author  Alan Kaminsky
 * @version 11-Dec-2007
 */
public class XYPlotSeries
	extends PlotSeries
	{

// Exported constructors.

	/**
	 * Construct a new XY plot series. This constructor is intended for use only
	 * by object deserialization.
	 */
	public XYPlotSeries()
		{
		super();
		}

	/**
	 * Construct a new XY plot series. The given data series contains the data
	 * points to be plotted. Dots are plotted on the data points using the given
	 * drawable object. Lines are not drawn between data points.
	 *
	 * @param  theDataSeries
	 *     Data series.
	 * @param  theDots
	 *     Drawable object for plotting data points. If null, dots are not
	 *     plotted on the data points.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDataSeries</TT> is null.
	 */
	public XYPlotSeries
		(XYSeries theDataSeries,
		 Dots theDots)
		{
		super (theDataSeries, theDots, null, null, false);
		}

	/**
	 * Construct a new XY plot series. The given data series contains the data
	 * points to be plotted. No dots are plotted on the data points. Lines are
	 * drawn between data points using the given stroke and color. The curve is
	 * not smooth.
	 *
	 * @param  theDataSeries
	 *     Data series.
	 * @param  theStroke
	 *     Stroke for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 * @param  theColor
	 *     Color for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDataSeries</TT> is null.
	 */
	public XYPlotSeries
		(XYSeries theDataSeries,
		 BasicStroke theStroke,
		 Color theColor)
		{
		super (theDataSeries, null, theStroke, theColor, false);
		}

	/**
	 * Construct a new XY plot series. The given data series contains the data
	 * points to be plotted. No dots are plotted on the data points. Lines are
	 * drawn between data points using the given stroke and color. The curve is
	 * smooth or not as specified by the <TT>isSmooth</TT> argument.
	 *
	 * @param  theDataSeries
	 *     Data series.
	 * @param  theStroke
	 *     Stroke for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 * @param  theColor
	 *     Color for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 * @param  isSmooth
	 *     True to plot a smooth curve, false to plot straight line segments.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDataSeries</TT> is null.
	 */
	public XYPlotSeries
		(XYSeries theDataSeries,
		 BasicStroke theStroke,
		 Color theColor,
		 boolean isSmooth)
		{
		super (theDataSeries, null, theStroke, theColor, isSmooth);
		}

	/**
	 * Construct a new XY plot series. The given data series contains the data
	 * points to be plotted. Dots are plotted on the data points using the given
	 * drawable object. Lines are drawn between data points using the given
	 * stroke and color. The curve is not smooth.
	 *
	 * @param  theDataSeries
	 *     Data series.
	 * @param  theDots
	 *     Drawable object for plotting data points. If null, dots are not
	 *     plotted on the data points.
	 * @param  theStroke
	 *     Stroke for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 * @param  theColor
	 *     Color for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDataSeries</TT> is null.
	 */
	public XYPlotSeries
		(XYSeries theDataSeries,
		 Dots theDots,
		 BasicStroke theStroke,
		 Color theColor)
		{
		super (theDataSeries, theDots, theStroke, theColor, false);
		}

	/**
	 * Construct a new XY plot series. The given data series contains the data
	 * points to be plotted. Dots are plotted on the data points using the given
	 * drawable object. Lines are drawn between data points using the given
	 * stroke and color. The curve is smooth or not as specified by the
	 * <TT>isSmooth</TT> argument.
	 *
	 * @param  theDataSeries
	 *     Data series.
	 * @param  theDots
	 *     Drawable object for plotting data points. If null, dots are not
	 *     plotted on the data points.
	 * @param  theStroke
	 *     Stroke for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 * @param  theColor
	 *     Color for drawing lines between data points. If null, lines are not
	 *     drawn between data points.
	 * @param  isSmooth
	 *     True to plot a smooth curve, false to plot straight line segments.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDataSeries</TT> is null.
	 */
	public XYPlotSeries
		(XYSeries theDataSeries,
		 Dots theDots,
		 BasicStroke theStroke,
		 Color theColor,
		 boolean isSmooth)
		{
		super (theDataSeries, theDots, theStroke, theColor, isSmooth);
		}

// Exported operations.

	/**
	 * Returns the number of dots to be plotted.
	 */
	public int getDotCount()
		{
		return myDataSeries.length();
		}

	/**
	 * Returns the X coordinate of the given dot.
	 *
	 * @param  i  Dot index.
	 *
	 * @return  X coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getDotCount()-1</TT>.
	 */
	public double getDotX
		(int i)
		{
		return myDataSeries.x (i);
		}

	/**
	 * Returns the Y coordinate of the given dot.
	 * <P>
	 * This method must be overridden in a subclass.
	 *
	 * @param  i  Dot index.
	 *
	 * @return  Y coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getDotCount()-1</TT>.
	 */
	public double getDotY
		(int i)
		{
		return myDataSeries.y (i);
		}

	/**
	 * Returns the number of line segments to be plotted.
	 */
	public int getLineCount()
		{
		return Math.max (0, myDataSeries.length() - 1);
		}

	/**
	 * Returns the starting X coordinate of the given line segment.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Starting X coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public double getLineX1
		(int i)
		{
		return myDataSeries.x (i);
		}

	/**
	 * Returns the starting Y coordinate of the given line segment.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Starting Y coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public double getLineY1
		(int i)
		{
		return myDataSeries.y (i);
		}

	/**
	 * Returns the ending X coordinate of the given line segment.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Ending X coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public double getLineX2
		(int i)
		{
		return myDataSeries.x (i+1);
		}

	/**
	 * Returns the ending Y coordinate of the given line segment.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Ending Y coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public double getLineY2
		(int i)
		{
		return myDataSeries.y (i+1);
		}

	}
