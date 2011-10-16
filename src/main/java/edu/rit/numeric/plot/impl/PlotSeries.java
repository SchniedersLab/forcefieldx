//******************************************************************************
//
// File:    PlotSeries.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.PlotSeries
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

import edu.rit.numeric.ListXYSeries;
import edu.rit.numeric.XYSeries;

import edu.rit.numeric.plot.Dots;
import edu.rit.numeric.plot.Strokes;

import java.awt.BasicStroke;
import java.awt.Color;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class PlotSeries is the abstract base class for a plot series for an
 * {@linkplain XYPlot}.
 * <P>
 * Class PlotSeries provides operations to tell where to plot data points as
 * well as which {@linkplain edu.rit.swing.Drawable Drawable} object to use to
 * plot the points. Class {@linkplain edu.rit.numeric.plot.Dots Dots} provides
 * several shapes of drawable objects for plotting data points.
 * <P>
 * Class PlotSeries provides operations to tell where to plot line segments as
 * well as which stroke and color objects to use to plot the line segments.
 * Class {@linkplain edu.rit.numeric.plot.Strokes Strokes} provides several
 * kinds of strokes for drawing lines.
 * <P>
 * Class PlotSeries provides an operation to tell whether to plot a smooth curve
 * instead of straight line segments.
 *
 * @author  Alan Kaminsky
 * @version 11-Dec-2007
 */
public abstract class PlotSeries
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 272975169581660883L;

	/**
	 * Data series.
	 */
	protected XYSeries myDataSeries;

	/**
	 * Drawable object for plotting data points.
	 */
	protected Dots myDots;

	/**
	 * Stroke for drawing line segments.
	 */
	protected BasicStroke myStroke;

	/**
	 * Color for drawing line segments.
	 */
	protected Color myColor;

	/**
	 * True to plot a smooth curve, false to plot straight line segments.
	 */
	protected boolean mySmooth;

// Exported constructors.

	/**
	 * Construct a new plot series. This constructor is intended for use only by
	 * object deserialization.
	 */
	public PlotSeries()
		{
		}

	/**
	 * Construct a new plot series. Dots are plotted on the data points using
	 * the given drawable object. Line segments are drawn using the given stroke
	 * and color.
	 *
	 * @param  theDataSeries
	 *     Data series.
	 * @param  theDots
	 *     Drawable object for plotting data points. If null, dots are not
	 *     plotted on the data points.
	 * @param  theStroke
	 *     Stroke for drawing line segments. If null, line segments are not
	 *     drawn.
	 * @param  theColor
	 *     Color for drawing line segments. If null, line segments are not
	 *     drawn.
	 * @param  isSmooth
	 *     True to plot a smooth curve, false to plot straight line segments.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDataSeries</TT> is null.
	 */
	public PlotSeries
		(XYSeries theDataSeries,
		 Dots theDots,
		 BasicStroke theStroke,
		 Color theColor,
		 boolean isSmooth)
		{
		if (theDataSeries == null)
			{
			throw new NullPointerException();
			}
		myDataSeries = theDataSeries;
		myDots = theDots;
		myStroke = theStroke;
		myColor = theColor;
		mySmooth = isSmooth;
		}

// Exported operations.

	/**
	 * Returns the data series for this plot series.
	 */
	public XYSeries getDataSeries()
		{
		return myDataSeries;
		}

	/**
	 * Returns the drawable object for plotting data points for this XY plot
	 * series. If dots are not to be plotted, null is returned.
	 */
	public Dots getDots()
		{
		return myDots;
		}

	/**
	 * Returns the stroke for drawing line segments for this plot series. If
	 * line segments are not to be drawn, null is returned.
	 */
	public BasicStroke getStroke()
		{
		return myStroke;
		}

	/**
	 * Returns the color for drawing line segments for this plot series. If line
	 * segments are not to be drawn, null is returned.
	 */
	public Color getColor()
		{
		return myColor;
		}

	/**
	 * Determine whether to plot a smooth curve or straight line segments.
	 *
	 * @return  True to plot a smooth curve, false to plot straight line
	 *          segments.
	 */
	public boolean isSmooth()
		{
		return mySmooth;
		}

	/**
	 * Returns the number of dots to be plotted.
	 * <P>
	 * This method must be overridden in a subclass.
	 */
	public abstract int getDotCount();

	/**
	 * Returns the X coordinate of the given dot.
	 * <P>
	 * This method must be overridden in a subclass.
	 *
	 * @param  i  Dot index.
	 *
	 * @return  X coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getDotCount()-1</TT>.
	 */
	public abstract double getDotX
		(int i);

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
	public abstract double getDotY
		(int i);

	/**
	 * Returns the number of line segments to be plotted.
	 * <P>
	 * This method must be overridden in a subclass.
	 */
	public abstract int getLineCount();

	/**
	 * Returns the starting X coordinate of the given line segment.
	 * <P>
	 * This method must be overridden in a subclass.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Starting X coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public abstract double getLineX1
		(int i);

	/**
	 * Returns the starting Y coordinate of the given line segment.
	 * <P>
	 * This method must be overridden in a subclass.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Starting Y coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public abstract double getLineY1
		(int i);

	/**
	 * Returns the ending X coordinate of the given line segment.
	 * <P>
	 * This method must be overridden in a subclass.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Ending X coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public abstract double getLineX2
		(int i);

	/**
	 * Returns the ending Y coordinate of the given line segment.
	 * <P>
	 * This method must be overridden in a subclass.
	 *
	 * @param  i  Line segment index.
	 *
	 * @return  Ending Y coordinate.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0 ..
	 *     <TT>getLineCount()-1</TT>.
	 */
	public abstract double getLineY2
		(int i);

	/**
	 * Write this plot series to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		out.writeObject (myDots);
		Strokes.writeExternal (myStroke, out);
		out.writeObject (myColor);
		out.writeBoolean (mySmooth);
		int n = myDataSeries.length();
		out.writeInt (n);
		for (int i = 0; i < n; ++ i)
			{
			out.writeDouble (myDataSeries.x(i));
			out.writeDouble (myDataSeries.y(i));
			}
		}

	/**
	 * Read this plot series from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if a class needed to deserialize this plot series was not
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		myDots = (Dots) in.readObject();
		myStroke = Strokes.readExternal (in);
		myColor = (Color) in.readObject();
		mySmooth = in.readBoolean();
		ListXYSeries series = new ListXYSeries();
		int n = in.readInt();
		for (int i = 0; i < n; ++ i)
			{
			series.add
				(in.readDouble(),
				 in.readDouble());
			}
		myDataSeries = series;
		}

	}
