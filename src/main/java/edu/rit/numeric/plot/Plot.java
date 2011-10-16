//******************************************************************************
//
// File:    Plot.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.Plot
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

package edu.rit.numeric.plot;

import edu.rit.numeric.ListXYSeries;
import edu.rit.numeric.XYSeries;

import edu.rit.numeric.plot.impl.Axis;
import edu.rit.numeric.plot.impl.Grid;
import edu.rit.numeric.plot.impl.Label;
import edu.rit.numeric.plot.impl.LinearAxis;
import edu.rit.numeric.plot.impl.LogarithmicAxis;
import edu.rit.numeric.plot.impl.NumericalAxis;
import edu.rit.numeric.plot.impl.PlotSeries;
import edu.rit.numeric.plot.impl.SegmentedPlotSeries;
import edu.rit.numeric.plot.impl.Ticks;
import edu.rit.numeric.plot.impl.XYPlot;
import edu.rit.numeric.plot.impl.XYPlotSeries;

import edu.rit.swing.Displayable;
import edu.rit.swing.DisplayableFrame;
import edu.rit.swing.Drawable;
import edu.rit.swing.Viewable;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;

import java.awt.geom.Rectangle2D;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.Externalizable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;

import java.text.DecimalFormat;

import java.util.ArrayList;
import java.util.Scanner;

/**
 * Class Plot is a class for creating and displaying a plot.
 * <P>
 * To create a plot, construct an instance of class Plot, and call the desired
 * methods to specify attributes of the plot and to provide the data to be
 * plotted.
 * <P>
 * To display a plot, call the <TT>getFrame()</TT> method to get a {@linkplain
 * DisplayableFrame} that displays the plot in a window, then call the frame's
 * <TT>setVisible(true)</TT> method. The frame has menus for zooming the
 * display, for changing the plot's attributes (title, margins, grid, X axis, Y
 * axis, etc.), and for saving the plot in a PNG or PostScript file.
 * <P>
 * If the plot's attributes are changed after displaying the frame, call the
 * frame's <TT>display()</TT> method, passing in a reference to the plot object.
 * If the data in the plot's X-Y series are changed after displaying the frame,
 * but the plot's attributes are not changed, call the frame's
 * <TT>repaint()</TT> method.
 * <P>
 * A Plot object can be written (serialized) into a file. The {@linkplain View}
 * program can be used to display the file.
 *
 * @author  Alan Kaminsky
 * @version 04-Aug-2008
 */
public class Plot
	implements Viewable, Externalizable
	{

// Exported enumerations.

	/**
	 * Enum Plot.AxisKind enumerates the possible kinds of axis on a {@linkplain
	 * Plot}.
	 *
	 * @author  Alan Kaminsky
	 * @version 11-Jun-2007
	 */
	public static enum AxisKind
		{
		/**
		 * Linear axis.
		 */
		LINEAR,

		/**
		 * Logarithmic axis.
		 */
		LOGARITHMIC,
		}

	/**
	 * Linear axis.
	 */
	public static final AxisKind LINEAR = AxisKind.LINEAR;

	/**
	 * Logarithmic axis.
	 */
	public static final AxisKind LOGARITHMIC = AxisKind.LOGARITHMIC;

// Exported constants.

	/**
	 * The label is positioned above and to the left of its location.
	 */
	public static final int ABOVE_LEFT = Label.ABOVE_LEFT;

	/**
	 * The label is positioned above its location.
	 */
	public static final int ABOVE = Label.ABOVE;

	/**
	 * The label is positioned above and to the right of its location.
	 */
	public static final int ABOVE_RIGHT = Label.ABOVE_RIGHT;

	/**
	 * The label is positioned to the left of its location.
	 */
	public static final int LEFT = Label.LEFT;

	/**
	 * The label is centered on its location.
	 */
	public static final int CENTER = Label.CENTER;

	/**
	 * The label is positioned to the right of its location.
	 */
	public static final int RIGHT = Label.RIGHT;

	/**
	 * The label is positioned below and to the left of its location.
	 */
	public static final int BELOW_LEFT = Label.BELOW_LEFT;

	/**
	 * The label is positioned below its location.
	 */
	public static final int BELOW = Label.BELOW;

	/**
	 * The label is positioned below and to the right of its location.
	 */
	public static final int BELOW_RIGHT = Label.BELOW_RIGHT;

	/**
	 * The label is rotated 90 degrees left (counterclockwise).
	 */
	public static final int ROTATE_LEFT = Label.ROTATE_LEFT;

	/**
	 * The label is rotated 90 degrees right (clockwise).
	 */
	public static final int ROTATE_RIGHT = Label.ROTATE_RIGHT;

	/**
	 * The label is rotated 180 degrees around.
	 */
	public static final int ROTATE_AROUND = Label.ROTATE_AROUND;

	/**
	 * The label's location is specified in pixel coordinates relative to the
	 * lower left corner of the plot area, rather than in plot coordinates.
	 */
	public static final int PIXEL_COORDINATES = Label.PIXEL_COORDINATES;

// Hidden data members.

	private static final long serialVersionUID = 5469967079348825097L;

	private static final Font thePlotFont =
		new Font ("SansSerif", Font.BOLD, 14);
	private static final Font theAxisFont =
		new Font ("SansSerif", Font.PLAIN, 12);

	private static final BasicStroke theAxisStroke =
		Strokes.solid (0.6f);
	private static final BasicStroke theMajorGridLineStroke =
		Strokes.solid (0.2f);
	private static final BasicStroke theMinorGridLineStroke =
		Strokes.solid (0.1f);

	// Plot attributes.
	private String frameTitle = null;
	private String plotTitle = null;
	private Font plotTitleFont = thePlotFont;
	private double plotTitleOffset = 9.0;
	private double leftMargin = 54.0;
	private double topMargin = 27.0;
	private double rightMargin = 27.0;
	private double bottomMargin = 54.0;
	private boolean majorGridLines = true;
	private boolean minorGridLines = false;

	// X axis attributes.
	private AxisKind xAxisKind = LINEAR;
	private double xAxisStart = Double.NaN; // Automatic
	private double xAxisEnd = Double.NaN; // Automatic
	private int xAxisMajorDivisions = 10;
	private int xAxisMinorDivisions = 1;
	private double xAxisCrossing = Double.NaN; // Automatic
	private double xAxisLength = 288.0; // 4 inches
	private DecimalFormat xAxisTickFormat = new DecimalFormat ("0");
	private double xAxisTickScale = 1.0;
	private Font xAxisTickFont = theAxisFont;
	private String xAxisTitle = null;
	private Font xAxisTitleFont = theAxisFont;
	private double xAxisTitleOffset = 30.0;

	// Y axis attributes.
	private AxisKind yAxisKind = LINEAR;
	private double yAxisStart = Double.NaN; // Automatic
	private double yAxisEnd = Double.NaN; // Automatic
	private int yAxisMajorDivisions = 10;
	private int yAxisMinorDivisions = 1;
	private double yAxisCrossing = Double.NaN; // Automatic
	private double yAxisLength = 288.0; // 4 inches
	private DecimalFormat yAxisTickFormat = new DecimalFormat ("0");
	private double yAxisTickScale = 1.0;
	private Font yAxisTickFont = theAxisFont;
	private String yAxisTitle = null;
	private Font yAxisTitleFont = theAxisFont;
	private double yAxisTitleOffset = 36.0;

	// Plot series attributes.
	private Dots seriesDots = Dots.circle (5);
	private BasicStroke seriesStroke = Strokes.solid (2);
	private Color seriesColor = Color.black;
	private boolean seriesSmooth = false;

	// Label attributes.
	private int labelPosition = CENTER;
	private double labelOffset = 0.0;
	private Font labelFont = theAxisFont;
	private Color labelColor = Color.black;
	private Color labelBackground = null;

	// Lists of plot series and labels.
	private ArrayList<PlotSeries> myPlotSeries = new ArrayList<PlotSeries>();
	private ArrayList<Label> myLabels = new ArrayList<Label>();

	// X-Y plot object to be displayed.
	private XYPlot myXYPlot;

// Exported constructors.

	/**
	 * Construct a new plot.
	 */
	public Plot()
		{
		}

// Exported operations.

	/**
	 * Set the title for the frame that displays the plot. The default is the
	 * same as the plot title.
	 *
	 * @param  theTitle  Title.
	 *
	 * @return  This plot object.
	 */
	public Plot frameTitle
		(String theTitle)
		{
		frameTitle = theTitle;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the title for the plot. The default is no title.
	 *
	 * @param  theTitle  Title.
	 *
	 * @return  This plot object.
	 */
	public Plot plotTitle
		(String theTitle)
		{
		plotTitle = theTitle;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the title font for the plot. The default is sans-serif, plain, 14
	 * point.
	 *
	 * @param  theTitleFont  Title font.
	 *
	 * @return  This plot object.
	 */
	public Plot plotTitleFont
		(Font theTitleFont)
		{
		plotTitleFont = theTitleFont;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the offset from the plot area to the title for the plot. The default
	 * is 9 points.
	 *
	 * @param  theTitleOffset  Title offset (points).
	 *
	 * @return  This plot object.
	 */
	public Plot plotTitleOffset
		(double theTitleOffset)
		{
		plotTitleOffset = theTitleOffset;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the margins for the plot. This is the distance from the left, top,
	 * right, and bottom sides of the plot area to the left, top, right, and
	 * bottom sides of the display.
	 *
	 * @param  theMargin  Margin (points).
	 *
	 * @return  This plot object.
	 */
	public Plot margins
		(double theMargin)
		{
		leftMargin = theMargin;
		topMargin = theMargin;
		rightMargin = theMargin;
		bottomMargin = theMargin;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the left margin for the plot. This is the distance from the left side
	 * of the plot area to the left side of the display. The default is 54
	 * points.
	 *
	 * @param  theMargin  Left margin (points).
	 *
	 * @return  This plot object.
	 */
	public Plot leftMargin
		(double theMargin)
		{
		leftMargin = theMargin;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the top margin for the plot. This is the distance from the top side
	 * of the plot area to the top side of the display. The default is 27
	 * points.
	 *
	 * @param  theMargin  Top margin (points).
	 *
	 * @return  This plot object.
	 */
	public Plot topMargin
		(double theMargin)
		{
		topMargin = theMargin;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the right margin for the plot. This is the distance from the right
	 * side of the plot area to the right side of the display. The default is 27
	 * points.
	 *
	 * @param  theMargin  Right margin (points).
	 *
	 * @return  This plot object.
	 */
	public Plot rightMargin
		(double theMargin)
		{
		rightMargin = theMargin;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the bottom margin for the plot. This is the distance from the bottom
	 * side of the plot area to the bottom side of the display. The default is
	 * 54 points.
	 *
	 * @param  theMargin  Bottom margin (points).
	 *
	 * @return  This plot object.
	 */
	public Plot bottomMargin
		(double theMargin)
		{
		bottomMargin = theMargin;
		myXYPlot = null;
		return this;
		}

	/**
	 * Specify whether major division grid lines are drawn on the plot. The
	 * default is to draw major division grid lines.
	 *
	 * @param  theGridLines  True to draw major division grid lines, false not
	 *                       to draw them.
	 *
	 * @return  This plot object.
	 */
	public Plot majorGridLines
		(boolean theGridLines)
		{
		majorGridLines = theGridLines;
		myXYPlot = null;
		return this;
		}

	/**
	 * Specify whether minor division grid lines are drawn on the plot. The
	 * default is not to draw minor division grid lines.
	 *
	 * @param  theGridLines  True to draw minor division grid lines, false not
	 *                       to draw them.
	 *
	 * @return  This plot object.
	 */
	public Plot minorGridLines
		(boolean theGridLines)
		{
		minorGridLines = theGridLines;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the kind of axis for the X axis. The default is a linear axis.
	 *
	 * @param  theKind  {@link #LINEAR} or {@link #LOGARITHMIC}.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisKind
		(AxisKind theKind)
		{
		xAxisKind = theKind;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the starting value for the X axis. For a linear axis, this is the
	 * actual starting value. For a logarithmic axis, this the base-10 logarithm
	 * of the starting value, and it is truncated to the nearest integer. The
	 * default is to determine the starting value automatically from the plotted
	 * data.
	 *
	 * @param  theStart  Starting value.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisStart
		(double theStart)
		{
		xAxisStart = theStart;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the ending value for the X axis. For a linear axis, this is the
	 * actual ending value. For a logarithmic axis, this the base-10 logarithm
	 * of the ending value, and it is truncated to the nearest integer. The
	 * default is to determine the ending value automatically from the plotted
	 * data.
	 *
	 * @param  theEnd  Ending value.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisEnd
		(double theEnd)
		{
		xAxisEnd = theEnd;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the number of major divisions for the X axis. The default is 10 major
	 * divisions.
	 * <P>
	 * For a logarithmic axis, <TT>theMajorDivisions</TT> is ignored; each
	 * factor of 10 between the starting and ending values is one major
	 * division.
	 *
	 * @param  theMajorDivisions  Number of major divisions.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisMajorDivisions
		(int theMajorDivisions)
		{
		xAxisMajorDivisions = theMajorDivisions;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the number of minor divisions <I>per major division</I> for the X
	 * axis. The default is 1 minor division per major division.
	 *
	 * @param  theMinorDivisions  Number of minor divisions per major division.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisMinorDivisions
		(int theMinorDivisions)
		{
		xAxisMinorDivisions = theMinorDivisions;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the perpendicular axis crossing value for the X axis. For a linear
	 * axis, this is the actual crossing value. For a logarithmic axis, this the
	 * base-10 logarithm of the crossing value, and it is truncated to the
	 * nearest integer. The default is for the perpendicular axis to cross at
	 * the X axis's starting value.
	 *
	 * @param  theCrossing  Perpendicular axis crossing value.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisCrossing
		(double theCrossing)
		{
		xAxisCrossing = theCrossing;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the length of the X axis. The default is 288 points.
	 *
	 * @param  theLength  Axis length (points).
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisLength
		(double theLength)
		{
		xAxisLength = theLength;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the tick label format for the X axis. If <TT>theTickFormat</TT> is
	 * null, tick labels will not be drawn on the X axis. The default is a
	 * decimal format of <TT>"0"</TT>.
	 *
	 * @param  theTickFormat  Tick format, or null.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisTickFormat
		(DecimalFormat theTickFormat)
		{
		xAxisTickFormat = theTickFormat;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the tick label scale factor for the X axis. Each tick mark's
	 * numerical value is divided by the scale factor before being displayed.
	 * The default is 1 (no scaling).
	 *
	 * @param  theTickScale  Tick scale factor.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisTickScale
		(double theTickScale)
		{
		xAxisTickScale = theTickScale;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the tick label font for the X axis. The default is sans-serif, plain,
	 * 12 point.
	 *
	 * @param  theTickFont  Tick font.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisTickFont
		(Font theTickFont)
		{
		xAxisTickFont = theTickFont;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the title for the X axis. The default is no title.
	 *
	 * @param  theTitle  Title.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisTitle
		(String theTitle)
		{
		xAxisTitle = theTitle;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the title font for the X axis. The default is sans-serif, plain, 12
	 * point.
	 *
	 * @param  theTitleFont  Title font.
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisTitleFont
		(Font theTitleFont)
		{
		xAxisTitleFont = theTitleFont;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the offset from the axis to the title for the X axis. The default is
	 * 30 points.
	 *
	 * @param  theTitleOffset  Title offset (points).
	 *
	 * @return  This plot object.
	 */
	public Plot xAxisTitleOffset
		(double theTitleOffset)
		{
		xAxisTitleOffset = theTitleOffset;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the kind of axis for the Y axis. The default is a linear axis.
	 *
	 * @param  theKind  {@link #LINEAR} or {@link #LOGARITHMIC}.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisKind
		(AxisKind theKind)
		{
		yAxisKind = theKind;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the starting value for the Y axis. For a linear axis, this is the
	 * actual starting value. For a logarithmic axis, this the base-10 logarithm
	 * of the starting value, and it is truncated to the nearest integer. The
	 * default is to determine the starting value automatically from the plotted
	 * data.
	 *
	 * @param  theStart  Starting value.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisStart
		(double theStart)
		{
		yAxisStart = theStart;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the ending value for the Y axis. For a linear axis, this is the
	 * actual ending value. For a logarithmic axis, this the base-10 logarithm
	 * of the ending value, and it is truncated to the nearest integer. The
	 * default is to determine the ending value automatically from the plotted
	 * data.
	 *
	 * @param  theEnd  Ending value.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisEnd
		(double theEnd)
		{
		yAxisEnd = theEnd;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the number of major divisions for the Y axis. The default is 10 major
	 * divisions.
	 * <P>
	 * For a logarithmic axis, <TT>theMajorDivisions</TT> is ignored; each
	 * factor of 10 between the starting and ending values is one major
	 * division.
	 *
	 * @param  theMajorDivisions  Number of major divisions.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisMajorDivisions
		(int theMajorDivisions)
		{
		yAxisMajorDivisions = theMajorDivisions;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the number of minor divisions <I>per major division</I> for the X
	 * axis. The default is 1 minor division per major division.
	 *
	 * @param  theMinorDivisions  Number of minor divisions per major division.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisMinorDivisions
		(int theMinorDivisions)
		{
		yAxisMinorDivisions = theMinorDivisions;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the perpendicular axis crossing value for the Y axis. For a linear
	 * axis, this is the actual crossing value. For a logarithmic axis, this the
	 * base-10 logarithm of the crossing value, and it is truncated to the
	 * nearest integer. The default is for the perpendicular axis to cross at
	 * the Y axis's starting value.
	 *
	 * @param  theCrossing  Perpendicular axis crossing value.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisCrossing
		(double theCrossing)
		{
		yAxisCrossing = theCrossing;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the length of the Y axis. The default is 288 points.
	 *
	 * @param  theLength  Axis length (points).
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisLength
		(double theLength)
		{
		yAxisLength = theLength;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the tick label format for the Y axis. If <TT>theTickFormat</TT> is
	 * null, tick labels will not be drawn on the Y axis. The default is a
	 * decimal format of <TT>"0"</TT>.
	 *
	 * @param  theTickFormat  Tick format, or null.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisTickFormat
		(DecimalFormat theTickFormat)
		{
		yAxisTickFormat = theTickFormat;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the tick label scale factor for the Y axis. Each tick mark's
	 * numerical value is divided by the scale factor before being displayed.
	 * The default is 1 (no scaling).
	 *
	 * @param  theTickScale  Tick scale factor.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisTickScale
		(double theTickScale)
		{
		yAxisTickScale = theTickScale;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the tick label font for the Y axis. The default is sans-serif, plain,
	 * 12 point.
	 *
	 * @param  theTickFont  Tick font.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisTickFont
		(Font theTickFont)
		{
		yAxisTickFont = theTickFont;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the title for the Y axis. The default is no title.
	 *
	 * @param  theTitle  Title.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisTitle
		(String theTitle)
		{
		yAxisTitle = theTitle;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the title font for the Y axis. The default is sans-serif, plain, 12
	 * point.
	 *
	 * @param  theTitleFont  Title font.
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisTitleFont
		(Font theTitleFont)
		{
		yAxisTitleFont = theTitleFont;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the offset from the axis to the title for the Y axis. The default is
	 * 36 points.
	 *
	 * @param  theTitleOffset  Title offset (points).
	 *
	 * @return  This plot object.
	 */
	public Plot yAxisTitleOffset
		(double theTitleOffset)
		{
		yAxisTitleOffset = theTitleOffset;
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the dots to use when plotting a data series. All data series added
	 * hereafter will use the given dots. The default is a solid black circle 5
	 * points in diameter. See class {@linkplain Dots} for other shapes.
	 *
	 * @param  theDots  Dots, or null not to plot dots.
	 *
	 * @return  This plot object.
	 */
	public Plot seriesDots
		(Dots theDots)
		{
		seriesDots = theDots;
		return this;
		}

	/**
	 * Set the stroke to use for line segments when plotting a data series. All
	 * data series added hereafter will use the given stroke. The default is a
	 * solid stroke 2 points wide. See class {@linkplain Strokes} for other
	 * strokes.
	 *
	 * @param  theStroke  Stroke, or null not to plot line segments.
	 *
	 * @return  This plot object.
	 */
	public Plot seriesStroke
		(BasicStroke theStroke)
		{
		seriesStroke = theStroke;
		return this;
		}

	/**
	 * Set the color to use for line segments when plotting a data series. All
	 * data series added hereafter will use the given color. The default is
	 * black.
	 *
	 * @param  theColor  Color, or null not to plot line segments.
	 *
	 * @return  This plot object.
	 */
	public Plot seriesColor
		(Color theColor)
		{
		seriesColor = theColor;
		return this;
		}

	/**
	 * Specify whether to draw a smooth curve when plotting a data series. All
	 * data series added hereafter will use the given smoothing. The default is
	 * false (draw straight line segments).
	 *
	 * @param  isSmooth  True to draw a smooth curve, false to draw straight
	 *                   line segments.
	 *
	 * @return  This plot object.
	 */
	public Plot seriesSmooth
		(boolean isSmooth)
		{
		seriesSmooth = isSmooth;
		return this;
		}

	/**
	 * Add an X-Y data series to this plot with data points read from the file
	 * with the given name. Double values are read from the file until there are
	 * no more. Each pair of double values becomes an (<I>x,y</I>) point in the
	 * data series. Any leftover double value is discarded.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theFileName  File name.
	 *
	 * @return  This plot object.
	 *
	 * @exception  FileNotFoundException
	 *     Thrown if the file could not be read.
	 */
	public Plot xySeries
		(String theFileName)
		throws FileNotFoundException
		{
		return xySeries (new Scanner (new File (theFileName)));
		}

	/**
	 * Add an X-Y data series to this plot with data points read from the given
	 * file. Double values are read from the file until there are no more. Each
	 * pair of double values becomes an (<I>x,y</I>) point in the data series.
	 * Any leftover double value is discarded.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theFile  File.
	 *
	 * @return  This plot object.
	 *
	 * @exception  FileNotFoundException
	 *     Thrown if the file could not be read.
	 */
	public Plot xySeries
		(File theFile)
		throws FileNotFoundException
		{
		return xySeries (new Scanner (theFile));
		}

	/**
	 * Add an X-Y data series to this plot with data points read from the given
	 * input stream. Double values are read from the input stream until there
	 * are no more. Each pair of double values becomes an (<I>x,y</I>) point in
	 * the data series. Any leftover double value is discarded. After adding the
	 * data series, the input stream is closed.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theInputStream  Input stream.
	 *
	 * @return  This plot object.
	 */
	public Plot xySeries
		(InputStream theInputStream)
		{
		return xySeries (new Scanner (theInputStream));
		}

	/**
	 * Add an X-Y data series to this plot with data points read from the given
	 * scanner. Double values are read from the scanner until there are no more.
	 * Each pair of double values becomes an (<I>x,y</I>) point in the data
	 * series. Any leftover double value is discarded. After constructing the
	 * data series, the scanner is closed.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theScanner  Scanner.
	 *
	 * @return  This plot object.
	 */
	public Plot xySeries
		(Scanner theScanner)
		{
		ListXYSeries series = new ListXYSeries().add (theScanner);
		theScanner.close();
		return xySeries (series);
		}

	/**
	 * Add an X-Y data series to this plot with data points taken from the
	 * arguments. The argument values are grouped into pairs, where the first
	 * argument value is the X value and the second argument value is the Y
	 * value. Each pair of argument values becomes one data point in the X-Y
	 * data series. It is assumed that there are an even number of argument
	 * values.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theData  Data point X and Y coordinates.
	 *
	 * @return  This plot object.
	 */
	public Plot xySeries
		(double... theData)
		{
		return xySeries (new ListXYSeries().add (theData));
		}

	/**
	 * Add an X-Y data series to this plot with data points taken from the two
	 * arrays. The first array contains the X coordinates, the second array
	 * contains the Y coordinates.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theXArray  Array of data point X coordinates.
	 * @param  theYArray  Array of data point Y coordinates.
	 *
	 * @return  This plot object.
	 */
	public Plot xySeries
		(double[] theXArray,
		 double[] theYArray)
		{
		return xySeries (new ListXYSeries().add (theXArray, theYArray));
		}

	/**
	 * Add an X-Y data series to this plot with data points taken from the given
	 * {@linkplain XYSeries}.
	 * <P>
	 * A dot is drawn at each data point as specified by the
	 * <TT>seriesDot()</TT> method. A line segment is drawn from each data point
	 * to the next as specified by the <TT>seriesStroke()</TT> and
	 * <TT>seriesColor()</TT> methods. A smooth curve is drawn, or not, as
	 * specified by the <TT>seriesSmooth()</TT> method.
	 *
	 * @param  theXYSeries  X-Y series.
	 *
	 * @return  This plot object.
	 */
	public Plot xySeries
		(XYSeries theXYSeries)
		{
		myPlotSeries.add
			(new XYPlotSeries
				(theXYSeries,
				 seriesDots,
				 seriesStroke,
				 seriesColor,
				 seriesSmooth));
		myXYPlot = null;
		return this;
		}

	/**
	 * Add a segmented data series to this plot with data points read from the
	 * file with the given name. Double values are read from the file until
	 * there are no more. Each pair of double values becomes an (<I>x,y</I>)
	 * point in the data series. Any leftover double value is discarded.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theFileName  File name.
	 *
	 * @return  This plot object.
	 *
	 * @exception  FileNotFoundException
	 *     Thrown if the file could not be read.
	 */
	public Plot segmentedSeries
		(String theFileName)
		throws FileNotFoundException
		{
		return segmentedSeries (new Scanner (new File (theFileName)));
		}

	/**
	 * Add a segmented data series to this plot with data points read from the
	 * given file. Double values are read from the file until there are no more.
	 * Each pair of double values becomes an (<I>x,y</I>) point in the data
	 * series. Any leftover double value is discarded.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theFile  File.
	 *
	 * @return  This plot object.
	 *
	 * @exception  FileNotFoundException
	 *     Thrown if the file could not be read.
	 */
	public Plot segmentedSeries
		(File theFile)
		throws FileNotFoundException
		{
		return segmentedSeries (new Scanner (theFile));
		}

	/**
	 * Add a segmented data series to this plot with data points read from the
	 * given input stream. Double values are read from the input stream until
	 * there are no more. Each pair of double values becomes an (<I>x,y</I>)
	 * point in the data series. Any leftover double value is discarded. After
	 * adding the data series, the input stream is closed.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theInputStream  Input stream.
	 *
	 * @return  This plot object.
	 */
	public Plot segmentedSeries
		(InputStream theInputStream)
		{
		return segmentedSeries (new Scanner (theInputStream));
		}

	/**
	 * Add a segmented data series to this plot with data points read from the
	 * given scanner. Double values are read from the scanner until there are no
	 * more. Each pair of double values becomes an (<I>x,y</I>) point in the
	 * data series. Any leftover double value is discarded. After constructing
	 * the data series, the scanner is closed.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theScanner  Scanner.
	 *
	 * @return  This plot object.
	 */
	public Plot segmentedSeries
		(Scanner theScanner)
		{
		ListXYSeries series = new ListXYSeries().add (theScanner);
		theScanner.close();
		return segmentedSeries (series);
		}

	/**
	 * Add a segmented data series to this plot with data points taken from the
	 * arguments. The argument values are grouped into pairs, where the first
	 * argument value is the X value and the second argument value is the Y
	 * value. Each pair of argument values becomes one data point in the
	 * segmented data series. It is assumed that there are an even number of
	 * argument values.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theData  Data point X and Y coordinates.
	 *
	 * @return  This plot object.
	 */
	public Plot segmentedSeries
		(double... theData)
		{
		return segmentedSeries (new ListXYSeries().add (theData));
		}

	/**
	 * Add a segmented data series to this plot with data points taken from the
	 * two arrays. The first array contains the X coordinates, the second array
	 * contains the Y coordinates.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theXArray  Array of data point X coordinates.
	 * @param  theYArray  Array of data point Y coordinates.
	 *
	 * @return  This plot object.
	 */
	public Plot segmentedSeries
		(double[] theXArray,
		 double[] theYArray)
		{
		return segmentedSeries (new ListXYSeries().add (theXArray, theYArray));
		}

	/**
	 * Add a segmented data series to this plot with data points taken from the
	 * given {@linkplain XYSeries}.
	 * <P>
	 * The data points are divided into pairs. Each pair of data points
	 * specifies the starting and ending coordinates of a line segment. A dot is
	 * drawn at each endpoint of each line segment. A line segment's starting
	 * coordinates need not be the same as the previous line segment's ending
	 * coordinates. If the data series has an odd number of data points, the
	 * last data point is ignored.
	 *
	 * @param  theXYSeries  X-Y series.
	 *
	 * @return  This plot object.
	 */
	public Plot segmentedSeries
		(XYSeries theXYSeries)
		{
		myPlotSeries.add
			(new SegmentedPlotSeries
				(theXYSeries,
				 seriesDots,
				 seriesStroke,
				 seriesColor));
		myXYPlot = null;
		return this;
		}

	/**
	 * Set the positioning to use for labels. All labels added hereafter will
	 * use the given positioning. The label positioning argument must be one of
	 * the following constants to specify how the label is positioned relative
	 * to the label's (<I>x,y</I>) location: {@link #ABOVE_LEFT}, {@link
	 * #ABOVE}, {@link #ABOVE_RIGHT}, {@link #LEFT}, {@link #CENTER}, {@link
	 * #RIGHT}, {@link #BELOW_LEFT}, {@link #BELOW}, or {@link #BELOW_RIGHT}.
	 * The default is {@link #CENTER}.
	 * <P>
	 * One of the following constants may be added to rotate the label around
	 * the label's (<I>x,y</I>) location: {@link #ROTATE_LEFT}, {@link
	 * #ROTATE_RIGHT}, or {@link #ROTATE_AROUND}. For example,
	 * <TT>labelPosition(Plot.LEFT+Plot.ROTATE_LEFT)</TT> will rotate each label
	 * 90 degrees left and position it to the left of its (<I>x,y</I>) location.
	 * The default is no rotation.
	 * <P>
	 * The following constant may also be added: {@link #PIXEL_COORDINATES}. If
	 * this constant is omitted, the label's (<I>x,y</I>) location is specified
	 * in plot coordinates as determined by the X and Y axes' attributes. If
	 * this constant is included, the label's (<I>x,y</I>) location is specified
	 * in pixel coordinates relative to the lower left corner of the plot area,
	 * regardless of the X and Y axes' attributes.
	 *
	 * @param  thePosition  Label positioning.
	 *
	 * @return  This plot object.
	 */
	public Plot labelPosition
		(int thePosition)
		{
		labelPosition = thePosition;
		return this;
		}

	/**
	 * Set the offset to use for labels. All labels added hereafter will use the
	 * given offset. The default is 0. The offset specifies, in pixel
	 * coordinates, how far away from the label's (<I>x,y</I>) location the
	 * label is positioned in the direction determined by the label positioning
	 * attribute.
	 *
	 * @param  theOffset  Label offset.
	 *
	 * @return  This plot object.
	 */
	public Plot labelOffset
		(double theOffset)
		{
		labelOffset = theOffset;
		return this;
		}

	/**
	 * Set the font to use for labels. All labels added hereafter will use the
	 * given font. The default is sans-serif, plain, 12 point.
	 *
	 * @param  theFont  Label offset.
	 *
	 * @return  This plot object.
	 */
	public Plot labelFont
		(Font theFont)
		{
		labelFont = theFont;
		return this;
		}

	/**
	 * Set the color to use for labels. All labels added hereafter will use the
	 * given color. The default is black.
	 *
	 * @param  theColor  Label color.
	 *
	 * @return  This plot object.
	 */
	public Plot labelColor
		(Color theColor)
		{
		labelColor = theColor;
		return this;
		}

	/**
	 * Set the color to use for label backgrounds. All labels added hereafter
	 * will use the given color for the background. The default is transparent.
	 *
	 * @param  theColor  Label background color, or null for transparent.
	 *
	 * @return  This plot object.
	 */
	public Plot labelBackground
		(Color theColor)
		{
		labelBackground = theColor;
		return this;
		}

	/**
	 * Add a label to this plot. The label's position is (<I>x,y</I>).
	 * <P>
	 * The label supports a subset of HTML for specifying attributes of
	 * characters. The following HTML tags may be embedded in a label's text:
	 * <UL>
	 * <LI><TT>&lt;B&gt;Bold&lt;/B&gt;</TT>
	 * <LI><TT>&lt;I&gt;Italic&lt;/I&gt;</TT>
	 * <LI><TT>&lt;SUP&gt;Superscript&lt;/SUP&gt;</TT>
	 * <LI><TT>&lt;SUB&gt;Subscript&lt;/SUB&gt;</TT>
	 * </UL>
	 * <P>
	 * To get literal less-than, greater-than, and ampersand characters in the
	 * text, use the following HTML character entities:
	 * <UL>
	 * <LI><TT>&amp;lt;</TT> -- Less-than (&lt;)
	 * <LI><TT>&amp;gt;</TT> -- Greater-than (&gt;)
	 * <LI><TT>&amp;amp;</TT> -- Ampersand (&amp;)
	 * </UL>
	 *
	 * @param  text  Label text.
	 * @param  x     X coordinate of label position.
	 * @param  y     Y coordinate of label position.
	 *
	 * @return  This plot object.
	 */
	public Plot label
		(String text,
		 double x,
		 double y)
		{
		myLabels.add
			(new Label
				(text,
				 x,
				 y,
				 labelPosition,
				 labelOffset,
				 labelFont,
				 labelColor,
				 labelBackground));
		myXYPlot = null;
		return this;
		}

// Exported operations implemented from interface Drawable.

	/**
	 * Draw this drawable object in the given graphics context. Upon return from
	 * this method, the given graphics context's state (color, font, transform,
	 * clip, and so on) is the same as it was upon entry to this method.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		getXYPlot().draw (g2d);
		}

// Exported operations implemented from interface Displayable.

	/**
	 * Returns this displayable object's bounding box. This is the smallest
	 * rectangle that encloses all of this displayable object.
	 */
	public Rectangle2D getBoundingBox()
		{
		return getXYPlot().getBoundingBox();
		}

	/**
	 * Returns this displayable object's background paint.
	 */
	public Paint getBackgroundPaint()
		{
		return getXYPlot().getBackgroundPaint();
		}

// Exported operations implemented from interface Viewable.

	/**
	 * Get a displayable frame in which to view this viewable object. Initially,
	 * the returned frame is displaying this viewable object.
	 *
	 * @return  Displayable frame.
	 */
	public DisplayableFrame getFrame()
		{
		return new PlotFrame (getTitle(), this);
		}

	/**
	 * Get the title for the frame used to view this viewable object. If the
	 * title is null, a default title is used.
	 *
	 * @return  Title.
	 */
	public String getTitle()
		{
		return frameTitle != null ? frameTitle : plotTitle;
		}

// Exported operations implemented from interface Externalizable.

	/**
	 * Write this plot to the given object output stream.
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
		// Plot attributes.
		out.writeObject (frameTitle);
		out.writeObject (plotTitle);
		out.writeObject (plotTitleFont);
		out.writeDouble (plotTitleOffset);
		out.writeDouble (leftMargin);
		out.writeDouble (topMargin);
		out.writeDouble (rightMargin);
		out.writeDouble (bottomMargin);
		out.writeBoolean (majorGridLines);
		out.writeBoolean (minorGridLines);

		// X axis attributes.
		out.writeObject (xAxisKind);
		out.writeDouble (xAxisStart);
		out.writeDouble (xAxisEnd);
		out.writeInt (xAxisMajorDivisions);
		out.writeInt (xAxisMinorDivisions);
		out.writeDouble (xAxisCrossing);
		out.writeDouble (xAxisLength);
		out.writeObject (xAxisTickFormat);
		out.writeDouble (xAxisTickScale);
		out.writeObject (xAxisTickFont);
		out.writeObject (xAxisTitle);
		out.writeObject (xAxisTitleFont);
		out.writeDouble (xAxisTitleOffset);

		// Y axis attributes.
		out.writeObject (yAxisKind);
		out.writeDouble (yAxisStart);
		out.writeDouble (yAxisEnd);
		out.writeInt (yAxisMajorDivisions);
		out.writeInt (yAxisMinorDivisions);
		out.writeDouble (yAxisCrossing);
		out.writeDouble (yAxisLength);
		out.writeObject (yAxisTickFormat);
		out.writeDouble (yAxisTickScale);
		out.writeObject (yAxisTickFont);
		out.writeObject (yAxisTitle);
		out.writeObject (yAxisTitleFont);
		out.writeDouble (yAxisTitleOffset);

		// Plot series attributes.
		out.writeObject (seriesDots);
		Strokes.writeExternal (seriesStroke, out);
		out.writeObject (seriesColor);
		out.writeBoolean (seriesSmooth);

		// Label attributes.
		out.writeInt (labelPosition);
		out.writeDouble (labelOffset);
		out.writeObject (labelFont);
		out.writeObject (labelColor);
		out.writeObject (labelBackground);

		// Lists of plot series and labels.
		out.writeInt (myPlotSeries.size());
		for (PlotSeries series : myPlotSeries)
			{
			out.writeObject (series);
			}
		out.writeInt (myLabels.size());
		for (Label label : myLabels)
			{
			out.writeObject (label);
			}
		}

	/**
	 * Read this plot from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if a class needed to deserialize this plot was not found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		// Plot attributes.
		frameTitle = (String) in.readObject();
		plotTitle = (String) in.readObject();
		plotTitleFont = (Font) in.readObject();
		plotTitleOffset = in.readDouble();
		leftMargin = in.readDouble();
		topMargin = in.readDouble();
		rightMargin = in.readDouble();
		bottomMargin = in.readDouble();
		majorGridLines = in.readBoolean();
		minorGridLines = in.readBoolean();

		// X axis attributes.
		xAxisKind = (AxisKind) in.readObject();
		xAxisStart = in.readDouble();
		xAxisEnd = in.readDouble();
		xAxisMajorDivisions = in.readInt();
		xAxisMinorDivisions = in.readInt();
		xAxisCrossing = in.readDouble();
		xAxisLength = in.readDouble();
		xAxisTickFormat = (DecimalFormat) in.readObject();
		xAxisTickScale = in.readDouble();
		xAxisTickFont = (Font) in.readObject();
		xAxisTitle = (String) in.readObject();
		xAxisTitleFont = (Font) in.readObject();
		xAxisTitleOffset = in.readDouble();

		// Y axis attributes.
		yAxisKind = (AxisKind) in.readObject();
		yAxisStart = in.readDouble();
		yAxisEnd = in.readDouble();
		yAxisMajorDivisions = in.readInt();
		yAxisMinorDivisions = in.readInt();
		yAxisCrossing = in.readDouble();
		yAxisLength = in.readDouble();
		yAxisTickFormat = (DecimalFormat) in.readObject();
		yAxisTickScale = in.readDouble();
		yAxisTickFont = (Font) in.readObject();
		yAxisTitle = (String) in.readObject();
		yAxisTitleFont = (Font) in.readObject();
		yAxisTitleOffset = in.readDouble();

		// Plot series attributes.
		seriesDots = (Dots) in.readObject();
		seriesStroke = Strokes.readExternal (in);
		seriesColor = (Color) in.readObject();
		seriesSmooth = in.readBoolean();

		// Label attributes.
		labelPosition = in.readInt();
		labelOffset = in.readDouble();
		labelFont = (Font) in.readObject();
		labelColor = (Color) in.readObject();
		labelBackground = (Color) in.readObject();

		// Lists of plot series and labels.
		myPlotSeries.clear();
		int n = in.readInt();
		for (int i = 0; i < n; ++ i)
			{
			myPlotSeries.add ((PlotSeries) in.readObject());
			}
		myLabels.clear();
		n = in.readInt();
		for (int i = 0; i < n; ++ i)
			{
			myLabels.add ((Label) in.readObject());
			}
		}

	/**
	 * Write the given plot to the file with the given name. <TT>thePlot</TT> is
	 * written in serialized form to the file.
	 *
	 * @param  thePlot      Plot.
	 * @param  theFileName  File name.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void write
		(Plot thePlot,
		 String theFileName)
		throws IOException
		{
		write (thePlot, new File (theFileName));
		}

	/**
	 * Write the given plot to the given file. <TT>thePlot</TT> is written in
	 * serialized form to <TT>theFile</TT>.
	 *
	 * @param  thePlot  Plot.
	 * @param  theFile  File.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public static void write
		(Plot thePlot,
		 File theFile)
		throws IOException
		{
		FileOutputStream fos = null;
		BufferedOutputStream bos = null;
		ObjectOutputStream oos = null;

		try
			{
			fos = new FileOutputStream (theFile);
			bos = new BufferedOutputStream (fos);
			oos = new ObjectOutputStream (bos);
			oos.writeObject (thePlot);
			oos.close();
			}

		catch (IOException exc)
			{
			if (fos != null)
				{
				try { fos.close(); } catch (IOException exc2) {}
				}
			throw exc;
			}
		}

	/**
	 * Read a plot from the file with the given name. The file must contain one
	 * instance of class Plot in serialized form; for example, as written by the
	 * static <TT>Plot.write()</TT> method.
	 *
	 * @param  theFileName  File name.
	 *
	 * @return  Plot.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize the plot cannot be found.
	 */
	public static Plot read
		(String theFileName)
		throws IOException, ClassNotFoundException
		{
		return read (new File (theFileName));
		}

	/**
	 * Read a plot from the given file. The file must contain one instance of
	 * class Plot in serialized form; for example, as written by the static
	 * <TT>Plot.write()</TT> method.
	 *
	 * @param  theFile  File.
	 *
	 * @return  Plot.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize the plot cannot be found.
	 */
	public static Plot read
		(File theFile)
		throws IOException, ClassNotFoundException
		{
		FileInputStream fis = null;
		BufferedInputStream bis = null;
		ObjectInputStream ois = null;
		Plot result = null;

		try
			{
			fis = new FileInputStream (theFile);
			bis = new BufferedInputStream (fis);
			ois = new ObjectInputStream (bis);
			result = (Plot) ois.readObject();
			ois.close();
			return result;
			}

		catch (IOException exc)
			{
			if (fis != null)
				{
				try { fis.close(); } catch (IOException exc2) {}
				}
			throw exc;
			}
		}

// Hidden operations.

	/**
	 * Get this plot's X-Y plot object.
	 *
	 * @return  X-Y plot object.
	 */
	private XYPlot getXYPlot()
		{
		if (myXYPlot != null) return myXYPlot;

		// Set up grid.
		Grid grid =
			new Grid
				(majorGridLines ? theMajorGridLineStroke : null,
				 Axis.DEFAULT_PAINT,
				 minorGridLines ? theMinorGridLineStroke : null,
				 Axis.DEFAULT_PAINT);

		// Find minimum and maximum X and Y values for all data points.
		double xmin = Double.POSITIVE_INFINITY;
		double xmax = Double.NEGATIVE_INFINITY;
		double ymin = Double.POSITIVE_INFINITY;
		double ymax = Double.NEGATIVE_INFINITY;
		for (PlotSeries plotseries : myPlotSeries)
			{
			XYSeries series = plotseries.getDataSeries();
			xmin = Math.min (xmin, series.minX());
			xmax = Math.max (xmax, series.maxX());
			ymin = Math.min (ymin, series.minY());
			ymax = Math.max (ymax, series.maxY());
			}

		// Set up X axis.
		NumericalAxis xAxis = null;
		double xcross = 0.0;
		Ticks xTicks =
			new Ticks
				(theAxisStroke,
				 Axis.DEFAULT_PAINT,
				 4.0,
				 xAxisTickFormat,
				 xAxisTickScale,
				 xAxisTickFont,
				 Axis.DEFAULT_PAINT,
				 4.0);
		switch (xAxisKind)
			{
			case LINEAR:
				if (Double.isNaN (xAxisStart))
					{
					xmin = xmin < 0.0 ? LinearAxis.autoscale (xmin) : 0.0;
					}
				else
					{
					xmin = xAxisStart;
					}
				if (Double.isNaN (xAxisEnd))
					{
					xmax = xmax > 0.0 ? LinearAxis.autoscale (xmax) : 0.0;
					}
				else
					{
					xmax = xAxisEnd;
					}
				if (Double.isNaN (xAxisCrossing))
					{
					xcross = xmin;
					}
				else
					{
					xcross = xAxisCrossing;
					}
				xAxis =
					new LinearAxis
						(xmin, xmax,
						 xAxisMajorDivisions, xAxisMinorDivisions,
						 xcross, xAxisLength,
						 theAxisStroke, Axis.DEFAULT_PAINT,
						 xTicks, (Ticks) null);
				break;
			case LOGARITHMIC:
				int[] xauto = null;
				if (Double.isNaN (xAxisStart) || Double.isNaN (xAxisEnd))
					{
					xauto = LogarithmicAxis.autoscale (xmin, xmax);
					}
				if (Double.isNaN (xAxisStart))
					{
					xmin = xauto[0];
					}
				else
					{
					xmin = xAxisStart;
					}
				if (Double.isNaN (xAxisEnd))
					{
					xmax = xauto[1];
					}
				else
					{
					xmax = xAxisEnd;
					}
				if (Double.isNaN (xAxisCrossing))
					{
					xcross = xmin;
					}
				else
					{
					xcross = xAxisCrossing;
					}
				xAxis =
					new LogarithmicAxis
						((int) xmin, (int) xmax,
						 xAxisMinorDivisions,
						 (int) xcross, xAxisLength,
						 theAxisStroke, Axis.DEFAULT_PAINT,
						 xTicks, (Ticks) null);
				break;
			}

		// Set up Y axis.
		NumericalAxis yAxis = null;
		double ycross = 0.0;
		Ticks yTicks =
			new Ticks
				(theAxisStroke,
				 Axis.DEFAULT_PAINT,
				 4.0,
				 yAxisTickFormat,
				 yAxisTickScale,
				 yAxisTickFont,
				 Axis.DEFAULT_PAINT,
				 4.0);
		switch (yAxisKind)
			{
			case LINEAR:
				if (Double.isNaN (yAxisStart))
					{
					ymin = ymin < 0.0 ? LinearAxis.autoscale (ymin) : 0.0;
					}
				else
					{
					ymin = yAxisStart;
					}
				if (Double.isNaN (yAxisEnd))
					{
					ymax = ymax > 0.0 ? LinearAxis.autoscale (ymax) : 0.0;
					}
				else
					{
					ymax = yAxisEnd;
					}
				if (Double.isNaN (yAxisCrossing))
					{
					ycross = ymin;
					}
				else
					{
					ycross = yAxisCrossing;
					}
				yAxis =
					new LinearAxis
						(ymin, ymax,
						 yAxisMajorDivisions, yAxisMinorDivisions,
						 ycross, yAxisLength,
						 theAxisStroke, Axis.DEFAULT_PAINT,
						 yTicks, (Ticks) null);
				break;
			case LOGARITHMIC:
				int[] yauto = null;
				if (Double.isNaN (yAxisStart) || Double.isNaN (yAxisEnd))
					{
					yauto = LogarithmicAxis.autoscale (ymin, ymax);
					}
				if (Double.isNaN (yAxisStart))
					{
					ymin = yauto[0];
					}
				else
					{
					ymin = yAxisStart;
					}
				if (Double.isNaN (yAxisEnd))
					{
					ymax = yauto[1];
					}
				else
					{
					ymax = yAxisEnd;
					}
				if (Double.isNaN (yAxisCrossing))
					{
					ycross = ymin;
					}
				else
					{
					ycross = yAxisCrossing;
					}
				yAxis =
					new LogarithmicAxis
						((int) ymin, (int) ymax,
						 yAxisMinorDivisions,
						 (int) ycross, yAxisLength,
						 theAxisStroke, Axis.DEFAULT_PAINT,
						 yTicks, (Ticks) null);
				break;
			}

		// Set up plot.
		XYPlot xyplot =
			new XYPlot
				(xAxis, yAxis, grid,
				 leftMargin, topMargin, rightMargin, bottomMargin);

		// Add plot series.
		for (PlotSeries plotseries : myPlotSeries)
			{
			xyplot.addPlotSeries (plotseries);
			}

		// Add labels.
		for (Label label : myLabels)
			{
			xyplot.addLabel (label);
			}

		// Add plot title.
		if (plotTitle != null)
			{
			xyplot.addLabel
				(new Label
					(plotTitle,
					 0.5*xAxisLength,
					 yAxisLength,
					 ABOVE+PIXEL_COORDINATES,
					 plotTitleOffset,
					 plotTitleFont));
			}

		// Add X axis title.
		if (xAxisTitle != null)
			{
			xyplot.addLabel
				(new Label
					(xAxisTitle,
					 0.5*xAxisLength,
					 0.0,
					 BELOW+PIXEL_COORDINATES,
					 xAxisTitleOffset,
					 xAxisTitleFont));
			}

		// Add Y axis title.
		if (yAxisTitle != null)
			{
			xyplot.addLabel
				(new Label
					(yAxisTitle,
					 0.0,
					 0.5*yAxisLength,
					 LEFT+ROTATE_LEFT+PIXEL_COORDINATES,
					 yAxisTitleOffset,
					 yAxisTitleFont));
			}

		myXYPlot = xyplot;
		return xyplot;
		}

	/**
	 * Format this plot's title interactively.
	 *
	 * @param  dialog  Plot title dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatPlotTitle
		(TitleDialog dialog)
		{
		dialog.setTitleText (plotTitle);
		dialog.setTitleFont (plotTitleFont);
		dialog.setTitleOffset (plotTitleOffset);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			plotTitle = dialog.getTitleText();
			plotTitleFont = dialog.getTitleFont();
			plotTitleOffset = dialog.getTitleOffset();
			myXYPlot = null;
			}
		return okay;
		}

	/**
	 * Format this plot's margins interactively.
	 *
	 * @param  dialog  Plot margins dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatPlotMargins
		(PlotMarginsDialog dialog)
		{
		dialog.setTopMargin (topMargin);
		dialog.setLeftMargin (leftMargin);
		dialog.setBottomMargin (bottomMargin);
		dialog.setRightMargin (rightMargin);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			topMargin = dialog.getTopMargin();
			leftMargin = dialog.getLeftMargin();
			bottomMargin = dialog.getBottomMargin();
			rightMargin = dialog.getRightMargin();
			myXYPlot = null;
			}
		return okay;
		}

	/**
	 * Format this plot's area interactively.
	 *
	 * @param  dialog  Plot area dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatPlotArea
		(PlotAreaDialog dialog)
		{
		dialog.setPlotAreaWidth (xAxisLength);
		dialog.setPlotAreaHeight (yAxisLength);
		dialog.setMajorGridLines (majorGridLines);
		dialog.setMinorGridLines (minorGridLines);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			xAxisLength = dialog.getPlotAreaWidth();
			yAxisLength = dialog.getPlotAreaHeight();
			majorGridLines = dialog.getMajorGridLines();
			minorGridLines = dialog.getMinorGridLines();
			myXYPlot = null;
			}
		return okay;
		}

	/**
	 * Format this plot's X axis interactively.
	 *
	 * @param  dialog  X axis dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatXAxis
		(AxisDialog dialog)
		{
		dialog.setAxisKind (xAxisKind);
		dialog.setAxisStart (xAxisStart);
		dialog.setAxisEnd (xAxisEnd);
		dialog.setAxisMajorDivisions (xAxisMajorDivisions);
		dialog.setAxisMinorDivisions (xAxisMinorDivisions);
		dialog.setAxisCrossing (xAxisCrossing);
		dialog.setTickFont (xAxisTickFont);
		dialog.setTickFormat (xAxisTickFormat);
		dialog.setTickScale (xAxisTickScale);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			xAxisKind = dialog.getAxisKind();
			xAxisStart = dialog.getAxisStart();
			xAxisEnd = dialog.getAxisEnd();
			xAxisMajorDivisions = dialog.getAxisMajorDivisions();
			xAxisMinorDivisions = dialog.getAxisMinorDivisions();
			xAxisCrossing = dialog.getAxisCrossing();
			xAxisTickFont = dialog.getTickFont();
			xAxisTickFormat = dialog.getTickFormat();
			xAxisTickScale = dialog.getTickScale();
			myXYPlot = null;
			}
		return okay;
		}

	/**
	 * Format this plot's X axis title interactively.
	 *
	 * @param  dialog  X axis Plot title dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatXAxisTitle
		(TitleDialog dialog)
		{
		dialog.setTitleText (xAxisTitle);
		dialog.setTitleFont (xAxisTitleFont);
		dialog.setTitleOffset (xAxisTitleOffset);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			xAxisTitle = dialog.getTitleText();
			xAxisTitleFont = dialog.getTitleFont();
			xAxisTitleOffset = dialog.getTitleOffset();
			myXYPlot = null;
			}
		return okay;
		}

	/**
	 * Format this plot's Y axis interactively.
	 *
	 * @param  dialog  Y axis dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatYAxis
		(AxisDialog dialog)
		{
		dialog.setAxisKind (yAxisKind);
		dialog.setAxisStart (yAxisStart);
		dialog.setAxisEnd (yAxisEnd);
		dialog.setAxisMajorDivisions (yAxisMajorDivisions);
		dialog.setAxisMinorDivisions (yAxisMinorDivisions);
		dialog.setAxisCrossing (yAxisCrossing);
		dialog.setTickFont (yAxisTickFont);
		dialog.setTickFormat (yAxisTickFormat);
		dialog.setTickScale (yAxisTickScale);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			yAxisKind = dialog.getAxisKind();
			yAxisStart = dialog.getAxisStart();
			yAxisEnd = dialog.getAxisEnd();
			yAxisMajorDivisions = dialog.getAxisMajorDivisions();
			yAxisMinorDivisions = dialog.getAxisMinorDivisions();
			yAxisCrossing = dialog.getAxisCrossing();
			yAxisTickFont = dialog.getTickFont();
			yAxisTickFormat = dialog.getTickFormat();
			yAxisTickScale = dialog.getTickScale();
			myXYPlot = null;
			}
		return okay;
		}

	/**
	 * Format this plot's Y axis title interactively.
	 *
	 * @param  dialog  Y axis Plot title dialog.
	 *
	 * @return  True if the user clicked OK, false if the user clicked Cancel.
	 */
	boolean formatYAxisTitle
		(TitleDialog dialog)
		{
		dialog.setTitleText (yAxisTitle);
		dialog.setTitleFont (yAxisTitleFont);
		dialog.setTitleOffset (yAxisTitleOffset);
		dialog.setVisible (true);
		boolean okay = dialog.isOkay();
		if (okay)
			{
			yAxisTitle = dialog.getTitleText();
			yAxisTitleFont = dialog.getTitleFont();
			yAxisTitleOffset = dialog.getTitleOffset();
			myXYPlot = null;
			}
		return okay;
		}

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		Plot plot1 = new Plot()
//			.plotTitle ("This Is A Test 1")
//			.xAxisTitle ("Number of Processors, K")
//			.yAxisTitle ("Running Time, T (sec)")
//			.xySeries (new double[] {0.15, 15}, new double[] {0.15, 15})
//			.labelPosition (RIGHT)
//			.labelOffset (6)
//			.label ("N = 1000", 15, 15)
//			.display();
//		Plot plot2 = new Plot()
//			.plotTitle ("This Is A Test 2")
//			.xAxisTitle ("Number of Processors, K")
//			.xAxisKind (LOGARITHMIC)
//			.xAxisMinorDivisions (10)
//			.xAxisTickFormat (new DecimalFormat ("0.0"))
//			.yAxisTitle ("Running Time, T (sec)")
//			.yAxisKind (LOGARITHMIC)
//			.yAxisMinorDivisions (10)
//			.yAxisTickFormat (new DecimalFormat ("0.0"))
//			.minorGridLines (true)
//			.xySeries (new double[] {0.15, 15}, new double[] {0.15, 15})
//			.display();
//		Plot plot3 = new Plot()
//			.plotTitle ("This Is A Test 3")
//			.xAxisTitle ("X")
//			.yAxisTitle ("Y")
//			.segmentedSeries
//				(new double[] {1, 3, 7, 9}, new double[] {1, 1, 9, 9})
//			.display();
//		}

	}
