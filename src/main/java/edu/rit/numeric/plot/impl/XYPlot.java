//******************************************************************************
//
// File:    XYPlot.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.XYPlot
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

import edu.rit.numeric.CurveSmoothing;

import edu.rit.swing.Displayable;
import edu.rit.swing.Drawable;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;

import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;

import java.util.ArrayList;

/**
 * Class XYPlot provides an XY plot. The plot has an X axis, a Y axis, and one
 * or more plot series (type {@linkplain PlotSeries}). The plot may optionally
 * have a {@linkplain Grid}. Labels (class {@linkplain Label}) may be added to a
 * plot and are drawn after drawing all the plot series.
 *
 * @author  Alan Kaminsky
 * @version 06-Oct-2009
 */
public class XYPlot
	implements Displayable
	{

// Exported constants.

	/**
	 * The default margin (54).
	 */
	public static final double DEFAULT_MARGIN = 54.0;

	/**
	 * The default background paint (white).
	 */
	public static final Paint DEFAULT_BACKGROUND = Color.white;

// Hidden data members.

	/**
	 * X axis.
	 */
	protected NumericalAxis myXAxis;

	/**
	 * Y axis.
	 */
	protected NumericalAxis myYAxis;

	/**
	 * Grid.
	 */
	protected Grid myGrid;

	/**
	 * List of plot series.
	 */
	protected ArrayList<PlotSeries> myPlotSeries = new ArrayList<PlotSeries>();

	/**
	 * List of labels.
	 */
	protected ArrayList<Label> myLabels = new ArrayList<Label>();

	/**
	 * Left margin.
	 */
	protected double myLeftMargin;

	/**
	 * Top margin.
	 */
	protected double myTopMargin;

	/**
	 * Right margin.
	 */
	protected double myRightMargin;

	/**
	 * Bottom margin.
	 */
	protected double myBottomMargin;

	/**
	 * Width.
	 */
	protected double myWidth;

	/**
	 * Height.
	 */
	protected double myHeight;

	/**
	 * Bounding box.
	 */
	protected Rectangle2D.Double myBoundingBox;

	/**
	 * Background paint.
	 */
	protected Paint myBackgroundPaint;

// Hidden constructors.

	/**
	 * Construct a new XY plot with the given axes. The plot has no grid. The
	 * plot uses the default margin (54) on all sides. The plot uses the default
	 * background (white).
	 *
	 * @param  theXAxis  X axis.
	 * @param  theYAxis  Y axis.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theXAxis</TT> is null or
	 *     <TT>theYAxis</TT> is null.
	 */
	public XYPlot
		(NumericalAxis theXAxis,
		 NumericalAxis theYAxis)
		{
		this
			(theXAxis, theYAxis, null,
			 DEFAULT_MARGIN, DEFAULT_MARGIN, DEFAULT_MARGIN, DEFAULT_MARGIN,
			 DEFAULT_BACKGROUND);
		}

	/**
	 * Construct a new XY plot with the given axes and grid. The plot uses the
	 * default margin (54) on all sides. The plot uses the default background
	 * (white).
	 *
	 * @param  theXAxis  X axis.
	 * @param  theYAxis  Y axis.
	 * @param  theGrid   Grid, or null for no grid.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theXAxis</TT> is null or
	 *     <TT>theYAxis</TT> is null.
	 */
	public XYPlot
		(NumericalAxis theXAxis,
		 NumericalAxis theYAxis,
		 Grid theGrid)
		{
		this
			(theXAxis, theYAxis, theGrid,
			 DEFAULT_MARGIN, DEFAULT_MARGIN, DEFAULT_MARGIN, DEFAULT_MARGIN,
			 DEFAULT_BACKGROUND);
		}

	/**
	 * Construct a new XY plot with the given axes, grid, and margin. The plot
	 * uses the same margin on all sides. The plot uses the default background
	 * (white).
	 *
	 * @param  theXAxis   X axis.
	 * @param  theYAxis   Y axis.
	 * @param  theGrid    Grid, or null for no grid.
	 * @param  theMargin  Margin.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theXAxis</TT> is null or
	 *     <TT>theYAxis</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theMargin</TT> is less than 0.
	 */
	public XYPlot
		(NumericalAxis theXAxis,
		 NumericalAxis theYAxis,
		 Grid theGrid,
		 double theMargin)
		{
		this
			(theXAxis, theYAxis, theGrid,
			 theMargin, theMargin, theMargin, theMargin,
			 DEFAULT_BACKGROUND);
		}

	/**
	 * Construct a new XY plot with the given axes, grid, and margins. The plot
	 * uses the default background (white).
	 *
	 * @param  theXAxis         X axis.
	 * @param  theYAxis         Y axis.
	 * @param  theGrid          Grid, or null for no grid.
	 * @param  theLeftMargin    Left margin.
	 * @param  theTopMargin     Top margin.
	 * @param  theRightMargin   Right margin.
	 * @param  theBottomMargin  Left margin.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theXAxis</TT> is null or
	 *     <TT>theYAxis</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any margin is less than 0.
	 */
	public XYPlot
		(NumericalAxis theXAxis,
		 NumericalAxis theYAxis,
		 Grid theGrid,
		 double theLeftMargin,
		 double theTopMargin,
		 double theRightMargin,
		 double theBottomMargin)
		{
		this
			(theXAxis, theYAxis, theGrid,
			 theLeftMargin, theTopMargin, theRightMargin, theBottomMargin,
			 DEFAULT_BACKGROUND);
		}

	/**
	 * Construct a new XY plot with the given axes, grid, margins, and
	 * background.
	 *
	 * @param  theXAxis            X axis.
	 * @param  theYAxis            Y axis.
	 * @param  theGrid             Grid, or null for no grid.
	 * @param  theLeftMargin       Left margin.
	 * @param  theTopMargin        Top margin.
	 * @param  theRightMargin      Right margin.
	 * @param  theBottomMargin     Left margin.
	 * @param  theBackgroundPaint  Background paint.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theXAxis</TT> is null,
	 *     <TT>theYAxis</TT> is null, or <TT>theBackgroundPaint</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if any margin is less than 0.
	 */
	public XYPlot
		(NumericalAxis theXAxis,
		 NumericalAxis theYAxis,
		 Grid theGrid,
		 double theLeftMargin,
		 double theTopMargin,
		 double theRightMargin,
		 double theBottomMargin,
		 Paint theBackgroundPaint)
		{
		if (theXAxis == null || theYAxis == null || theBackgroundPaint == null)
			{
			throw new NullPointerException();
			}
		if (theLeftMargin < 0.0 || theTopMargin < 0.0 || theRightMargin < 0.0 ||
					theBottomMargin < 0.0)
			{
			throw new IllegalArgumentException();
			}
		myXAxis = theXAxis;
		myYAxis = theYAxis;
		myGrid = theGrid;
		myLeftMargin = theLeftMargin;
		myTopMargin = theTopMargin;
		myRightMargin = theRightMargin;
		myBottomMargin = theBottomMargin;
		myWidth = myLeftMargin + myXAxis.getLength() + myRightMargin;
		myHeight = myTopMargin + myYAxis.getLength() + myBottomMargin;
		myBoundingBox = new Rectangle2D.Double (0, 0, myWidth, myHeight);
		myBackgroundPaint = theBackgroundPaint;
		}

// Exported operations.

	/**
	 * Add the given plot series to this XY plot. When this XY plot is drawn,
	 * the plot series are drawn in the order they were added; that is, last
	 * added plot series on top.
	 *
	 * @param  thePlotSeries  Plot series.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlotSeries</TT> is null.
	 */
	public void addPlotSeries
		(PlotSeries thePlotSeries)
		{
		if (thePlotSeries == null)
			{
			throw new NullPointerException();
			}
		myPlotSeries.add (thePlotSeries);
		}

	/**
	 * Add the given label to this XY plot. When this XY plot is drawn, the
	 * labels are drawn after the plot series, and the labels are drawn in the
	 * order they were added; that is, last added label on top.
	 *
	 * @param  theLabel  Label.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theLabel</TT> is null.
	 */
	public void addLabel
		(Label theLabel)
		{
		if (theLabel == null)
			{
			throw new NullPointerException();
			}
		myLabels.add (theLabel);
		}

	/**
	 * Draw this drawable object in the given graphics context. Upon return from
	 * this method, the given graphics context's state (color, font, transform,
	 * clip, and so on) is the same as it was upon entry to this method.
	 * <P>
	 * This XY plot is drawn so that the upper left corner (including the
	 * margins) is located at the display coordinates (0, 0) in the graphics
	 * context.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		double sx1, sy1, sx2, sy2;
		double[] x1, y1, x2, y2;
		double dx1, dy1, dx2, dy2;

		// Save graphics context.
		AffineTransform oldTransform = g2d.getTransform();
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		Shape oldClip = g2d.getClip();

		// Compute origin of plotting area.
		double originx = myLeftMargin;
		double originy = myHeight - myBottomMargin;
		
		// Compute clip rectangle for plotting area.
		double plotWidth = myWidth - myLeftMargin - myRightMargin;
		double plotHeight = myHeight - myTopMargin - myBottomMargin;
		Shape plottingArea =
			new Rectangle2D.Double
				(/*x*/ myLeftMargin - 3,
				 /*y*/ myTopMargin - 3,
				 /*w*/ plotWidth + 6,
				 /*h*/ plotHeight + 6);

		// Draw grid if any.
		if (myGrid != null)
			{
			g2d.translate (originx, originy);
			myGrid.drawHorizontalGridlines (g2d, myYAxis, plotWidth);
			myGrid.drawVerticalGridlines (g2d, myXAxis, plotHeight);
			g2d.setTransform (oldTransform);
			}

		// Draw X axis.
		double xaxisy = myYAxis.getDisplayDistance (myYAxis.getCrossing());
		g2d.translate (originx, originy - xaxisy);
		myXAxis.drawHorizontal (g2d);
		g2d.setTransform (oldTransform);

		// Draw Y axis.
		double yaxisx =  myXAxis.getDisplayDistance (myXAxis.getCrossing());
		g2d.translate (originx + yaxisx, originy);
		myYAxis.drawVertical (g2d);
		g2d.setTransform (oldTransform);

		// Draw each plot series.
		for (PlotSeries thePlotSeries : myPlotSeries)
			{
			Drawable theDots = thePlotSeries.getDots();
			Stroke theStroke = thePlotSeries.getStroke();
			Paint thePaint = thePlotSeries.getColor();

			// Draw line segments or smooth Bezier curves, clipped to the
			// plotting area.
			if (theStroke != null && thePaint != null)
				{
				// If there are no lines, draw nothing.
				int n = thePlotSeries.getLineCount();
				if (n > 0)
					{
					g2d.clip (plottingArea);
					g2d.translate (originx, originy);
					g2d.setStroke (theStroke);
					g2d.setPaint (thePaint);

					// Get starting and ending coordinates of each segment.
					x1 = new double [n];
					y1 = new double [n];
					x2 = new double [n];
					y2 = new double [n];
					for (int i = 0; i < n; ++ i)
						{
						sx1 = thePlotSeries.getLineX1 (i);
						sy1 = thePlotSeries.getLineY1 (i);
						sx2 = thePlotSeries.getLineX2 (i);
						sy2 = thePlotSeries.getLineY2 (i);
						x1[i] = myXAxis.getDisplayDistance (sx1);
						y1[i] = -myYAxis.getDisplayDistance (sy1);
						x2[i] = myXAxis.getDisplayDistance (sx2);
						y2[i] = -myYAxis.getDisplayDistance (sy2);
						}

					// Get and draw a general path object.
					g2d.draw
						(thePlotSeries.isSmooth() ?
							getCurvedPath (x1, y1, x2, y2, n) :
							getStraightPath (x1, y1, x2, y2, n));

					g2d.setTransform (oldTransform);
					g2d.setClip (oldClip);
					}
				}

			// Draw each dot that falls within the plotting area.
			if (theDots != null)
				{
				int n = thePlotSeries.getDotCount();
				for (int i = 0; i < n; ++ i)
					{
					// Get point location.
					sx1 = thePlotSeries.getDotX (i);
					sy1 = thePlotSeries.getDotY (i);
					if (myXAxis.includesValue (sx1) &&
								myYAxis.includesValue (sy1))
						{
						// Translate the display origin to the point's location.
						dx1 = myXAxis.getDisplayDistance (sx1);
						dy1 = myYAxis.getDisplayDistance (sy1);
						g2d.translate (originx + dx1, originy - dy1);

						// Draw the dot.
						theDots.draw (g2d);
						g2d.setTransform (oldTransform);
						}
					}
				}
			}
		
		// Draw each label.
		for (Label theLabel : myLabels)
			{
			// Translate the display origin to the point's location.
			if (theLabel.getPositioning() / Label.PIXEL_COORDINATES == 1)
				{
				dx1 = theLabel.getLocationX();
				dy1 = theLabel.getLocationY();
				}
			else
				{
				dx1 = myXAxis.getDisplayDistance (theLabel.getLocationX());
				dy1 = myYAxis.getDisplayDistance (theLabel.getLocationY());
				}
			g2d.translate (originx + dx1, originy - dy1);

			// Draw the label.
			theLabel.draw (g2d);
			g2d.setTransform (oldTransform);
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Returns a path consisting of straight line segments.
	 */
	private static GeneralPath getStraightPath
		(double[] x1,
		 double[] y1,
		 double[] x2,
		 double[] y2,
		 int n)
		{
		GeneralPath path = new GeneralPath();
		path.moveTo ((float) x1[0], (float) y1[0]);
		for (int i = 0; i < n; ++ i)
			{
			path.lineTo ((float) x2[i], (float) y2[i]);
			if (i < n-1 && (x2[i] != x1[i+1] || y2[i] != y1[i+1]))
				{
				path.moveTo ((float) x1[i+1], (float) y1[i+1]);
				}
			}
		return path;
		}

	/**
	 * Returns a path consisting of Bezier curve segments.
	 */
	private static GeneralPath getCurvedPath
		(double[] x1,
		 double[] y1,
		 double[] x2,
		 double[] y2,
		 int n)
		{
		GeneralPath path = new GeneralPath();

		// Working storage for Bezier control points.
		double[] xu = new double [n+1];
		double[] xa = new double [n+1];
		double[] xc = new double [n+1];
		double[] yu = new double [n+1];
		double[] ya = new double [n+1];
		double[] yc = new double [n+1];

		// Iterate over all sequences of contiguous segments. first = index of
		// first segment in sequence. last = 1 + index of last segment in
		// sequence.
		int first = 0;
		int last;
		while (first < n)
			{
			// Scan until we hit a noncontiguous segment.
			xu[first] = x1[first];
			yu[first] = y1[first];
			last = first + 1;
			while (last < n && x2[last-1] == x1[last] && y2[last-1] == y1[last])
				{
				xu[last] = x1[last];
				yu[last] = y1[last];
				++ last;
				}
			xu[last] = x2[last-1];
			yu[last] = y2[last-1];

			// Get number of segments in sequence.
			int len = last - first;
			if (len == 1)
				{
				// Only one segment in sequence. Draw a straight line.
				path.moveTo ((float) xu[first], (float) yu[first]);
				path.lineTo ((float) xu[last], (float) yu[last]);
				}
			else
				{
				// More than one segment in sequence. Compute Bezier control
				// points.
				CurveSmoothing.computeBezierOpen (xu, xa, xc, first, len+1);
				CurveSmoothing.computeBezierOpen (yu, ya, yc, first, len+1);

				// Draw Bezier curves.
				path.moveTo ((float) xu[first], (float) yu[first]);
				for (int i = first; i < last; ++ i)
					{
					path.curveTo
						((float) xa[i], (float) ya[i],
						 (float) xc[i], (float) yc[i],
						 (float) xu[i+1], (float) yu[i+1]);
					}
				}

			// Go to the next sequence of contiguous segments.
			first = last;
			}

		return path;
		}

	/**
	 * Returns this displayable object's bounding box. This is the smallest
	 * rectangle that encloses all of this displayable object.
	 */
	public Rectangle2D getBoundingBox()
		{
		return myBoundingBox;
		}

	/**
	 * Returns this displayable object's background paint.
	 */
	public Paint getBackgroundPaint()
		{
		return myBackgroundPaint;
		}

	}
