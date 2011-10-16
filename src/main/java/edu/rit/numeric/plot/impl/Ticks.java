//******************************************************************************
//
// File:    Ticks.java
// Package: edu.rit.numeric.plot.impl
// Unit:    Class edu.rit.numeric.plot.impl.Ticks
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
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;

import java.awt.font.FontRenderContext;
import java.awt.font.LineMetrics;

import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

import java.text.DecimalFormat;

/**
 * Class Ticks provides information about tick marks to be drawn on a
 * {@linkplain NumericalAxis}.
 *
 * @author  Alan Kaminsky
 * @version 11-Jun-2007
 */
public class Ticks
	{

// Exported constants.

	/**
	 * The default tick mark stroke (solid, width=1).
	 */
	public static final Stroke DEFAULT_STROKE = Strokes.solid (1);

	/**
	 * The default tick mark paint (black).
	 */
	public static final Paint DEFAULT_PAINT = Color.black;

	/**
	 * The default tick mark length (5).
	 */
	public static final double DEFAULT_LENGTH = 5.0;

	/**
	 * The default tick label font (sanserif, plain, 12).
	 */
	public static final Font DEFAULT_LABEL_FONT =
		new Font ("sanserif", Font.PLAIN, 12);

	/**
	 * The default tick label paint (black).
	 */
	public static final Paint DEFAULT_LABEL_PAINT = Color.black;

	/**
	 * The default offset from the tick mark to the tick label (5).
	 */
	public static final double DEFAULT_LABEL_OFFSET = 5.0;

// Hidden data members.

	/**
	 * Stroke for drawing tick marks on the display.
	 */
	protected Stroke myStroke;

	/**
	 * Paint for drawing tick marks on the display.
	 */
	protected Paint myPaint;

	/**
	 * Length of the tick marks on the display.
	 */
	protected double myLength;

	/**
	 * Format for the tick labels, or null for no tick labels.
	 */
	protected DecimalFormat myLabelFormat;

	/**
	 * Scale factor for the tick labels.
	 */
	protected double myLabelScale;

	/**
	 * Font for drawing tick labels on the display, if any.
	 */
	protected Font myLabelFont;

	/**
	 * Paint for drawing tick labels on the display, if any.
	 */
	protected Paint myLabelPaint;

	/**
	 * Offset from the tick marks to the tick labels, if any.
	 */
	protected double myLabelOffset;

// Exported constructors.

	/**
	 * Construct a new ticks information object. The tick marks have the default
	 * stroke, the default paint, and the default length. There are no tick
	 * labels.
	 */
	public Ticks()
		{
		this
			(DEFAULT_STROKE, DEFAULT_PAINT, DEFAULT_LENGTH,
			 null, 1.0, DEFAULT_LABEL_FONT,
			 DEFAULT_LABEL_PAINT, DEFAULT_LABEL_OFFSET);
		}

	/**
	 * Construct a new ticks information object. The tick marks have the default
	 * stroke, the default paint, and the default length. The tick labels have
	 * the given decimal format, no scaling, the default font, the default
	 * paint, and the default offset from the tick marks.
	 *
	 * @param  theLabelFormat
	 *     Format for the tick labels. If null, no tick labels are drawn.
	 */
	public Ticks
		(DecimalFormat theLabelFormat)
		{
		this
			(DEFAULT_STROKE, DEFAULT_PAINT, DEFAULT_LENGTH,
			 theLabelFormat, 1.0, DEFAULT_LABEL_FONT,
			 DEFAULT_LABEL_PAINT, DEFAULT_LABEL_OFFSET);
		}

	/**
	 * Construct a new ticks information object. The tick marks have the default
	 * stroke, the default paint, and the default length. The tick labels have
	 * the given decimal format, the given scaling, the default font, the
	 * default paint, and the default offset from the tick marks.
	 *
	 * @param  theLabelFormat
	 *     Format for the tick labels. If null, no tick labels are drawn.
	 * @param  theLabelScale
	 *     Scale factor for the tick labels. Each tick label's numerical value
	 *     is divided by <TT>theLabelScale</TT> before being displayed. Use a
	 *     value of 1.0 for no scaling.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theLabelFormat</TT> is non-null
	 *     and <TT>theLabelScale</TT> is equal to 0.
	 */
	public Ticks
		(DecimalFormat theLabelFormat,
		 double theLabelScale)
		{
		this
			(DEFAULT_STROKE, DEFAULT_PAINT, DEFAULT_LENGTH,
			 theLabelFormat, theLabelScale, DEFAULT_LABEL_FONT,
			 DEFAULT_LABEL_PAINT, DEFAULT_LABEL_OFFSET);
		}

	/**
	 * Construct a new ticks information object. The tick marks have the given
	 * stroke, the given paint, and the given length. The tick labels have the
	 * given decimal format, the given scaling, the given font, the given paint,
	 * and the given offset from the tick marks.
	 *
	 * @param  theStroke
	 *     Stroke for drawing tick marks on the display.
	 * @param  thePaint
	 *     Paint for drawing tick marks on the display.
	 * @param  theLength
	 *     Length of the tick marks on the display.
	 * @param  theLabelFormat
	 *     Format for the tick labels. If null, no tick labels are drawn.
	 * @param  theLabelScale
	 *     Scale factor for the tick labels. Each tick label's numerical value
	 *     is divided by <TT>theLabelScale</TT> before being displayed. Use a
	 *     value of 1.0 for no scaling.
	 * @param  theLabelFont
	 *     Font for drawing tick labels on the display, if any. Ignored if no
	 *     tick labels are to be drawn.
	 * @param  theLabelPaint
	 *     Paint for drawing tick labels on the display, if any. Ignored if no
	 *     tick labels are to be drawn.
	 * @param  theLabelOffset
	 *     Offset from the tick marks to the tick labels, if any. Ignored if no
	 *     tick labels are to be drawn.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStroke</TT> is null,
	 *     <TT>thePaint</TT> is null, <TT>theLabelFormat</TT> is non-null and
	 *     <TT>theLabelFont</TT> is null, or <TT>theLabelFormat</TT> is non-null
	 *     and <TT>theLabelPaint</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theLength</TT> is less than or
	 *     equal to 0, <TT>theLabelFormat</TT> is non-null and
	 *     <TT>theLabelScale</TT> is equal to 0, or <TT>theLabelFormat</TT> is
	 *     non-null and <TT>theLabelOffset</TT> is less than or equal to 0.
	 */
	public Ticks
		(Stroke theStroke,
		 Paint thePaint,
		 double theLength,
		 DecimalFormat theLabelFormat,
		 double theLabelScale,
		 Font theLabelFont,
		 Paint theLabelPaint,
		 double theLabelOffset)
		{
		if
			(theStroke == null ||
			 thePaint == null ||
			 (theLabelFormat != null && theLabelFont == null) ||
			 (theLabelFormat != null && theLabelPaint == null))
			{
			throw new NullPointerException();
			}
		if
			(theLength <= 0.0 ||
			 (theLabelFormat != null && theLabelScale == 0.0) ||
			 (theLabelFormat != null && theLabelOffset <= 0.0))
			{
			throw new IllegalArgumentException();
			}
		myStroke = theStroke;
		myPaint = thePaint;
		myLength = theLength;
		myLabelFormat = theLabelFormat;
		if (theLabelFormat == null)
			{
			myLabelScale = 1.0;
			myLabelFont = null;
			myLabelPaint = null;
			myLabelOffset = 0.0;
			}
		else
			{
			myLabelScale = theLabelScale;
			myLabelFont = theLabelFont;
			myLabelPaint = theLabelPaint;
			myLabelOffset = theLabelOffset;
			}
		}

// Exported operations.

	/**
	 * Returns the stroke for drawing tick marks on the display.
	 */
	public Stroke getStroke()
		{
		return myStroke;
		}

	/**
	 * Returns the paint for drawing tick marks on the display.
	 */
	public Paint getPaint()
		{
		return myPaint;
		}

	/**
	 * Returns the length of the tick marks on the display.
	 */
	public double getLength()
		{
		return myLength;
		}

	/**
	 * Returns the format for the tick labels. If no tick labels are to be
	 * drawn, null is returned.
	 */
	public DecimalFormat getLabelFormat()
		{
		return myLabelFormat;
		}

	/**
	 * Returns the scale factor for the tick labels. Each tick label's numerical
	 * value is divided by the scale factor before being displayed. If no tick
	 * labels are to be drawn, 1.0 is returned.
	 */
	public double getLabelScale()
		{
		return myLabelScale;
		}

	/**
	 * Returns the font for drawing tick labels on the display. If no tick
	 * labels are to be drawn, null is returned.
	 */
	public Font getLabelFont()
		{
		return myLabelFont;
		}

	/**
	 * Returns the paint for drawing tick labels on the display. If no tick
	 * labels are to be drawn, null is returned.
	 */
	public Paint getLabelPaint()
		{
		return myLabelPaint;
		}

	/**
	 * Returns the offset from the tick marks to the tick labels. If no tick
	 * labels are to be drawn, 0 is returned.
	 */
	public double getLabelOffset()
		{
		return myLabelOffset;
		}

	/**
	 * Draw a tick to the left. The tick mark is drawn in the given graphics
	 * context starting at coordinates (0, 0) and going to the left. There is no
	 * tick label.
	 *
	 * @param  g2d            2-D graphics context.
	 */
	public void drawLeft
		(Graphics2D g2d)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, -myLength, 0));

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Draw a tick to the left. The tick mark is drawn in the given graphics
	 * context starting at coordinates (0, 0) and going to the left. The tick
	 * label, if any, is drawn to the left of the tick mark with the given
	 * numerical value.
	 *
	 * @param  g2d            2-D graphics context.
	 * @param  theLabelValue  Numerical value for the tick label, if any.
	 */
	public void drawLeft
		(Graphics2D g2d,
		 double theLabelValue)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		Font oldFont = g2d.getFont();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, -myLength, 0));

		// Draw tick label if any.
		if (myLabelFormat != null)
			{
			g2d.setFont (myLabelFont);
			g2d.setPaint (myLabelPaint);
			String s = myLabelFormat.format (theLabelValue / myLabelScale);
			FontRenderContext frc = g2d.getFontRenderContext();
			Rectangle2D bounds = myLabelFont.getStringBounds (s, frc);
			LineMetrics metrics = myLabelFont.getLineMetrics (s, frc);
			g2d.drawString
				(s,
				 (float) (- myLength - myLabelOffset - bounds.getWidth()),
				 (float) (metrics.getAscent() / 2));
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		g2d.setFont (oldFont);
		}

	/**
	 * Draw a tick to the right. The tick mark is drawn in the given graphics
	 * context starting at coordinates (0, 0) and going to the right. There is
	 * no tick label.
	 *
	 * @param  g2d            2-D graphics context.
	 */
	public void drawRight
		(Graphics2D g2d)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, myLength, 0));

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Draw a tick to the right. The tick mark is drawn in the given graphics
	 * context starting at coordinates (0, 0) and going to the right. The tick
	 * label, if any, is drawn to the right of the tick mark with the given
	 * numerical value.
	 *
	 * @param  g2d            2-D graphics context.
	 * @param  theLabelValue  Numerical value for the tick label, if any.
	 */
	public void drawRight
		(Graphics2D g2d,
		 double theLabelValue)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		Font oldFont = g2d.getFont();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, myLength, 0));

		// Draw tick label if any.
		if (myLabelFormat != null)
			{
			g2d.setFont (myLabelFont);
			g2d.setPaint (myLabelPaint);
			String s = myLabelFormat.format (theLabelValue / myLabelScale);
			FontRenderContext frc = g2d.getFontRenderContext();
			Rectangle2D bounds = myLabelFont.getStringBounds (s, frc);
			LineMetrics metrics = myLabelFont.getLineMetrics (s, frc);
			g2d.drawString
				(s,
				 (float) (myLength + myLabelOffset),
				 (float) (metrics.getAscent() / 2));
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		g2d.setFont (oldFont);
		}

	/**
	 * Draw a tick below. The tick mark is drawn in the given graphics context
	 * starting at coordinates (0, 0) and going downwards. There is no tick
	 * label.
	 *
	 * @param  g2d            2-D graphics context.
	 */
	public void drawBelow
		(Graphics2D g2d)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, 0, myLength));

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Draw a tick below. The tick mark is drawn in the given graphics context
	 * starting at coordinates (0, 0) and going downwards. The tick label, if
	 * any, is drawn below the tick mark with the given numerical value.
	 *
	 * @param  g2d            2-D graphics context.
	 * @param  theLabelValue  Numerical value for the tick label, if any.
	 */
	public void drawBelow
		(Graphics2D g2d,
		 double theLabelValue)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		Font oldFont = g2d.getFont();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, 0, myLength));

		// Draw tick label if any.
		if (myLabelFormat != null)
			{
			g2d.setFont (myLabelFont);
			g2d.setPaint (myLabelPaint);
			String s = myLabelFormat.format (theLabelValue / myLabelScale);
			FontRenderContext frc = g2d.getFontRenderContext();
			Rectangle2D bounds = myLabelFont.getStringBounds (s, frc);
			LineMetrics metrics = myLabelFont.getLineMetrics (s, frc);
			g2d.drawString
				(s,
				 (float) (- bounds.getWidth() / 2),
				 (float) (myLength + myLabelOffset + metrics.getAscent()));
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		g2d.setFont (oldFont);
		}

	/**
	 * Draw a tick above. The tick mark is drawn in the given graphics context
	 * starting at coordinates (0, 0) and going upwards. There is no tick label.
	 *
	 * @param  g2d            2-D graphics context.
	 */
	public void drawAbove
		(Graphics2D g2d)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, 0, - myLength));

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		}

	/**
	 * Draw a tick above. The tick mark is drawn in the given graphics context
	 * starting at coordinates (0, 0) and going upwards. The tick label, if any,
	 * is drawn above the tick mark with the given numerical value.
	 *
	 * @param  g2d            2-D graphics context.
	 * @param  theLabelValue  Numerical value for the tick label, if any.
	 */
	public void drawAbove
		(Graphics2D g2d,
		 double theLabelValue)
		{
		// Save graphics context.
		Stroke oldStroke = g2d.getStroke();
		Paint oldPaint = g2d.getPaint();
		Font oldFont = g2d.getFont();

		// Draw tick mark.
		g2d.setStroke (myStroke);
		g2d.setPaint (myPaint);
		g2d.draw (new Line2D.Double (0, 0, 0, - myLength));

		// Draw tick label if any.
		if (myLabelFormat != null)
			{
			g2d.setFont (myLabelFont);
			g2d.setPaint (myLabelPaint);
			String s = myLabelFormat.format (theLabelValue / myLabelScale);
			FontRenderContext frc = g2d.getFontRenderContext();
			Rectangle2D bounds = myLabelFont.getStringBounds (s, frc);
			LineMetrics metrics = myLabelFont.getLineMetrics (s, frc);
			g2d.drawString
				(s,
				 (float) (- bounds.getWidth() / 2),
				 (float) (- myLength - myLabelOffset - metrics.getDescent()));
			}

		// Restore graphics context.
		g2d.setStroke (oldStroke);
		g2d.setPaint (oldPaint);
		g2d.setFont (oldFont);
		}

	}
