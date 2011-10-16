//******************************************************************************
//
// File:    Dots.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.Dots
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

import edu.rit.swing.Drawable;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;

import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Dots provides operations for creating several common shapes of dots for
 * plotting points.
 *
 * @author  Alan Kaminsky
 * @version 11-Dec-2007
 */
public abstract class Dots
	implements Drawable, Externalizable
	{

// Hidden helper classes.

	private static class Ellipse
		extends Dots
		{
		private static final long serialVersionUID = -3905475094865863992L;

		private Color myFillColor;
		private BasicStroke myOutlineStroke;
		private Color myOutlineColor;
		private Ellipse2D.Double myEllipse;

		public Ellipse()
			{
			}

		public Ellipse
			(Color theFillColor,
			 BasicStroke theOutlineStroke,
			 Color theOutlineColor,
			 double theWidth,
			 double theHeight)
			{
			myFillColor = theFillColor;
			myOutlineStroke = theOutlineStroke;
			myOutlineColor = theOutlineColor;
			myEllipse = new Ellipse2D.Double
				(- theWidth/2, - theHeight/2, theWidth, theHeight);
			}

		public void draw
			(Graphics2D g2d)
			{
			Paint oldPaint = g2d.getPaint();
			if (myFillColor != null)
				{
				g2d.setPaint (myFillColor);
				g2d.fill (myEllipse);
				}
			if (myOutlineStroke != null && myOutlineColor != null)
				{
				Stroke oldStroke = g2d.getStroke();
				g2d.setStroke (myOutlineStroke);
				g2d.setPaint (myOutlineColor);
				g2d.draw (myEllipse);
				g2d.setStroke (oldStroke);
				}
			g2d.setPaint (oldPaint);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeObject (myFillColor);
			out.writeObject (myOutlineColor);
			Strokes.writeExternal (myOutlineStroke, out);
			out.writeDouble (myEllipse.x);
			out.writeDouble (myEllipse.y);
			out.writeDouble (myEllipse.width);
			out.writeDouble (myEllipse.height);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException, ClassNotFoundException
			{
			myFillColor = (Color) in.readObject();
			myOutlineColor = (Color) in.readObject();
			myOutlineStroke = Strokes.readExternal (in);
			myEllipse = new Ellipse2D.Double
				(in.readDouble(),
				 in.readDouble(),
				 in.readDouble(),
				 in.readDouble());
			}
		}

	private static class Rectangle
		extends Dots
		{
		private static final long serialVersionUID = 3000182552146442565L;

		private Color myFillColor;
		private BasicStroke myOutlineStroke;
		private Color myOutlineColor;
		private Rectangle2D.Double myRectangle;
		private double myAngle;

		public Rectangle()
			{
			}

		public Rectangle
			(Color theFillColor,
			 BasicStroke theOutlineStroke,
			 Color theOutlineColor,
			 double theWidth,
			 double theHeight,
			 double theAngle)
			{
			myFillColor = theFillColor;
			myOutlineStroke = theOutlineStroke;
			myOutlineColor = theOutlineColor;
			myRectangle = new Rectangle2D.Double
				(- theWidth/2, - theHeight/2, theWidth, theHeight);
			myAngle = theAngle;
			}

		public void draw
			(Graphics2D g2d)
			{
			Paint oldPaint = g2d.getPaint();
			AffineTransform oldTransform = g2d.getTransform();
			g2d.rotate (myAngle);
			if (myFillColor != null)
				{
				g2d.setPaint (myFillColor);
				g2d.fill (myRectangle);
				}
			if (myOutlineStroke != null && myOutlineColor != null)
				{
				Stroke oldStroke = g2d.getStroke();
				g2d.setStroke (myOutlineStroke);
				g2d.setPaint (myOutlineColor);
				g2d.draw (myRectangle);
				g2d.setStroke (oldStroke);
				}
			g2d.setPaint (oldPaint);
			g2d.setTransform (oldTransform);
			}

		public void writeExternal
			(ObjectOutput out)
			throws IOException
			{
			out.writeObject (myFillColor);
			out.writeObject (myOutlineColor);
			Strokes.writeExternal (myOutlineStroke, out);
			out.writeDouble (myRectangle.x);
			out.writeDouble (myRectangle.y);
			out.writeDouble (myRectangle.width);
			out.writeDouble (myRectangle.height);
			out.writeDouble (myAngle);
			}

		public void readExternal
			(ObjectInput in)
			throws IOException, ClassNotFoundException
			{
			myFillColor = (Color) in.readObject();
			myOutlineColor = (Color) in.readObject();
			myOutlineStroke = Strokes.readExternal (in);
			myRectangle = new Rectangle2D.Double
				(in.readDouble(),
				 in.readDouble(),
				 in.readDouble(),
				 in.readDouble());
			myAngle = in.readDouble();
			}
		}

// Exported constructors.

	/**
	 * Construct a new Dots object.
	 */
	public Dots()
		{
		}

// Exported operations.

	/**
	 * Returns an object that draws a circle. The circle uses the default fill
	 * color (black), the default outline stroke (none), the default outline
	 * color (none), and the default diameter (5).
	 */
	public static Dots circle()
		{
		return new Ellipse (Color.black, null, null, 5, 5);
		}

	/**
	 * Returns an object that draws a circle with the given diameter. The circle
	 * uses the default fill color (black), the default outline stroke (none),
	 * and the default outline color (none).
	 *
	 * @param  theDiameter
	 *     Circle's diameter.
	 */
	public static Dots circle
		(double theDiameter)
		{
		return new Ellipse (Color.black, null, null, theDiameter, theDiameter);
		}

	/**
	 * Returns an object that draws a circle with the given fill color. The
	 * circle uses the default outline stroke (none), the default outline color
	 * (none), and the default diameter (5).
	 *
	 * @param  theFillColor
	 *     Color to fill the circle's interior, or null not to fill the circle's
	 *     interior.
	 */
	public static Dots circle
		(Color theFillColor)
		{
		return new Ellipse (theFillColor, null, null, 5, 5);
		}

	/**
	 * Returns an object that draws a circle with the given fill color, the
	 * given outline stroke, and the given outline color. The circle uses the
	 * default diameter (5).
	 *
	 * @param  theFillColor
	 *     Color to fill the circle's interior, or null not to fill the circle's
	 *     interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the circle's outline, or null not to draw the circle's
	 *     outline.
	 * @param  theOutlineColor
	 *     Color to draw the circle's outline, or null not to draw the circle's
	 *     outline.
	 */
	public static Dots circle
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor)
		{
		return new Ellipse (theFillColor, theOutlineStroke, theOutlineColor,
			5, 5);
		}

	/**
	 * Returns an object that draws a circle with the given fill color, the
	 * given outline stroke, the given outline color, and the given diameter.
	 *
	 * @param  theFillColor
	 *     Color to fill the circle's interior, or null not to fill the circle's
	 *     interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the circle's outline, or null not to draw the circle's
	 *     outline.
	 * @param  theOutlineColor
	 *     Color to draw the circle's outline, or null not to draw the circle's
	 *     outline.
	 * @param  theDiameter
	 *     Circle's diameter.
	 */
	public static Dots circle
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor,
		 double theDiameter)
		{
		return new Ellipse (theFillColor, theOutlineStroke, theOutlineColor,
			theDiameter, theDiameter);
		}

	/**
	 * Returns an object that draws an ellipse with the given dimensions. The
	 * ellipse uses the default fill color (black), the default outline stroke
	 * (none), and the default outline color (none).
	 *
	 * @param  theWidth
	 *     Ellipse's width.
	 * @param  theHeight
	 *     Ellipse's height.
	 */
	public static Dots ellipse
		(double theWidth,
		 double theHeight)
		{
		return new Ellipse (Color.black, null, null, theWidth, theHeight);
		}

	/**
	 * Returns an object that draws an ellipse with the given fill color and the
	 * given dimensions. The ellipse uses the default outline stroke (none) and
	 * the default outline color (none).
	 *
	 * @param  theFillColor
	 *     Color to fill the ellipse's interior, or null not to fill the
	 *     ellipse's interior.
	 * @param  theWidth
	 *     Ellipse's width.
	 * @param  theHeight
	 *     Ellipse's height.
	 */
	public static Dots ellipse
		(Color theFillColor,
		 double theWidth,
		 double theHeight)
		{
		return new Ellipse (theFillColor, null, null, theWidth, theHeight);
		}

	/**
	 * Returns an object that draws an ellipse with the given fill color, the
	 * given outline stroke, the given outline color, and the given dimensions.
	 *
	 * @param  theFillColor
	 *     Color to fill the ellipse's interior, or null not to fill the
	 *     ellipse's interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the ellipse's outline, or null not to draw the
	 *     ellipse's outline.
	 * @param  theOutlineColor
	 *     Color to draw the ellipse's outline, or null not to draw the
	 *     ellipse's outline.
	 * @param  theWidth
	 *     Ellipse's width.
	 * @param  theHeight
	 *     Ellipse's height.
	 */
	public static Dots ellipse
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor,
		 double theWidth,
		 double theHeight)
		{
		return new Ellipse (theFillColor, theOutlineStroke, theOutlineColor,
			theWidth, theHeight);
		}

	/**
	 * Returns an object that draws a diamond. The diamond uses the default fill
	 * color (black), the default outline stroke (none), the default outline
	 * color (none), and the default side (5).
	 */
	public static Dots diamond()
		{
		return new Rectangle (Color.black, null, null, 5, 5, Math.PI/4);
		}

	/**
	 * Returns an object that draws a diamond with the given side. The diamond
	 * uses the default fill color (black), the default outline stroke (none),
	 * and the default outline color (none).
	 *
	 * @param  theSide
	 *     Diamond's side.
	 */
	public static Dots diamond
		(double theSide)
		{
		return new Rectangle
			(Color.black, null, null, theSide, theSide, Math.PI/4);
		}

	/**
	 * Returns an object that draws a diamond with the given fill color. The
	 * diamond uses the default outline stroke (none), the default outline color
	 * (none), and the default side (5).
	 *
	 * @param  theFillColor
	 *     Color to fill the diamond's interior, or null not to fill the
	 *     diamond's interior.
	 */
	public static Dots diamond
		(Color theFillColor)
		{
		return new Rectangle (theFillColor, null, null, 5, 5, Math.PI/4);
		}

	/**
	 * Returns an object that draws a diamond with the given fill color, the
	 * given outline stroke, and the given outline color. The diamond uses the
	 * default side (5).
	 *
	 * @param  theFillColor
	 *     Color to fill the diamond's interior, or null not to fill the
	 *     diamond's interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the diamond's outline, or null not to draw the
	 *     diamond's outline.
	 * @param  theOutlineColor
	 *     Color to draw the diamond's outline, or null not to draw the
	 *     diamond's outline.
	 */
	public static Dots diamond
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor)
		{
		return new Rectangle
			(theFillColor, theOutlineStroke, theOutlineColor, 5, 5, Math.PI/4);
		}

	/**
	 * Returns an object that draws a diamond with the given fill color, the
	 * given outline stroke, the given outline color, and the given side.
	 *
	 * @param  theFillColor
	 *     Color to fill the diamond's interior, or null not to fill the
	 *     diamond's interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the diamond's outline, or null not to draw the
	 *     diamond's outline.
	 * @param  theOutlineColor
	 *     Color to draw the diamond's outline, or null not to draw the
	 *     diamond's outline.
	 * @param  theSide
	 *     Diamond's side.
	 */
	public static Dots diamond
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor,
		 double theSide)
		{
		return new Rectangle (theFillColor, theOutlineStroke, theOutlineColor,
			theSide, theSide, Math.PI/4);
		}

	/**
	 * Returns an object that draws a square. The square uses the default fill
	 * color (black), the default outline stroke (none), the default outline
	 * color (none), and the default side (5).
	 */
	public static Dots square()
		{
		return new Rectangle (Color.black, null, null, 5, 5, 0);
		}

	/**
	 * Returns an object that draws a square with the given side. The square
	 * uses the default fill color (black), the default outline stroke (none),
	 * and the default outline color (none).
	 *
	 * @param  theSide
	 *     Square's side.
	 */
	public static Dots square
		(double theSide)
		{
		return new Rectangle
			(Color.black, null, null, theSide, theSide, 0);
		}

	/**
	 * Returns an object that draws a square with the given fill color. The
	 * square uses the default outline stroke (none), the default outline color
	 * (none), and the default side (5).
	 *
	 * @param  theFillColor
	 *     Color to fill the square's interior, or null not to fill the square's
	 *     interior.
	 */
	public static Dots square
		(Color theFillColor)
		{
		return new Rectangle (theFillColor, null, null, 5, 5, 0);
		}

	/**
	 * Returns an object that draws a square with the given fill color, the
	 * given outline stroke, and the given outline color. The square uses the
	 * default side (5).
	 *
	 * @param  theFillColor
	 *     Color to fill the square's interior, or null not to fill the square's
	 *     interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the square's outline, or null not to draw the square's
	 *     outline.
	 * @param  theOutlineColor
	 *     Color to draw the square's outline, or null not to draw the square's
	 *     outline.
	 */
	public static Dots square
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor)
		{
		return new Rectangle
			(theFillColor, theOutlineStroke, theOutlineColor, 5, 5, 0);
		}

	/**
	 * Returns an object that draws a square with the given fill color, the
	 * given outline stroke, the given outline color, and the given side.
	 *
	 * @param  theFillColor
	 *     Color to fill the square's interior, or null not to fill the square's
	 *     interior.
	 * @param  theOutlineStroke
	 *     Stroke to draw the square's outline, or null not to draw the square's
	 *     outline.
	 * @param  theOutlineColor
	 *     Color to draw the square's outline, or null not to draw the square's
	 *     outline.
	 * @param  theSide
	 *     Diamond's side.
	 */
	public static Dots square
		(Color theFillColor,
		 BasicStroke theOutlineStroke,
		 Color theOutlineColor,
		 double theSide)
		{
		return new Rectangle (theFillColor, theOutlineStroke, theOutlineColor,
			theSide, theSide, 0);
		}

	}
