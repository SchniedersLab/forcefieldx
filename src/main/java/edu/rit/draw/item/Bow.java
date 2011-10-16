//******************************************************************************
//
// File:    Bow.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Bow
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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

package edu.rit.draw.item;

import edu.rit.draw.Drawing;

import java.awt.BasicStroke;
import java.awt.Graphics2D;

import java.awt.geom.GeneralPath;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Bow provides a {@linkplain DrawingItem} that joins two points with a
 * curved (bowed) line. The curve is specified by a "curve factor." If the curve
 * factor is positive, the line curves to the right when going from the starting
 * point to the ending point. If the curve factor is negative, the line curves
 * to the left. The larger the magnitude of the curve factor, the farther the
 * line curves. (The curve factor is the tangent of the angle at which the line
 * leaves the starting point and approaches the ending point.) The line may have
 * an {@linkplain Arrow} at either or both ends. The line has an outline but no
 * filled interior.
 * <P>
 * The static <TT>defaultCurve()</TT>, <TT>defaultStartArrow()</TT>, and
 * <TT>defaultEndArrow()</TT> methods are provided to set the default curve
 * factor, starting arrow, and ending arrow. If the curve factor, starting
 * arrow, or ending arrow is not specified, the current default curve factor,
 * starting arrow, or ending arrow is used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Bow
	extends OutlinedItem
	implements Externalizable
	{

// Exported constants.

	/**
	 * The normal starting arrow: None.
	 */
	public static final Arrow NORMAL_START_ARROW = Arrow.NONE;

	/**
	 * The normal ending arrow: None.
	 */
	public static final Arrow NORMAL_END_ARROW = Arrow.NONE;

	/**
	 * The normal curve factor (0.5).
	 */
	public static final double NORMAL_CURVE = 0.5;

// Hidden data members.

	private static final long serialVersionUID = 7585879104799512917L;

	// Line endpoints.
	Point myStartPoint;
	Point myEndPoint;

	// Attributes.
	Arrow myStartArrow = theDefaultStartArrow;
	Arrow myEndArrow = theDefaultEndArrow;
	double myCurve = theDefaultCurve;

	// The default attributes.
	static Arrow theDefaultStartArrow = NORMAL_START_ARROW;
	static Arrow theDefaultEndArrow = NORMAL_END_ARROW;
	static double theDefaultCurve = NORMAL_CURVE;

	// For computing the bounding box.
	Point myNw;
	Point mySe;
	Size mySize;

	// For drawing the line and arrows.
	double x1, y1, x2, y2, x3, y3, xb, yb;

// Exported constructors.

	/**
	 * Construct a new bowed line. The line has no points initially. The default
	 * outline, starting arrow, ending arrow, and curve factor are used.
	 */
	public Bow()
		{
		super();
		}

	/**
	 * Construct a bowed line with the same points, outline, starting arrow,
	 * ending arrow, and curve factor as the given line.
	 *
	 * @param  theBow  Bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theBow</TT> is null.
	 */
	public Bow
		(Bow theBow)
		{
		super (theBow);
		myStartPoint = theBow.myStartPoint;
		myEndPoint = theBow.myEndPoint;
		myStartArrow = theBow.myStartArrow;
		myEndArrow = theBow.myEndArrow;
		myCurve = theBow.myCurve;
		}

// Exported operations.

	/**
	 * Returns the default starting arrow.
	 *
	 * @return  Default starting arrow.
	 */
	public static Arrow defaultStartArrow()
		{
		return theDefaultStartArrow;
		}

	/**
	 * Set the default starting arrow. Before calling this method the first
	 * time, the default starting arrow is <TT>Arrow.NONE</TT>.
	 *
	 * @param  theArrow  Default starting arrow.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArrow</TT> is null.
	 */
	public static void defaultStartArrow
		(Arrow theArrow)
		{
		if (theArrow == null) throw new NullPointerException();
		theDefaultStartArrow = theArrow;
		}

	/**
	 * Returns the default ending arrow.
	 *
	 * @return  Default ending arrow.
	 */
	public static Arrow defaultEndArrow()
		{
		return theDefaultEndArrow;
		}

	/**
	 * Set the default ending arrow. Before calling this method the first
	 * time, the default ending arrow is <TT>Arrow.NONE</TT>.
	 *
	 * @param  theArrow  Default ending arrow.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArrow</TT> is null.
	 */
	public static void defaultEndArrow
		(Arrow theArrow)
		{
		if (theArrow == null) throw new NullPointerException();
		theDefaultEndArrow = theArrow;
		}

	/**
	 * Returns the default curve factor.
	 *
	 * @return  Default curve factor.
	 */
	public static double defaultCurve()
		{
		return theDefaultCurve;
		}

	/**
	 * Set the default curve factor. Before calling this method the first time,
	 * the default curve factor is 0.5.
	 *
	 * @param  theCurve  Default curve factor.
	 */
	public static void defaultCurve
		(double theCurve)
		{
		theDefaultCurve = theCurve;
		}

	/**
	 * Returns the size of this drawing item's bounding box.
	 *
	 * @return  Size.
	 */
	public Size size()
		{
		computeBoundingBox();
		return mySize;
		}

	/**
	 * Returns the width of this drawing item's bounding box.
	 *
	 * @return  Width.
	 */
	public double width()
		{
		computeBoundingBox();
		return mySize.width;
		}

	/**
	 * Returns the height of this drawing item's bounding box.
	 *
	 * @return  Height.
	 */
	public double height()
		{
		computeBoundingBox();
		return mySize.height;
		}

	/**
	 * Returns the northwest corner point of this drawing item's bounding box.
	 *
	 * @return  Northwest corner point.
	 */
	public Point nw()
		{
		computeBoundingBox();
		return myNw;
		}

	/**
	 * Set this bowed line's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This bowed line.
	 */
	public Bow outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Returns this bowed line's starting point.
	 *
	 * @return  Starting point.
	 */
	public Point start()
		{
		return myStartPoint;
		}

	/**
	 * Set this bowed line's starting point.
	 *
	 * @param  x  Starting point X coordinate.
	 * @param  y  Starting point Y coordinate.
	 *
	 * @return  This bowed line.
	 */
	public Bow start
		(double x,
		 double y)
		{
		myStartPoint = new Point (x, y);
		return this;
		}

	/**
	 * Set this bowed line's starting point.
	 *
	 * @param  thePoint  Starting point.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Bow start
		(Point thePoint)
		{
		if (thePoint == null) throw new NullPointerException();
		myStartPoint = thePoint;
		return this;
		}

	/**
	 * Returns this bowed line's ending point.
	 *
	 * @return  Ending point.
	 */
	public Point end()
		{
		return myEndPoint;
		}

	/**
	 * Set this bowed line's ending point.
	 *
	 * @param  x  Ending point X coordinate.
	 * @param  y  Ending point Y coordinate.
	 *
	 * @return  This bowed line.
	 */
	public Bow end
		(double x,
		 double y)
		{
		myEndPoint = new Point (x, y);
		return this;
		}

	/**
	 * Set this bowed line's ending point.
	 *
	 * @param  thePoint  Ending point.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Bow end
		(Point thePoint)
		{
		if (thePoint == null) throw new NullPointerException();
		myEndPoint = thePoint;
		return this;
		}

	/**
	 * Returns this bowed line's starting arrow.
	 *
	 * @return  Starting arrow.
	 */
	public Arrow startArrow()
		{
		return myStartArrow;
		}

	/**
	 * Set this bowed line's starting arrow.
	 *
	 * @param  theArrow  Starting arrow.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArrow</TT> is null.
	 */
	public Bow startArrow
		(Arrow theArrow)
		{
		if (theArrow == null) throw new NullPointerException();
		myStartArrow = theArrow;
		return this;
		}

	/**
	 * Returns this bowed line's ending arrow.
	 *
	 * @return  Ending arrow.
	 */
	public Arrow endArrow()
		{
		return myEndArrow;
		}

	/**
	 * Set this bowed line's ending arrow.
	 *
	 * @param  theArrow  Ending arrow.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArrow</TT> is null.
	 */
	public Bow endArrow
		(Arrow theArrow)
		{
		if (theArrow == null) throw new NullPointerException();
		myEndArrow = theArrow;
		return this;
		}

	/**
	 * Returns this bowed line's curve factor.
	 *
	 * @return  Curve factor.
	 */
	public double curve()
		{
		return myCurve;
		}

	/**
	 * Set this bowed line's curve factor.
	 *
	 * @param  theCurve  Curve factor.
	 *
	 * @return  This bowed line.
	 */
	public Bow curve
		(double theCurve)
		{
		myCurve = theCurve;
		return this;
		}

	/**
	 * Flip this bowed line from curving right to curving left or vice versa.
	 *
	 * @return  This bowed line.
	 */
	public Bow flip()
		{
		myCurve = -myCurve;
		return this;
		}

	/**
	 * Add this bowed line to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Bow add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this bowed line to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Bow add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this bowed line to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Bow addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this bowed line to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This bowed line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Bow addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this bowed line to the given object output stream.
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
		super.writeExternal (out);
		out.writeObject (myStartPoint);
		out.writeObject (myEndPoint);
		out.writeObject (myStartArrow);
		out.writeObject (myEndArrow);
		out.writeDouble (myCurve);
		}

	/**
	 * Read this bowed line from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this line cannot be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		myStartPoint = (Point) in.readObject();
		myEndPoint = (Point) in.readObject();
		myStartArrow = (Arrow) in.readObject();
		myEndArrow = (Arrow) in.readObject();
		myCurve = in.readDouble();
		myNw = null;
		mySe = null;
		mySize = null;
		}

	/**
	 * Draw this drawing item in the given graphics context. This method is
	 * allowed to change the graphics context's paint, stroke, and transform,
	 * and it doesn't have to change them back.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void draw
		(Graphics2D g2d)
		{
		super.draw (g2d);

		// Draw line.
		GeneralPath path = computePath();
		myOutline.setGraphicsContext (g2d);
		g2d.draw (path);

		// Return if there are no arrows.
		if (myStartArrow == null && myEndArrow == null) return;

		// Prepare to draw arrows. Set stroke to a solid stroke the same width
		// as the current outline.
		float width = myOutline.getStrokeWidth();
		g2d.setStroke (new BasicStroke (width));
		double dx;
		double dy;
		double phi;

		if (myStartArrow != null)
			{
			dx = x1 - x2;
			dy = y1 - y2;
			phi =
				dx == 0.0 && dy == 0.0 ?
					0.0 :
					Math.atan2 (dy, dx);
			myStartArrow.draw (g2d, width, x1, y1, phi);
			}

		if (myEndArrow != null)
			{
			dx = x3 - x2;
			dy = y3 - y2;
			phi =
				dx == 0.0 && dy == 0.0 ?
					0.0 :
					Math.atan2 (dy, dx);
			myEndArrow.draw (g2d, width, x3, y3, phi);
			}
		}

// Hidden operations.

	/**
	 * Compute this line's three control points and store their (x,y)
	 * coordinates in x1, y1, x2, y2, x3, and y3.
	 */
	private void computeControlPoints()
		{
		// Get starting and ending coordinates.
		x1 = myStartPoint.x();
		y1 = myStartPoint.y();
		x3 = myEndPoint.x();
		y3 = myEndPoint.y();

		// Get midpoint coordinates.
		double xm = 0.5 * (x1 + x3);
		double ym = 0.5 * (y1 + y3);

		// Get unit vector from starting point to ending point.
		double dx = x3 - x1;
		double dy = y3 - y1;
		double d = Math.sqrt (dx*dx + dy*dy);
		double ux = dx / d;
		double uy = dy / d;

		// Get unit vector perpendicular to the above vector.
		double vx = -uy;
		double vy =  ux;

		// Get second control point.
		x2 = xm + 0.5 * d * myCurve * vx;
		y2 = ym + 0.5 * d * myCurve * vy;

		// Get point through which curve passes, for bounding box.
		xb = xm + 0.25 * d * myCurve * vx;
		yb = ym + 0.25 * d * myCurve * vy;
		}

	/**
	 * Compute this line's bounding box from its control points. The results are
	 * stored in myNw, mySe, and mySize.
	 */
	private void computeBoundingBox()
		{
		if (myNw == null)
			{
			computeControlPoints();
			myNw = new Point
				(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
			mySe = new Point
				(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
			Point p = new Point (x1, y1);
			myNw = myNw.min (p);
			mySe = mySe.max (p);
			p = new Point (xb, yb);
			myNw = myNw.min (p);
			mySe = mySe.max (p);
			p = new Point (x3, y3);
			myNw = myNw.min (p);
			mySe = mySe.max (p);
			mySize = mySe.difference (myNw);
			}
		}

	/**
	 * Compute this line's path from its control points.
	 */
	private GeneralPath computePath()
		{
		computeControlPoints();
		GeneralPath path = new GeneralPath();
		path.moveTo ((float) x1, (float) y1);
		path.quadTo ((float) x2, (float) y2, (float) x3, (float) y3);
		return path;
		}

	}
