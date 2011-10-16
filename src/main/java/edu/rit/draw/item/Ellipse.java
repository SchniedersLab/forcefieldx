//******************************************************************************
//
// File:    Ellipse.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Ellipse
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

import java.awt.Shape;

import java.awt.geom.Ellipse2D;

/**
 * Class Ellipse provides an ellipse {@linkplain DrawingItem}. A circle is a
 * special case of an ellipse with equal width and height.
 * <P>
 * The static <TT>defaultSize()</TT>, <TT>defaultWidth()</TT>, and
 * <TT>defaultHeight()</TT> methods are provided to set the default size, width,
 * and height. If an ellipse's size, width, or height is not specified, the
 * current default size, width, or height is used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Ellipse
	extends RectangularItem
	{

// Exported constants.

	/**
	 * The normal width for ellipses (72). <I>Note:</I> 72 points = 1 inch.
	 */
	public static final double NORMAL_WIDTH = 72.0;

	/**
	 * The normal height for ellipses (72). <I>Note:</I> 72 points = 1 inch.
	 */
	public static final double NORMAL_HEIGHT = 72.0;

// Hidden data members.

	private static final long serialVersionUID = 8085149752960383250L;

	private static double theDefaultWidth = NORMAL_WIDTH;
	private static double theDefaultHeight = NORMAL_HEIGHT;

// Exported constructors.

	/**
	 * Construct a new ellipse. The ellipse's northwest corner is located at
	 * (0,0). The ellipse's size is the default size.
	 */
	public Ellipse()
		{
		super();
		this.width = theDefaultWidth;
		this.height = theDefaultHeight;
		}

	/**
	 * Construct a new ellipse with the same outline, fill paint, location,
	 * and size as the given ellipse.
	 *
	 * @param  theItem  Ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Ellipse
		(Ellipse theItem)
		{
		super (theItem);
		}

// Exported operations.

	/**
	 * Returns the default size for ellipses.
	 *
	 * @return  Default size.
	 */
	public static Size defaultSize()
		{
		return new Size (theDefaultWidth, theDefaultHeight);
		}

	/**
	 * Set the default size for ellipses.
	 *
	 * @param  theSize  Default size.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	public static void defaultSize
		(Size theSize)
		{
		if (theSize.width < 0.0 || theSize.height < 0.0)
			{
			throw new IllegalArgumentException();
			}
		theDefaultWidth = theSize.width;
		theDefaultHeight = theSize.height;
		}

	/**
	 * Returns the default width for ellipses.
	 *
	 * @return  Default width.
	 */
	public static double defaultWidth()
		{
		return theDefaultWidth;
		}

	/**
	 * Set the default width for ellipses.
	 *
	 * @param  theWidth  Default width.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	public static void defaultWidth
		(double theWidth)
		{
		if (theWidth < 0.0) throw new IllegalArgumentException();
		theDefaultWidth = theWidth;
		}

	/**
	 * Returns the default height for ellipses.
	 *
	 * @return  Default height.
	 */
	public static double defaultHeight()
		{
		return theDefaultHeight;
		}

	/**
	 * Set the default height for ellipses.
	 *
	 * @param  theHeight  Default height.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	public static void defaultHeight
		(double theHeight)
		{
		if (theHeight < 0.0) throw new IllegalArgumentException();
		theDefaultHeight = theHeight;
		}

	/**
	 * Set this ellipse's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Set this ellipse's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Set the size of this ellipse's bounding box.
	 *
	 * @param  theSize  Size.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	public Ellipse size
		(Size theSize)
		{
		doSize (theSize);
		return this;
		}

	/**
	 * Set the width of this ellipse's bounding box.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	public Ellipse width
		(double theWidth)
		{
		doWidth (theWidth);
		return this;
		}

	/**
	 * Set the height of this ellipse's bounding box.
	 *
	 * @param  theHeight  Height.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	public Ellipse height
		(double theHeight)
		{
		doHeight (theHeight);
		return this;
		}

	/**
	 * Set the width and height of this ellipse's bounding box to make this
	 * ellipse a circle.
	 *
	 * @param  theDiameter  Diameter.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theDiameter</TT> is less than 0.
	 */
	public Ellipse diameter
		(double theDiameter)
		{
		doWidth (theDiameter);
		doHeight (theDiameter);
		return this;
		}

	/**
	 * Set the northwest corner point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse nw
		(double x,
		 double y)
		{
		doNw (x, y);
		return this;
		}

	/**
	 * Set the northwest corner point of this ellipse's bounding box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the north middle point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse n
		(double x,
		 double y)
		{
		doN (x, y);
		return this;
		}

	/**
	 * Set the north middle point of this ellipse's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the northeast corner point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse ne
		(double x,
		 double y)
		{
		doNe (x, y);
		return this;
		}

	/**
	 * Set the northeast corner point of this ellipse's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the west middle point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse w
		(double x,
		 double y)
		{
		doW (x, y);
		return this;
		}

	/**
	 * Set the west middle point of this ellipse's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the center point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse c
		(double x,
		 double y)
		{
		doC (x, y);
		return this;
		}

	/**
	 * Set the center point of this ellipse's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the east middle point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse e
		(double x,
		 double y)
		{
		doE (x, y);
		return this;
		}

	/**
	 * Set the east middle point of this ellipse's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southwest corner point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse sw
		(double x,
		 double y)
		{
		doSw (x, y);
		return this;
		}

	/**
	 * Set the southwest corner point of this ellipse's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the south middle point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse s
		(double x,
		 double y)
		{
		doS (x, y);
		return this;
		}

	/**
	 * Set the south middle point of this ellipse's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southeast corner point of this ellipse's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This ellipse.
	 */
	public Ellipse se
		(double x,
		 double y)
		{
		doSe (x, y);
		return this;
		}

	/**
	 * Set the southeast corner point of this ellipse's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Ellipse se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Add this ellipse to the end of the default drawing's sequence of drawing
	 * items.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Ellipse add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this ellipse to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Ellipse add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this ellipse to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Ellipse addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this ellipse to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This ellipse.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Ellipse addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Returns the point on the circumference of this ellipse at the given angle
	 * relative to this ellipse's center. The point on the ellipse's
	 * circumference is determined by <TT>angle</TT>, measured in radians. An
	 * angle of 0 gives the east point; an angle of pi/2 gives the south point;
	 * an angle of pi gives the west point; an angle of 3*pi/2 gives the north
	 * point.
	 *
	 * @param  angle  Angle (radians).
	 *
	 * @return  Point on this ellipse's circumference.
	 */
	public Point circumferencePoint
		(double angle)
		{
		Point c = c();
		double xc = c.x;
		double yc = c.y;
		double xr = width() / 2;
		double yr = height() / 2;
		double x1 = xc + xr * Math.cos (angle);
		double y1 = yc + yr * Math.sin (angle);
		return new Point (x1, y1);
		}

	/**
	 * Returns the point on the circumference of this ellipse at the given angle
	 * relative to this ellipse's center. The point on the ellipse's
	 * circumference is determined by <TT>angle</TT>, measured in degrees. An
	 * angle of 0 gives the east point; an angle of 90 gives the south point;
	 * an angle of 180 gives the west point; an angle of 270 gives the north
	 * point.
	 *
	 * @param  angle  Angle (degrees).
	 *
	 * @return  Point on this ellipse's circumference.
	 */
	public Point circumferencePointDegrees
		(double angle)
		{
		return circumferencePoint (angle / 180.0 * Math.PI);
		}

	/**
	 * Returns the point on the circumference of this ellipse at the given angle
	 * relative to this ellipse's center, extended outward by the given length.
	 * First, a point on this ellipse's circumference is determined by
	 * <TT>angle</TT>, measured in radians. An angle of 0 gives the east point;
	 * an angle of pi/2 gives the south point; an angle of pi gives the west
	 * point; an angle of 3*pi/2 gives the north point. A line is projected from
	 * the center of this ellipse, through the point on the circumference, for a
	 * further distance of <TT>length</TT>; the endpoint of this line is
	 * returned.
	 *
	 * @param  angle   Angle (radians).
	 * @param  length  Length to extend from circumference.
	 *
	 * @return  Point extended outward from this ellipse's circumference.
	 */
	public Point circumferencePoint
		(double angle,
		 double length)
		{
		Point c = c();
		double xc = c.x;
		double yc = c.y;
		double xr = width() / 2;
		double yr = height() / 2;
		double dx = xr * Math.cos (angle);
		double dy = yr * Math.sin (angle);
		double d = Math.sqrt (dx*dx + dy*dy);
		double x1 = xc + dx + dx / d * length;
		double y1 = yc + dy + dy / d * length;
		return new Point (x1, y1);
		}

	/**
	 * Returns the point on the circumference of this ellipse at the given angle
	 * relative to this ellipse's center, extended outward by the given length.
	 * First, a point on this ellipse's circumference is determined by
	 * <TT>angle</TT>, measured in degrees. An angle of 0 gives the east point;
	 * an angle of 90 gives the south point; an angle of 180 gives the west
	 * point; an angle of 270 gives the north point. A line is projected from
	 * the center of this ellipse, through the point on the circumference, for a
	 * further distance of <TT>length</TT>; the endpoint of this line is
	 * returned.
	 *
	 * @param  angle   Angle (degrees).
	 * @param  length  Length to extend from circumference.
	 *
	 * @return  Point extended outward from this ellipse's circumference.
	 */
	public Point circumferencePointDegrees
		(double angle,
		 double length)
		{
		return circumferencePoint (angle / 180.0 * Math.PI, length);
		}

	/**
	 * Returns the point of intersection between (1) this ellipse's perimeter
	 * and (2) the line passing through this ellipse's center and the given
	 * point.
	 *
	 * @param  point  One point on the line.
	 *
	 * @return  Point of intersection between perimeter and line.
	 */
	public Point perimeterIntersection
		(Point point)
		{
		return perimeterIntersection (point.x, point.y);
		}

	/**
	 * Returns the point of intersection between (1) this ellipse's perimeter
	 * and (2) the line passing through this ellipse's center and the given
	 * point.
	 *
	 * @param  x  X coordinate of one point on the line.
	 * @param  y  Y coordinate of one point on the line.
	 *
	 * @return  Point of intersection between perimeter and line.
	 */
	public Point perimeterIntersection
		(double x,
		 double y)
		{
		Point center = c();
		double cx = center.x;
		double cy = center.y;
		double w = width()/2;
		double h = height()/2;
		double dx = x - cx;
		double dy = y - cy;
		double d = Math.sqrt (dx*dx + dy*dy);
		return new Point (cx + w*dx/d, cy + h*dy/d);
		}

// Hidden operations.

	/**
	 * Determine the 2-D graphics shape that this shape object will draw.
	 *
	 * @return  Shape.
	 */
	Shape getShape()
		{
		Point nw = nw();
		return new Ellipse2D.Double (nw.x, nw.y, width, height);
		}

	}
