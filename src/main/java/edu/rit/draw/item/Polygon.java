//******************************************************************************
//
// File:    Polygon.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Polygon
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

import java.awt.geom.GeneralPath;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.util.Iterator;
import java.util.LinkedList;

/**
 * Class Polygon provides a polygon {@linkplain DrawingItem} whose boundary
 * consists of one or more straight line segments. The polygon has an outline
 * and a filled interior.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Polygon
	extends ShapeItem
	{

// Hidden data members.

	private static final long serialVersionUID = -8680229678155189048L;

	// Line segment endpoints.
	LinkedList<Point> myPoints = new LinkedList<Point>();

	// For adding line segments.
	Point myLastPoint;

	// For computing the bounding box.
	Point myNw;
	Point mySe;
	Size mySize;

// Exported constructors.

	/**
	 * Construct a new polygon. The polygon has no points initially. The default
	 * outline and fill paint are used.
	 */
	public Polygon()
		{
		super();
		}

	/**
	 * Construct a new polygon with the same points, outline, and fill paint as
	 * the given polygon.
	 *
	 * @param  thePolygon  Polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePolygon</TT> is null.
	 */
	public Polygon
		(Polygon thePolygon)
		{
		super (thePolygon);
		myPoints.clear();
		myPoints.addAll (thePolygon.myPoints);
		if (! myPoints.isEmpty())
			{
			myLastPoint = myPoints.get (myPoints.size()-1);
			}
		}

// Exported operations.

	/**
	 * Returns the size of this polygon's bounding box.
	 *
	 * @return  Size.
	 */
	public Size size()
		{
		computeBoundingBox();
		return mySize;
		}

	/**
	 * Returns the width of this polygon's bounding box.
	 *
	 * @return  Width.
	 */
	public double width()
		{
		computeBoundingBox();
		return mySize.width;
		}

	/**
	 * Returns the height of this polygon's bounding box.
	 *
	 * @return  Height.
	 */
	public double height()
		{
		computeBoundingBox();
		return mySize.height;
		}

	/**
	 * Returns the northwest corner point of this polygon's bounding box.
	 *
	 * @return  Northwest corner point.
	 */
	public Point nw()
		{
		computeBoundingBox();
		return myNw;
		}

	/**
	 * Set this polygon's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This polygon.
	 */
	public Polygon outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Set this polygon's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This polygon.
	 */
	public Polygon fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Clear this polygon. All of this polygon's line segments are removed.
	 *
	 * @return  This polygon.
	 */
	public Polygon clear()
		{
		myPoints.clear();
		myLastPoint = null;
		myNw = null;
		return this;
		}

	/**
	 * Add a line segment to this polygon. The new endpoint is (<TT>x</TT>,
	 * <TT>y</TT>). If this polygon has no line segments, the <TT>to()</TT>
	 * method specifies this polygon's starting endpoint.
	 *
	 * @param  x  Endpoint's X coordinate.
	 * @param  y  Endpoint's Y coordinate.
	 *
	 * @return  This polygon.
	 */
	public Polygon to
		(double x,
		 double y)
		{
		return addEndpoint (new Point (x, y));
		}

	/**
	 * Add a line segment to this polygon. The new endpoint is
	 * <TT>thePoint</TT>. If this polygon has no line segments, the
	 * <TT>to()</TT> method specifies this polygon's starting endpoint.
	 *
	 * @param  thePoint  Endpoint.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Polygon to
		(Point thePoint)
		{
		if (thePoint == null) throw new NullPointerException();
		return addEndpoint (thePoint);
		}

	/**
	 * Add a horizontal line segment to this polygon. The new endpoint is
	 * (<TT>x</TT>, previous endpoint's Y coordinate).
	 *
	 * @param  x  Endpoint's X coordinate.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Polygon hto
		(double x)
		{
		return addEndpoint (new Point (x, myLastPoint.y));
		}

	/**
	 * Add a horizontal line segment to this polygon. The new endpoint is
	 * (<TT>thePoint</TT>'s X coordinate, previous endpoint's Y coordinate).
	 *
	 * @param  thePoint  Endpoint.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Polygon hto
		(Point thePoint)
		{
		return addEndpoint (new Point (thePoint.x, myLastPoint.y));
		}

	/**
	 * Add a vertical line segment to this polygon. The new endpoint is
	 * (previous endpoint's X coordinate, <TT>y</TT>).
	 *
	 * @param  y  Endpoint's Y coordinate.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Polygon vto
		(double y)
		{
		return addEndpoint (new Point (myLastPoint.x, y));
		}

	/**
	 * Add a vertical line segment to this polygon. The new endpoint is
	 * (previous endpoint's X coordinate, <TT>thePoint</TT>'s Y coordinate).
	 *
	 * @param  thePoint  Endpoint.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Polygon vto
		(Point thePoint)
		{
		return addEndpoint (new Point (myLastPoint.x, thePoint.y));
		}

	/**
	 * Add a line segment to this polygon. The new endpoint is (previous
	 * endpoint's X coordinate + <TT>dx</TT>, previous endpoint's Y coordinate +
	 * <TT>dy</TT>).
	 *
	 * @param  dx  X distance.
	 * @param  dy  Y distance.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Polygon by
		(double dx,
		 double dy)
		{
		return addEndpoint (new Point (myLastPoint.x+dx, myLastPoint.y+dy));
		}

	/**
	 * Add a line segment to this polygon. The new endpoint is (previous
	 * endpoint's X coordinate + <TT>theSize.width()</TT>, previous endpoint's Y
	 * coordinate + <TT>theSize.height()</TT>).
	 *
	 * @param  theSize  Distance.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Polygon by
		(Size theSize)
		{
		return addEndpoint
			(new Point
				(myLastPoint.x+theSize.width, myLastPoint.y+theSize.height));
		}

	/**
	 * Add a horizontal line segment to this polygon. The new endpoint is
	 * (previous endpoint's X coordinate + <TT>dx</TT>, previous endpoint's Y
	 * coordinate).
	 *
	 * @param  dx  X distance.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Polygon hby
		(double dx)
		{
		return addEndpoint (new Point (myLastPoint.x+dx, myLastPoint.y));
		}

	/**
	 * Add a horizontal line segment to this polygon. The new endpoint is
	 * (previous endpoint's X coordinate + <TT>theSize.width()</TT>, previous
	 * endpoint's Y coordinate).
	 *
	 * @param  theSize  Distance.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Polygon hby
		(Size theSize)
		{
		return addEndpoint
			(new Point (myLastPoint.x+theSize.width, myLastPoint.y));
		}

	/**
	 * Add a vertical line segment to this polygon. The new endpoint is
	 * (previous endpoint's X coordinate, previous endpoint's Y coordinate +
	 * <TT>dy</TT>).
	 *
	 * @param  dy  Y distance.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Polygon vby
		(double dy)
		{
		return addEndpoint (new Point (myLastPoint.x, myLastPoint.y+dy));
		}

	/**
	 * Add a vertical line segment to this polygon. The new endpoint is
	 * (previous endpoint's X coordinate, previous endpoint's Y coordinate +
	 * <TT>theSize.height()</TT>).
	 *
	 * @param  theSize  Distance.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Polygon vby
		(Size theSize)
		{
		return addEndpoint
			(new Point (myLastPoint.x, myLastPoint.y+theSize.height));
		}

	/**
	 * Add this polygon to the end of the default drawing's sequence of drawing
	 * items.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Polygon add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this polygon to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Polygon add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this polygon to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Polygon addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this polygon to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This polygon.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Polygon addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this polygon to the given object output stream.
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
		out.writeInt (myPoints.size());
		for (Point point : myPoints)
			{
			out.writeObject (point);
			}
		}

	/**
	 * Read this polygon from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this polygon cannot be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		int n = in.readInt();
		myPoints.clear();
		for (int i = 0; i < n; ++ i)
			{
			myPoints.add ((Point) in.readObject());
			}
		if (! myPoints.isEmpty())
			{
			myLastPoint = myPoints.get (myPoints.size()-1);
			}
		myNw = null;
		mySe = null;
		mySize = null;
		}

// Hidden operations.

	/**
	 * Add the given endpoint to this polygon.
	 */
	private Polygon addEndpoint
		(Point p)
		{
		myLastPoint = p;
		myPoints.add (p);
		myNw = null;
		return this;
		}

	/**
	 * Compute this polygon's bounding box from myPoints. The results are stored
	 * in myNw, mySe, and mySize.
	 */
	private void computeBoundingBox()
		{
		if (myNw == null)
			{
			myNw = new Point
				(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
			mySe = new Point
				(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
			for (Point p : myPoints)
				{
				myNw = myNw.min (p);
				mySe = mySe.max (p);
				}
			mySize = mySe.difference (myNw);
			}
		}

	/**
	 * Determine the 2-D graphics shape that this shape object will draw.
	 *
	 * @return  Shape.
	 */
	Shape getShape()
		{
		GeneralPath path = new GeneralPath();
		int n = myPoints.size();
		Iterator<Point> iter = myPoints.iterator();
		Point point = null;

		if (n > 0)
			{
			point = iter.next();
			path.moveTo ((float) point.x, (float) point.y);
			while (iter.hasNext())
				{
				point = iter.next();
				path.lineTo ((float) point.x, (float) point.y);
				}
			path.closePath();
			}

		return path;
		}

	}
