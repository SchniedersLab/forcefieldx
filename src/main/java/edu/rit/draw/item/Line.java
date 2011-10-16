//******************************************************************************
//
// File:    Line.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Line
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

import java.util.Iterator;
import java.util.LinkedList;

/**
 * Class Line provides a {@linkplain DrawingItem} that consists of one or more
 * straight line segments. The line may have an {@linkplain Arrow} at either or
 * both ends. The line may have sharp corners or round corners. The line has an
 * outline but no filled interior.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0 WIDTH=100%>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * If the line has round corners, the round corner distance gives the distance
 * from where the (sharp) corner would normally be to the starting and ending
 * points of the round corner arc, as shown in the adjoining diagram. (The round
 * corner distance is <I>d</I>.) A round corner distance of 0 results in a sharp
 * corner. In a particular line, the same round corner distance is used for
 * every corner; that is, either all the corners are sharp or all the corners
 * are rounded the same amount. Different lines can have different round corner
 * distances.
 * </TD>
 * <TD WIDTH=10> </TD>
 * <TD ALIGN="left" VALIGN="top" WIDTH=99><IMG SRC="doc-files/LineFig01.png"></TD>
 * </TR>
 * </TABLE>
 * <P>
 * The static <TT>defaultStartArrow()</TT>, <TT>defaultEndArrow()</TT>, and
 * <TT>defaultRound()</TT> methods are provided to set the default starting
 * arrow, ending arrow, and round corner distance. If the starting arrow, ending
 * arrow, or round corner distance is not specified, the current default
 * starting arrow, ending arrow, or round corner distance is used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Line
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
	 * The normal round corner distance (0), signifying sharp corners.
	 */
	public static final double NORMAL_ROUND = 0.0;

// Hidden data members.

	private static final long serialVersionUID = -4094248146984537516L;

	// Line segment endpoints.
	LinkedList<Point> myPoints = new LinkedList<Point>();

	// Attributes.
	Arrow myStartArrow = theDefaultStartArrow;
	Arrow myEndArrow = theDefaultEndArrow;
	double myRound = theDefaultRound;

	// The default attributes.
	static Arrow theDefaultStartArrow = NORMAL_START_ARROW;
	static Arrow theDefaultEndArrow = NORMAL_END_ARROW;
	static double theDefaultRound = NORMAL_ROUND;

	// For adding line segments.
	Point myLastPoint;

	// For computing the bounding box.
	Point myNw;
	Point mySe;
	Size mySize;

// Exported constructors.

	/**
	 * Construct a new line. The line has no points initially. The default
	 * outline, starting arrow, ending arrow, and round corner distance are
	 * used.
	 */
	public Line()
		{
		super();
		}

	/**
	 * Construct a new line with the same points, outline, starting arrow,
	 * ending arrow, and round corner distance as the given line.
	 *
	 * @param  theLine  Line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theLine</TT> is null.
	 */
	public Line
		(Line theLine)
		{
		super (theLine);
		myPoints.clear();
		myPoints.addAll (theLine.myPoints);
		myStartArrow = theLine.myStartArrow;
		myEndArrow = theLine.myEndArrow;
		myRound = theLine.myRound;
		if (! myPoints.isEmpty())
			{
			myLastPoint = myPoints.get (myPoints.size()-1);
			}
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
	 * Returns the default round corner distance. A value of 0 signifies sharp
	 * corners.
	 *
	 * @return  Default round corner distance.
	 */
	public static double defaultRound()
		{
		return theDefaultRound;
		}

	/**
	 * Set the default round corner distance. A value of 0 signifies sharp
	 * corners. Before calling this method the first time, the default round
	 * corner distance is 0.
	 *
	 * @param  theRound  Default round corner distance.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theRound</TT> is less than 0.
	 */
	public static void defaultRound
		(double theRound)
		{
		if (theRound < 0.0) throw new IllegalArgumentException();
		theDefaultRound = theRound;
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
	 * Set this line's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This line.
	 */
	public Line outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Clear this line. All of this point's line segments are removed.
	 *
	 * @return  This line.
	 */
	public Line clear()
		{
		myPoints.clear();
		myLastPoint = null;
		myNw = null;
		return this;
		}

	/**
	 * Add a line segment to this line. The new endpoint is (<TT>x</TT>,
	 * <TT>y</TT>). If this line has no line segments, the <TT>to()</TT> method
	 * specifies this line's starting endpoint.
	 *
	 * @param  x  Endpoint's X coordinate.
	 * @param  y  Endpoint's Y coordinate.
	 *
	 * @return  This line.
	 */
	public Line to
		(double x,
		 double y)
		{
		return addEndpoint (new Point (x, y));
		}

	/**
	 * Add a line segment to this line. The new endpoint is <TT>thePoint</TT>.
	 * If this line has no line segments, the <TT>to()</TT> method specifies
	 * this line's starting endpoint.
	 *
	 * @param  thePoint  Endpoint.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Line to
		(Point thePoint)
		{
		if (thePoint == null) throw new NullPointerException();
		return addEndpoint (thePoint);
		}

	/**
	 * Add a horizontal line segment to this line. The new endpoint is
	 * (<TT>x</TT>, previous endpoint's Y coordinate).
	 *
	 * @param  x  Endpoint's X coordinate.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Line hto
		(double x)
		{
		return addEndpoint (new Point (x, myLastPoint.y));
		}

	/**
	 * Add a horizontal line segment to this line. The new endpoint is
	 * (<TT>thePoint</TT>'s X coordinate, previous endpoint's Y coordinate).
	 *
	 * @param  thePoint  Endpoint.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Line hto
		(Point thePoint)
		{
		return addEndpoint (new Point (thePoint.x, myLastPoint.y));
		}

	/**
	 * Add a vertical line segment to this line. The new endpoint is (previous
	 * endpoint's X coordinate, <TT>y</TT>).
	 *
	 * @param  y  Endpoint's Y coordinate.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Line vto
		(double y)
		{
		return addEndpoint (new Point (myLastPoint.x, y));
		}

	/**
	 * Add a vertical line segment to this line. The new endpoint is
	 * (previous endpoint's X coordinate, <TT>thePoint</TT>'s Y coordinate).
	 *
	 * @param  thePoint  Endpoint.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Line vto
		(Point thePoint)
		{
		return addEndpoint (new Point (myLastPoint.x, thePoint.y));
		}

	/**
	 * Add a line segment to this line. The new endpoint is (previous endpoint's
	 * X coordinate + <TT>dx</TT>, previous endpoint's Y coordinate +
	 * <TT>dy</TT>).
	 *
	 * @param  dx  X distance.
	 * @param  dy  Y distance.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Line by
		(double dx,
		 double dy)
		{
		return addEndpoint (new Point (myLastPoint.x+dx, myLastPoint.y+dy));
		}

	/**
	 * Add a line segment to this line. The new endpoint is (previous endpoint's
	 * X coordinate + <TT>theSize.width()</TT>, previous endpoint's Y coordinate
	 * + <TT>theSize.height()</TT>).
	 *
	 * @param  theSize  Distance.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Line by
		(Size theSize)
		{
		return addEndpoint
			(new Point
				(myLastPoint.x+theSize.width, myLastPoint.y+theSize.height));
		}

	/**
	 * Add a horizontal line segment to this line. The new endpoint is (previous
	 * endpoint's X coordinate + <TT>dx</TT>, previous endpoint's Y coordinate).
	 *
	 * @param  dx  X distance.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Line hby
		(double dx)
		{
		return addEndpoint (new Point (myLastPoint.x+dx, myLastPoint.y));
		}

	/**
	 * Add a horizontal line segment to this line. The new endpoint is (previous
	 * endpoint's X coordinate + <TT>theSize.width()</TT>, previous endpoint's Y
	 * coordinate).
	 *
	 * @param  theSize  Distance.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Line hby
		(Size theSize)
		{
		return addEndpoint
			(new Point (myLastPoint.x+theSize.width, myLastPoint.y));
		}

	/**
	 * Add a vertical line segment to this line. The new endpoint is (previous
	 * endpoint's X coordinate, previous endpoint's Y coordinate + <TT>dy</TT>).
	 *
	 * @param  dy  Y distance.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no previous endpoint.
	 */
	public Line vby
		(double dy)
		{
		return addEndpoint (new Point (myLastPoint.x, myLastPoint.y+dy));
		}

	/**
	 * Add a vertical line segment to this line. The new endpoint is (previous
	 * endpoint's X coordinate, previous endpoint's Y coordinate +
	 * <TT>theSize.height()</TT>).
	 *
	 * @param  theSize  Distance.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSize</TT> is null. Thrown if
	 *     there is no previous endpoint.
	 */
	public Line vby
		(Size theSize)
		{
		return addEndpoint
			(new Point (myLastPoint.x, myLastPoint.y+theSize.height));
		}

	/**
	 * Returns this line's starting arrow.
	 *
	 * @return  Starting arrow.
	 */
	public Arrow startArrow()
		{
		return myStartArrow;
		}

	/**
	 * Set this line's starting arrow.
	 *
	 * @param  theArrow  Starting arrow.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArrow</TT> is null.
	 */
	public Line startArrow
		(Arrow theArrow)
		{
		if (theArrow == null) throw new NullPointerException();
		myStartArrow = theArrow;
		return this;
		}

	/**
	 * Returns this line's ending arrow.
	 *
	 * @return  Ending arrow.
	 */
	public Arrow endArrow()
		{
		return myEndArrow;
		}

	/**
	 * Set this line's ending arrow.
	 *
	 * @param  theArrow  Ending arrow.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArrow</TT> is null.
	 */
	public Line endArrow
		(Arrow theArrow)
		{
		if (theArrow == null) throw new NullPointerException();
		myEndArrow = theArrow;
		return this;
		}

	/**
	 * Returns this line's round corner distance. A value of 0 signifies sharp
	 * corners.
	 *
	 * @return  Round corner distance.
	 */
	public double round()
		{
		return myRound;
		}

	/**
	 * Set this line's round corner distance. A value of 0 signifies sharp
	 * corners.
	 *
	 * @param  theRound  Round corner distance.
	 *
	 * @return  This line.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theRound</TT> is less than 0.
	 */
	public Line round
		(double theRound)
		{
		if (theRound < 0.0) throw new IllegalArgumentException();
		myRound = theRound;
		return this;
		}

	/**
	 * Add this line to the end of the default drawing's sequence of drawing
	 * items.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Line add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this line to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Line add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this line to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Line addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this line to the beginning of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This line.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Line addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this line to the given object output stream.
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
		out.writeObject (myStartArrow);
		out.writeObject (myEndArrow);
		out.writeDouble (myRound);
		}

	/**
	 * Read this line from the given object input stream.
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
		int n = in.readInt();
		myPoints.clear();
		for (int i = 0; i < n; ++ i)
			{
			myPoints.add ((Point) in.readObject());
			}
		myStartArrow = (Arrow) in.readObject();
		myEndArrow = (Arrow) in.readObject();
		myRound = in.readDouble();
		if (! myPoints.isEmpty())
			{
			myLastPoint = myPoints.get (myPoints.size()-1);
			}
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
		int n = myPoints.size();
		if (n < 2) return;
		if (myStartArrow == null && myEndArrow == null) return;

		// Prepare to draw arrows. Set stroke to a solid stroke the same width
		// as the current outline.
		float width = myOutline.getStrokeWidth();
		g2d.setStroke (new BasicStroke (width));
		double x;
		double y;
		double dx;
		double dy;
		double phi;

		if (myStartArrow != null)
			{
			Point p0 = myPoints.get (0);
			Point p1 = myPoints.get (1);
			x = p0.x;
			y = p0.y;
			dx = x - p1.x;
			dy = y - p1.y;
			phi =
				dx == 0.0 && dy == 0.0 ?
					0.0 :
					Math.atan2 (dy, dx);
			myStartArrow.draw (g2d, width, x, y, phi);
			}

		if (myEndArrow != null)
			{
			Point pnm1 = myPoints.get (n-1);
			Point pnm2 = myPoints.get (n-2);
			x = pnm1.x;
			y = pnm1.y;
			dx = x - pnm2.x;
			dy = y - pnm2.y;
			phi =
				dx == 0.0 && dy == 0.0 ?
					0.0 :
					Math.atan2 (dy, dx);
			myEndArrow.draw (g2d, width, x, y, phi);
			}
		}

// Hidden operations.

	/**
	 * Add the given endpoint to this line.
	 */
	private Line addEndpoint
		(Point p)
		{
		myLastPoint = p;
		myPoints.add (p);
		myNw = null;
		return this;
		}

	/**
	 * Compute this line's bounding box from myPoints. The results are stored in
	 * myNw, mySe, and mySize.
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
	 * Compute this line's path from myPoints.
	 */
	private GeneralPath computePath()
		{
		GeneralPath path = new GeneralPath();
		int n = myPoints.size();
		Iterator<Point> iter = myPoints.iterator();
		Point point = null;

		if (n == 0)
			{
			// No points. Empty path.
			}

		else if (n == 1)
			{
			// One point.
			point = iter.next();
			path.moveTo ((float) point.x, (float) point.y);
			}

		else if (n == 2)
			{
			// Two points.
			point = iter.next();
			path.moveTo ((float) point.x, (float) point.y);
			point = iter.next();
			path.lineTo ((float) point.x, (float) point.y);
			}

		else if (myRound == 0.0)
			{
			// Three or more points, sharp corners.
			point = iter.next();
			path.moveTo ((float) point.x, (float) point.y);
			while (iter.hasNext())
				{
				point = iter.next();
				path.lineTo ((float) point.x, (float) point.y);
				}
			}

		else
			{
			// Three or more points, round corners.
			double a_x, a_y;
			double b_x, b_y;
			double c_x, c_y;
			double d_x, d_y;
			double e_x, e_y;
			double f_x, f_y;
			double g_x, g_y;
			double current_x, current_y;
			double delta_x, delta_y, dist, p, q;
			double startLimit, endLimit;

			point = iter.next();
			path.moveTo ((float) point.x, (float) point.y);
			a_x = point.x; a_y = point.y;
			current_x = a_x; current_y = a_y;
			point = iter.next();
			b_x = point.x; b_y = point.y;
			startLimit = 1.0;
			endLimit = 0.5;
			while (iter.hasNext())
				{
				point = iter.next();
				c_x = point.x; c_y = point.y;
				if (! iter.hasNext()) endLimit = 1.0;

				// Points a, b, and c define the two line segments that join at
				// a corner.

				// Point d is the cubic Bezier curve starting point.
				delta_x = b_x - a_x;
				delta_y = b_y - a_y;
				dist = Math.sqrt (delta_x * delta_x + delta_y * delta_y);
				p = Math.min (myRound / dist, startLimit);
				q = 1.0 - p;
				d_x = p * a_x + q * b_x;
				d_y = p * a_y + q * b_y;

				// Point g is the cubic Bezier curve ending point.
				delta_x = b_x - c_x;
				delta_y = b_y - c_y;
				dist = Math.sqrt (delta_x * delta_x + delta_y * delta_y);
				p = Math.min (myRound / dist, endLimit);
				q = 1.0 - p;
				g_x = p * c_x + q * b_x;
				g_y = p * c_y + q * b_y;

				// Point e is the cubic Bezier curve first control point.
				p = 0.4477152501692067;
				q = 0.5522847498307933;
				e_x = p * d_x + q * b_x;
				e_y = p * d_y + q * b_y;

				// Point f is the cubic Bezier curve first control point.
				f_x = p * g_x + q * b_x;
				f_y = p * g_y + q * b_y;

				// Line from current point to d if necessary.
				if (current_x != d_x || current_y != d_y)
					{
					path.lineTo ((float) d_x, (float) d_y);
					}

				// Cubic Bezier curve from d to g.
				path.curveTo
					((float) e_x, (float) e_y,
					 (float) f_x, (float) f_y,
					 (float) g_x, (float) g_y);
				current_x = g_x; current_y = g_y;

				// Advance to the next corner.
				a_x = b_x; a_y = b_y;
				b_x = c_x; b_y = c_y;
				startLimit = 0.5;
				}

			// Line to final endpoint if necessary.
			if (current_x != b_x || current_y != b_y)
				{
				path.lineTo ((float) b_x, (float) b_y);
				}
			}

		return path;
		}

	}
