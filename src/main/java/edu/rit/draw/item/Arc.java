//******************************************************************************
//
// File:    Arc.java
// Package: edu.rit.draw.Arc
// Unit:    Class edu.rit.draw.item.Arc
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

import java.awt.Graphics2D;
import java.awt.Shape;

import java.awt.geom.Arc2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Arc provides a {@linkplain DrawingItem} that is a circular arc. The arc
 * has an outline but no filled interior. The arc is defined by the following
 * attributes:
 * <UL>
 * <LI>
 * <TT>center()</TT> -- The arc's center point. The default is (0,0).
 * <BR>&nbsp;
 * <LI>
 * <TT>radius()</TT> -- The arc's radius. The default is 1.
 * <BR>&nbsp;
 * <LI>
 * <TT>start()</TT> -- The arc's starting angle in radians. This is the angle
 * from the positive X axis to the arc's starting point. The default is 0.
 * <BR>&nbsp;
 * <LI>
 * <TT>extent()</TT> -- The arc's angular extent in radians. This is the angle
 * from the arc's starting point to the arc's ending point. The default is
 * &pi;/2.
 * </UL>
 * <P>
 * Angles increase from the positive X axis towards the positive Y axis. Since
 * the positive X axis points to the right and the positive Y axis points down,
 * angles increase in the clockwise direction. (This is the opposite of the
 * usual mathematical convention, where the positive Y axis points upwards.)
 * <P>
 * The static <TT>defaultCenter()</TT>, <TT>defaultRadius()</TT>,
 * <TT>defaultStart()</TT>, and <TT>defaultExtent()</TT> methods are provided to
 * set the default center point, radius, starting angle, and angular extent. If
 * the center point, radius, starting angle, or angular extent is not specified,
 * the current default center point, radius, starting angle, or angular extent
 * is used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Arc
	extends OutlinedItem
	implements Externalizable
	{

// Exported constants.

	/**
	 * The normal center point: (0,0).
	 */
	public static final Point NORMAL_CENTER = new Point (0, 0);

	/**
	 * The normal radius: 1.
	 */
	public static final double NORMAL_RADIUS = 1.0;

	/**
	 * The normal starting angle: 0.
	 */
	public static final double NORMAL_START = 0.0;

	/**
	 * The normal angular extent: &pi;/2.
	 */
	public static final double NORMAL_EXTENT = Math.PI / 2.0;

// Hidden data members.

	private static final long serialVersionUID = 4876674535249335766L;

	private static final double DEGREES_PER_RADIAN = 180.0 / Math.PI;
	private static final double A_PI_2 = Math.PI / 2.0;
	private static final double A_PI = Math.PI;
	private static final double A_3_PI_2 = 3.0 * Math.PI / 2.0;
	private static final double A_2_PI = 2.0 * Math.PI;
	private static final double A_5_PI_2 = 5.0 * Math.PI / 2.0;
	private static final double A_3_PI = 3.0 * Math.PI;
	private static final double A_7_PI_2 = 7.0 * Math.PI / 2.0;

	// Attributes.
	Point myCenter = theDefaultCenter;
	double myRadius = theDefaultRadius;
	double myStart = theDefaultStart;
	double myExtent = theDefaultExtent;

	// The default attributes.
	static Point theDefaultCenter = NORMAL_CENTER;
	static double theDefaultRadius = NORMAL_RADIUS;
	static double theDefaultStart = NORMAL_START;
	static double theDefaultExtent = NORMAL_EXTENT;

	// For computing the bounding box.
	Point myNw;
	Point mySe;
	Size mySize;

// Exported constructors.

	/**
	 * Construct a new arc. The default center point, radius, starting angle,
	 * and angular extent are used.
	 */
	public Arc()
		{
		super();
		}

	/**
	 * Construct a new arc with the same center point, radius, starting angle,
	 * and angular extent as the given arc.
	 *
	 * @param  theArc  Arc.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArc</TT> is null.
	 */
	public Arc
		(Arc theArc)
		{
		super (theArc);
		myCenter = theArc.myCenter;
		myRadius = theArc.myRadius;
		myStart = theArc.myStart;
		myExtent = theArc.myExtent;
		}

// Exported operations.

	/**
	 * Returns the default center point.
	 *
	 * @return  Default center point.
	 */
	public static Point defaultCenter()
		{
		return theDefaultCenter;
		}

	/**
	 * Set the default center point. Before calling this method the first time,
	 * the default center point is (0,0).
	 *
	 * @param  theCenter  Default center point.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theCenter</TT> is null.
	 */
	public static void defaultCenter
		(Point theCenter)
		{
		if (theCenter == null) throw new NullPointerException();
		theDefaultCenter = theCenter;
		}

	/**
	 * Returns the default radius.
	 *
	 * @return  Default radius.
	 */
	public static double defaultRadius()
		{
		return theDefaultRadius;
		}

	/**
	 * Set the default radius. Before calling this method the first time, the
	 * default radius is 1.
	 *
	 * @param  theRadius  Default radius.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theRadius</TT> &lt;= 0.
	 */
	public static void defaultRadius
		(double theRadius)
		{
		if (theRadius <= 0.0) throw new IllegalArgumentException();
		theDefaultRadius = theRadius;
		}

	/**
	 * Returns the default starting angle.
	 *
	 * @return  Default starting angle (radians).
	 */
	public static double defaultStart()
		{
		return theDefaultStart;
		}

	/**
	 * Set the default starting angle. Before calling this method the first
	 * time, the default starting angle is 0.
	 *
	 * @param  theStart  Default starting angle (radians).
	 */
	public static void defaultStart
		(double theStart)
		{
		theDefaultStart = theStart;
		}

	/**
	 * Returns the default angular extent.
	 *
	 * @return  Default angular extent (radians).
	 */
	public static double defaultExtent()
		{
		return theDefaultExtent;
		}

	/**
	 * Set the default angular extent. Before calling this method the first
	 * time, the default angular extent is &pi;/2.
	 *
	 * @param  theExtent  Default angular extent (radians).
	 */
	public static void defaultExtent
		(double theExtent)
		{
		theDefaultExtent = theExtent;
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
	 * Set this arc's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This arc.
	 */
	public Arc outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Returns this arc's center point.
	 *
	 * @return  Center point.
	 */
	public Point center()
		{
		return myCenter;
		}

	/**
	 * Set this arc's center point.
	 *
	 * @param  theCenter  Center point.
	 *
	 * @return  This arc.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theCenter</TT> is null.
	 */
	public Arc center
		(Point theCenter)
		{
		if (theCenter == null) throw new NullPointerException();
		myCenter = theCenter;
		return this;
		}

	/**
	 * Returns this arc's radius.
	 *
	 * @return  Radius.
	 */
	public double radius()
		{
		return myRadius;
		}

	/**
	 * Set this arc's radius.
	 *
	 * @param  theRadius  Radius.
	 *
	 * @return  This arc.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theRadius</TT> &lt;= 0.
	 */
	public Arc radius
		(double theRadius)
		{
		if (theRadius <= 0.0) throw new IllegalArgumentException();
		myRadius = theRadius;
		return this;
		}

	/**
	 * Returns this arc's starting angle.
	 *
	 * @return  Starting angle (radians).
	 */
	public double start()
		{
		return myStart;
		}

	/**
	 * Set this arc's starting angle.
	 *
	 * @param  theStart  Starting angle (radians).
	 *
	 * @return  This arc.
	 */
	public Arc start
		(double theStart)
		{
		myStart = theStart;
		return this;
		}

	/**
	 * Returns this arc's angular extent.
	 *
	 * @return  Angular extent (radians).
	 */
	public double extent()
		{
		return myExtent;
		}

	/**
	 * Set this arc's angular extent.
	 *
	 * @param  theExtent  Angular extent (radians).
	 *
	 * @return  This arc.
	 */
	public Arc extent
		(double theExtent)
		{
		myExtent = theExtent;
		return this;
		}

	/**
	 * Add this arc to the end of the default drawing's sequence of drawing
	 * items.
	 *
	 * @return  This arc.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Arc add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this arc to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This arc.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Arc add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this arc to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This arc.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Arc addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this arc to the beginning of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This arc.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Arc addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this arc to the given object output stream.
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
		out.writeObject (myCenter);
		out.writeDouble (myRadius);
		out.writeDouble (myStart);
		out.writeDouble (myExtent);
		}

	/**
	 * Read this arc from the given object input stream.
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
		myCenter = (Point) in.readObject();
		myRadius = in.readDouble();
		myStart = in.readDouble();
		myExtent = in.readDouble();
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

		// Draw arc.
		Shape shape = computeShape();
		myOutline.setGraphicsContext (g2d);
		g2d.draw (shape);
		}

// Hidden operations.

	/**
	 * Compute this arc's bounding box. The results are stored in
	 * myNw, mySe, and mySize.
	 */
	private void computeBoundingBox()
		{
		if (myNw == null)
			{
			double extent, start1, start2, sx, sy, ex, ey, x1, y1, x2, y2;

			// Make the angular extent positive if necessary.
			if (myExtent >= 0.0)
				{
				extent = myExtent;
				start1 = myStart;
				}
			else
				{
				extent = - myExtent;
				start1 = myStart - extent;
				}
				
			// Compute starting angle (mod 2 pi).
			start2 = Math.abs (start1);
			start2 = start2 - Math.floor (start2 / A_2_PI) * A_2_PI;
			if (start1 < 0.0) start2 = A_2_PI - start2;

			// Compute starting point.
			sx = myCenter.x() + myRadius * Math.cos (start2);
			sy = myCenter.y() + myRadius * Math.sin (start2);

			// Compute ending point.
			ex = myCenter.x() + myRadius * Math.cos (start2 + extent);
			ey = myCenter.y() + myRadius * Math.sin (start2 + extent);

			// Compute initial bounding box from starting and ending points.
			x1 = Math.min (sx, ex);
			y1 = Math.min (sy, ey);
			x2 = Math.max (sx, ex);
			y2 = Math.max (sy, ey);

			// Adjust bounding box if arc crosses the positive X axis.
			if (start2 < A_2_PI && start2 + extent > A_2_PI)
				{
				x2 = myCenter.x() + myRadius;
				}

			// Adjust bounding box if arc crosses the positive Y axis.
			if
				((start2 < A_PI_2 && start2 + extent > A_PI_2) ||
				 (start2 < A_5_PI_2 && start2 + extent > A_5_PI_2))
				{
				y2 = myCenter.y() + myRadius;
				}

			// Adjust bounding box if arc crosses the negative X axis.
			if
				((start2 < A_PI && start2 + extent > A_PI) ||
				 (start2 < A_3_PI && start2 + extent > A_3_PI))
				{
				x1 = myCenter.x() - myRadius;
				}

			// Adjust bounding box if arc crosses the negative Y axis.
			if
				((start2 < A_3_PI_2 && start2 + extent > A_3_PI_2) ||
				 (start2 < A_7_PI_2 && start2 + extent > A_7_PI_2))
				{
				y1 = myCenter.y() - myRadius;
				}

			// Compute final bounding box.
			myNw = new Point (x1, y1);
			mySe = new Point (x2, y2);
			mySize = mySe.difference (myNw);
			}
		}

	/**
	 * Compute this arc's shape from its attributes.
	 */
	private Shape computeShape()
		{
		return
			new Arc2D.Double
				(/*x     */ myCenter.x() - myRadius,
				 /*y     */ myCenter.y() - myRadius,
				 /*w     */ 2 * myRadius,
				 /*h     */ 2 * myRadius,
				 /*start */ -myStart * DEGREES_PER_RADIAN,
				 /*extent*/ -myExtent * DEGREES_PER_RADIAN,
				 /*type  */ Arc2D.OPEN);
		}

	}
