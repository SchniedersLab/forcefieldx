//******************************************************************************
//
// File:    Group.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Group
//
// This Java source file is copyright (C) 2011 by Alan Kaminsky. All rights
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

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.util.ArrayList;

/**
 * Class Group provides a {@linkplain DrawingItem} that consists of a group of
 * other {@linkplain DrawingItem}s. To add drawing items to a group, call the
 * group's <TT>append()</TT> and <TT>prepend()</TT> methods.
 * <P>
 * The group can be scaled as a unit by calling the <TT>xScale()</TT> and
 * <TT>yScale()</TT> methods. The group can be sheared as a unit by calling the
 * <TT>xShear()</TT> and <TT>yShear()</TT> methods. The group can be rotated as
 * a unit by calling the <TT>rotationAngle()</TT> method; the group is always
 * rotated around its center point. The order in which these transforms are
 * applied is first scaling, then shearing, then rotation.
 * <P>
 * The group's "bounding box" is the smallest rectangle that encloses all the
 * contained drawing items' bounding boxes, <I>after</I> applying the scale,
 * shear, and rotation transforms if any. The <TT>nw()</TT>, <TT>n()</TT>,
 * <TT>ne()</TT>, <TT>w()</TT>, <TT>c()</TT>, <TT>e()</TT>, <TT>sw()</TT>,
 * <TT>s()</TT>, and <TT>se()</TT> methods refer to points on the group's
 * bounding box. The setter versions of these methods translate the group so the
 * designated point on the group's bounding box coincides with the point
 * specified as the argument. The getter versions of these methods return the
 * designated point on the group's bounding box.
 * <P>
 * The group's "original bounding box" is the smallest rectangle that encloses
 * all the contained drawing items' bounding boxes, <I>before</I> applying the
 * scale, shear, and rotation transforms if any. The <TT>orig_nw()</TT>,
 * <TT>orig_n()</TT>, <TT>orig_ne()</TT>, <TT>orig_w()</TT>, <TT>orig_c()</TT>,
 * <TT>orig_e()</TT>, <TT>orig_sw()</TT>, <TT>orig_s()</TT>, and
 * <TT>orig_se()</TT> methods refer to points on the group's original bounding
 * box. The setter versions of these methods translate the group so the
 * designated point on the group's original bounding box coincides with the
 * point specified as the argument. The getter versions of these methods return
 * the designated point on the group's original bounding box.
 * <P>
 * The <TT>transform()</TT> method takes a point, applies the group's scaling,
 * shearing, rotation, and translation transforms to the point, and returns the
 * transformed point. You can use this method, for example, to find out where a
 * point defined before applying transforms to the group would end up after
 * applying transforms to the group.
 *
 * @author  Alan Kaminsky
 * @version 29-Apr-2011
 */
public class Group
	extends DrawingItem
	implements Externalizable
	{

// Exported constants.

	/**
	 * The normal X scale factor (1).
	 */
	public static final double NORMAL_X_SCALE = 1.0;

	/**
	 * The normal Y scale factor (1).
	 */
	public static final double NORMAL_Y_SCALE = 1.0;

	/**
	 * The normal X shear factor (0).
	 */
	public static final double NORMAL_X_SHEAR = 0.0;

	/**
	 * The normal Y shear factor (0).
	 */
	public static final double NORMAL_Y_SHEAR = 0.0;

	/**
	 * The normal rotation angle (0).
	 */
	public static final double NORMAL_ROTATION_ANGLE = 0.0;

// Hidden data members.

	private static final long serialVersionUID = -2445041343028527776L;

	// Default attributes.
	static double theDefaultXScale = NORMAL_X_SCALE;
	static double theDefaultYScale = NORMAL_Y_SCALE;
	static double theDefaultXShear = NORMAL_X_SHEAR;
	static double theDefaultYShear = NORMAL_Y_SHEAR;
	static double theDefaultRotationAngle = NORMAL_ROTATION_ANGLE;

	// List of drawing items in this group.
	ArrayList<DrawingItem> myItemList;

	// Coordinates of most recently specified corner.
	double x;
	double y;

	// Factors for going from specified corner to northwest corner.
	double xFactor;
	double yFactor;

	// True if corner was specified with respect to original (untransformed)
	// bounding box. False if corner was specified with respect to (transformed)
	// bounding box.
	boolean cornerIsOriginal;

	// Transform parameters.
	double myXScale = theDefaultXScale;
	double myYScale = theDefaultYScale;
	double myXShear = theDefaultXShear;
	double myYShear = theDefaultYShear;
	double myRotationAngle = theDefaultRotationAngle;

	// Transform. If myTransform is null, the transform must be recomputed.
	AffineTransform myTransform;

	// Northwest point and size of original bounding box. If myOriginalNw is
	// null, these objects must be recomputed.
	Point myOriginalNw;
	Size myOriginalSize;

	// Northwest point and size of transformed bounding box. If myTransformedNw
	// is null, these objects must be recomputed.
	Point myTransformedNw;
	Size myTransformedSize;

// Exported constructors.

	/**
	 * Construct a new empty group. The group's northwest corner is located at
	 * (0,0).
	 */
	public Group()
		{
		super();
		this.myItemList = new ArrayList<DrawingItem>();
		}

	/**
	 * Construct a new group that contains the same drawing items, has the same
	 * transforms, and is at the same location, as the given group.
	 *
	 * @param  theGroup  Group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theGroup</TT> is null.
	 */
	public Group
		(Group theGroup)
		{
		super (theGroup);
		this.myItemList = new ArrayList<DrawingItem> (theGroup.myItemList);
		this.x = theGroup.x;
		this.y = theGroup.y;
		this.xFactor = theGroup.xFactor;
		this.yFactor = theGroup.yFactor;
		this.cornerIsOriginal = theGroup.cornerIsOriginal;
		this.myXScale = theGroup.myXScale;
		this.myYScale = theGroup.myYScale;
		this.myXShear = theGroup.myXShear;
		this.myYShear = theGroup.myYShear;
		this.myRotationAngle = theGroup.myRotationAngle;
		}

// Exported operations.

	/**
	 * Returns the default X scale factor.
	 *
	 * @return  Default X scale factor.
	 */
	public static double defaultXScale()
		{
		return theDefaultXScale;
		}

	/**
	 * Set the default X scale factor. Before calling this method the first
	 * time, the default X scale factor is 1.
	 *
	 * @param  theXScale  Default X scale factor.
	 */
	public static void defaultXScale
		(double theXScale)
		{
		theDefaultXScale = theXScale;
		}

	/**
	 * Returns the default Y scale factor.
	 *
	 * @return  Default Y scale factor.
	 */
	public static double defaultYScale()
		{
		return theDefaultYScale;
		}

	/**
	 * Set the default Y scale factor. Before calling this method the first
	 * time, the default Y scale factor is 1.
	 *
	 * @param  theYScale  Default Y scale factor.
	 */
	public static void defaultYScale
		(double theYScale)
		{
		theDefaultYScale = theYScale;
		}

	/**
	 * Returns the default X shear factor.
	 *
	 * @return  Default X shear factor.
	 */
	public static double defaultXShear()
		{
		return theDefaultXShear;
		}

	/**
	 * Set the default X shear factor. Before calling this method the first
	 * time, the default X shear factor is 0.
	 *
	 * @param  theXShear  Default X shear factor.
	 */
	public static void defaultXShear
		(double theXShear)
		{
		theDefaultXShear = theXShear;
		}

	/**
	 * Returns the default Y shear factor.
	 *
	 * @return  Default Y shear factor.
	 */
	public static double defaultYShear()
		{
		return theDefaultYShear;
		}

	/**
	 * Set the default Y shear factor. Before calling this method the first
	 * time, the default Y shear factor is 0.
	 *
	 * @param  theYShear  Default Y shear factor.
	 */
	public static void defaultYShear
		(double theYShear)
		{
		theDefaultYShear = theYShear;
		}

	/**
	 * Returns the default rotation angle.
	 *
	 * @return  Default rotation angle.
	 */
	public static double defaultRotationAngle()
		{
		return theDefaultRotationAngle;
		}

	/**
	 * Set the default rotation angle. Before calling this method the first
	 * time, the default rotation angle is 0.
	 *
	 * @param  theAngle  Default rotation angle.
	 */
	public static void defaultRotationAngle
		(double theAngle)
		{
		theDefaultRotationAngle = theAngle;
		}

	/**
	 * Returns this group's X scale factor.
	 *
	 * @return  X scale factor.
	 */
	public double xScale()
		{
		return myXScale;
		}

	/**
	 * Set this group's X scale factor.
	 *
	 * @param  theXScale  X scale factor.
	 *
	 * @return  This group.
	 */
	public Group xScale
		(double theXScale)
		{
		myXScale = theXScale;
		myTransform = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Returns this group's Y scale factor.
	 *
	 * @return  Y scale factor.
	 */
	public double yScale()
		{
		return myYScale;
		}

	/**
	 * Set this group's Y scale factor.
	 *
	 * @param  theYScale  Y scale factor.
	 *
	 * @return  This group.
	 */
	public Group yScale
		(double theYScale)
		{
		myYScale = theYScale;
		myTransform = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Returns this group's X shear factor.
	 *
	 * @return  X shear factor.
	 */
	public double xShear()
		{
		return myXShear;
		}

	/**
	 * Set this group's X shear factor.
	 *
	 * @param  theXShear  X shear factor.
	 *
	 * @return  This group.
	 */
	public Group xShear
		(double theXShear)
		{
		myXShear = theXShear;
		myTransform = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Returns this group's Y shear factor.
	 *
	 * @return  Y shear factor.
	 */
	public double yShear()
		{
		return myYShear;
		}

	/**
	 * Set this group's Y shear factor.
	 *
	 * @param  theYShear  Y shear factor.
	 *
	 * @return  This group.
	 */
	public Group yShear
		(double theYShear)
		{
		myYShear = theYShear;
		myTransform = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Returns this group's rotation angle.
	 *
	 * @return  Rotation angle.
	 */
	public double rotationAngle()
		{
		return myRotationAngle;
		}

	/**
	 * Set this group's rotation angle.
	 *
	 * @param  theAngle  Rotation angle.
	 *
	 * @return  This group.
	 */
	public Group rotationAngle
		(double theAngle)
		{
		myRotationAngle = theAngle;
		myTransform = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Clear this group. All drawing items in this group are removed.
	 *
	 * @return  This group.
	 */
	public Group clear()
		{
		myItemList.clear();
		myOriginalNw = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Add the given drawing item to the end of this group's list of drawing
	 * items.
	 *
	 * @param  theItem  Drawing item.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Group append
		(DrawingItem theItem)
		{
		if (theItem == null) throw new NullPointerException();
		myItemList.add (theItem);
		myOriginalNw = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Add the given drawing item to the beginning of this group's list of
	 * drawing items.
	 *
	 * @param  theItem  Drawing item.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Group prepend
		(DrawingItem theItem)
		{
		if (theItem == null) throw new NullPointerException();
		myItemList.add (0, theItem);
		myOriginalNw = null;
		myTransformedNw = null;
		return this;
		}

	/**
	 * Returns the size of this group's bounding box.
	 *
	 * @return  Size.
	 */
	public Size size()
		{
		computeTransformedBoundingBox();
		return myTransformedSize;
		}

	/**
	 * Returns the width of this group's bounding box.
	 *
	 * @return  Width.
	 */
	public double width()
		{
		computeTransformedBoundingBox();
		return myTransformedSize.width();
		}

	/**
	 * Returns the height of this group's bounding box.
	 *
	 * @return  Height.
	 */
	public double height()
		{
		computeTransformedBoundingBox();
		return myTransformedSize.height();
		}

	/**
	 * Returns the northwest corner point of this group's bounding box.
	 *
	 * @return  Northwest corner point.
	 */
	public Point nw()
		{
		computeTransformedBoundingBox();
		return myTransformedNw;
		}

	/**
	 * Set the northwest corner point of this group's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This group.
	 */
	public Group nw
		(double x,
		 double y)
		{
		doNw (x, y, false);
		return this;
		}

	/**
	 * Set the northwest corner point of this group's bounding box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the north middle point of this group's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This group.
	 */
	public Group n
		(double x,
		 double y)
		{
		doN (x, y, false);
		return this;
		}

	/**
	 * Set the north middle point of this group's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the northeast corner point of this group's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This group.
	 */
	public Group ne
		(double x,
		 double y)
		{
		doNe (x, y, false);
		return this;
		}

	/**
	 * Set the northeast corner point of this group's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the west middle point of this group's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This group.
	 */
	public Group w
		(double x,
		 double y)
		{
		doW (x, y, false);
		return this;
		}

	/**
	 * Set the west middle point of this group's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the center point of this group's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This group.
	 */
	public Group c
		(double x,
		 double y)
		{
		doC (x, y, false);
		return this;
		}

	/**
	 * Set the center point of this group's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the east middle point of this group's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This group.
	 */
	public Group e
		(double x,
		 double y)
		{
		doE (x, y, false);
		return this;
		}

	/**
	 * Set the east middle point of this group's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the southwest corner point of this group's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This group.
	 */
	public Group sw
		(double x,
		 double y)
		{
		doSw (x, y, false);
		return this;
		}

	/**
	 * Set the southwest corner point of this group's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the south middle point of this group's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This group.
	 */
	public Group s
		(double x,
		 double y)
		{
		doS (x, y, false);
		return this;
		}

	/**
	 * Set the south middle point of this group's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Set the southeast corner point of this group's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This group.
	 */
	public Group se
		(double x,
		 double y)
		{
		doSe (x, y, false);
		return this;
		}

	/**
	 * Set the southeast corner point of this group's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y, false);
		return this;
		}

	/**
	 * Returns the northwest corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  Northwest corner point.
	 */
	public Point orig_nw()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x,
				 myOriginalNw.y);
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the northwest corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This group.
	 */
	public Group orig_nw
		(double x,
		 double y)
		{
		doNw (x, y, true);
		return this;
		}

	/**
	 * Set the northwest corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the north center point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  North center point.
	 */
	public Point orig_n()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x + 0.5*myOriginalSize.width,
				 myOriginalNw.y);
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the north center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of north center point.
	 * @param  y  Y coordinate of north center point.
	 *
	 * @return  This group.
	 */
	public Group orig_n
		(double x,
		 double y)
		{
		doN (x, y, true);
		return this;
		}

	/**
	 * Set the north center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  thePoint  North center point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the northeast corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  Northeast corner point.
	 */
	public Point orig_ne()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x + myOriginalSize.width(),
				 myOriginalNw.y);
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the northeast corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This group.
	 */
	public Group orig_ne
		(double x,
		 double y)
		{
		doNe (x, y, true);
		return this;
		}

	/**
	 * Set the northeast corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the west center point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  West center point.
	 */
	public Point orig_w()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x,
				 myOriginalNw.y + 0.5*myOriginalSize.height());
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the west center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of west center point.
	 * @param  y  Y coordinate of west center point.
	 *
	 * @return  This group.
	 */
	public Group orig_w
		(double x,
		 double y)
		{
		doW (x, y, true);
		return this;
		}

	/**
	 * Set the west center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  thePoint  West center point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @return  Center point.
	 */
	public Point orig_c()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x + 0.5*myOriginalSize.width(),
				 myOriginalNw.y + 0.5*myOriginalSize.height());
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This group.
	 */
	public Group orig_c
		(double x,
		 double y)
		{
		doC (x, y, true);
		return this;
		}

	/**
	 * Set the center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the east center point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  East center point.
	 */
	public Point orig_e()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x + myOriginalSize.width(),
				 myOriginalNw.y + 0.5*myOriginalSize.height());
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the east center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of east center point.
	 * @param  y  Y coordinate of east center point.
	 *
	 * @return  This group.
	 */
	public Group orig_e
		(double x,
		 double y)
		{
		doE (x, y, true);
		return this;
		}

	/**
	 * Set the east center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  thePoint  East center point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the southwest corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  Southwest corner point.
	 */
	public Point orig_sw()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x,
				 myOriginalNw.y + myOriginalSize.height());
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the southwest corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This group.
	 */
	public Group orig_sw
		(double x,
		 double y)
		{
		doSw (x, y, true);
		return this;
		}

	/**
	 * Set the southwest corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the south center point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  South center point.
	 */
	public Point orig_s()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x + 0.5*myOriginalSize.width(),
				 myOriginalNw.y + myOriginalSize.height());
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the south center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of south center point.
	 * @param  y  Y coordinate of south center point.
	 *
	 * @return  This group.
	 */
	public Group orig_s
		(double x,
		 double y)
		{
		doS (x, y, true);
		return this;
		}

	/**
	 * Set the south center point of this group's original bounding box, after
	 * applying the group's transformations if any.
	 *
	 * @param  thePoint  South center point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the southeast corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @return  Southeast corner point.
	 */
	public Point orig_se()
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d =
			new Point2D.Double
				(myOriginalNw.x + myOriginalSize.width(),
				 myOriginalNw.y + myOriginalSize.height());
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

	/**
	 * Set the southeast corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This group.
	 */
	public Group orig_se
		(double x,
		 double y)
		{
		doSe (x, y, true);
		return this;
		}

	/**
	 * Set the southeast corner point of this group's original bounding box,
	 * after applying the group's transformations if any.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Group orig_se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y, true);
		return this;
		}

	/**
	 * Returns the northwest corner point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  Northwest corner point.
	 */
	public Point content_nw()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x,
			 myOriginalNw.y);
		}

	/**
	 * Returns the north center point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  North center point.
	 */
	public Point content_n()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x + 0.5*myOriginalSize.width,
			 myOriginalNw.y);
		}

	/**
	 * Returns the northeast corner point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  Northeast corner point.
	 */
	public Point content_ne()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x + myOriginalSize.width,
			 myOriginalNw.y);
		}

	/**
	 * Returns the west center point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  West center point.
	 */
	public Point content_w()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x,
			 myOriginalNw.y + 0.5*myOriginalSize.height);
		}

	/**
	 * Returns the center point of this group's original bounding box, without
	 * applying the group's transformations if any.
	 *
	 * @return  Center point.
	 */
	public Point content_c()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x + 0.5*myOriginalSize.width,
			 myOriginalNw.y + 0.5*myOriginalSize.height);
		}

	/**
	 * Returns the east center point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  East center point.
	 */
	public Point content_e()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x + myOriginalSize.width,
			 myOriginalNw.y + 0.5*myOriginalSize.height);
		}

	/**
	 * Returns the southwest corner point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  Southwest corner point.
	 */
	public Point content_sw()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x,
			 myOriginalNw.y + myOriginalSize.height);
		}

	/**
	 * Returns the south center point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  South center point.
	 */
	public Point content_s()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x + 0.5*myOriginalSize.width,
			 myOriginalNw.y + myOriginalSize.height);
		}

	/**
	 * Returns the southeast corner point of this group's original bounding box,
	 * without applying the group's transformations if any.
	 *
	 * @return  Southeast corner point.
	 */
	public Point content_se()
		{
		computeOriginalBoundingBox();
		return new Point
			(myOriginalNw.x + myOriginalSize.width,
			 myOriginalNw.y + myOriginalSize.height);
		}

	/**
	 * Add this group to the end of the default drawing's sequence of drawing
	 * items.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Group add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this group to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Group add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this group to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Group addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this group to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This group.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Group addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this group to the given object output stream.
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
		int n = myItemList.size();
		out.writeInt (n);
		for (DrawingItem item : myItemList)
			{
			out.writeObject (item);
			}
		out.writeDouble (x);
		out.writeDouble (y);
		out.writeDouble (xFactor);
		out.writeDouble (yFactor);
		out.writeBoolean (cornerIsOriginal);
		out.writeDouble (myXScale);
		out.writeDouble (myYScale);
		out.writeDouble (myXShear);
		out.writeDouble (myYShear);
		out.writeDouble (myRotationAngle);
		}

	/**
	 * Read this group from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this group cannot be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		myItemList.clear();
		int n = in.readInt();
		for (int i = 0; i < n; ++ i)
			{
			myItemList.add ((DrawingItem) in.readObject());
			}
		x = in.readDouble();
		y = in.readDouble();
		xFactor = in.readDouble();
		yFactor = in.readDouble();
		cornerIsOriginal = in.readBoolean();
		myXScale = in.readDouble();
		myYScale = in.readDouble();
		myXShear = in.readDouble();
		myYShear = in.readDouble();
		myRotationAngle = in.readDouble();
		myTransform = null;
		myOriginalNw = null;
		myTransformedNw = null;
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

		// Apply the scaling, shearing, rotation, and translation transforms.
		computeTransformedBoundingBox();
		AffineTransform oldTransform = g2d.getTransform();
		g2d.transform (myTransform);

		// Draw each item.
		for (DrawingItem item : myItemList)
			{
			item.draw (g2d);
			}

		// Restore transform. This is necessary to draw nested groups.
		g2d.setTransform (oldTransform);
		}

	/**
	 * Transform the given point by applying this group's scaling, shearing,
	 * rotation, and translation transforms.
	 *
	 * @param  thePoint  Point before applying transforms.
	 *
	 * @return  Point after applying transforms.
	 */
	public Point transform
		(Point thePoint)
		{
		computeTransformedBoundingBox();
		Point2D.Double p2d = new Point2D.Double (thePoint.x, thePoint.y);
		myTransform.transform (p2d, p2d);
		return new Point (p2d.x, p2d.y);
		}

// Hidden operations.

	/**
	 * Compute the original bounding box of this group. Store the results in the
	 * myOriginalNw and myOriginalSize fields.
	 */
	private void computeOriginalBoundingBox()
		{
		if (myOriginalNw == null)
			{
			Point p;
			double minX = Double.MAX_VALUE;
			double minY = Double.MAX_VALUE;
			double maxX = Double.MIN_VALUE;
			double maxY = Double.MIN_VALUE;
			for (DrawingItem item : myItemList)
				{
				p = item.nw();
				minX = Math.min (minX, p.x);
				minY = Math.min (minY, p.y);
				p = item.se();
				maxX = Math.max (maxX, p.x);
				maxY = Math.max (maxY, p.y);
				}
			myOriginalNw = new Point (minX, minY);
			myOriginalSize = new Size (maxX - minX, maxY - minY);
			}
		}

	/**
	 * Compute the transformed bounding box of this group. Store the results in
	 * the myTransformedNw and myTransformedSize fields.
	 */
	private void computeTransformedBoundingBox()
		{
		if (myTransformedNw == null)
			{
			Point p;
			Point2D.Double p2d = new Point2D.Double();
			AffineTransform t = new AffineTransform();
			AffineTransform t2 = new AffineTransform();

			// Compute original bounding box, to find its center.
			computeOriginalBoundingBox();

			// Set up transform with just scale and shear.
			t.shear
				(/*shx*/ myXShear,
				 /*shy*/ myYShear);
			t.scale
				(/*sx*/ myXScale,
				 /*sy*/ myYScale);

			// Scale and shear the original bounding box's center.
			p2d.x = myOriginalNw.x() + myOriginalSize.width()/2;
			p2d.y = myOriginalNw.y() + myOriginalSize.height()/2;
			t.transform (p2d, p2d);

			// Prepend rotation to transform.
			t2.setToRotation
				(/*theta*/ myRotationAngle,
				 /*x    */ p2d.x,
				 /*y    */ p2d.y);
			t.preConcatenate (t2);

			// Compute bounding box with just scale, shear, and rotation.
			double minX = Double.MAX_VALUE;
			double minY = Double.MAX_VALUE;
			double maxX = Double.MIN_VALUE;
			double maxY = Double.MIN_VALUE;
			for (DrawingItem item : myItemList)
				{
				p = item.nw(); p2d.x = p.x; p2d.y = p.y;
				t.transform (p2d, p2d);
				minX = Math.min (minX, p2d.x);
				minY = Math.min (minY, p2d.y);
				maxX = Math.max (maxX, p2d.x);
				maxY = Math.max (maxY, p2d.y);
				p = item.ne(); p2d.x = p.x; p2d.y = p.y;
				t.transform (p2d, p2d);
				minX = Math.min (minX, p2d.x);
				minY = Math.min (minY, p2d.y);
				maxX = Math.max (maxX, p2d.x);
				maxY = Math.max (maxY, p2d.y);
				p = item.sw(); p2d.x = p.x; p2d.y = p.y;
				t.transform (p2d, p2d);
				minX = Math.min (minX, p2d.x);
				minY = Math.min (minY, p2d.y);
				maxX = Math.max (maxX, p2d.x);
				maxY = Math.max (maxY, p2d.y);
				p = item.se(); p2d.x = p.x; p2d.y = p.y;
				t.transform (p2d, p2d);
				minX = Math.min (minX, p2d.x);
				minY = Math.min (minY, p2d.y);
				maxX = Math.max (maxX, p2d.x);
				maxY = Math.max (maxY, p2d.y);
				}

			// Record transformed bounding box's size.
			myTransformedSize = new Size (maxX - minX, maxY - minY);

			// Corner was specified with respect to original bounding box.
			if (cornerIsOriginal)
				{
				// Compute original corner point based on (xFactor,yFactor).
				p2d.x = myOriginalNw.x - xFactor*myOriginalSize.width();
				p2d.y = myOriginalNw.y - yFactor*myOriginalSize.height();

				// Transform it by scale, shear, and rotation.
				t.transform (p2d, p2d);

				// Translate it to the point (x,y).
				t2.setToTranslation
					(/*tx*/ x - p2d.x,
					 /*ty*/ y - p2d.y);
				t.preConcatenate (t2);

				// Compute transformed NW point.
				myTransformedNw =
					new Point (minX + x - p2d.x, minY + y - p2d.y);
				}

			// Corner was specified with respect to transformed bounding box.
			else
				{
				// Compute transformed NW point based on (x,y) and
				// (xFactor,yFactor).
				myTransformedNw =
					new Point
						(x + xFactor*myTransformedSize.width(),
						 y + yFactor*myTransformedSize.height());

				// Translate the point (minX,minY) to the point myTransformedNw.
				t2.setToTranslation
					(/*tx*/ myTransformedNw.x - minX,
					 /*ty*/ myTransformedNw.y - minY);
				t.preConcatenate (t2);
				}

			// Save final transformation.
			myTransform = t;
			}
		}

	/**
	 * Set the northwest corner point of this group's bounding
	 * box.
	 *
	 * @param  x     X coordinate of northwest corner point.
	 * @param  y     Y coordinate of northwest corner point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doNw
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = 0.0;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the north middle point of this group's bounding box.
	 *
	 * @param  x     X coordinate of north middle point.
	 * @param  y     Y coordinate of north middle point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doN
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = 0.0;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the northeast corner point of this group's bounding box.
	 *
	 * @param  x     X coordinate of northeast corner point.
	 * @param  y     Y coordinate of northeast corner point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doNe
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = 0.0;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the west middle point of this group's bounding box.
	 *
	 * @param  x     X coordinate of west middle point.
	 * @param  y     Y coordinate of west middle point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doW
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = -0.5;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the center point of this group's bounding box.
	 *
	 * @param  x     X coordinate of center point.
	 * @param  y     Y coordinate of center point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doC
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = -0.5;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the east middle point of this group's bounding box.
	 *
	 * @param  x     X coordinate of east middle point.
	 * @param  y     Y coordinate of east middle point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doE
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = -0.5;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the southwest corner point of this group's bounding box.
	 *
	 * @param  x     X coordinate of southwest corner point.
	 * @param     y  Y coordinate of southwest corner point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doSw
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = -1.0;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the south middle point of this group's bounding box.
	 *
	 * @param  x     X coordinate of south middle point.
	 * @param  y     Y coordinate of south middle point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doS
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = -1.0;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	/**
	 * Set the southeast corner point of this group's bounding box.
	 *
	 * @param  x     X coordinate of southeast corner point.
	 * @param  y     Y coordinate of southeast corner point.
	 * @param  orig  True if original corner, false if transformed corner.
	 */
	void doSe
		(double x,
		 double y,
		 boolean orig)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = -1.0;
		this.cornerIsOriginal = orig;
		this.myOriginalNw = null;
		this.myTransformedNw = null;
		}

	}
