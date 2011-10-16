//******************************************************************************
//
// File:    Oval.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.Oval
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

import java.awt.geom.RoundRectangle2D;

/**
 * Class Oval provides an oval {@linkplain DrawingItem}. An oval is a rectangle
 * whose narrower sides have been rounded into half-circles:
 * <P>
 * <CENTER><IMG SRC="doc-files/OvalFig01.png"></CENTER>
 * <P>
 * The static <TT>defaultSize()</TT>, <TT>defaultWidth()</TT>, and
 * <TT>defaultHeight()</TT> methods are provided to set the default size, width,
 * and height. If an oval's size, width, or height is not specified, the current
 * default size, width, or height is used.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public class Oval
	extends RectangularItem
	{

// Exported constants.

	/**
	 * The normal width for ovals (72). <I>Note:</I> 72 points = 1 inch.
	 */
	public static final double NORMAL_WIDTH = 72.0;

	/**
	 * The normal height for ovals (72). <I>Note:</I> 72 points = 1 inch.
	 */
	public static final double NORMAL_HEIGHT = 72.0;

// Hidden data members.

	private static final long serialVersionUID = 2920893107679792184L;

	private static double theDefaultWidth = NORMAL_WIDTH;
	private static double theDefaultHeight = NORMAL_HEIGHT;

// Exported constructors.

	/**
	 * Construct a new oval. The oval's northwest corner is located at (0,0).
	 * The oval's size is the default size.
	 */
	public Oval()
		{
		super();
		this.width = theDefaultWidth;
		this.height = theDefaultHeight;
		}

	/**
	 * Construct a new oval with the same outline, fill paint, location, and
	 * size as the given oval.
	 *
	 * @param  theItem  Oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public Oval
		(Oval theItem)
		{
		super (theItem);
		}

// Exported operations.

	/**
	 * Returns the default size for ovals.
	 *
	 * @return  Default size.
	 */
	public static Size defaultSize()
		{
		return new Size (theDefaultWidth, theDefaultHeight);
		}

	/**
	 * Set the default size for ovals.
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
	 * Returns the default width for ovals.
	 *
	 * @return  Default width.
	 */
	public static double defaultWidth()
		{
		return theDefaultWidth;
		}

	/**
	 * Set the default width for ovals.
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
	 * Returns the default height for ovals.
	 *
	 * @return  Default height.
	 */
	public static double defaultHeight()
		{
		return theDefaultHeight;
		}

	/**
	 * Set the default height for ovals.
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
	 * Set this oval's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This oval.
	 */
	public Oval outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Set this oval's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This oval.
	 */
	public Oval fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Set the size of this oval's bounding box.
	 *
	 * @param  theSize  Size.
	 *
	 * @return  This oval.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	public Oval size
		(Size theSize)
		{
		doSize (theSize);
		return this;
		}

	/**
	 * Set the width of this oval's bounding box.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This oval.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	public Oval width
		(double theWidth)
		{
		doWidth (theWidth);
		return this;
		}

	/**
	 * Set the height of this oval's bounding box.
	 *
	 * @param  theHeight  Height.
	 *
	 * @return  This oval.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	public Oval height
		(double theHeight)
		{
		doHeight (theHeight);
		return this;
		}

	/**
	 * Set the northwest corner point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This oval.
	 */
	public Oval nw
		(double x,
		 double y)
		{
		doNw (x, y);
		return this;
		}

	/**
	 * Set the northwest corner point of this oval's bounding box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the north middle point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This oval.
	 */
	public Oval n
		(double x,
		 double y)
		{
		doN (x, y);
		return this;
		}

	/**
	 * Set the north middle point of this oval's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the northeast corner point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This oval.
	 */
	public Oval ne
		(double x,
		 double y)
		{
		doNe (x, y);
		return this;
		}

	/**
	 * Set the northeast corner point of this oval's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the west middle point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This oval.
	 */
	public Oval w
		(double x,
		 double y)
		{
		doW (x, y);
		return this;
		}

	/**
	 * Set the west middle point of this oval's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the center point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This oval.
	 */
	public Oval c
		(double x,
		 double y)
		{
		doC (x, y);
		return this;
		}

	/**
	 * Set the center point of this oval's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the east middle point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This oval.
	 */
	public Oval e
		(double x,
		 double y)
		{
		doE (x, y);
		return this;
		}

	/**
	 * Set the east middle point of this oval's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southwest corner point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This oval.
	 */
	public Oval sw
		(double x,
		 double y)
		{
		doSw (x, y);
		return this;
		}

	/**
	 * Set the southwest corner point of this oval's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the south middle point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This oval.
	 */
	public Oval s
		(double x,
		 double y)
		{
		doS (x, y);
		return this;
		}

	/**
	 * Set the south middle point of this oval's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southeast corner point of this oval's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This oval.
	 */
	public Oval se
		(double x,
		 double y)
		{
		doSe (x, y);
		return this;
		}

	/**
	 * Set the southeast corner point of this oval's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public Oval se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Add this oval to the end of the default drawing's sequence of drawing
	 * items.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Oval add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this oval to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Oval add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this oval to the beginning of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public Oval addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this oval to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This oval.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public Oval addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
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
		double d = Math.min (width, height);
		return new RoundRectangle2D.Double (nw.x, nw.y, width, height, d, d);
		}

	}
