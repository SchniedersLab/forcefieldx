//******************************************************************************
//
// File:    RectangularItem.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.RectangularItem
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

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class RectangularItem is the abstract base class for a {@linkplain
 * DrawingItem} that has an outline, is filled with a paint, and occupies a
 * rectangular area of a certain size at a certain location.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public abstract class RectangularItem
	extends ShapeItem
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = -6695492802147352641L;

	// Coordinates of most recently specified corner.
	double x;
	double y;

	// Factors for going from specified corner to northwest corner.
	double xFactor;
	double yFactor;

	// Size.
	double width;
	double height;

// Exported constructors.

	/**
	 * Construct a new rectangular item. The rectangular item's northwest corner
	 * is located at (0,0). The rectangular item's size is the default size
	 * (determined by the subclass).
	 */
	public RectangularItem()
		{
		super();
		}

	/**
	 * Construct a new rectangular item with the same outline, fill paint,
	 * location, and size as the given rectangular item.
	 *
	 * @param  theItem  Rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public RectangularItem
		(RectangularItem theItem)
		{
		super (theItem);
		this.x = theItem.x;
		this.y = theItem.y;
		this.xFactor = theItem.xFactor;
		this.yFactor = theItem.yFactor;
		this.width = theItem.width;
		this.height = theItem.height;
		}

// Exported operations.

	/**
	 * Set this rectangular item's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Set this rectangular item's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Set the size of this rectangular item's bounding box.
	 *
	 * @param  theSize  Size.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	public RectangularItem size
		(Size theSize)
		{
		doSize (theSize);
		return this;
		}

	/**
	 * Returns the width of this rectangular item's bounding box.
	 *
	 * @return  Width.
	 */
	public double width()
		{
		return this.width;
		}

	/**
	 * Set the width of this rectangular item's bounding box.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	public RectangularItem width
		(double theWidth)
		{
		doWidth (theWidth);
		return this;
		}

	/**
	 * Returns the height of this rectangular item's bounding box.
	 *
	 * @return  Height.
	 */
	public double height()
		{
		return this.height;
		}

	/**
	 * Set the height of this rectangular item's bounding box.
	 *
	 * @param  theHeight  Height.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	public RectangularItem height
		(double theHeight)
		{
		doHeight (theHeight);
		return this;
		}

	/**
	 * Returns the northwest corner point of this rectangular item's bounding
	 * box.
	 *
	 * @return  Northwest corner point.
	 */
	public Point nw()
		{
		return new Point (x + xFactor * width, y + yFactor * height);
		}

	/**
	 * Set the northwest corner point of this rectangular item's bounding
	 * box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem nw
		(double x,
		 double y)
		{
		doNw (x, y);
		return this;
		}

	/**
	 * Set the northwest corner point of this rectangular item's bounding
	 * box.
	 *
	 * @param  thePoint  Northwest corner point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem nw
		(Point thePoint)
		{
		doNw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the north middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem n
		(double x,
		 double y)
		{
		doN (x, y);
		return this;
		}

	/**
	 * Set the north middle point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  North middle point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem n
		(Point thePoint)
		{
		doN (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the northeast corner point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem ne
		(double x,
		 double y)
		{
		doNe (x, y);
		return this;
		}

	/**
	 * Set the northeast corner point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  Northeast corner point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem ne
		(Point thePoint)
		{
		doNe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the west middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem w
		(double x,
		 double y)
		{
		doW (x, y);
		return this;
		}

	/**
	 * Set the west middle point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  West middle point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem w
		(Point thePoint)
		{
		doW (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the center point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem c
		(double x,
		 double y)
		{
		doC (x, y);
		return this;
		}

	/**
	 * Set the center point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  Center point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem c
		(Point thePoint)
		{
		doC (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the east middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem e
		(double x,
		 double y)
		{
		doE (x, y);
		return this;
		}

	/**
	 * Set the east middle point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  East middle point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem e
		(Point thePoint)
		{
		doE (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southwest corner point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem sw
		(double x,
		 double y)
		{
		doSw (x, y);
		return this;
		}

	/**
	 * Set the southwest corner point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  Southwest corner point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem sw
		(Point thePoint)
		{
		doSw (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the south middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem s
		(double x,
		 double y)
		{
		doS (x, y);
		return this;
		}

	/**
	 * Set the south middle point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  South middle point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem s
		(Point thePoint)
		{
		doS (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Set the southeast corner point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 *
	 * @return  This rectangular item.
	 */
	public RectangularItem se
		(double x,
		 double y)
		{
		doSe (x, y);
		return this;
		}

	/**
	 * Set the southeast corner point of this rectangular item's bounding box.
	 *
	 * @param  thePoint  Southeast corner point.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePoint</TT> is null.
	 */
	public RectangularItem se
		(Point thePoint)
		{
		doSe (thePoint.x, thePoint.y);
		return this;
		}

	/**
	 * Add this rectangular item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public RectangularItem add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this rectangular item to the end of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public RectangularItem add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this rectangular item to the beginning of the default drawing's
	 * sequence of drawing items.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public RectangularItem addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this rectangular item to the beginning of the given drawing's
	 * sequence of drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This rectangular item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public RectangularItem addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this rectangular item to the given object output stream.
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
		out.writeDouble (x);
		out.writeDouble (y);
		out.writeDouble (xFactor);
		out.writeDouble (yFactor);
		out.writeDouble (width);
		out.writeDouble (height);
		}

	/**
	 * Read this rectangular item from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this rectangular item
	 *     cannot be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		super.readExternal (in);
		x = in.readDouble();
		y = in.readDouble();
		xFactor = in.readDouble();
		yFactor = in.readDouble();
		width = in.readDouble();
		height = in.readDouble();
		}

// Hidden operations.

	/**
	 * Set the size of this rectangular item's bounding box.
	 *
	 * @param  theSize  Size.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the width or the height of
	 *     <TT>theSize</TT> is less than 0.
	 */
	void doSize
		(Size theSize)
		{
		if (theSize.width < 0.0 || theSize.height < 0.0)
			{
			throw new IllegalArgumentException();
			}
		this.width = theSize.width;
		this.height = theSize.height;
		}

	/**
	 * Set the width of this rectangular item's bounding box.
	 *
	 * @param  theWidth  Width.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than 0.
	 */
	void doWidth
		(double theWidth)
		{
		if (theWidth < 0.0) throw new IllegalArgumentException();
		this.width = theWidth;
		}

	/**
	 * Set the height of this rectangular item's bounding box.
	 *
	 * @param  theHeight  Height.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> is less than 0.
	 */
	void doHeight
		(double theHeight)
		{
		if (theHeight < 0.0) throw new IllegalArgumentException();
		this.height = theHeight;
		}

	/**
	 * Set the northwest corner point of this rectangular item's bounding
	 * box.
	 *
	 * @param  x  X coordinate of northwest corner point.
	 * @param  y  Y coordinate of northwest corner point.
	 */
	void doNw
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = 0.0;
		}

	/**
	 * Set the north middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of north middle point.
	 * @param  y  Y coordinate of north middle point.
	 */
	void doN
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = 0.0;
		}

	/**
	 * Set the northeast corner point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of northeast corner point.
	 * @param  y  Y coordinate of northeast corner point.
	 */
	void doNe
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = 0.0;
		}

	/**
	 * Set the west middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of west middle point.
	 * @param  y  Y coordinate of west middle point.
	 */
	void doW
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = -0.5;
		}

	/**
	 * Set the center point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of center point.
	 * @param  y  Y coordinate of center point.
	 */
	void doC
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = -0.5;
		}

	/**
	 * Set the east middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of east middle point.
	 * @param  y  Y coordinate of east middle point.
	 */
	void doE
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = -0.5;
		}

	/**
	 * Set the southwest corner point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of southwest corner point.
	 * @param  y  Y coordinate of southwest corner point.
	 */
	void doSw
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = 0.0;
		this.yFactor = -1.0;
		}

	/**
	 * Set the south middle point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of south middle point.
	 * @param  y  Y coordinate of south middle point.
	 */
	void doS
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -0.5;
		this.yFactor = -1.0;
		}

	/**
	 * Set the southeast corner point of this rectangular item's bounding box.
	 *
	 * @param  x  X coordinate of southeast corner point.
	 * @param  y  Y coordinate of southeast corner point.
	 */
	void doSe
		(double x,
		 double y)
		{
		this.x = x;
		this.y = y;
		this.xFactor = -1.0;
		this.yFactor = -1.0;
		}

	}
