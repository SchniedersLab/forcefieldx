//******************************************************************************
//
// File:    DrawingItem.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.DrawingItem
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

import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class DrawingItem is the abstract base class for each item on a {@linkplain
 * edu.rit.draw.Drawing}. Subclasses provide different kinds of drawing items,
 * such as text, lines, and shapes.
 * <P>
 * Each drawing item is contained within a rectangular <B>bounding box</B>.
 * Class DrawingItem's methods return the location and dimensions of the
 * bounding box.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public abstract class DrawingItem
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = -7001285340255383829L;

// Exported constructors.

	/**
	 * Construct a new drawing item.
	 */
	public DrawingItem()
		{
		}

	/**
	 * Construct a new drawing item that is the same as the given drawing item.
	 *
	 * @param  theItem  Drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public DrawingItem
		(DrawingItem theItem)
		{
		}

// Exported operations.

	/**
	 * Returns the size of this drawing item's bounding box.
	 *
	 * @return  Size.
	 */
	public Size size()
		{
		return new Size (width(), height());
		}

	/**
	 * Returns the width of this drawing item's bounding box.
	 *
	 * @return  Width.
	 */
	public abstract double width();

	/**
	 * Returns the height of this drawing item's bounding box.
	 *
	 * @return  Height.
	 */
	public abstract double height();

	/**
	 * Returns the northwest corner point of this drawing item's bounding box.
	 *
	 * @return  Northwest corner point.
	 */
	public abstract Point nw();

	/**
	 * Returns the north middle point of this drawing item's bounding box.
	 *
	 * @return  North middle point.
	 */
	public Point n()
		{
		return nw().add (0.5 * width(), 0.0);
		}

	/**
	 * Returns the northeast corner point of this drawing item's bounding box.
	 *
	 * @return  Northeast corner point.
	 */
	public Point ne()
		{
		return nw().add (width(), 0.0);
		}

	/**
	 * Returns the west middle point of this drawing item's bounding box.
	 *
	 * @return  West middle point.
	 */
	public Point w()
		{
		return nw().add (0.0, 0.5 * height());
		}

	/**
	 * Returns the center point of this drawing item's bounding box.
	 *
	 * @return  Center point.
	 */
	public Point c()
		{
		return nw().add (0.5 * width(), 0.5 * height());
		}

	/**
	 * Returns the east middle point of this drawing item's bounding box.
	 *
	 * @return  East middle point.
	 */
	public Point e()
		{
		return nw().add (width(), 0.5 * height());
		}

	/**
	 * Returns the southwest corner point of this drawing item's bounding box.
	 *
	 * @return  Southwest corner point.
	 */
	public Point sw()
		{
		return nw().add (0.0, height());
		}

	/**
	 * Returns the south middle point of this drawing item's bounding box.
	 *
	 * @return  South middle point.
	 */
	public Point s()
		{
		return nw().add (0.5 * width(), height());
		}

	/**
	 * Returns the southeast corner point of this drawing item's bounding box.
	 *
	 * @return  Southeast corner point.
	 */
	public Point se()
		{
		return nw().add (width(), height());
		}

	/**
	 * Add this drawing item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public DrawingItem add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this drawing item to the end of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public DrawingItem add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this drawing item to the beginning of the default drawing's sequence
	 * of drawing items.
	 *
	 * @return  This drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public DrawingItem addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this drawing item to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This drawing item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public DrawingItem addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
		}

	/**
	 * Write this drawing item to the given object output stream.
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
		}

	/**
	 * Read this drawing item from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this drawing item cannot be
	 *     found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
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
		}

	/**
	 * Returns the rectangular region this drawing item occupies.
	 *
	 * @return  Bounding box.
	 */
	public Rectangle2D boundingBox()
		{
		Point nw = nw();
		return new Rectangle2D.Double (nw.x, nw.y, width(), height());
		}

// Hidden operations.

	/**
	 * Add this drawing item to the end of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	void doAdd
		(Drawing theDrawing)
		{
		theDrawing.add (this);
		}

	/**
	 * Add this drawing item to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	void doAddFirst
		(Drawing theDrawing)
		{
		theDrawing.addFirst (this);
		}

	}
