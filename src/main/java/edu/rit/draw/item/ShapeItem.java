//******************************************************************************
//
// File:    ShapeItem.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.ShapeItem
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

/**
 * Class ShapeItem is the abstract base class for a {@linkplain DrawingItem}
 * that consists of a single 2-D graphics shape with an outline and an interior.
 *
 * @author  Alan Kaminsky
 * @version 23-Jun-2008
 */
public abstract class ShapeItem
	extends FilledItem
	{

// Hidden data members.

	private static final long serialVersionUID = -6592516179733916611L;

// Exported constructors.

	/**
	 * Construct a new shape item.
	 */
	public ShapeItem()
		{
		super();
		}

	/**
	 * Construct a new shape item with the same outline and fill paint as the
	 * given shape item.
	 *
	 * @param  theItem  Shape item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theItem</TT> is null.
	 */
	public ShapeItem
		(ShapeItem theItem)
		{
		super (theItem);
		}

// Exported operations.

	/**
	 * Set this shape item's outline.
	 *
	 * @param  theOutline  Outline, or <TT>Outline.NONE</TT>.
	 *
	 * @return  This shape item.
	 */
	public ShapeItem outline
		(Outline theOutline)
		{
		doOutline (theOutline);
		return this;
		}

	/**
	 * Set this shape item's fill paint.
	 *
	 * @param  theFill  Fill paint, or <TT>Fill.NONE</TT>.
	 *
	 * @return  This shape item.
	 */
	public ShapeItem fill
		(Fill theFill)
		{
		doFill (theFill);
		return this;
		}

	/**
	 * Add this shape item to the end of the default drawing's sequence of
	 * drawing items.
	 *
	 * @return  This shape item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public ShapeItem add()
		{
		doAdd (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this shape item to the end of the given drawing's sequence of drawing
	 * items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This shape item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public ShapeItem add
		(Drawing theDrawing)
		{
		doAdd (theDrawing);
		return this;
		}

	/**
	 * Add this shape item to the beginning of the default drawing's sequence
	 * of drawing items.
	 *
	 * @return  This shape item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if there is no default drawing.
	 *
	 * @see  Drawing#defaultDrawing()
	 */
	public ShapeItem addFirst()
		{
		doAddFirst (Drawing.defaultDrawing());
		return this;
		}

	/**
	 * Add this shape item to the beginning of the given drawing's sequence of
	 * drawing items.
	 *
	 * @param  theDrawing  Drawing.
	 *
	 * @return  This shape item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public ShapeItem addFirst
		(Drawing theDrawing)
		{
		doAddFirst (theDrawing);
		return this;
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

		Shape shape = getShape();
		if (myFill != null)
			{
			myFill.setGraphicsContext (g2d);
			g2d.fill (shape);
			}
		if (myOutline != null)
			{
			myOutline.setGraphicsContext (g2d);
			g2d.draw (shape);
			}
		}

// Hidden operations.

	/**
	 * Determine the 2-D graphics shape that this shape object will draw.
	 *
	 * @return  Shape.
	 */
	abstract Shape getShape();

	}
