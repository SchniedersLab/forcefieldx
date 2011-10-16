//******************************************************************************
//
// File:    DottedOutline.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.DottedOutline
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

import java.awt.BasicStroke;
import java.awt.Graphics2D;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class DottedOutline provides an object that outlines an area in a {@linkplain
 * DrawingItem} with a round dotted stroke in a solid color.
 *
 * @author  Alan Kaminsky
 * @version 12-Jul-2006
 */
public class DottedOutline
	implements Outline
	{

// Exported constants.

	/**
	 * The normal dotted outline width (1).
	 */
	public static final float NORMAL_WIDTH = 1.0f;

	/**
	 * The normal dotted outline fill paint (black).
	 */
	public static final Fill NORMAL_FILL = ColorFill.BLACK;

	/**
	 * The normal dotted outline (width = 1, fill paint = black).
	 */
	public static final DottedOutline NORMAL_OUTLINE = new DottedOutline();

// Hidden data members.

	private static final long serialVersionUID = -5814995045228731216L;

	private float myWidth;
	private Fill myFill;
	private transient BasicStroke myStroke;

// Exported constructors.

	/**
	 * Construct a new dotted outline object with the normal width (1) and the
	 * normal fill paint (black).
	 */
	public DottedOutline()
		{
		myWidth = NORMAL_WIDTH;
		myFill = NORMAL_FILL;
		computeStroke();
		}

	/**
	 * Construct a new dotted outline object with the same width and fill paint
	 * as the given dotted outline object.
	 *
	 * @param  theOutline  Dotted outline object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theOutline</TT> is null.
	 */
	public DottedOutline
		(DottedOutline theOutline)
		{
		myWidth = theOutline.myWidth;
		myFill = theOutline.myFill;
		computeStroke();
		}

// Exported operations.

	/**
	 * Returns this dotted outline object's width.
	 *
	 * @return  Width.
	 */
	public float width()
		{
		return myWidth;
		}

	/**
	 * Set this dotted outline object's width.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This dotted outline object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than or
	 *     equal to 0.
	 */
	public DottedOutline width
		(float theWidth)
		{
		if (theWidth <= 0.0)
			{
			throw new IllegalArgumentException();
			}
		myWidth = theWidth;
		computeStroke();
		return this;
		}

	/**
	 * Returns this dotted outline object's fill paint.
	 *
	 * @return  Fill paint.
	 */
	public Fill fill()
		{
		return myFill;
		}

	/**
	 * Set this dotted outline object's fill paint.
	 *
	 * @param  theFill  Fill paint.
	 *
	 * @return  This dotted outline object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFill</TT> is null.
	 */
	public DottedOutline fill
		(Fill theFill)
		{
		if (theFill == null)
			{
			throw new NullPointerException();
			}
		myFill = theFill;
		computeStroke();
		return this;
		}

	/**
	 * Returns the stroke width of this outline.
	 */
	public float getStrokeWidth()
		{
		return myWidth;
		}

	/**
	 * Set the given graphics context's stroke and paint attributes as specified
	 * by this outline object.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void setGraphicsContext
		(Graphics2D g2d)
		{
		g2d.setStroke (myStroke);
		myFill.setGraphicsContext (g2d);
		}

	/**
	 * Write this dotted outline object to the given object output stream.
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
		out.writeFloat (myWidth);
		out.writeObject (myFill);
		}

	/**
	 * Read this dotted outline object from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this dotted outline object
	 *     cannot be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		myWidth = in.readFloat();
		myFill = (Fill) in.readObject();
		computeStroke();
		}

// Hidden operations.

	private void computeStroke()
		{
		myStroke = 
			new BasicStroke
				(myWidth,
				 BasicStroke.CAP_ROUND,
				 BasicStroke.JOIN_MITER,
				 10.0f,
				 new float[] {0.0f, 1.5f * myWidth},
				 0);
		}

	}
