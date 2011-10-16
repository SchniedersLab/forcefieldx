//******************************************************************************
//
// File:    SquareDottedOutline.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.SquareDottedOutline
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
 * Class SquareDottedOutline provides an object that outlines an area in a
 * {@linkplain DrawingItem} with a square-cornered dotted stroke in a solid
 * color.
 *
 * @author  Alan Kaminsky
 * @version 12-Jul-2006
 */
public class SquareDottedOutline
	implements Outline
	{

// Exported constants.

	/**
	 * The normal square dotted outline width (1).
	 */
	public static final float NORMAL_WIDTH = 1.0f;

	/**
	 * The normal square dotted outline fill paint (black).
	 */
	public static final Fill NORMAL_FILL = ColorFill.BLACK;

	/**
	 * The normal square dotted outline (width = 1, fill paint = black).
	 */
	public static final SquareDottedOutline NORMAL_OUTLINE =
		new SquareDottedOutline();

// Hidden data members.

	private static final long serialVersionUID = 7878713408282698254L;

	private float myWidth;
	private Fill myFill;
	private transient BasicStroke myStroke;

// Exported constructors.

	/**
	 * Construct a new square dotted outline object with the normal width (1)
	 * and the normal fill paint (black).
	 */
	public SquareDottedOutline()
		{
		myWidth = NORMAL_WIDTH;
		myFill = NORMAL_FILL;
		computeStroke();
		}

	/**
	 * Construct a new square dotted outline object with the same width and fill
	 * paint as the given dotted outline object.
	 *
	 * @param  theOutline  Square dotted outline object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theOutline</TT> is null.
	 */
	public SquareDottedOutline
		(SquareDottedOutline theOutline)
		{
		myWidth = theOutline.myWidth;
		myFill = theOutline.myFill;
		computeStroke();
		}

// Exported operations.

	/**
	 * Returns this square dotted outline object's width.
	 *
	 * @return  Width.
	 */
	public float width()
		{
		return myWidth;
		}

	/**
	 * Set this square dotted outline object's width.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This square dotted outline object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than or
	 *     equal to 0.
	 */
	public SquareDottedOutline width
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
	 * Returns this square dotted outline object's fill paint.
	 *
	 * @return  Fill paint.
	 */
	public Fill fill()
		{
		return myFill;
		}

	/**
	 * Set this square dotted outline object's fill paint.
	 *
	 * @param  theFill  Fill paint.
	 *
	 * @return  This square dotted outline object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFill</TT> is null.
	 */
	public SquareDottedOutline fill
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
	 * Write this square dotted outline object to the given object output
	 * stream.
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
	 * Read this square dotted outline object from the given object input
	 * stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this square dotted outline
	 *     object cannot be found.
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
				 BasicStroke.CAP_SQUARE,
				 BasicStroke.JOIN_MITER,
				 10.0f,
				 new float[] {0.0f, 2.0f * myWidth},
				 0);
		}

	}
