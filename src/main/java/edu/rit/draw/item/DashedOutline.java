//******************************************************************************
//
// File:    DashedOutline.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.DashedOutline
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
 * Class DashedOutline provides an object that outlines an area in a {@linkplain
 * DrawingItem} with a square-cornered dashed stroke in a solid color.
 *
 * @author  Alan Kaminsky
 * @version 12-Jul-2006
 */
public class DashedOutline
	implements Outline
	{

// Exported constants.

	/**
	 * The normal dashed outline width (1).
	 */
	public static final float NORMAL_WIDTH = 1.0f;

	/**
	 * The normal dashed outline fill paint (black).
	 */
	public static final Fill NORMAL_FILL = ColorFill.BLACK;

	/**
	 * The normal dashed outline down factor (10).
	 */
	public static final float NORMAL_DOWN_FACTOR = 10.0f;

	/**
	 * The normal dashed outline up factor (6).
	 */
	public static final float NORMAL_UP_FACTOR = 6.0f;

	/**
	 * The normal dashed outline (width = 1, fill paint = black, down factor =
	 * 5, up factor = 3).
	 */
	public static final DashedOutline NORMAL_OUTLINE = new DashedOutline();

// Hidden data members.

	private static final long serialVersionUID = 8450835114943482027L;

	private float myWidth;
	private Fill myFill;
	private float myDownFactor;
	private float myUpFactor;
	private transient BasicStroke myStroke;

// Exported constructors.

	/**
	 * Construct a new dashed outline object with the normal width (1), the
	 * normal fill paint (black), the normal down factor (10), and the normal up
	 * factor (6).
	 */
	public DashedOutline()
		{
		myWidth = NORMAL_WIDTH;
		myFill = NORMAL_FILL;
		myDownFactor = NORMAL_DOWN_FACTOR;
		myUpFactor = NORMAL_UP_FACTOR;
		computeStroke();
		}

	/**
	 * Construct a new dashed outline object with the same width, fill paint,
	 * down factor, and up factor as the given dashed outline object.
	 *
	 * @param  theOutline  Dashed outline object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theOutline</TT> is null.
	 */
	public DashedOutline
		(DashedOutline theOutline)
		{
		myWidth = theOutline.myWidth;
		myFill = theOutline.myFill;
		myDownFactor = theOutline.myDownFactor;
		myUpFactor = theOutline.myUpFactor;
		computeStroke();
		}

// Exported operations.

	/**
	 * Returns this dashed outline object's width.
	 *
	 * @return  Width.
	 */
	public float width()
		{
		return myWidth;
		}

	/**
	 * Set this dashed outline object's width.
	 *
	 * @param  theWidth  Width.
	 *
	 * @return  This dashed outline object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theWidth</TT> is less than or
	 *     equal to 0.
	 */
	public DashedOutline width
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
	 * Returns this dashed outline object's fill paint.
	 *
	 * @return  Fill paint.
	 */
	public Fill fill()
		{
		return myFill;
		}

	/**
	 * Set this dashed outline object's fill paint.
	 *
	 * @param  theFill  Fill paint.
	 *
	 * @return  This dashed outline object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFill</TT> is null.
	 */
	public DashedOutline fill
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
	 * Returns this dashed outline object's down factor. The down factor is the
	 * length, relative to this dashed outline object's width, that the "pen" is
	 * "down" drawing a dash.
	 *
	 * @return  Down factor.
	 */
	public float down()
		{
		return myDownFactor;
		}

	/**
	 * Set this dashed outline object's down factor. The down factor is the
	 * length, relative to this dashed outline object's width, that the "pen" is
	 * "down" drawing a dash.
	 *
	 * @param  theDownFactor  Down factor.
	 *
	 * @return  This dashed outline object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theDownFactor</TT> is less than
	 *     or equal to 0.
	 */
	public DashedOutline down
		(float theDownFactor)
		{
		if (theDownFactor <= 0.0f)
			{
			throw new IllegalArgumentException();
			}
		myDownFactor = theDownFactor;
		computeStroke();
		return this;
		}

	/**
	 * Returns this dashed outline object's up factor. The up factor is the
	 * length, relative to this dashed outline object's width, that the "pen" is
	 * "up" between dashes.
	 *
	 * @return  Up factor.
	 */
	public float up()
		{
		return myUpFactor;
		}

	/**
	 * Set this dashed outline object's up factor. The down factor is the
	 * length, relative to this dashed outline object's width, that the "pen" is
	 * "up" between dashes.
	 *
	 * @param  theUpFactor  Up factor.
	 *
	 * @return  This dashed outline object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theUpFactor</TT> is less than
	 *     or equal to 0.
	 */
	public DashedOutline up
		(float theUpFactor)
		{
		if (theUpFactor <= 0.0f)
			{
			throw new IllegalArgumentException();
			}
		myUpFactor = theUpFactor;
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
	 * Write this dashed outline object to the given object output stream.
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
		out.writeFloat (myDownFactor);
		out.writeFloat (myUpFactor);
		out.writeObject (myFill);
		}

	/**
	 * Read this dashed outline object from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  ClassNotFoundException
	 *     Thrown if any class needed to deserialize this dashed outline object
	 *     cannot be found.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException, ClassNotFoundException
		{
		myWidth = in.readFloat();
		myDownFactor = in.readFloat();
		myUpFactor = in.readFloat();
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
				 new float[] {myDownFactor * myWidth, myUpFactor * myWidth},
				 0);
		}

	}
