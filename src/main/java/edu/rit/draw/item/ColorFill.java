//******************************************************************************
//
// File:    ColorFill.java
// Package: edu.rit.draw.item
// Unit:    Class edu.rit.draw.item.ColorFill
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

import java.awt.Color;
import java.awt.Graphics2D;

import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class ColorFill provides an object that fills an area in a {@linkplain
 * DrawingItem} with a solid color.
 *
 * @author  Alan Kaminsky
 * @version 10-Jul-2006
 */
public class ColorFill
	implements Fill
	{

// Exported constants.

	/**
	 * Color fill object for the color white.
	 */
	public static final ColorFill WHITE = new ColorFill();

	/**
	 * Color fill object for the color light gray.
	 */
	public static final ColorFill LIGHT_GRAY =
		new ColorFill().color (Color.LIGHT_GRAY);

	/**
	 * Color fill object for the color gray.
	 */
	public static final ColorFill GRAY =
		new ColorFill().color (Color.GRAY);

	/**
	 * Color fill object for the color dark gray.
	 */
	public static final ColorFill DARK_GRAY =
		new ColorFill().color (Color.DARK_GRAY);

	/**
	 * Color fill object for the color black.
	 */
	public static final ColorFill BLACK =
		new ColorFill().color (Color.BLACK);

	/**
	 * Color fill object for the color red.
	 */
	public static final ColorFill RED =
		new ColorFill().color (Color.RED);

	/**
	 * Color fill object for the color pink.
	 */
	public static final ColorFill PINK =
		new ColorFill().color (Color.PINK);

	/**
	 * Color fill object for the color orange.
	 */
	public static final ColorFill ORANGE =
		new ColorFill().color (Color.ORANGE);

	/**
	 * Color fill object for the color yellow.
	 */
	public static final ColorFill YELLOW =
		new ColorFill().color (Color.YELLOW);

	/**
	 * Color fill object for the color green.
	 */
	public static final ColorFill GREEN =
		new ColorFill().color (Color.GREEN);

	/**
	 * Color fill object for the color magenta.
	 */
	public static final ColorFill MAGENTA =
		new ColorFill().color (Color.MAGENTA);

	/**
	 * Color fill object for the color cyan.
	 */
	public static final ColorFill CYAN =
		new ColorFill().color (Color.CYAN);

	/**
	 * Color fill object for the color blue.
	 */
	public static final ColorFill BLUE =
		new ColorFill().color (Color.BLUE);

	/**
	 * The normal color fill object (white).
	 */
	public static final ColorFill NORMAL_FILL = WHITE;

// Hidden data members.

	private static final long serialVersionUID = 1702092537732335278L;

	private transient Color myColor;

// Exported constructors.

	/**
	 * Construct a new color fill object with the normal fill color (white).
	 */
	public ColorFill()
		{
		myColor = Color.white;
		}

	/**
	 * Construct a new color fill object with the same color as the given color
	 * fill object.
	 *
	 * @param  theColorFill  Color fill object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theColorFill</TT> is null.
	 */
	public ColorFill
		(ColorFill theColorFill)
		{
		this.myColor = theColorFill.myColor;
		}

// Exported operations.

	/**
	 * Returns the color of this color fill object.
	 *
	 * @return  Color.
	 */
	public Color color()
		{
		return myColor;
		}

	/**
	 * Set this color fill object to the given shade of gray.
	 *
	 * @param  grayLevel  Gray level in the range 0 (black) through 255 (white).
	 *
	 * @return  This color fill object.
	 */
	public ColorFill gray
		(int grayLevel)
		{
		myColor = new Color (grayLevel, grayLevel, grayLevel);
		return this;
		}

	/**
	 * Set this color fill object to the given shade of gray.
	 *
	 * @param  grayLevel  Gray level in the range 0.0f (black) through 1.0f
	 *                    (white).
	 *
	 * @return  This color fill object.
	 */
	public ColorFill gray
		(float grayLevel)
		{
		myColor = new Color (grayLevel, grayLevel, grayLevel);
		return this;
		}

	/**
	 * Set this color fill object to the color with the given red, green, and
	 * blue components.
	 *
	 * @param  red    Red component in the range 0 through 255.
	 * @param  green  Green component in the range 0 through 255.
	 * @param  blue   Blue component in the range 0 through 255.
	 *
	 * @return  This color fill object.
	 */
	public ColorFill rgb
		(int red,
		 int green,
		 int blue)
		{
		myColor = new Color (red, green, blue);
		return this;
		}

	/**
	 * Set this color fill object to the color with the given red, green, and
	 * blue components.
	 *
	 * @param  red    Red component in the range 0.0f through 1.0f.
	 * @param  green  Green component in the range 0.0f through 1.0f.
	 * @param  blue   Blue component in the range 0.0f through 1.0f.
	 *
	 * @return  This color fill object.
	 */
	public ColorFill rgb
		(float red,
		 float green,
		 float blue)
		{
		myColor = new Color (red, green, blue);
		return this;
		}

	/**
	 * Set this color fill object to the given color.
	 *
	 * @param  theColor  Color.
	 *
	 * @return  This color fill object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theColor</TT> is null.
	 */
	public ColorFill color
		(Color theColor)
		{
		if (theColor == null)
			{
			throw new NullPointerException();
			}
		myColor = theColor;
		return this;
		}

	/**
	 * Set the given graphics context's paint attribute as specified by this
	 * fill object.
	 *
	 * @param  g2d  2-D graphics context.
	 */
	public void setGraphicsContext
		(Graphics2D g2d)
		{
		g2d.setPaint (myColor);
		}

	/**
	 * Write this color fill object to the given object output stream.
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
		out.writeInt (myColor.getRGB());
		}

	/**
	 * Read this color fill object from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		myColor = new Color (in.readInt(), true);
		}

	}
