//******************************************************************************
//
// File:    RGB.java
// Package: edu.rit.color
// Unit:    Class edu.rit.color.RGB
//
// This Java source file is copyright (C) 2005 by Alan Kaminsky. All rights
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

package edu.rit.color;

/**
 * Class RGB provides a color represented as floating point red, green, and blue
 * components.
 * <P>
 * Class RGB includes methods for packing and unpacking a floating point RGB
 * color object into and from an integer. The packed color representation uses 8
 * bits for each component, with the red component in bits 23-16, the green
 * component in bits 15-8, and the blue component in bits 7-0.
 *
 * @author  Alan Kaminsky
 * @version 12-Dec-2005
 */
public class RGB
	{

// Exported data members.

	/**
	 * The red component in the range 0.0 through 1.0.
	 */
	public float red;

	/**
	 * The green component in the range 0.0 through 1.0.
	 */
	public float green;

	/**
	 * The blue component in the range 0.0 through 1.0.
	 */
	public float blue;

// Exported constructors.

	/**
	 * Construct a new floating point RGB color. The red, green, and blue
	 * components are all 0.
	 */
	public RGB()
		{
		}

	/**
	 * Construct a new floating point RGB color with the given red, green, and
	 * blue components.
	 *
	 * @param  red    Red component.
	 * @param  green  Green component.
	 * @param  blue   Blue component.
	 */
	public RGB
		(float red,
		 float green,
		 float blue)
		{
		this.red = red;
		this.green = green;
		this.blue = blue;
		}

	/**
	 * Construct a new floating point RGB color from the given packed color.
	 *
	 * @param  color  Packed color.
	 */
	public RGB
		(int color)
		{
		this.unpack (color);
		}

// Exported operations.

	/**
	 * Pack this floating point RGB color into a packed color.
	 *
	 * @return  Packed color.
	 */
	public int pack()
		{
		return pack (this.red, this.green, this.blue);
		}

	/**
	 * Pack the given red, green, and blue components into a packed color.
	 *
	 * @param  red    Red component.
	 * @param  green  Green component.
	 * @param  blue   Blue component.
	 *
	 * @return  Packed color.
	 */
	public static int pack
		(float red,
		 float green,
		 float blue)
		{
		int r = (int) (red * 256.0f);
		if (r < 0) r = 0;
		else if (r > 255) r = 255;
		int g = (int) (green * 256.0f);
		if (g < 0) g = 0;
		else if (g > 255) g = 255;
		int b = (int) (blue * 256.0f);
		if (b < 0) b = 0;
		else if (b > 255) b = 255;
		return (r << 16) | (g << 8) | b;
		}

	/**
	 * Unpack this floating point RGB color from a packed color.
	 *
	 * @param  color  Packed color.
	 */
	public void unpack
		(int color)
		{
		this.red   = ((color >> 16) & 0xFF) / 256.0f;
		this.green = ((color >>  8) & 0xFF) / 256.0f;
		this.blue  = ((color      ) & 0xFF) / 256.0f;
		}

	}
