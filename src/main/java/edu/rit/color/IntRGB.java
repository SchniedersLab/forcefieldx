//******************************************************************************
//
// File:    IntRGB.java
// Package: edu.rit.color
// Unit:    Class edu.rit.color.IntRGB
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
 * Class IntRGB provides a color represented as integer red, green, and blue
 * components.
 * <P>
 * Class IntRGB includes methods for packing and unpacking an integer RGB color
 * object into and from an integer. The packed color representation uses 8 bits
 * for each component, with the red component in bits 23-16, the green component
 * in bits 15-8, and the blue component in bits 7-0.
 *
 * @author  Alan Kaminsky
 * @version 12-Dec-2005
 */
public class IntRGB
	{

// Exported data members.

	/**
	 * The red component in the range 0 through 255.
	 */
	public int red;

	/**
	 * The green component in the range 0 through 255.
	 */
	public int green;

	/**
	 * The blue component in the range 0 through 255.
	 */
	public int blue;

// Exported constructors.

	/**
	 * Construct a new integer RGB color. The red, green, and blue components
	 * are all 0.
	 */
	public IntRGB()
		{
		}

	/**
	 * Construct a new integer RGB color with the given red, green, and blue
	 * components.
	 *
	 * @param  red    Red component.
	 * @param  green  Green component.
	 * @param  blue   Blue component.
	 */
	public IntRGB
		(int red,
		 int green,
		 int blue)
		{
		this.red = red;
		this.green = green;
		this.blue = blue;
		}

	/**
	 * Construct a new integer RGB color from the given packed color.
	 *
	 * @param  color  Packed color.
	 */
	public IntRGB
		(int color)
		{
		this.unpack (color);
		}

// Exported operations.

	/**
	 * Pack this integer RGB color into a packed color.
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
		(int red,
		 int green,
		 int blue)
		{
		if (red < 0) red = 0;
		else if (red > 255) red = 255;
		if (green < 0) green = 0;
		else if (green > 255) green = 255;
		if (blue < 0) blue = 0;
		else if (blue > 255) blue = 255;
		return (red << 16) | (green << 8) | blue;
		}

	/**
	 * Unpack this integer RGB color from a packed color.
	 *
	 * @param  color  Packed color.
	 */
	public void unpack
		(int color)
		{
		this.red   = (color >> 16) & 0xFF;
		this.green = (color >>  8) & 0xFF;
		this.blue  = (color      ) & 0xFF;
		}

	}
