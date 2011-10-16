//******************************************************************************
//
// File:    HSB.java
// Package: edu.rit.color
// Unit:    Class edu.rit.color.HSB
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
 * Class HSB provides a color represented as floating point hue, saturation, and
 * brightness components.
 * <P>
 * The hue component gives the basic color. A hue of 0 = red; 1/6 = yellow; 2/6
 * = green; 3/6 = cyan; 4/6 = blue; 5/6 = magenta; 1 = red again. Intermediate
 * hue values yield intermediate colors.
 * <P>
 * The saturation component specifies how gray or colored the color is. A
 * saturation of 0 yields fully gray; a saturation of 1 yields fully colored.
 * Intermediate saturation values yield mixtures of gray and colored.
 * <P>
 * The brightness component specifies how dark or light the color is. A
 * brightness of 0 yields fully dark (black); a brightness of 1 yields fully
 * light (somewhere between white and colored depending on the saturation).
 * Intermediate brightness values yield somewhere between a gray shade and a
 * darkened color (depending on the saturation).
 * <P>
 * Class HSB includes methods for packing and unpacking a floating point HSB
 * color object into and from an integer. The packed color representation uses 8
 * bits for each component, with the red component in bits 23-16, the green
 * component in bits 15-8, and the blue component in bits 7-0.
 *
 * @author  Alan Kaminsky
 * @version 22-Apr-2010
 */
public class HSB
	{

// Exported data members.

	/**
	 * The hue component in the range 0.0 through 1.0.
	 */
	public float hue;

	/**
	 * The saturation component in the range 0.0 through 1.0.
	 */
	public float sat;

	/**
	 * The brightness component in the range 0.0 through 1.0.
	 */
	public float bri;

// Exported constructors.

	/**
	 * Construct a new floating point HSB color. The hue, saturation, and
	 * brightness components are all 0.
	 */
	public HSB()
		{
		}

	/**
	 * Construct a new floating point HSB color with the given hue, saturation,
	 * and brightness components.
	 *
	 * @param  hue  Hue component.
	 * @param  sat  Saturation component.
	 * @param  bri  Brightness component.
	 */
	public HSB
		(float hue,
		 float sat,
		 float bri)
		{
		this.hue = hue;
		this.sat = sat;
		this.bri = bri;
		}

	/**
	 * Construct a new floating point HSB color from the given packed color.
	 *
	 * @param  color  Packed color.
	 */
	public HSB
		(int color)
		{
		this.unpack (color);
		}

// Exported operations.

	/**
	 * Pack this floating point HSB color into a packed color.
	 *
	 * @return  Packed color.
	 */
	public int pack()
		{
		return pack (this.hue, this.sat, this.bri);
		}

	/**
	 * Pack the given hue, saturation, and brightness components into a packed
	 * color.
	 *
	 * @param  hue  Hue component.
	 * @param  sat  Saturation component.
	 * @param  bri  Brightness component.
	 *
	 * @return  Packed color.
	 */
	public static int pack
		(float hue,
		 float sat,
		 float bri)
		{
		if (hue < 0.0f) hue = 0.0f;
		else if (hue > 1.0f) hue = 1.0f;
		if (sat < 0.0f) sat = 0.0f;
		else if (sat > 1.0f) sat = 1.0f;
		if (bri < 0.0f) bri = 0.0f;
		else if (bri > 1.0f) bri = 1.0f;

		bri *= 256.0f;

		int red, green, blue;

		if (sat == 0.0f)
			{
			red = green = blue = (int)(bri);
			}
		else
			{
			hue = hue * 6.0f;
			int huecase = (int) hue;
			hue = hue - huecase;
			switch (huecase)
				{
				case 0:
				case 6:
					red   = (int)(bri);
					green = (int)(bri*(1.0f - (sat*(1.0f - hue))));
					blue  = (int)(bri*(1.0f - sat));
					// red >= green >= blue
					break;
				case 1:
					red   = (int)(bri*(1.0f - sat*hue));
					green = (int)(bri);
					blue  = (int)(bri*(1.0f - sat));
					// green >= red >= blue
					break;
				case 2:
					red   = (int)(bri*(1.0f - sat));
					green = (int)(bri);
					blue  = (int)(bri*(1.0f - (sat*(1.0f - hue))));
					// green >= blue >= red
					break;
				case 3:
					red   = (int)(bri*(1.0f - sat));
					green = (int)(bri*(1.0f - sat*hue));
					blue  = (int)(bri);
					// blue >= green >= red
					break;
				case 4:
					red   = (int)(bri*(1.0f - (sat*(1.0f - hue))));
					green = (int)(bri*(1.0f - sat));
					blue  = (int)(bri);
					// blue >= red >= green
					break;
				case 5:
					red   = (int)(bri);
					green = (int)(bri*(1.0f - sat));
					blue  = (int)(bri*(1.0f - sat*hue));
					// red >= blue >= green
					break;
				default:
					red = green = blue = (int)(bri);
					// red == green == blue
					break;
				}
			}

		if (red > 255) red = 255;
		if (green > 255) green = 255;
		if (blue > 255) blue = 255;

		return (red << 16) | (green << 8) | blue;
		}

	/**
	 * Unpack this floating point HSB color from a packed color.
	 *
	 * @param  color  Packed color.
	 */
	public void unpack
		(int color)
		{
		float red   = ((color >> 16) & 0xFF) / 256.0f;
		float green = ((color >>  8) & 0xFF) / 256.0f;
		float blue  = ((color      ) & 0xFF) / 256.0f;
		if (red == green && green == blue)
			{
			this.bri = red;
			this.sat = 0.0f;
			this.hue = 0.0f;
			}
		else if (red >= green && green >= blue)
			{
			// huecase = 0
			bri = red;
			sat = 1.0f - blue/bri;
			hue = (1.0f - (1.0f - green/bri)/sat)/6.0f;
			}
		else if (green >= red && red >= blue)
			{
			// huecase = 1
			bri = green;
			sat = 1.0f - blue/bri;
			hue = (1.0f + (1.0f - red/bri)/sat)/6.0f;
			}
		else if (green >= blue && blue >= red)
			{
			// huecase = 2
			bri = green;
			sat = 1.0f - red/bri;
			hue = (3.0f - (1.0f - blue/bri)/sat)/6.0f;
			}
		else if (blue >= green && green >= red)
			{
			// huecase = 3
			bri = blue;
			sat = 1.0f - red/bri;
			hue = (3.0f + (1.0f - green/bri)/sat)/6.0f;
			}
		else if (blue >= red && red >= green)
			{
			// huecase = 4
			bri = blue;
			sat = 1.0f - green/bri;
			hue = (5.0f - (1.0f - red/bri)/sat)/6.0f;
			}
		else // if (red >= blue && blue >= green)
			{
			// huecase = 5
			bri = red;
			sat = 1.0f - green/bri;
			hue = (5.0f + (1.0f - blue/bri)/sat)/6.0f;
			}
		}

	}
