//******************************************************************************
//
// File:    ColorImageRow.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.ColorImageRow
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

package edu.rit.image;

import edu.rit.color.HSB;
import edu.rit.color.IntRGB;
import edu.rit.color.RGB;

import java.awt.Color;

import java.util.Arrays;

/**
 * Class ColorImageRow provides one row of a color image. The image row is
 * layered on top of an integer array (type <TT>int[]</TT>). Use the
 * <TT>setArray()</TT> method to specify which array this is. The width of the
 * image row is equal to the number of elements in the underlying array. The
 * image row's underlying array usually is (but does not have to be) one row of
 * an image's underlying matrix.
 * <P>
 * To read and write the pixels of an image row, use the <TT>getPixel()</TT>,
 * <TT>getPixelColor()</TT>, <TT>setPixel()</TT>, <TT>setPixelColor()</TT>, and
 * <TT>setPixelHSB()</TT> methods.
 * <P>
 * Changing the contents of the underlying array directly will also change the
 * image. The color information is stored in an array element as follows:
 * <UL>
 * <LI>Bits 31 .. 24 -- Unused, must be 0
 * <LI>Bits 23 .. 16 -- Red component in the range 0 .. 255
 * <LI>Bits 15 .. 8 -- Green component in the range 0 .. 255
 * <LI>Bits 7 .. 0 -- Blue component in the range 0 .. 255
 * </UL>
 * <P>
 * <B>HSB Colors.</B> A color may be specified using hue, saturation, and
 * brightness components instead of red, green, and blue components.
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
 * 
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public class ColorImageRow
	{

// Hidden data members.

	int[] myArray;

// Exported constructors.

	/**
	 * Construct a new color image row on top of the given array.
	 *
	 * @param  theArray  Underlying array.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public ColorImageRow
		(int[] theArray)
		{
		setArray (theArray);
		}

// Exported operations.

	/**
	 * Obtain this image row's underlying array.
	 *
	 * @return  Underlying array.
	 */
	public int[] getArray()
		{
		return myArray;
		}

	/**
	 * Set this image row's underlying array.
	 *
	 * @param  theArray  Underlying array.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public void setArray
		(int[] theArray)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("ColorImageRow.setArray(): theArray is null");
			}
		myArray = theArray;
		}

	/**
	 * Obtain the pixel at the given column in this image row.
	 *
	 * @param  c      Column index.
	 * @param  color  Floating point RGB color object in which to store the
	 *                pixel's color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void getPixel
		(int c,
		 RGB color)
		{
		color.unpack (myArray[c]);
		}

	/**
	 * Obtain the pixel at the given column in this image row.
	 *
	 * @param  c      Column index.
	 * @param  color  Integer RGB color object in which to store the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void getPixel
		(int c,
		 IntRGB color)
		{
		color.unpack (myArray[c]);
		}

	/**
	 * Obtain the pixel at the given column in this image row.
	 *
	 * @param  c      Column index.
	 * @param  color  Floating point HSB color object in which to store the
	 *                pixel's color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void getPixel
		(int c,
		 HSB color)
		{
		color.unpack (myArray[c]);
		}

	/**
	 * Obtain the pixel at the given column in this image row.
	 *
	 * @param  c  Column index.
	 *
	 * @return  AWT Color of pixel.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public Color getPixelColor
		(int c)
		{
		return new Color (myArray[c] & 0x00FFFFFF);
		}

	/**
	 * Set the pixel at the given column in this image row. If any component
	 * of <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  c      Column index.
	 * @param  color  Floating point RGB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixel
		(int c,
		 RGB color)
		{
		myArray[c] = color.pack();
		}

	/**
	 * Set the pixel at the given column in this image row. If any component
	 * of <TT>color</TT> is outside the range 0 .. 255, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  c      Column index.
	 * @param  color  Integer RGB color object containing the pixel's color
	 *                components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixel
		(int c,
		 IntRGB color)
		{
		myArray[c] = color.pack();
		}

	/**
	 * Set the pixel at the given column in this image row. If any component
	 * of <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  c      Column index.
	 * @param  color  Floating point HSB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixel
		(int c,
		 HSB color)
		{
		myArray[c] = color.pack();
		}

	/**
	 * Set the pixel at the given column in this image row.
	 *
	 * @param  c      Column index.
	 * @param  color  The pixel's AWT color.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixelColor
		(int c,
		 Color color)
		{
		myArray[c] = color.getRGB();
		}

	/**
	 * Set the pixel at the given column in this image row. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0.0
	 * .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  c      Column index in the range <TT>matrix().colRange()</TT>.
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixel
		(int c,
		 float red,
		 float green,
		 float blue)
		{
		myArray[c] = RGB.pack (red, green, blue);
		}

	/**
	 * Set the pixel at the given column in this image row. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0 ..
	 * 255, it is pinned to the appropriate boundary.
	 *
	 * @param  c      Column index in the range <TT>matrix().colRange()</TT>.
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixel
		(int c,
		 int red,
		 int green,
		 int blue)
		{
		myArray[c] = IntRGB.pack (red, green, blue);
		}

	/**
	 * Set the pixel at the given column in this image row. If any value
	 * <TT>hue</TT>, <TT>sat</TT>, or <TT>bri</TT> is outside the range 0.0
	 * .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  c    Column index in the range <TT>matrix().colRange()</TT>.
	 * @param  hue  Pixel's hue component.
	 * @param  sat  Pixel's saturation component.
	 * @param  bri  Pixel's brightness component.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixelHSB
		(int c,
		 float hue,
		 float sat,
		 float bri)
		{
		myArray[c] = HSB.pack (hue, sat, bri);
		}

	/**
	 * Set all pixels in this image row to the given color. If any component of
	 * <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  color  Floating point RGB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(RGB color)
		{
		Arrays.fill (myArray, color.pack());
		}

	/**
	 * Set all pixels in this image row to the given color. If any component of
	 * <TT>color</TT> is outside the range 0 .. 255, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  color  Integer RGB color object containing the pixel's color
	 *                components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(IntRGB color)
		{
		Arrays.fill (myArray, color.pack());
		}

	/**
	 * Set all pixels in this image row to the given color. If any component of
	 * <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  color  Floating point HSB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(HSB color)
		{
		Arrays.fill (myArray, color.pack());
		}

	/**
	 * Set all pixels in this image row to the given color.
	 *
	 * @param  color  The pixel's AWT color.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified. Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(Color color)
		{
		Arrays.fill (myArray, color.getRGB());
		}

	/**
	 * Set all pixels in this image row to the given color. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0.0
	 * .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 */
	public void fill
		(float red,
		 float green,
		 float blue)
		{
		Arrays.fill (myArray, RGB.pack (red, green, blue));
		}

	/**
	 * Set all pixels in this image row to the given color. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0 ..
	 * 255, it is pinned to the appropriate boundary.
	 *
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 */
	public void fill
		(int red,
		 int green,
		 int blue)
		{
		Arrays.fill (myArray, IntRGB.pack (red, green, blue));
		}

	/**
	 * Set all pixels in this image row to the given color. If any value
	 * <TT>hue</TT>, <TT>sat</TT>, or <TT>bri</TT> is outside the range 0.0 ..
	 * 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  hue  Pixel's hue component.
	 * @param  sat  Pixel's saturation component.
	 * @param  bri  Pixel's brightness component.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 */
	public void fillHSB
		(float hue,
		 float sat,
		 float bri)
		{
		Arrays.fill (myArray, HSB.pack (hue, sat, bri));
		}

	}
