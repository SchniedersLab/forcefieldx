//******************************************************************************
//
// File:    BaseColorImage.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.BaseColorImage
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

package edu.rit.image;

import edu.rit.color.HSB;
import edu.rit.color.IntRGB;
import edu.rit.color.RGB;

import edu.rit.swing.Displayable;

import java.awt.Color;

import java.awt.image.BufferedImage;

import java.util.Arrays;

/**
 * Class BaseColorImage is the abstract superclass for a color image file in
 * Parallel Java Graphics (PJG) format. The image is layered on top of an
 * integer matrix (type <TT>int[][]</TT>). The height and width of the image are
 * equal to the number of rows and columns in the underlying matrix.
 * <P>
 * To get and set the image's pixel data, use the <TT>getPixel()</TT>,
 * <TT>getPixelColor()</TT>, <TT>setPixel()</TT>, <TT>setPixelColor()</TT>, and
 * <TT>setPixelHSB()</TT> methods. You only need to allocate storage in the
 * pixel data matrix for the portions of the image you are actually accessing;
 * the complete matrix need not be allocated. Class {@linkplain
 * edu.rit.util.Arrays} has static methods for allocating portions of a matrix.
 * <P>
 * Changing the contents of the underlying matrix directly will also change the
 * image. The color information is stored in a matrix element as follows:
 * <UL>
 * <LI>Bits 31 .. 24 -- Unused, must be 0
 * <LI>Bits 23 .. 16 -- Red component in the range 0 .. 255
 * <LI>Bits 15 .. 8 -- Green component in the range 0 .. 255
 * <LI>Bits 7 .. 0 -- Blue component in the range 0 .. 255
 * </UL>
 * <P>
 * A color may be specified using hue, saturation, and brightness components
 * instead of red, green, and blue components.
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
 * To write a BaseColorImage object to a PJG image file, call the
 * <TT>prepareToWrite()</TT> method, specifying the output stream to write. The
 * <TT>prepareToWrite()</TT> method returns an instance of class {@linkplain
 * PJGImage.Writer}. Call the methods of the PJG image writer object to write
 * the pixel data, or sections of the pixel data, to the output stream. When
 * finished, close the PJG image writer.
 * <P>
 * To read a BaseColorImage object from a PJG image file, call the
 * <TT>prepareToRead()</TT> method, specifying the input stream to read. The
 * <TT>prepareToRead()</TT> method returns an instance of class {@linkplain
 * PJGImage.Reader}. Call the methods of the PJG image reader object to read
 * the pixel data, or sections of the pixel data, from the input stream. When
 * finished, close the PJG image reader.
 * <P>
 * To get a BufferedImage object that uses the same underlying pixel data matrix
 * as the BaseColorImage object, call the <TT>getBufferedImage()</TT> method.
 * You can then do all the following with the BufferedImage: display it on the
 * screen, draw into it using a graphics context, copy another BufferedImage
 * into it, read it from or write it to a file using package javax.imageio
 * (which typically supports PNG, JPG, and GIF formats). The rows and columns of
 * the underlying matrix need not all be allocated when accessing the
 * BufferedImage. If you get a pixel from the BufferedImage in an unallocated
 * row or column, a pixel value of 0 (black) is returned. If you set a pixel in
 * the BufferedImage in an unallocated row or column, the pixel value is
 * discarded.
 * <P>
 * <I>Note:</I> Class BaseColorImage is not multiple thread safe.
 * 
 * @author  Alan Kaminsky
 * @version 08-Apr-2008
 */
public abstract class BaseColorImage
	extends PJGImage
	{

// Hidden data members.

	int[][] myMatrix;

// Hidden constructors.

	/**
	 * Construct a new base color image. The image's height and width are
	 * uninitialized. Before accessing the image's pixels, specify the height
	 * and width by calling a subclass method or by reading the image from an
	 * input stream.
	 *
	 * @param  theImageType  Image type.
	 */
	BaseColorImage
		(int theImageType)
		{
		super (theImageType);
		}

// Exported operations.

	/**
	 * Obtain this image's underlying matrix.
	 *
	 * @return  Underlying matrix, or null if none.
	 */
	public int[][] getMatrix()
		{
		return myMatrix;
		}

	/**
	 * Set this image's height, width, and underlying matrix.
	 *
	 * @param  theHeight  Image height in pixels.
	 * @param  theWidth   Image width in pixels.
	 * @param  theMatrix  Underlying matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> &lt;= 0. Thrown if
	 *     <TT>theWidth</TT> &lt;= 0. Thrown if <TT>theMatrix.length</TT> does
	 *     not equal <TT>theHeight</TT>.
	 */
	public void setMatrix
		(int theHeight,
		 int theWidth,
		 int[][] theMatrix)
		{
		setHeightAndWidth (theHeight, theWidth);
		if (theMatrix.length != myHeight)
			{
			throw new IllegalArgumentException
				(getClass().getName() + ".setMatrix(): theMatrix.length (= " +
				 theMatrix.length + ") does not equal image height (= " +
				 myHeight + ")");
			}
		myMatrix = theMatrix;
		}

	/**
	 * Obtain the pixel at the given row and column in this image.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  Floating point RGB color object in which to store the
	 *                pixel's color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void getPixel
		(int r,
		 int c,
		 RGB color)
		{
		color.unpack (myMatrix[r][c]);
		}

	/**
	 * Obtain the pixel at the given row and column in this image.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  Integer RGB color object in which to store the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void getPixel
		(int r,
		 int c,
		 IntRGB color)
		{
		color.unpack (myMatrix[r][c]);
		}

	/**
	 * Obtain the pixel at the given row and column in this image.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  Floating point HSB color object in which to store the
	 *                pixel's color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void getPixel
		(int r,
		 int c,
		 HSB color)
		{
		color.unpack (myMatrix[r][c]);
		}

	/**
	 * Obtain the pixel at the given row and column in this image.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 *
	 * @return  AWT color of pixel.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public Color getPixelColor
		(int r,
		 int c)
		{
		return new Color (myMatrix[r][c] & 0x00FFFFFF);
		}

	/**
	 * Set the pixel at the given row and column in this image. If any component
	 * of <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  Floating point RGB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixel
		(int r,
		 int c,
		 RGB color)
		{
		myMatrix[r][c] = color.pack();
		}

	/**
	 * Set the pixel at the given row and column in this image. If any component
	 * of <TT>color</TT> is outside the range 0 .. 255, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  Integer RGB color object containing the pixel's color
	 *                components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixel
		(int r,
		 int c,
		 IntRGB color)
		{
		myMatrix[r][c] = color.pack();
		}

	/**
	 * Set the pixel at the given row and column in this image. If any component
	 * of <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  Floating point HSB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixel
		(int r,
		 int c,
		 HSB color)
		{
		myMatrix[r][c] = color.pack();
		}

	/**
	 * Set the pixel at the given row and column in this image.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  color  The pixel's AWT color.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixelColor
		(int r,
		 int c,
		 Color color)
		{
		myMatrix[r][c] = color.getRGB();
		}

	/**
	 * Set the pixel at the given row and column in this image. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0.0
	 * .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  r      Row index.
	 * @param  c      Column index.
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixel
		(int r,
		 int c,
		 float red,
		 float green,
		 float blue)
		{
		myMatrix[r][c] = RGB.pack (red, green, blue);
		}

	/**
	 * Set the pixel at the given row and column in this image. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0 ..
	 * 255, it is pinned to the appropriate boundary.
	 *
	 * @param  r      Row index in the range <TT>matrix().rowRange()</TT>.
	 * @param  c      Column index in the range <TT>matrix().colRange()</TT>.
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixel
		(int r,
		 int c,
		 int red,
		 int green,
		 int blue)
		{
		myMatrix[r][c] = IntRGB.pack (red, green, blue);
		}

	/**
	 * Set the pixel at the given row and column in this image. If any value
	 * <TT>hue</TT>, <TT>sat</TT>, or <TT>bri</TT> is outside the range 0.0
	 * .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  r    Row index in the range <TT>matrix().rowRange()</TT>.
	 * @param  c    Column index in the range <TT>matrix().colRange()</TT>.
	 * @param  hue  Pixel's hue component.
	 * @param  sat  Pixel's saturation component.
	 * @param  bri  Pixel's brightness component.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixelHSB
		(int r,
		 int c,
		 float hue,
		 float sat,
		 float bri)
		{
		myMatrix[r][c] = HSB.pack (hue, sat, bri);
		}

	/**
	 * Set all pixels in this image to the given color. If any component of
	 * <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  color  Floating point RGB color object containing the pixel's
	 *                color components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(RGB color)
		{
		int val = color.pack();
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Set all pixels in this image to the given color. If any component of
	 * <TT>color</TT> is outside the range 0 .. 255, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  color  Integer RGB color object containing the pixel's color
	 *                components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(IntRGB color)
		{
		int val = color.pack();
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Set all pixels in this image to the given color. If any component of
	 * <TT>color</TT> is outside the range 0.0 .. 1.0, it is pinned to the
	 * appropriate boundary.
	 *
	 * @param  color  Floating point HSB object containing the pixel's color
	 *                components.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(HSB color)
		{
		int val = color.pack();
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Set all pixels in this image to the given color.
	 *
	 * @param  color  The pixel's AWT color.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>color</TT> is null.
	 */
	public void fill
		(Color color)
		{
		int val = color.getRGB();
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Set all pixels in this image to the given color. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0.0
	 * .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 */
	public void fill
		(float red,
		 float green,
		 float blue)
		{
		int val = RGB.pack (red, green, blue);
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Set all pixels in this image to the given color. If any value
	 * <TT>red</TT>, <TT>green</TT>, or <TT>blue</TT> is outside the range 0 ..
	 * 255, it is pinned to the appropriate boundary.
	 *
	 * @param  red    Pixel's red component.
	 * @param  green  Pixel's green component.
	 * @param  blue   Pixel's blue component.
	 */
	public void fill
		(int red,
		 int green,
		 int blue)
		{
		int val = IntRGB.pack (red, green, blue);
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Set all pixels in this image to the given color. If any value
	 * <TT>hue</TT>, <TT>sat</TT>, or <TT>bri</TT> is outside the range 0.0 ..
	 * 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  hue  Pixel's hue component.
	 * @param  sat  Pixel's saturation component.
	 * @param  bri  Pixel's brightness component.
	 */
	public void fillHSB
		(float hue,
		 float sat,
		 float bri)
		{
		int val = HSB.pack (hue, sat, bri);
		int rows = myMatrix.length;
		for (int r = 0; r < rows; ++ r)
			{
			Arrays.fill (myMatrix[r], val);
			}
		}

	/**
	 * Obtain a BufferedImage whose pixel data comes from this image's
	 * underlying matrix.
	 *
	 * @return  BufferedImage.
	 */
	public BufferedImage getBufferedImage()
		{
		return new ColorBufferedImage (myHeight, myWidth, myMatrix);
		}

	/**
	 * Obtain a Displayable object with which to display this image in a Swing
	 * UI.
	 *
	 * @return  Displayable object.
	 */
	public Displayable getDisplayable()
		{
		return new ColorDisplayable (myHeight, myWidth, myMatrix);
		}

	}
