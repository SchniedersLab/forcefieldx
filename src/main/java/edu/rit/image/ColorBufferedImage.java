//******************************************************************************
//
// File:    ColorBufferedImage.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.ColorBufferedImage
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

import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.Raster;
import java.awt.image.WritableRaster;

/**
 * Class ColorBufferedImage provides a color image. The image is layered on top
 * of an integer matrix (type <TT>int[][]</TT>). The height and width of the
 * image are equal to the number of rows and columns in the underlying matrix.
 * The color information is stored in a matrix element as follows:
 * <UL>
 * <LI>Bits 31 .. 24 -- Unused, must be 0
 * <LI>Bits 23 .. 16 -- Red component in the range 0 .. 255
 * <LI>Bits 15 .. 8 -- Green component in the range 0 .. 255
 * <LI>Bits 7 .. 0 -- Blue component in the range 0 .. 255
 * </UL>
 * <P>
 * The rows and columns of the underlying matrix need not all be allocated. If
 * you get a pixel in an unallocated row or column, a pixel value of 0 (black)
 * is returned. If you set a pixel in an unallocated row or column, the pixel
 * value is discarded.
 * <P>
 * <B>BufferedImage.</B> Since class ColorBufferedImage is a subclass of class
 * {@linkplain java.awt.image.BufferedImage}, you can do all the following with
 * a color image just as you would with a BufferedImage: display it on the
 * screen, draw into it using a graphics context, copy another BufferedImage
 * into it, read it from or write it to a file using package {@linkplain
 * javax.imageio}.
 * 
 * @author  Alan Kaminsky
 * @version 01-Nov-2007
 */
class ColorBufferedImage
	extends BufferedImage
	{

// Exported constructors.

	/**
	 * Construct a new color image.
	 *
	 * @param  theHeight  Image height in pixels.
	 * @param  theWidth   Image width in pixels.
	 * @param  theMatrix  Underlying matrix.
	 */
	public ColorBufferedImage
		(int theHeight,
		 int theWidth,
		 int[][] theMatrix)
		{
		this
			(theMatrix,
			 new IntegerColorModel(),
			 makeRaster (theHeight, theWidth, theMatrix));
		}

	/**
	 * Construct a new color image with the given color model and raster.
	 *
	 * @param  theMatrix      Underlying matrix.
	 * @param  theColorModel  Color model.
	 * @param  theRaster      Writable raster.
	 */
	private ColorBufferedImage
		(int[][] theMatrix,
		 ColorModel theColorModel,
		 WritableRaster theRaster)
		{
		// Construct superclass BufferedImage.
		super
			(/*ColorModel cm                */ theColorModel,
			 /*WritableRaster raster        */ theRaster,
			 /*boolean isRasterPremultiplied*/ false,
			 /*Hashtable<?,?> properties    */ null);
		}

	private static WritableRaster makeRaster
		(int theHeight,
		 int theWidth,
		 int[][] theMatrix)
		{
		return Raster.createWritableRaster
			(/*SampleModel sm*/
				new SampleModelIntegerMatrix (theWidth, theHeight),
			 /*DataBuffer db */
				new DataBufferIntegerMatrix (theMatrix),
			 /*Point location*/
				null);
		}

	}
