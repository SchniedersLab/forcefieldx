//******************************************************************************
//
// File:    GrayBufferedImage.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.GrayBufferedImage
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
 * Class GrayBufferedImage provides a grayscale image. The image is layered on
 * top of a byte matrix (type <TT>byte[][]</TT>). The height and width of the
 * image are equal to the number of rows and columns in the underlying matrix.
 * The gray information is stored in a matrix element as a number from 0 (black)
 * to 255 (white).
 * <P>
 * The rows and columns of the underlying matrix need not all be allocated. If
 * you get a pixel in an unallocated row or column, a pixel value of 0 (black)
 * is returned. If you set a pixel in an unallocated row or column, the pixel
 * value is discarded.
 * <P>
 * Since class GrayBufferedImage is a subclass of class {@linkplain
 * java.awt.image.BufferedImage}, you can do all the following with a grayscale
 * image just as you would with a BufferedImage: display it on the screen, draw
 * into it using a graphics context, copy another grayscale BufferedImage into
 * it, read it from or write it to a file using package {@linkplain
 * javax.imageio}.
 * 
 * @author  Alan Kaminsky
 * @version 10-Nov-2007
 */
class GrayBufferedImage
	extends BufferedImage
	{

// Exported constructors.

	/**
	 * Construct a new grayscale image.
	 *
	 * @param  theHeight  Image height in pixels.
	 * @param  theWidth   Image width in pixels.
	 * @param  theMatrix  Underlying byte matrix.
	 */
	public GrayBufferedImage
		(int theHeight,
		 int theWidth,
		 byte[][] theMatrix)
		{
		this
			(theMatrix,
			 new ByteGrayColorModel(),
			 makeRaster (theHeight, theWidth, theMatrix));
		}

	/**
	 * Construct a new grayscale image with the given color model and raster.
	 *
	 * @param  theMatrix      Underlying byte matrix.
	 * @param  theColorModel  Color model.
	 * @param  theRaster      Raster.
	 */
	private GrayBufferedImage
		(byte[][] theMatrix,
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
		 byte[][] theMatrix)
		{
		return Raster.createWritableRaster
			(/*SampleModel sm*/
				new SampleModelByteMatrix (theWidth, theHeight),
			 /*DataBuffer db */
				new DataBufferByteMatrix (theMatrix),
			 /*Point location*/
				null);
		}

	}
