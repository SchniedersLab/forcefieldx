//******************************************************************************
//
// File:    ByteGrayColorModel.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.ByteGrayColorModel
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

import java.awt.Transparency;

import java.awt.color.ColorSpace;

import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;

/**
 * Class ByteGrayColorModel provides a ColorModel for a grayscale image that
 * obtains the color components for a pixel from a single byte.
 *
 * @author  Alan Kaminsky
 * @version 09-Jan-2007
 */
class ByteGrayColorModel
	extends ColorModel
	{

// Exported constructors.

	/**
	 * Construct a new byte gray color model.
	 */
	public ByteGrayColorModel()
		{
		super
			(/*int pixel_bits*/
				8,
			 /*int[] bits*/
				new int[] {8},
			 /*ColorSpace cspace*/
				ColorSpace.getInstance (ColorSpace.CS_GRAY),
			 /*boolean hasAlpha*/
				false,
			 /*boolean isAlphaPremultiplied*/
				false,
			 /*int transparency*/
				Transparency.OPAQUE,
			 /*int transferType*/
				DataBuffer.TYPE_BYTE);
		}

// Exported operations.

	/**
	 * Obtain the red color component for the given pixel.
	 *
	 * @param  pixel  Pixel value.
	 *
	 * @return  Red color component.
	 */
	public int getRed
		(int pixel)
		{
		return pixel & 0xFF;
		}

	/**
	 * Obtain the green color component for the given pixel.
	 *
	 * @param  pixel  Pixel value.
	 *
	 * @return  Green color component.
	 */
	public int getGreen
		(int pixel)
		{
		return pixel & 0xFF;
		}

	/**
	 * Obtain the blue color component for the given pixel.
	 *
	 * @param  pixel  Pixel value.
	 *
	 * @return  Blue color component.
	 */
	public int getBlue
		(int pixel)
		{
		return pixel & 0xFF;
		}

	/**
	 * Obtain the alpha component for the given pixel.
	 *
	 * @param  pixel  Pixel value.
	 *
	 * @return  Alpha component.
	 */
	public int getAlpha
		(int pixel)
		{
		return 255;
		}

	/**
	 * Obtain the red, green, blue, and alpha components for the given pixel.
	 *
	 * @param  pixel  Pixel value.
	 *
	 * @return  ARGB components.
	 */
	public int getRGB
		(int pixel)
		{
		int b = pixel & 0xFF;
		return 0xFF000000 | (b << 16) | (b << 8) | b;
		}

	/**
	 * Determine if the given raster is compatible with this color model.
	 *
	 * @param  raster  Raster.
	 *
	 * @return  True if <TT>raster</TT> is compatible with this color model,
	 *          false otherwise.
	 */
	public boolean isCompatibleRaster
		(Raster raster)
		{
		return
			raster.getSampleModel() instanceof SampleModelByteMatrix &&
			raster.getDataBuffer() instanceof DataBufferByteMatrix;
		}

	/**
	 * Create a writable raster compatible with this color model.
	 *
	 * @param  w  Width of raster.
	 * @param  h  Height of raster.
	 *
	 * @return  Writable raster.
	 */
	public WritableRaster createCompatibleWritableRaster
		(int w,
		 int h)
		{
		byte[][] matrix = new byte [h] [w];
		return Raster.createWritableRaster
			(/*SampleModel sm*/ new SampleModelByteMatrix (w, h),
			 /*DataBuffer db */ new DataBufferByteMatrix (matrix),
			 /*Point location*/ null);
		}

	/**
	 * Create a sample model compatible with this color model.
	 *
	 * @param  w  Width of sample model.
	 * @param  h  Height of sample model.
	 *
	 * @return  Sample model.
	 */
	public SampleModel createCompatibleSampleModel
		(int w,
		 int h)
		{
		return new SampleModelByteMatrix (w, h);
		}

	/**
	 * Determine if the given sample model is compatible with this color model.
	 *
	 * @param  sm  Sample model.
	 *
	 * @return  True if <TT>sm</TT> is compatible with this color model, false
	 *          otherwise.
	 */
	public boolean isCompatibleSampleModel
		(SampleModel sm)
		{
		return sm instanceof SampleModelByteMatrix;
		}

	}
