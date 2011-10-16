//******************************************************************************
//
// File:    IntegerColorModel.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.IntegerColorModel
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

package edu.rit.image;

import java.awt.Transparency;

import java.awt.color.ColorSpace;

import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;

/**
 * Class IntegerColorModel provides a ColorModel that obtains the color
 * components for a pixel from a single integer.
 *
 * @author  Alan Kaminsky
 * @version 12-Feb-2006
 */
class IntegerColorModel
	extends ColorModel
	{

// Exported constructors.

	/**
	 * Construct a new integer color model.
	 */
	public IntegerColorModel()
		{
		super
			(/*int pixel_bits*/
				24,
			 /*int[] bits*/
				new int[] {8, 8, 8},
			 /*ColorSpace cspace*/
				ColorSpace.getInstance (ColorSpace.CS_sRGB),
			 /*boolean hasAlpha*/
				false,
			 /*boolean isAlphaPremultiplied*/
				false,
			 /*int transparency*/
				Transparency.OPAQUE,
			 /*int transferType*/
				DataBuffer.TYPE_INT);
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
		return (pixel >> 16) & 0xFF;
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
		return (pixel >> 8) & 0xFF;
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
		return 0xFF000000 | pixel;
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
			raster.getSampleModel() instanceof SampleModelIntegerMatrix &&
			raster.getDataBuffer() instanceof DataBufferIntegerMatrix;
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
		int[][] matrix = new int [h] [w];
		return Raster.createWritableRaster
			(/*SampleModel sm*/ new SampleModelIntegerMatrix (w, h),
			 /*DataBuffer db */ new DataBufferIntegerMatrix (matrix),
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
		return new SampleModelIntegerMatrix (w, h);
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
		return sm instanceof SampleModelIntegerMatrix;
		}

	}
