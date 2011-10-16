//******************************************************************************
//
// File:    SampleModelByteMatrix.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.SampleModelByteMatrix
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

import java.awt.image.DataBuffer;
import java.awt.image.SampleModel;

/**
 * Class SampleModelByteMatrix provides a SampleModel that stores the value of a
 * pixel in an element of an underlying byte matrix.
 *
 * @author  Alan Kaminsky
 * @version 10-Nov-2007
 */
class SampleModelByteMatrix
	extends SampleModel
	{

// Exported constructors.

	/**
	 * Construct a new byte matrix sample model of the given size.
	 *
	 * @param  w  Width.
	 * @param  h  Height.
	 */
	public SampleModelByteMatrix
		(int w,
		 int h)
		{
		super (DataBuffer.TYPE_BYTE, w, h, 1);
		}

// Exported operations.

	/**
	 * Obtain the number of data elements needed to transfer a pixel via the
	 * <TT>getDataElements()</TT> and <TT>setDataElements()</TT> methods.
	 *
	 * @return  Number of data elements.
	 */
	public int getNumDataElements()
		{
		return 1;
		}

	/**
	 * Obtain the data for the given pixel in a primitive array of the transfer
	 * type.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  obj   Array in which to store the pixel data. If null, a new
	 *               array is created.
	 * @param  data  Data buffer.
	 *
	 * @return  Array containing the pixel data.
	 */
	public Object getDataElements
		(int x,
		 int y,
		 Object obj,
		 DataBuffer data)
		{
		byte[] bArray = (byte[]) obj;
		if (bArray == null) bArray = new byte [1];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			bArray[0] = matrix_y[x];
			}
		else
			{
			bArray[0] = (byte) 0;
			}
		return bArray;
		}

	/**
	 * Obtain the data for the given block of pixels in a primitive array of the
	 * transfer type.
	 *
	 * @param  x     Upper left pixel X coordinate.
	 * @param  y     Upper left pixel Y coordinate.
	 * @param  w     Block width.
	 * @param  h     Block height.
	 * @param  obj   Array in which to store the pixel data. If null, a new
	 *               array is created.
	 * @param  data  Data buffer.
	 *
	 * @return  Array containing the pixel data.
	 */
	public Object getDataElements
		(int x,
		 int y,
		 int w,
		 int h,
		 Object obj,
		 DataBuffer data)
		{
		byte[] bArray = (byte[]) obj;
		if (bArray == null) bArray = new byte [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					bArray[p++] = (byte) 0;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					bArray[p++] = matrix_r[c];
					}
				for (int c = max; c < colub; ++ c)
					{
					bArray[p++] = (byte) 0;
					}
				}
			}
		return bArray;
		}

	/**
	 * Set the data for the given pixel from a primitive array of the transfer
	 * type.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  obj   Array containing the pixel data.
	 * @param  data  Data buffer.
	 */
	public void setDataElements
		(int x,
		 int y,
		 Object obj,
		 DataBuffer data)
		{
		byte[] bArray = (byte[]) obj;
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = bArray[0];
			}
		}

	/**
	 * Set the data for the given block of pixels from a primitive array of the
	 * transfer type.
	 *
	 * @param  x     Upper left pixel X coordinate.
	 * @param  y     Upper left pixel Y coordinate.
	 * @param  w     Block width.
	 * @param  h     Block height.
	 * @param  obj   Array containing the pixel data.
	 * @param  data  Data buffer.
	 */
	public void setDataElements
		(int x,
		 int y,
		 int w,
		 int h,
		 Object obj,
		 DataBuffer data)
		{
		byte[] bArray = (byte[]) obj;
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = bArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Obtain the samples for the given pixel in an <TT>int</TT> array, one
	 * sample per array element.
	 *
	 * @param  x       Pixel X coordinate.
	 * @param  y       Pixel Y coordinate.
	 * @param  iArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public int[] getPixel
		(int x,
		 int y,
		 int[] iArray,
		 DataBuffer data)
		{
		if (iArray == null) iArray = new int [1];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			iArray[0] = matrix_y[x] & 0xFF;
			}
		else
			{
			iArray[0] = 0;
			}
		return iArray;
		}

	/**
	 * Obtain the samples for the given pixel in a <TT>float</TT> array, one
	 * sample per array element.
	 *
	 * @param  x       Pixel X coordinate.
	 * @param  y       Pixel Y coordinate.
	 * @param  fArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public float[] getPixel
		(int x,
		 int y,
		 float[] fArray,
		 DataBuffer data)
		{
		if (fArray == null) fArray = new float [1];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			fArray[0] = matrix_y[x] & 0xFF;
			}
		else
			{
			fArray[0] = 0.0f;
			}
		return fArray;
		}

	/**
	 * Obtain the samples for the given pixel in a <TT>double</TT> array, one
	 * sample per array element.
	 *
	 * @param  x       Pixel X coordinate.
	 * @param  y       Pixel Y coordinate.
	 * @param  dArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public double[] getPixel
		(int x,
		 int y,
		 double[] dArray,
		 DataBuffer data)
		{
		if (dArray == null) dArray = new double [1];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			dArray[0] = matrix_y[x] & 0xFF;
			}
		else
			{
			dArray[0] = 0.0;
			}
		return dArray;
		}

	/**
	 * Obtain the samples for the given block of pixels in an <TT>int</TT>
	 * array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  iArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public int[] getPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 int[] iArray,
		 DataBuffer data)
		{
		if (iArray == null) iArray = new int [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					iArray[p++] = 0;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					iArray[p++] = matrix_r[c] & 0xFF;
					}
				for (int c = max; c < colub; ++ c)
					{
					iArray[p++] = 0;
					}
				}
			}
		return iArray;
		}

	/**
	 * Obtain the samples for the given block of pixels in a <TT>float</TT>
	 * array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  fArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public float[] getPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 float[] fArray,
		 DataBuffer data)
		{
		if (fArray == null) fArray = new float [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					fArray[p++] = 0.0f;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					fArray[p++] = matrix_r[c] & 0xFF;
					}
				for (int c = max; c < colub; ++ c)
					{
					fArray[p++] = 0.0f;
					}
				}
			}
		return fArray;
		}

	/**
	 * Obtain the samples for the given block of pixels in a <TT>double</TT>
	 * array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  fArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public double[] getPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 double[] dArray,
		 DataBuffer data)
		{
		if (dArray == null) dArray = new double [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					dArray[p++] = 0.0;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					dArray[p++] = matrix_r[c] & 0xFF;
					}
				for (int c = max; c < colub; ++ c)
					{
					dArray[p++] = 0.0;
					}
				}
			}
		return dArray;
		}

	/**
	 * Obtain the sample in the given band for the given pixel as an
	 * <TT>int</TT>.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  b     Band.
	 * @param  data  Data buffer.
	 *
	 * @return  Sample value.
	 */
	public int getSample
		(int x,
		 int y,
		 int b,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			return matrix_y[x] & 0xFF;
			}
		else
			{
			return 0;
			}
		}

	/**
	 * Obtain the sample in the given band for the given pixel as a
	 * <TT>float</TT>.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  b     Band.
	 * @param  data  Data buffer.
	 *
	 * @return  Sample value.
	 */
	public float getSampleFloat
		(int x,
		 int y,
		 int b,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			return matrix_y[x] & 0xFF;
			}
		else
			{
			return 0.0f;
			}
		}

	/**
	 * Obtain the sample in the given band for the given pixel as a
	 * <TT>double</TT>.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  b     Band.
	 * @param  data  Data buffer.
	 *
	 * @return  Sample value.
	 */
	public double getSampleDouble
		(int x,
		 int y,
		 int b,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			return matrix_y[x] & 0xFF;
			}
		else
			{
			return 0.0;
			}
		}

	/**
	 * Obtain the samples in the given band for the given block of pixels in an
	 * <TT>int</TT> array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  b       Band.
	 * @param  iArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public int[] getSamples
		(int x,
		 int y,
		 int w,
		 int h,
		 int b,
		 int[] iArray,
		 DataBuffer data)
		{
		if (iArray == null) iArray = new int [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					iArray[p++] = 0;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					iArray[p++] = matrix_r[c] & 0xFF;
					}
				for (int c = max; c < colub; ++ c)
					{
					iArray[p++] = 0;
					}
				}
			}
		return iArray;
		}

	/**
	 * Obtain the samples in the given band for the given block of pixels in a
	 * <TT>float</TT> array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  b       Band.
	 * @param  fArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public float[] getSamples
		(int x,
		 int y,
		 int w,
		 int h,
		 int b,
		 float[] fArray,
		 DataBuffer data)
		{
		if (fArray == null) fArray = new float [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					fArray[p++] = 0.0f;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					fArray[p++] = matrix_r[c] & 0xFF;
					}
				for (int c = max; c < colub; ++ c)
					{
					fArray[p++] = 0.0f;
					}
				}
			}
		return fArray;
		}

	/**
	 * Obtain the samples in the given band for the given block of pixels in a
	 * <TT>double</TT> array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  b       Band.
	 * @param  fArray  Array in which to store the pixel samples. If null, a new
	 *                 array is created.
	 * @param  data    Data buffer.
	 *
	 * @return  Array containing the pixel samples.
	 */
	public double[] getSamples
		(int x,
		 int y,
		 int w,
		 int h,
		 int b,
		 double[] dArray,
		 DataBuffer data)
		{
		if (dArray == null) dArray = new double [w*h];
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				for (int c = x; c < colub; ++ c)
					{
					dArray[p++] = 0.0;
					}
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					dArray[p++] = matrix_r[c] & 0xFF;
					}
				for (int c = max; c < colub; ++ c)
					{
					dArray[p++] = 0.0;
					}
				}
			}
		return dArray;
		}

	/**
	 * Set the samples for the given pixel from an <TT>int</TT> array, one
	 * sample per array element.
	 *
	 * @param  x       Pixel X coordinate.
	 * @param  y       Pixel Y coordinate.
	 * @param  iArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixel
		(int x,
		 int y,
		 int[] iArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = (byte) iArray[0];
			}
		}

	/**
	 * Set the samples for the given pixel from a <TT>float</TT> array, one
	 * sample per array element.
	 *
	 * @param  x       Pixel X coordinate.
	 * @param  y       Pixel Y coordinate.
	 * @param  fArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixel
		(int x,
		 int y,
		 float[] fArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = (byte) fArray[0];
			}
		}

	/**
	 * Set the samples for the given pixel from a <TT>double</TT> array, one
	 * sample per array element.
	 *
	 * @param  x       Pixel X coordinate.
	 * @param  y       Pixel Y coordinate.
	 * @param  dArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixel
		(int x,
		 int y,
		 double[] dArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = (byte) dArray[0];
			}
		}

	/**
	 * Set the samples for the given block of pixels from an <TT>int</TT>
	 * array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  iArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 int[] iArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = (byte) iArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Set the samples for the given block of pixels from a <TT>float</TT>
	 * array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  fArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 float[] fArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = (byte) fArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Set the samples for the given block of pixels from a <TT>double</TT>
	 * array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  dArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 double[] dArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = (byte) dArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Set the sample in the given band for the given pixel from an
	 * <TT>int</TT>.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  b     Band.
	 * @param  s     Sample value.
	 * @param  data  Data buffer.
	 */
	public void setSample
		(int x,
		 int y,
		 int b,
		 int s,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = (byte) s;
			}
		}

	/**
	 * Set the sample in the given band for the given pixel from a
	 * <TT>float</TT>.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  b     Band.
	 * @param  s     Sample value.
	 * @param  data  Data buffer.
	 */
	public void setSample
		(int x,
		 int y,
		 int b,
		 float s,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = (byte) s;
			}
		}

	/**
	 * Set the sample in the given band for the given pixel from a
	 * <TT>double</TT>.
	 *
	 * @param  x     Pixel X coordinate.
	 * @param  y     Pixel Y coordinate.
	 * @param  b     Band.
	 * @param  s     Sample value.
	 * @param  data  Data buffer.
	 */
	public void setSample
		(int x,
		 int y,
		 int b,
		 double s,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		byte[] matrix_y = matrix[y];
		if (matrix_y != null && x < matrix_y.length)
			{
			matrix_y[x] = (byte) s;
			}
		}

	/**
	 * Set the samples in the given band for the given block of pixels from an
	 * <TT>int</TT> array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  b       Band.
	 * @param  iArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 int b,
		 int[] iArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = (byte) iArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Set the samples in the given band for the given block of pixels from a
	 * <TT>float</TT> array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  b       Band.
	 * @param  fArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 int b,
		 float[] fArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = (byte) fArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Set the samples in the given band for the given block of pixels from a
	 * <TT>double</TT> array, one sample per array element.
	 *
	 * @param  x       Upper left pixel X coordinate.
	 * @param  y       Upper left pixel Y coordinate.
	 * @param  w       Block width.
	 * @param  h       Block height.
	 * @param  b       Band.
	 * @param  dArray  Array containing the pixel samples.
	 * @param  data    Data buffer.
	 */
	public void setPixels
		(int x,
		 int y,
		 int w,
		 int h,
		 int b,
		 double[] dArray,
		 DataBuffer data)
		{
		byte[][] matrix = ((DataBufferByteMatrix) data).myMatrix;
		int rowub = y + h;
		int colub = x + w;
		int p = 0;
		for (int r = y; r < rowub; ++ r)
			{
			byte[] matrix_r = matrix[r];
			if (matrix_r == null)
				{
				p += w;
				}
			else
				{
				int max = Math.min (colub, matrix_r.length);
				for (int c = x; c < max; ++ c)
					{
					matrix_r[c] = (byte) dArray[p++];
					}
				if (max < colub)
					{
					p += colub - max;
					}
				}
			}
		}

	/**
	 * Create a sample model compatible with this sample model, but with the
	 * given width and height.
	 *
	 * @param  w  Width.
	 * @param  h  Height.
	 *
	 * @return  New sample model.
	 */
	public SampleModel createCompatibleSampleModel
		(int w,
		 int h)
		{
		return new SampleModelByteMatrix (w, h);
		}

	/**
	 * Create a new sample model with a subset of the bands of this sample
	 * model. <I>Note:</I> This method ignores <TT>bands</TT> and just returns a
	 * compatible sample model with the same width and height as this one.
	 *
	 * @param  bands  Subset of the bands of this sample model.
	 *
	 * @return  New sample model.
	 */
	public SampleModel createSubsetSampleModel
		(int[] bands)
		{
		return new SampleModelByteMatrix (getWidth(), getHeight());
		}

	/**
	 * Create a data buffer that corresponds to this sample model. The data
	 * buffer's width and height are the same as this sample model's.
	 *
	 * @return  Data buffer.
	 */
	public DataBuffer createDataBuffer()
		{
		byte[][] matrix = new byte [getHeight()] [getWidth()];
		return new DataBufferByteMatrix (matrix);
		}

	/**
	 * Obtain the size in bits of samples for all bands.
	 *
	 * @return  Array of band bit sizes.
	 */
	public int[] getSampleSize()
		{
		return new int[] {8};
		}

	/**
	 * Obtain the size in bits of samples for the given band.
	 *
	 * @param  band  Band.
	 *
	 * @return  Band bit size.
	 */
	public int getSampleSize
		(int band)
		{
		return 8;
		}

	}
