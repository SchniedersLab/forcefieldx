//******************************************************************************
//
// File:    GrayImageRow.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.GrayImageRow
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

import java.util.Arrays;

/**
 * Class GrayImageRow provides one row of a grayscale image. The image row is
 * layered on top of a byte array (type <TT>byte[]</TT>). Use the
 * <TT>setArray()</TT> method to specify which array this is. The width of the
 * image row is equal to the number of elements in the underlying array. The
 * image row's underlying array usually is (but does not have to be) one row of
 * an image's underlying matrix.
 * <P>
 * To read and write the pixels of an image row, use the <TT>getIntPixel()</TT>,
 * <TT>getPixel()</TT>, <TT>setIntPixel()</TT>, and <TT>setPixel()</TT> methods.
 * These methods represent a pixel value as an integer in the range 0 .. 255 or
 * as a floating point number in the range 0.0 .. 1.0. The image row's
 * <I>interpretation</I> attribute specifies how to map between the numerical
 * pixel value and the shade of gray. If the interpretation is
 * <TT>ZERO_IS_WHITE</TT>, then a value of 0 or 0.0 represents white and 255 or
 * 1.0 represents black. If the interpretation is <TT>ZERO_IS_BLACK</TT>, then a
 * value of 0 or 0.0 represents black and 255 or 1.0 represents white. The
 * default interpretation is <TT>ZERO_IS_BLACK</TT>. To change the
 * interpretation, call the <TT>setInterpretation()</TT> method.
 * <P>
 * Changing the contents of the underlying matrix directly will also change the
 * image row. In the underlying matrix, a value of 0 represents black and a
 * value of 255 represents white.
 * 
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public class GrayImageRow
	{

// Exported enumerations.

	/**
	 * A pixel value of 0 or 0.0 is white, 255 or 1.0 is black.
	 */
	public static final PJGGrayImage.Interpretation ZERO_IS_WHITE =
		PJGGrayImage.Interpretation.ZERO_IS_WHITE;

	/**
	 * A pixel value of 0 or 0.0 is black, 255 or 1.0 is white.
	 */
	public static final PJGGrayImage.Interpretation ZERO_IS_BLACK =
		PJGGrayImage.Interpretation.ZERO_IS_BLACK;

// Hidden data members.

	PJGGrayImage.Transformation myTransformation =
		PJGGrayImage.ZERO_IS_BLACK_TRANSFORMATION;
	byte[] myArray;

// Exported constructors.

	/**
	 * Construct a new grayscale image row on top of the given array. The image
	 * row's interpretation is {@link #ZERO_IS_BLACK}.
	 *
	 * @param  theArray  Underlying array.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public GrayImageRow
		(byte[] theArray)
		{
		setArray (theArray);
		}

// Exported operations.

	/**
	 * Obtain this image row's underlying array.
	 *
	 * @return  Underlying array.
	 */
	public byte[] getArray()
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
		(byte[] theArray)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("GrayImageRow.setArray(): theArray is null");
			}
		myArray = theArray;
		}

	/**
	 * Get this image row's interpretation.
	 *
	 * @return  {@link #ZERO_IS_WHITE} or {@link #ZERO_IS_BLACK}.
	 */
	public PJGGrayImage.Interpretation getInterpretation()
		{
		if (myTransformation == PJGGrayImage.ZERO_IS_WHITE_TRANSFORMATION)
			{
			return ZERO_IS_WHITE;
			}
		else
			{
			return ZERO_IS_BLACK;
			}
		}

	/**
	 * Set this image row's interpretation. If this method is not called, this
	 * image row's default interpretation is {@link #ZERO_IS_BLACK}.
	 *
	 * @param  theInterpretation  {@link #ZERO_IS_WHITE} or
	 *                            {@link #ZERO_IS_BLACK}.
	 */
	public void setInterpretation
		(PJGGrayImage.Interpretation theInterpretation)
		{
		switch (theInterpretation)
			{
			case ZERO_IS_WHITE:
				myTransformation = PJGGrayImage.ZERO_IS_WHITE_TRANSFORMATION;
				break;
			case ZERO_IS_BLACK:
				myTransformation = PJGGrayImage.ZERO_IS_BLACK_TRANSFORMATION;
				break;
			}
		}

	/**
	 * Obtain the pixel at the given column in this image row.
	 *
	 * @param  c  Column index.
	 *
	 * @return  Pixel[<I>c</I>] as an integer in the range 0 .. 255.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public int getIntPixel
		(int c)
		{
		return myTransformation.inverseTransformInt (myArray[c]);
		}

	/**
	 * Obtain the pixel at the given column in this image row.
	 *
	 * @param  c  Column index.
	 *
	 * @return  Pixel[<I>c</I>] as a floating point number in the range 0.0 ..
	 *          1.0.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public float getPixel
		(int c)
		{
		return myTransformation.inverseTransformFloat (myArray[c]);
		}

	/**
	 * Set the pixel at the given column in this image row. If <TT>val</TT> is
	 * outside the range 0 .. 255, it is pinned to the appropriate boundary.
	 *
	 * @param  c    Column index.
	 * @param  val  New pixel value as an integer in the range 0 .. 255.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setIntPixel
		(int c,
		 int val)
		{
		myArray[c] = myTransformation.transformInt (val);
		}

	/**
	 * Set the pixel at the given column in this image row. If <TT>val</TT> is
	 * outside the range 0.0 .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  c    Column index.
	 * @param  val  New pixel value as a floating point number in the range 0.0
	 *              .. 1.0.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>c</TT> is out of bounds.
	 */
	public void setPixel
		(int c,
		 float val)
		{
		myArray[c] = myTransformation.transformFloat (val);
		}

	/**
	 * Set all pixels in this image row to the given value. If <TT>val</TT> is
	 * outside the range 0 .. 255, it is pinned to the appropriate boundary.
	 *
	 * @param  val  New pixel value as an integer in the range 0 .. 255.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 */
	public void fill
		(int val)
		{
		Arrays.fill (myArray, myTransformation.transformInt (val));
		}

	/**
	 * Set all pixels in this image row to the given value. If <TT>val</TT> is
	 * outside the range 0.0 .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  val  New pixel value as a floating point number in the range 0.0
	 *              .. 1.0.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix row has not
	 *     been specified.
	 */
	public void fill
		(float val)
		{
		Arrays.fill (myArray, myTransformation.transformFloat (val));
		}

	}
