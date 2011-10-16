//******************************************************************************
//
// File:    DataBufferIntegerMatrix.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.DataBufferIntegerMatrix
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

import java.awt.image.DataBuffer;

/**
 * Class DataBufferIntegerMatrix provides a DataBuffer that stores its data in
 * an underlying integer matrix.
 *
 * @author  Alan Kaminsky
 * @version 12-Feb-2006
 */
class DataBufferIntegerMatrix
	extends DataBuffer
	{

// Hidden data members.

	int[][] myMatrix;
	private int myWidth;

// Exported constructors.

	/**
	 * Construct a new data buffer for the given integer matrix.
	 *
	 * @param  theMatrix  Underlying matrix.
	 */
	public DataBufferIntegerMatrix
		(int[][] theMatrix)
		{
		super (DataBuffer.TYPE_INT, height (theMatrix) * width (theMatrix));
		myMatrix = theMatrix;
		myWidth = width (theMatrix);
		}

	static int height
		(int[][] theMatrix)
		{
		return theMatrix.length;
		}

	static int width
		(int[][] theMatrix)
		{
		return theMatrix.length == 0 ? 0 : theMatrix[0].length;
		}

// Exported operations.

	/**
	 * Obtain the element at the given index in this data buffer.
	 *
	 * @param  i  Index.
	 *
	 * @return  Element value.
	 */
	public int getElem
		(int i)
		{
		return myMatrix [i / myWidth] [i % myWidth];
		}

	/**
	 * Obtain the element at the given bank and index in this data buffer.
	 *
	 * @param  bank  Bank.
	 * @param  i     Index.
	 *
	 * @return  Element value.
	 */
	public int getElem
		(int bank,
		 int i)
		{
		return myMatrix [i / myWidth] [i % myWidth];
		}

	/**
	 * Set the element at the given index in this data buffer.
	 *
	 * @param  i    Index.
	 * @param  val  Element value.
	 */
	public void setElem
		(int i,
		 int val)
		{
		myMatrix [i / myWidth] [i % myWidth] = val;
		}

	/**
	 * Set the element at the given bank and index in this data buffer.
	 *
	 * @param  bank  Bank.
	 * @param  i     Index.
	 * @param  val   Element value.
	 */
	public void setElem
		(int bank,
		 int i,
		 int val)
		{
		myMatrix [i / myWidth] [i % myWidth] = val;
		}

	}
