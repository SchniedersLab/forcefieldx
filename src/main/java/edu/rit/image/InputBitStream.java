//******************************************************************************
//
// File:    InputBitStream.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.InputBitStream
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

import java.io.DataInputStream;
import java.io.IOException;

/**
 * Class InputBitStream provides an object for reading a stream of bits from an
 * input stream.
 *
 * @author  Alan Kaminsky
 * @version 08-Apr-2008
 */
class InputBitStream
	{

// Hidden data members.

	private DataInputStream dis;
	private int inData;
	private int inBits;

// Exported constructors.

	/**
	 * Construct a new input bit stream on top of the given data input
	 * stream.
	 *
	 * @param  dis  Data input stream.
	 */
	public InputBitStream
		(DataInputStream dis)
		{
		this.dis = dis;
		}

// Exported operations.

	/**
	 * Peek ahead at the given number of bits, but do not consume them.
	 *
	 * @param  bits  Number of bits to peek.
	 *
	 * @return  Bits peeked, at the least significant end of the int.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public int peekBits
		(int bits)
		throws IOException
		{
		// Make sure we've read enough bits.
		while (inBits < bits)
			{
			inData = (inData << 8) | dis.readUnsignedByte();
			inBits += 8;
			}

		// Shift the given number of bits to the right end and mask them in.
		return (inData >> (inBits - bits)) & (~ (0xFFFFFFFF << bits));
		}

	/**
	 * Read and consume the given number of bits.
	 *
	 * @param  bits  Number of bits to read.
	 *
	 * @return  Bits read, at the least significant end of the int.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public int readBits
		(int bits)
		throws IOException
		{
		// Make sure we've read enough bits.
		while (inBits < bits)
			{
			inData = (inData << 8) | dis.readUnsignedByte();
			inBits += 8;
			}

		// Shift the given number of bits to the right end and mask them in.
		int b = (inData >> (inBits - bits)) & (~ (0xFFFFFFFF << bits));

		// Consume the given number of bits.
		inBits -= bits;
		return b;
		}

	/**
	 * Skip over and consume the given number of bits.
	 *
	 * @param  bits  Number of bits to skip.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void skipBits
		(int bits)
		throws IOException
		{
		// Make sure we've read enough bits.
		while (inBits < bits)
			{
			inData = (inData << 8) | dis.readUnsignedByte();
			inBits += 8;
			}

		// Consume the given number of bits.
		inBits -= bits;
		}

	}
