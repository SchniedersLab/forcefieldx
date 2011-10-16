//******************************************************************************
//
// File:    OutputBitStream.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.OutputBitStream
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

import java.io.DataOutputStream;
import java.io.IOException;

/**
 * Class OutputBitStream provides an object for writing a stream of bits to an
 * output stream.
 *
 * @author  Alan Kaminsky
 * @version 08-Apr-2008
 */
class OutputBitStream
	{

// Hidden data members.

	private DataOutputStream dos;
	private int outData;
	private int outBits;

// Exported constructors.

	/**
	 * Construct a new output bit stream on top of the given data output
	 * stream.
	 *
	 * @param  dos  Data output stream.
	 */
	public OutputBitStream
		(DataOutputStream dos)
		{
		this.dos = dos;
		}

// Exported operations.

	/**
	 * Write the given number of least significant bits of the given value.
	 *
	 * @param  data  Data value.
	 * @param  bits  Number of bits.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeBits
		(int data,
		 int bits)
		throws IOException
		{
		outData = (outData << bits) | data;
		outBits += bits;
		while (outBits >= 8)
			{
			dos.writeByte (outData >> (outBits - 8));
			outBits -= 8;
			}
		}

	/**
	 * Flush any remaining accumulated bits.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void flush()
		throws IOException
		{
		while (outBits >= 8)
			{
			dos.writeByte (outData >> (outBits - 8));
			outBits -= 8;
			}
		if (outBits > 0)
			{
			dos.writeByte (outData << (8 - outBits));
			}
		}

	}
