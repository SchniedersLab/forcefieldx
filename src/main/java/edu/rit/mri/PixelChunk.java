//******************************************************************************
//
// File:    PixelChunk.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.PixelChunk
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

package edu.rit.mri;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class PixelChunk encapsulates a chunk of pixels on which to perform a spin
 * relaxometry analysis.
 *
 * @author  Alan Kaminsky
 * @version 25-Jun-2008
 */
public class PixelChunk
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 1506538883148760846L;

	private int fileIndex;
	private int pixelIndex;
	private int pixelCount;

// Exported constructors.

	/**
	 * Construct a new, uninitialized pixel chunk object. This constructor is
	 * for use only by object deserialization.
	 */
	public PixelChunk()
		{
		}

	/**
	 * Construct a new pixel chunk object.
	 *
	 * @param  fileIndex   File index.
	 * @param  pixelIndex  Pixel index of first pixel to analyze.
	 * @param  pixelCount  Number of pixels to analyze.
	 */
	public PixelChunk
		(int fileIndex,
		 int pixelIndex,
		 int pixelCount)
		{
		this.fileIndex = fileIndex;
		this.pixelIndex = pixelIndex;
		this.pixelCount = pixelCount;
		}

// Exported operations.

	/**
	 * Get the file index.
	 *
	 * @return  File index.
	 */
	public int fileIndex()
		{
		return fileIndex;
		}

	/**
	 * Get the pixel index of the first pixel to analyze.
	 *
	 * @return  Pixel index.
	 */
	public int pixelIndex()
		{
		return pixelIndex;
		}

	/**
	 * Get the number of pixels to analyze.
	 *
	 * @return  Pixel count.
	 */
	public int pixelCount()
		{
		return pixelCount;
		}

	/**
	 * Write this pixel chunk object to the given object output stream.
	 *
	 * @param  out  Object output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void writeExternal
		(ObjectOutput out)
		throws IOException
		{
		out.writeInt (fileIndex);
		out.writeInt (pixelIndex);
		out.writeInt (pixelCount);
		}

	/**
	 * Read this pixel chunk object from the given object input stream.
	 *
	 * @param  in  Object input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void readExternal
		(ObjectInput in)
		throws IOException
		{
		fileIndex = in.readInt();
		pixelIndex = in.readInt();
		pixelCount = in.readInt();
		}

	}
