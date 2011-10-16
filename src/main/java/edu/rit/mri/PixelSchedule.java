//******************************************************************************
//
// File:    PixelSchedule.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.PixelSchedule
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

import java.io.File;
import java.io.IOException;

/**
 * Class PixelSchedule provides a load balancing schedule object used by class
 * {@linkplain SpinRelaxometryClu} for partitioning the pixels of the input spin
 * signal data sets among the workers. The data sets are divided into chunks of
 * a fixed size.
 *
 * @author  Alan Kaminsky
 * @version 25-Jun-2008
 */
public class PixelSchedule
	{

// Hidden data members.

	// Chunk size.
	private int chunkSize;

	// Spin signal data set file names.
	private String[] signalFileName;

	// Index of current file.
	private int fileIndex;

	// Index of current pixel.
	private int pixelIndex;

	// Number of pixels remaining in current file.
	private int pixelsRemaining;

// Exported constructors.

	/**
	 * Construct a new pixel schedule object.
	 *
	 * @param  chunkSize       Maximum number of pixels in each chunk.
	 * @param  signalFileName  Array of spin signal data set file names.
	 */
	public PixelSchedule
		(int chunkSize,
		 String[] signalFileName)
		{
		this.chunkSize = chunkSize;
		this.signalFileName = signalFileName;
		this.fileIndex = -1;
		this.pixelIndex = 0;
		this.pixelsRemaining = 0;
		}

// Exported operations.

	/**
	 * Get the next chunk of pixels to analyze. If there are no more chunks,
	 * null is returned.
	 *
	 * @return  Chunk, or null.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public PixelChunk next()
		throws IOException
		{
		// Finished with current file?
		if (pixelsRemaining == 0)
			{
			// Yes. Go to next file.
			fileIndex = Math.min (fileIndex+1, signalFileName.length);

			// More files?
			if (fileIndex == signalFileName.length)
				{
				// No more files. 
				return null;
				}

			// Get number of pixels in current file.
			SignalDataSetReader reader =
				new SignalDataSetReader
					(new File (signalFileName[fileIndex]));
			pixelIndex = 0;
			pixelsRemaining = reader.getPixelCount();
			reader.close();
			}

		// Carve next chunk out of current file.
		int n = Math.min (chunkSize, pixelsRemaining);
		PixelChunk chunk = new PixelChunk (fileIndex, pixelIndex, n);
		pixelIndex += n;
		pixelsRemaining -= n;
		return chunk;
		}

	}
