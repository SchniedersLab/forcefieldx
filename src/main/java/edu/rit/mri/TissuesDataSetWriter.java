//******************************************************************************
//
// File:    TissuesDataSetWriter.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.TissuesDataSetWriter
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
import java.io.RandomAccessFile;

/**
 * Class TissuesDataSetWriter provides an object that writes a magnetic
 * resonance image tissues data set to a file. The data set includes the results
 * of a spin relaxometry analysis for each pixel in the image, including number
 * of tissues, spin density of each tissue, and spin-lattice relaxation rate of
 * each tissue.
 * <P>
 * The tissues data set file format is as follows. Primitive types are written
 * as with java.io.DataOutput.
 * <UL>
 * <LI>
 * Image height <I>H</I> (int, 4 bytes).
 * <LI>
 * Image width <I>W</I> (int, 4 bytes). Number of pixels <I>P</I> =
 * <I>H</I>&times;<I>W</I>.
 * <LI>
 * <I>P</I> offsets (long, 8 bytes each). Each is the offset to the start of the
 * pixel tissues data for the corresponding pixel index. An offset of 0 means no
 * tissues data for that pixel index.
 * <LI>
 * Pixel tissues data for each pixel. Written using
 * <TT>PixelTissues.write()</TT>.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class TissuesDataSetWriter
	{

// Hidden data members.

	private RandomAccessFile myFile;
	private int P;
	private int M;

// Exported constructors.

	/**
	 * Construct a new tissues data set writer.
	 *
	 * @param  theFile   File to write.
	 * @param  H         Height (number of rows).
	 * @param  W         Width (number of columns).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFile</TT> is null. Thrown if
	 *     <TT>t_series</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>H</TT> &lt; 0 or <TT>W</TT> &lt;
	 *     0.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public TissuesDataSetWriter
		(File theFile,
		 int H,
		 int W)
		throws IOException
		{
		// Verify preconditions.
		if (theFile == null)
			{
			throw new NullPointerException
				("TissuesDataSetWriter(): theFile is null");
			}
		if (H < 0)
			{
			throw new IllegalArgumentException
				("TissuesDataSetWriter(): Height (= "+H+") < 0, illegal");
			}
		if (W < 0)
			{
			throw new IllegalArgumentException
				("TissuesDataSetWriter(): Width (= "+W+") < 0, illegal");
			}

		// Open and truncate file.
		myFile = new RandomAccessFile (theFile, "rw");
		myFile.setLength (0L);

		// Write height and width.
		myFile.writeInt (H);
		myFile.writeInt (W);

		// Clear pixel signal data offsets.
		P = H*W;
		for (int i = 0; i < P; ++ i)
			{
			myFile.writeLong (0L);
			}
		}

// Exported operations.

	/**
	 * Add the given pixel tissues data to this data set. The pixel tissues data
	 * consists of the pixel index, a list of spin densities, and a list of
	 * spin-lattice relaxation rates. The length of the lists gives the number
	 * of tissues.
	 *
	 * @param  tissues  {@linkplain PixelTissues} object, with a pixel index
	 *                  <TT>i</TT> in the range 0 &le; <TT>i</TT> &le;
	 *                  <I>P</I>&minus;1.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>tissues</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void addPixelTissues
		(PixelTissues tissues)
		throws IOException
		{
		// Verify preconditions.
		if (tissues == null)
			{
			throw new NullPointerException
				("TissuesDataSetWriter.addPixelTissues(): tissues is null");
			}
		int index = tissues.pixelIndex();
		if (0 > index || index >= P)
			{
			throw new IndexOutOfBoundsException
				("TissuesDataSetWriter.addPixelTissues(): Pixel index (= "+
				 index+") out of bounds");
			}

		// Record offset to pixel tissues data.
		long offset = myFile.length();
		myFile.seek (8L + index*8L);
		myFile.writeLong (offset);

		// Write pixel tissues data.
		myFile.seek (offset);
		tissues.write (myFile);
		}

	/**
	 * Close this data set.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void close()
		throws IOException
		{
		myFile.close();
		}

	/**
	 * Finalize this data set.
	 */
	protected void finalize()
		{
		try { close(); } catch (IOException exc) {}
		}

	}
