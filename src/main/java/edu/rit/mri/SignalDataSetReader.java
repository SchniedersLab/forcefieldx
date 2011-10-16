//******************************************************************************
//
// File:    SignalDataSetReader.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SignalDataSetReader
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

import edu.rit.numeric.ArraySeries;
import edu.rit.numeric.Series;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * Class SignalDataSetReader provides an object that reads a magnetic resonance
 * image spin signal data set from a file. The data set includes the measured
 * spin signals for each pixel in the image.
 * <P>
 * The spin signal data set file format is defined in class {@linkplain
 * SignalDataSetWriter}.
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class SignalDataSetReader
	{

// Hidden data members.

	private RandomAccessFile myFile;
	private int H;
	private int W;
	private int P;
	private int M;
	private Series t_series;

// Exported constructors.

	/**
	 * Construct a new spin signal data set reader.
	 *
	 * @param  theFile   File to read.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFile</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public SignalDataSetReader
		(File theFile)
		throws IOException
		{
		// Verify preconditions.
		if (theFile == null)
			{
			throw new NullPointerException
				("SignalDataSetReader(): theFile is null");
			}

		// Open file for read-only.
		myFile = new RandomAccessFile (theFile, "r");

		// Read height and width.
		H = myFile.readInt();
		W = myFile.readInt();

		// Skip pixel signal data offsets.
		P = H*W;
		myFile.seek (8L + P*8L);

		// Read time series.
		M = myFile.readShort();
		double[] t_array = new double [M];
		for (int i = 0; i < M; ++ i)
			{
			t_array[i] = myFile.readDouble();
			}
		t_series = new ArraySeries (t_array);
		}

// Exported operations.

	/**
	 * Get the height <I>H</I> of this data set's magnetic resonance image.
	 *
	 * @return  Height (number of rows).
	 */
	public int getHeight()
		{
		return H;
		}

	/**
	 * Get the width <I>W</I> of this data set's magnetic resonance image.
	 *
	 * @return  Width (number of columns).
	 */
	public int getWidth()
		{
		return W;
		}

	/**
	 * Get the number of pixels <I>P</I> in this data set's magnetic resonance
	 * image.
	 *
	 * @return  Number of pixels.
	 */
	public int getPixelCount()
		{
		return P;
		}

	/**
	 * Get the pixel index for the given row and column indexes.
	 *
	 * @param  r  Row index, 0 &le; <TT>r</TT> &le; <I>H</I>&minus;1.
	 * @param  c  Column index, 0 &le; <TT>c</TT> &le; <I>W</I>&minus;1.
	 *
	 * @return  Pixel index.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public int indexFor
		(int r,
		 int c)
		{
		if (0 > r || r >= H)
			{
			throw new IndexOutOfBoundsException
				("SignalDataSetReader.indexFor(): r (= "+r+") out of bounds");
			}
		if (0 > c || c >= W)
			{
			throw new IndexOutOfBoundsException
				("SignalDataSetReader.indexFor(): c (= "+c+") out of bounds");
			}
		return r*W + c;
		}

	/**
	 * Get the row index corresponding to the given pixel index.
	 *
	 * @param  i  Pixel index, 0 &le; <TT>i</TT> &le; <I>P</I>&minus;1.
	 *
	 * @return  Row index.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int rowFor
		(int i)
		{
		if (0 > i || i >= P)
			{
			throw new IndexOutOfBoundsException
				("SignalDataSetReader.rowFor(): i (= "+i+") out of bounds");
			}
		return i/W;
		}

	/**
	 * Get the column index corresponding to the given pixel index.
	 *
	 * @param  i  Pixel index, 0 &le; <TT>i</TT> &le; <I>P</I>&minus;1.
	 *
	 * @return  Column index.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public int columnFor
		(int i)
		{
		if (0 > i || i >= P)
			{
			throw new IndexOutOfBoundsException
				("SignalDataSetReader.columnFor(): i (= "+i+") out of bounds");
			}
		return i%W;
		}

	/**
	 * Get the time series for this data set.
	 *
	 * @return  Time series.
	 */
	public Series getTimeSeries()
		{
		return t_series;
		}

	/**
	 * Get the pixel signal data for the given pixel in this data set. If the
	 * given pixel has no associated signal data, null is returned.
	 *
	 * @param  i  Pixel index, 0 &le; <TT>i</TT> &le; <I>P</I>&minus;1.
	 *
	 * @return  {@linkplain PixelSignal} object for pixel <TT>i</TT>, or null.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public PixelSignal getPixelSignal
		(int i)
		throws IOException
		{
		// Verify preconditions.
		if (0 > i || i >= P)
			{
			throw new IndexOutOfBoundsException
				("SignalDataSetReader.getPixelSignal(): i (= "+i+
				 ") out of bounds");
			}

		// Read offset to pixel signal data.
		myFile.seek (8L + i*8L);
		long offset = myFile.readLong();
		if (offset == 0L) return null;

		// Read pixel signal data.
		myFile.seek (offset);
		PixelSignal signal = new PixelSignal();
		signal.pixelIndex (i);
		signal.read (myFile);

		return signal;
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
