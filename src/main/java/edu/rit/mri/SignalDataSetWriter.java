//******************************************************************************
//
// File:    SignalDataSetWriter.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.SignalDataSetWriter
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

import edu.rit.numeric.Series;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * Class SignalDataSetWriter provides an object that writes a magnetic resonance
 * image spin signal data set to a file. The data set includes the measured spin
 * signals for each pixel in the image.
 * <P>
 * The spin signal data set file format is as follows. Primitive types are
 * written as with java.io.DataOutput.
 * <UL>
 * <LI>
 * Image height <I>H</I> (int, 4 bytes).
 * <LI>
 * Image width <I>W</I> (int, 4 bytes). Number of pixels <I>P</I> =
 * <I>H</I>&times;<I>W</I>.
 * <LI>
 * <I>P</I> offsets (long, 8 bytes each). Each is the offset to the start of the
 * pixel signal data for the corresponding pixel index. An offset of 0 means no
 * signal data for that pixel index.
 * <LI>
 * Length of time series <I>M</I> (short, 2 bytes).
 * <LI>
 * <I>M</I> time values (double, 8 bytes each).
 * <LI>
 * Pixel signal data for each pixel. Written using <TT>PixelSignal.write()</TT>.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class SignalDataSetWriter
	{

// Hidden data members.

	private RandomAccessFile myFile;
	private int P;
	private int M;

// Exported constructors.

	/**
	 * Construct a new spin signal data set writer.
	 *
	 * @param  theFile   File to write.
	 * @param  H         Height (number of rows).
	 * @param  W         Width (number of columns).
	 * @param  t_series  Time series.
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
	public SignalDataSetWriter
		(File theFile,
		 int H,
		 int W,
		 Series t_series)
		throws IOException
		{
		// Verify preconditions.
		if (theFile == null)
			{
			throw new NullPointerException
				("SignalDataSetWriter(): theFile is null");
			}
		if (H < 0)
			{
			throw new IllegalArgumentException
				("SignalDataSetWriter(): Height (= "+H+") < 0, illegal");
			}
		if (W < 0)
			{
			throw new IllegalArgumentException
				("SignalDataSetWriter(): Width (= "+W+") < 0, illegal");
			}
		if (t_series == null)
			{
			throw new NullPointerException
				("SignalDataSetWriter(): t_series is null");
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

		// Write time series.
		M = t_series.length();
		myFile.writeShort ((short) M);
		for (int i = 0; i < M; ++ i)
			{
			myFile.writeDouble (t_series.x(i));
			}
		}

// Exported operations.

	/**
	 * Add the given pixel signal data to this data set. The pixel signal data
	 * consists of the pixel index and the signal series. The signal series is a
	 * series of length <I>M</I> giving the spin signal measurements
	 * <I>S</I>(<I>t</I><SUB><I>i</I></SUB>), where the
	 * <I>t</I><SUB><I>i</I></SUB> are the elements of the time series.
	 *
	 * @param  signal  {@linkplain PixelSignal} object, with a pixel index
	 *                 <TT>i</TT> in the range 0 &le; <TT>i</TT> &le;
	 *                 <I>P</I>&minus;1.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>signal</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the length of the signal series &ne;
	 *     <I>M</I>, the length of the time series.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void addPixelSignal
		(PixelSignal signal)
		throws IOException
		{
		// Verify preconditions.
		if (signal == null)
			{
			throw new NullPointerException
				("SpinSignalDataSet.addPixelSignal(): signal is null");
			}
		int index = signal.pixelIndex();
		if (0 > index || index >= P)
			{
			throw new IndexOutOfBoundsException
				("SpinSignalDataSet.addPixelSignal(): Pixel index (= "+index+
				 ") out of bounds");
			}
		if (signal.length() != M)
			{
			throw new IllegalArgumentException
				("SpinSignalDataSet.addPixelSignal(): Signal series length (= "+
				 signal.length()+") != time series length (= "+M+"), illegal");
			}

		// Record offset to pixel signal data.
		long offset = myFile.length();
		myFile.seek (8L + index*8L);
		myFile.writeLong (offset);

		// Write pixel signal data.
		myFile.seek (offset);
		signal.write (myFile);
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
