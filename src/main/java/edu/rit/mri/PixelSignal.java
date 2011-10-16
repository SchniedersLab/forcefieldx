//******************************************************************************
//
// File:    PixelSignal.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.PixelSignal
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

import java.io.DataInput;
import java.io.DataOutput;
import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class PixelSignal encapsulates the measured spin signal on one pixel of a
 * magnetic resonance image. This includes:
 * <UL>
 * <LI>
 * File index.
 * <LI>
 * Pixel index.
 * <LI>
 * Pixel's measured spin signal values, <I>S</I><SUB><I>i</I></SUB>.
 * </UL>
 * <P>
 * Each spin signal value is stored as a value of type <TT>short</TT> in the
 * range &minus;32768..+32767.
 * <P>
 * Operations are provided to read a pixel signal object from a DataInputStream
 * and write a pixel signal object to a DataOutputStream. The format is:
 * <UL>
 * <LI>
 * Number of spin signal values, <I>M</I> (2-byte short).
 * <LI>
 * <I>M</I> spin signal values, each a 2-byte short.
 * </UL>
 * <P>
 * Operations are provided to read a pixel signal object from an
 * ObjectInputStream and write a pixel signal object to an ObjectOutputStream.
 * The format is:
 * <UL>
 * <LI>
 * File index (4-byte int).
 * <LI>
 * Pixel index (4-byte int).
 * <LI>
 * Number of spin signal values, <I>M</I> (2-byte short).
 * <LI>
 * <I>M</I> spin signal values, each a 2-byte short.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 25-Jun-2008
 */
public class PixelSignal
	implements Externalizable
	{

// Hidden data members.

	private static final long serialVersionUID = 3137091490421051491L;

	// File index.
	private int fileIndex;

	// Pixel index.
	private int pixelIndex;

	// Array of measured spin signal values.
	private short[] S;

// Exported constructors.

	/**
	 * Construct a new, uninitialized pixel signal object. This constructor is
	 * for use only by object deserialization.
	 */
	public PixelSignal()
		{
		}

	/**
	 * Construct a new pixel signal object. If there was no measured spin signal
	 * for this pixel, <TT>S_array</TT> is null.
	 *
	 * @param  fileIndex   File index.
	 * @param  pixelIndex  Pixel index.
	 * @param  S_array     Array containing measured spin signal values, or
	 *                     null.
	 */
	public PixelSignal
		(int fileIndex,
		 int pixelIndex,
		 short[] S_array)
		{
		fileIndex (fileIndex);
		pixelIndex (pixelIndex);
		S_measured (S_array);
		}

	/**
	 * Construct a new pixel signal object. If there was no measured spin signal
	 * for this pixel, <TT>S_series</TT> is null. Each value in
	 * <TT>S_series</TT> is truncated to type <TT>short</TT> before being
	 * stored; values outside the range &minus;32768..+32767 are pinned to the
	 * nearest boundary of that range.
	 *
	 * @param  fileIndex   File index.
	 * @param  pixelIndex  Pixel index.
	 * @param  S_series    Series containing measured spin signal values, or
	 *                     null.
	 */
	public PixelSignal
		(int fileIndex,
		 int pixelIndex,
		 Series S_series)
		{
		fileIndex (fileIndex);
		pixelIndex (pixelIndex);
		S_measured (S_series);
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
	 * Specify the file index.
	 *
	 * @param  index  File index.
	 */
	public void fileIndex
		(int index)
		{
		fileIndex = index;
		}

	/**
	 * Get the pixel index.
	 *
	 * @return  Pixel index.
	 */
	public int pixelIndex()
		{
		return pixelIndex;
		}

	/**
	 * Specify the pixel index.
	 *
	 * @param  index  Pixel index.
	 */
	public void pixelIndex
		(int index)
		{
		pixelIndex = index;
		}

	/**
	 * Get the number of measured spin signal values for this pixel. If there
	 * was no measured spin signal for this pixel, 0 is returned.
	 *
	 * @return  Number of measured spin signal values, or 0.
	 */
	public int length()
		{
		return S == null ? 0 : S.length;
		}

	/**
	 * Get a series containing this pixel's measured spin signal values. If
	 * there was no measured spin signal for this pixel, null is returned.
	 *
	 * @return  Series containing measured spin signal values, or null.
	 */
	public Series S_measured()
		{
		return S == null ? null : new Series()
			{
			public int length()
				{
				return S.length;
				}
			public double x (int i)
				{
				return S[i];
				}
			};
		}

	/**
	 * Specify an array containing this pixel's measured spin signal values. If
	 * there was no measured spin signal for this pixel, <TT>S_array</TT> is
	 * null.
	 *
	 * @param  S_array  Array containing measured spin signal values, or null.
	 */
	public void S_measured
		(short[] S_array)
		{
		if (S_array != null)
			{
			this.S = (short[]) S_array.clone();
			}
		else
			{
			this.S = null;
			}
		}

	/**
	 * Specify a series containing this pixel's measured spin signal values. If
	 * there was no measured spin signal for this pixel, <TT>S_series</TT> is
	 * null. Each value in <TT>S_series</TT> is truncated to type <TT>short</TT>
	 * before being stored; values outside the range &minus;32768..+32767 are
	 * pinned to the nearest boundary of that range.
	 *
	 * @param  S_series  Series containing measured spin signal values, or null.
	 */
	public void S_measured
		(Series S_series)
		{
		if (S_series != null)
			{
			int M = S_series.length();
			this.S = new short [M];
			for (int i = 0; i < M; ++ i)
				{
				this.S[i] = (short)
					Math.max (-32768.0, Math.min (S_series.x(i), 32767.0));
				}
			}
		else
			{
			this.S = null;
			}
		}

	/**
	 * Write this pixel signal object to the given data output stream.
	 * <P>
	 * <I>Note:</I> The file index and pixel index are not written.
	 *
	 * @param  out  Data output stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void write
		(DataOutput out)
		throws IOException
		{
		int M = S == null ? 0 : S.length;
		out.writeShort ((short) M);
		for (int i = 0; i < M; ++ i)
			{
			out.writeShort (S[i]);
			}
		}

	/**
	 * Read this pixel signal object from the given data input stream.
	 * <P>
	 * <I>Note:</I> The file index and pixel index are not read.
	 *
	 * @param  in  Data input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void read
		(DataInput in)
		throws IOException
		{
		int M = in.readShort();
		if (M == 0)
			{
			S = null;
			}
		else
			{
			S = new short [M];
			for (int i = 0; i < M; ++ i)
				{
				S[i] = in.readShort();
				}
			}
		}

	/**
	 * Write this pixel signal object to the given object output stream.
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
		write (out);
		}

	/**
	 * Read this pixel signal object from the given object input stream.
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
		read (in);
		}

	}
