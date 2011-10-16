//******************************************************************************
//
// File:    PJGColorImage.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.PJGColorImage
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

import edu.rit.util.Range;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Class PJGColorImage provides a color image that is read from or written to a
 * file in Parallel Java Graphics (PJG) format; class PJGColorImage is best
 * suited when the image has a small number of discrete colors. For further
 * information about the PJG format, see class {@linkplain PJGImage}. The image
 * is layered on top of an integer matrix (type <TT>int[][]</TT>). The height
 * and width of the image are equal to the number of rows and columns in the
 * underlying matrix.
 * <P>
 * To get and set the image's pixel data, use the <TT>getPixel()</TT>,
 * <TT>getPixelColor()</TT>, <TT>setPixel()</TT>, <TT>setPixelColor()</TT>, and
 * <TT>setPixelHSB()</TT> methods. You only need to allocate storage in the
 * pixel data matrix for the portions of the image you are actually accessing;
 * the complete matrix need not be allocated. Class {@linkplain
 * edu.rit.util.Arrays} has static methods for allocating portions of a matrix.
 * <P>
 * Changing the contents of the underlying matrix directly will also change the
 * image. The color information is stored in a matrix element as follows:
 * <UL>
 * <LI>Bits 31 .. 24 -- Unused, must be 0
 * <LI>Bits 23 .. 16 -- Red component in the range 0 .. 255
 * <LI>Bits 15 .. 8 -- Green component in the range 0 .. 255
 * <LI>Bits 7 .. 0 -- Blue component in the range 0 .. 255
 * </UL>
 * <P>
 * A color may be specified using hue, saturation, and brightness components
 * instead of red, green, and blue components.
 * <P>
 * The hue component gives the basic color. A hue of 0 = red; 1/6 = yellow; 2/6
 * = green; 3/6 = cyan; 4/6 = blue; 5/6 = magenta; 1 = red again. Intermediate
 * hue values yield intermediate colors.
 * <P>
 * The saturation component specifies how gray or colored the color is. A
 * saturation of 0 yields fully gray; a saturation of 1 yields fully colored.
 * Intermediate saturation values yield mixtures of gray and colored.
 * <P>
 * The brightness component specifies how dark or light the color is. A
 * brightness of 0 yields fully dark (black); a brightness of 1 yields fully
 * light (somewhere between white and colored depending on the saturation).
 * Intermediate brightness values yield somewhere between a gray shade and a
 * darkened color (depending on the saturation).
 * <P>
 * To write a PJGColorImage object to a PJG image file, call the
 * <TT>prepareToWrite()</TT> method, specifying the output stream to write. The
 * <TT>prepareToWrite()</TT> method returns an instance of class {@linkplain
 * PJGImage.Writer}. Call the methods of the PJG image writer object to write
 * the pixel data, or sections of the pixel data, to the output stream. When
 * finished, close the PJG image writer.
 * <P>
 * When writing an image to a PJG file, class PJGColorImage uses run length
 * encoded 24-bit color pixel data segments (see class {@linkplain PJGImage}) to
 * compress the pixel data. This technique will work well if the image has many
 * runs of consecutive identical pixels and the image has no more than 16,383
 * unique colors; it will work even better if the image has no more than 127
 * unique colors.
 * <P>
 * To read a PJGColorImage object from a PJG image file, call the
 * <TT>prepareToRead()</TT> method, specifying the input stream to read. The
 * <TT>prepareToRead()</TT> method returns an instance of class {@linkplain
 * PJGImage.Reader}. Call the methods of the PJG image reader object to read
 * the pixel data, or sections of the pixel data, from the input stream. When
 * finished, close the PJG image reader.
 * <P>
 * To get a BufferedImage object that uses the same underlying pixel data matrix
 * as the PJGColorImage object, call the <TT>getBufferedImage()</TT> method. You
 * can then do all the following with the BufferedImage: display it on the
 * screen, draw into it using a graphics context, copy another BufferedImage
 * into it, read it from or write it to a file using package javax.imageio
 * (which typically supports PNG, JPG, and GIF formats). The rows and columns of
 * the underlying matrix need not all be allocated when accessing the
 * BufferedImage. If you get a pixel from the BufferedImage in an unallocated
 * row or column, a pixel value of 0 (black) is returned. If you set a pixel in
 * the BufferedImage in an unallocated row or column, the pixel value is
 * discarded.
 * <P>
 * <I>Note:</I> Class PJGColorImage is not multiple thread safe.
 * 
 * @author  Alan Kaminsky
 * @version 08-Apr-2008
 */
public class PJGColorImage
	extends BaseColorImage
	{

// Exported constructors.

	/**
	 * Construct a new PJG color image. The image's height and width are
	 * uninitialized. Before accessing the image's pixels, specify the height
	 * and width by calling the <TT>setMatrix()</TT> method or by reading the
	 * image from an input stream.
	 */
	public PJGColorImage()
		{
		super (IMAGE_TYPE_24_BIT_COLOR);
		}

	/**
	 * Construct a new PJG color image with the given height, width, and
	 * underlying matrix.
	 *
	 * @param  theHeight  Image height in pixels.
	 * @param  theWidth   Image width in pixels.
	 * @param  theMatrix  Underlying matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theHeight</TT> &lt;= 0. Thrown if
	 *     <TT>theWidth</TT> &lt;= 0. Thrown if <TT>theMatrix.length</TT> does
	 *     not equal <TT>theHeight</TT>.
	 */
	public PJGColorImage
		(int theHeight,
		 int theWidth,
		 int[][] theMatrix)
		{
		super (IMAGE_TYPE_24_BIT_COLOR);
		setMatrix (theHeight, theWidth, theMatrix);
		}

// Exported operations.

	/**
	 * Prepare to write this image to the given output stream. Certain header
	 * information is written to the output stream at this time. To write this
	 * image's pixel data, call methods on the returned PJG image writer, then
	 * close the PJG image writer.
	 * <P>
	 * For improved performance, specify an output stream with buffering, such
	 * as an instance of class java.io.BufferedOutputStream.
	 *
	 * @param  theStream  Output stream.
	 *
	 * @return  PJG image writer object with which to write this image.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public PJGImage.Writer prepareToWrite
		(OutputStream theStream)
		throws IOException
		{
		return new ColorWriter (theStream);
		}

	/**
	 * Prepare to read this image from the given input stream. Certain header
	 * information is read from the input stream at this time. To read this
	 * image's pixel data, call methods on the returned PJG image reader, then
	 * close the PJG image reader.
	 * <P>
	 * For improved performance, specify an input stream with buffering, such
	 * as an instance of class java.io.BufferedInputStream.
	 *
	 * @param  theStream  Input stream.
	 *
	 * @return  PJG image reader object with which to read this image.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public PJGImage.Reader prepareToRead
		(InputStream theStream)
		throws IOException
		{
		return new ColorReader (theStream);
		}

// Hidden operations.

	/**
	 * Write a pixel data segment to the given data output stream.
	 *
	 * @param  dos         Data output stream.
	 * @param  dictionary  Dictionary.
	 * @param  row         Row index.
	 * @param  col         Column index.
	 * @param  rowlen      Number of rows.
	 * @param  collen      Number of columns.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void writePixelData
		(DataOutputStream dos,
		 HashMap<Integer,Integer> dictionary,
		 int row,
		 int col,
		 int rowlen,
		 int collen)
		throws IOException
		{
		// Write segment type, row index, column index, number of rows, number
		// of columns.
		dos.writeByte (SEGMENT_PIXEL_DATA_RLE_24_BIT_COLOR);
		dos.writeInt (row);
		dos.writeInt (col);
		dos.writeInt (rowlen);
		dos.writeInt (collen);

		// Scan rows.
		int rowub = row + rowlen;
		int colub = col + collen;
		for (int r = row; r < rowub; ++ r)
			{
			int[] matrix_r = myMatrix[r];

			// Scan columns and write runs of consecutive identical pixels.
			int runLength = 1;
			int runPixel = matrix_r[col];
			for (int c = col+1; c < colub; ++ c)
				{
				int pixel = matrix_r[c];
				if (pixel == runPixel)
					{
					++ runLength;
					}
				else
					{
					writeRun (dos, dictionary, runLength, runPixel);
					runLength = 1;
					runPixel = pixel;
					}
				}
			writeRun (dos, dictionary, runLength, runPixel);
			}
		}

	/**
	 * Write a run of consecutive identical pixels to the given data output
	 * stream.
	 *
	 * @param  dos         Data output stream.
	 * @param  dictionary  Dictionary.
	 * @param  runLength   Run length.
	 * @param  runPixel    Pixel value.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void writeRun
		(DataOutputStream dos,
		 HashMap<Integer,Integer> dictionary,
		 int runLength,
		 int runPixel)
		throws IOException
		{
		// Write one or more runs of run length at most 32768.
		while (runLength > 0)
			{
			// Write encoded run length.
			int n = Math.min (runLength-1, 32767);
			if (n <= 127)
				{
				dos.writeByte (n);
				}
			else
				{
				dos.writeByte (0x80 | (n >> 8));
				dos.writeByte (n);
				}

			// Write encoded pixel value.
			Integer index = dictionary.get (runPixel);
			if (index == null)
				{
				dos.writeByte (0xE0);
				dos.writeByte (runPixel >> 16);
				dos.writeByte (runPixel >>  8);
				dos.writeByte (runPixel);
				int x = dictionary.size();
				if (x < 2097152)
					{
					dictionary.put (runPixel, x);
					}
				}
			else
				{
				int x = index.intValue();
				if (x <= 127)
					{
					dos.writeByte (x);
					}
				else if (x <= 16383)
					{
					dos.writeByte (0x80 | (x >> 8));
					dos.writeByte (x);
					}
				else
					{
					dos.writeByte (0xC0 | (x >> 16));
					dos.writeByte (x >> 8);
					dos.writeByte (x);
					}
				}

			// Write further runs if necessary.
			runLength -= 32768;
			}
		}

	/**
	 * Read a pixel data segment from the given data input stream. Assumes the
	 * segment type, row index, column index, row count, and column count have
	 * already been read.
	 *
	 * @param  dis      Data input stream.
	 * @param  revdict  Reverse dictionary.
	 * @param  row      Row index.
	 * @param  col      Column index.
	 * @param  rowlen   Number of rows.
	 * @param  collen   Number of columns.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void readPixelData
		(DataInputStream dis,
		 ArrayList<Integer> revdict,
		 int row,
		 int col,
		 int rowlen,
		 int collen)
		throws IOException
		{
		// Allocate storage for matrix if necessary.
		if (myMatrix == null) myMatrix = new int [myHeight] [];

		// Scan rows.
		int rowub = row + rowlen;
		int colub = col + collen;
		for (int r = row; r < rowub; ++ r)
			{
			// Allocate storage for row if necessary.
			int[] matrix_r = myMatrix[r];
			if (matrix_r == null)
				{
				matrix_r = new int [colub];
				myMatrix[r] = matrix_r;
				}
			else if (matrix_r.length < colub)
				{
				int[] new_matrix_r = new int [colub];
				System.arraycopy
					(matrix_r, 0, new_matrix_r, 0, matrix_r.length);
				matrix_r = new_matrix_r;
				myMatrix[r] = new_matrix_r;
				}

			// Scan columns.
			int c = col;
			while (c < colub)
				{
				// Read and verify run length.
				int runLength = dis.readUnsignedByte();
				if (runLength >= 0x80)
					{
					runLength =
						((runLength & 0x7F)     << 8) |
						(dis.readUnsignedByte()     );
					}
				runLength += 1;
				if (c + runLength > colub)
					{
					throw new PJGImageFileFormatException
						("Invalid PJG pixel data run length (= " + runLength +
						 ")");
					}

				// Read and verify pixel value.
				int runPixel = dis.readUnsignedByte();
				if (runPixel > 0xE0)
					{
					throw new PJGImageFileFormatException
						("Invalid PJG pixel data value (= " + runPixel + ")");
					}
				else if (runPixel == 0xE0)
					{
					runPixel =
						(dis.readUnsignedByte() << 16) |
						(dis.readUnsignedByte() <<  8) |
						(dis.readUnsignedByte()      );
					if (revdict.size() < 2097152)
						{
						revdict.add (runPixel);
						}
					}
				else if (runPixel >= 0xC0)
					{
					runPixel =
						((runPixel & 0x1F)      << 16) |
						(dis.readUnsignedByte() <<  8) |
						(dis.readUnsignedByte()      );
					if (runPixel < revdict.size())
						{
						runPixel = revdict.get (runPixel);
						}
					else
						{
						throw new PJGImageFileFormatException
							("Invalid PJG pixel data index (= " + runPixel
							 + ")");
						}
					}
				else if (runPixel >= 0x80)
					{
					runPixel =
						((runPixel & 0x3F)      << 8) |
						(dis.readUnsignedByte()     );
					if (runPixel < revdict.size())
						{
						runPixel = revdict.get (runPixel);
						}
					else
						{
						throw new PJGImageFileFormatException
							("Invalid PJG pixel data index (= " + runPixel
							 + ")");
						}
					}
				else
					{
					if (runPixel < revdict.size())
						{
						runPixel = revdict.get (runPixel);
						}
					else
						{
						throw new PJGImageFileFormatException
							("Invalid PJG pixel data index (= " + runPixel
							 + ")");
						}
					}

				// Store pixel run.
				Arrays.fill (matrix_r, c, c + runLength, runPixel);
				c += runLength;
				}
			}
		}

// Hidden helper classes.

	/**
	 * Class PJGColorImage.ColorWriter provides an object with which to write a
	 * {@linkplain PJGColorImage} to a file.
	 *
	 * @author  Alan Kaminsky
	 * @version 02-Nov-2007
	 */
	private class ColorWriter
		extends PJGImage.Writer
		{

	// Hidden data members.

		private HashMap<Integer,Integer> myDictionary =
			new HashMap<Integer,Integer>();

	// Hidden constructors.

		/**
		 * Construct a new PJG color image writer.
		 *
		 * @param  theStream  Output stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private ColorWriter
			(OutputStream theStream)
			throws IOException
			{
			super (theStream);
			}

	// Exported operations.

		/**
		 * Write all rows and columns of the image to the output stream.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void write()
			throws IOException
			{
			writePixelData (myDos, myDictionary, 0, 0, myHeight, myWidth);
			}

		/**
		 * Write the given row slice of the image to the output stream. Pixels
		 * in the given range of rows and in all columns are written.
		 * <P>
		 * <I>Note:</I> <TT>theRowRange</TT>'s stride must be 1.
		 *
		 * @param  theRowRange  Range of pixel rows.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT> is null.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theRowRange</TT>
		 *     is outside the range 0 .. image height - 1.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writeRowSlice
			(Range theRowRange)
			throws IOException
			{
			int rowlb = theRowRange.lb();
			int rowub = theRowRange.ub();
			int rowlen = theRowRange.length();
			int rowstride = theRowRange.stride();
			if (0 > rowlb || rowub >= myHeight)
				{
				throw new IndexOutOfBoundsException
					("PJGImage.Writer.writeRowSlice(): Image row range = 0.." +
					 (myHeight-1) + ", theRowRange = " + theRowRange);
				}
			if (rowstride > 1)
				{
				throw new IllegalArgumentException
					("PJGImage.Writer.writeRowSlice(): theRowRange stride = " +
					 rowstride + " illegal");
				}
			writePixelData (myDos, myDictionary, rowlb, 0, rowlen, myWidth);
			}

		/**
		 * Write the given column slice of the image to the output stream.
		 * Pixels in all rows and in the given range of columns are written.
		 * <P>
		 * <I>Note:</I> <TT>theColRange</TT>'s stride must be 1.
		 *
		 * @param  theColRange  Range of pixel columns.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theColRange</TT>
		 *     is outside the range 0 .. image width - 1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writeColSlice
			(Range theColRange)
			throws IOException
			{
			int collb = theColRange.lb();
			int colub = theColRange.ub();
			int collen = theColRange.length();
			int colstride = theColRange.stride();
			if (0 > collb || colub >= myWidth)
				{
				throw new IndexOutOfBoundsException
					("PJGImage.Writer.writeColSlice(): Image column range = 0.." +
					 (myWidth-1) + ", theColRange = " + theColRange);
				}
			if (colstride > 1)
				{
				throw new IllegalArgumentException
					("PJGImage.Writer.writeColSlice(): theColRange stride = " +
					 colstride + " illegal");
				}
			writePixelData (myDos, myDictionary, 0, collb, myHeight, collen);
			}

		/**
		 * Write the given patch of the image to the output stream. Pixels in
		 * the given range of rows and in the given range of columns are
		 * written.
		 * <P>
		 * <I>Note:</I> <TT>theRowRange</TT>'s stride must be 1.
		 * <TT>theColRange</TT>'s stride must be 1.
		 *
		 * @param  theRowRange  Range of pixel rows.
		 * @param  theColRange  Range of pixel columns.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT> is null.
		 *     Thrown if <TT>theColRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT>'s stride is
		 *     greater than 1. Thrown if <TT>theColRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theRowRange</TT>
		 *     is outside the range 0 .. image height - 1. Thrown if any index
		 *     in <TT>theColRange</TT> is outside the range 0 .. image width -
		 *     1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writePatch
			(Range theRowRange,
			 Range theColRange)
			throws IOException
			{
			int rowlb = theRowRange.lb();
			int rowub = theRowRange.ub();
			int rowlen = theRowRange.length();
			int rowstride = theRowRange.stride();
			if (0 > rowlb || rowub >= myHeight)
				{
				throw new IndexOutOfBoundsException
					("PJGImage.Writer.writePatch(): Image row range = 0.." +
					 (myHeight-1) + ", theRowRange = " + theRowRange);
				}
			if (rowstride > 1)
				{
				throw new IllegalArgumentException
					("PJGImage.Writer.writePatch(): theRowRange stride = " +
					 rowstride + " illegal");
				}
			int collb = theColRange.lb();
			int colub = theColRange.ub();
			int collen = theColRange.length();
			int colstride = theColRange.stride();
			if (0 > collb || colub >= myWidth)
				{
				throw new IndexOutOfBoundsException
					("PJGImage.Writer.writePatch(): Image column range = 0.." +
					 (myWidth-1) + ", theColRange = " + theColRange);
				}
			if (colstride > 1)
				{
				throw new IllegalArgumentException
					("PJGImage.Writer.writePatch(): theColRange stride = " +
					 colstride + " illegal");
				}
			writePixelData (myDos, myDictionary, rowlb, collb, rowlen, collen);
			}

		}

	/**
	 * Class PJGColorImage.ColorReader provides an object with which to read a
	 * {@linkplain PJGColorImage} from an input stream.
	 * <P>
	 * When a PJG image reader is created, the header, image type, height,
	 * width, creation time, and comment segments are read from the input
	 * stream.
	 * <P>
	 * To read the pixel data segments one at a time, first call the
	 * <TT>getRowRange()</TT> and <TT>getColRange()</TT> methods to obtain the
	 * range of rows and columns in the next pixel data segment. At this point,
	 * allocate storage for the rows and columns in the underlying matrix if
	 * necessary. Then call the <TT>readSegment()</TT> method to read the actual
	 * pixel data. Repeat these steps if there are additional pixel data
	 * segments.
	 * <P>
	 * To read all the pixel data segments (or all the remaining pixel data
	 * segments), call the <TT>read()</TT> method.
	 * <P>
	 * When finished, call the <TT>close()</TT> method.
	 * <P>
	 * <I>Note:</I> Class PJGImage.Reader is not multiple thread safe.
	 *
	 * @author  Alan Kaminsky
	 * @version 02-Nov-2007
	 */
	private class ColorReader
		extends PJGImage.Reader
		{

	// Hidden data members.

		ArrayList<Integer> myReverseDictionary = new ArrayList<Integer>();
		int myRow;
		int myCol;
		int myRowLen;
		int myColLen;

	// Hidden constructors.

		/**
		 * Construct a new PJG color image reader.
		 *
		 * @param  theStream  Input stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		ColorReader
			(InputStream theStream)
			throws IOException
			{
			super (theStream);
			getPixelDataSegmentParameters();
			}

	// Exported operations.

		/**
		 * Read all pixel data segments from the input stream. If some pixel
		 * data segments have already been read, the <TT>read()</TT> method
		 * reads all remaining pixel data segments. If there are no more pixel
		 * data segments, the <TT>read()</TT> method does nothing. If storage is
		 * not already allocated in the underlying matrix for the pixel rows and
		 * columns in a pixel data segment, the <TT>read()</TT> method allocates
		 * the necessary storage.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void read()
			throws IOException
			{
			while (myNextSegmentType != -1) readSegment();
			}

		/**
		 * Obtain the row range of the next pixel data segment in the input
		 * stream. If there are no more pixel data segments, null is returned.
		 *
		 * @return  Row range, or null.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public Range getRowRange()
			{
			return
				myNextSegmentType == -1 ?
					null :
					new Range (myRow, myRow + myRowLen - 1);
			}

		/**
		 * Obtain the column range of the next pixel data segment in the input
		 * stream. If there are no more pixel data segments, null is returned.
		 *
		 * @return  Column range, or null.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public Range getColRange()
			{
			return
				myNextSegmentType == -1 ?
					null :
					new Range (myCol, myCol + myColLen - 1);
			}

		/**
		 * Read the next pixel data segment from the input stream. If there are
		 * no more pixel data segments, the <TT>readSegment()</TT> method does
		 * nothing. If storage is not already allocated in the underlying matrix
		 * for the pixel rows and columns in the pixel data segment, the
		 * <TT>readSegment()</TT> method allocates the necessary storage.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void readSegment()
			throws IOException
			{
			// Early return if end of stream.
			if (myNextSegmentType == -1) return;

			// Read pixel data.
			readPixelData
				(myDis, myReverseDictionary,
				 myRow, myCol, myRowLen, myColLen);

			// Prepare for next segment.
			myNextSegmentType = myDis.read();
			getPixelDataSegmentParameters();
			}

	// Hidden operations.

		/**
		 * Get the row index, column index, row count, and column count for the
		 * next pixel data segment. Store the results in myRowRange and
		 * myColRange. Assumes the segment type has been read into
		 * myNextSegmentType.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private void getPixelDataSegmentParameters()
			throws IOException
			{
			// Early return if end of stream.
			if (myNextSegmentType == -1) return;

			// Verify segment type.
			if (myNextSegmentType != SEGMENT_PIXEL_DATA_RLE_24_BIT_COLOR)
				{
				throw new PJGImageFileFormatException
					("Invalid PJG pixel data segment type (= " +
					 myNextSegmentType + ")");
				}

			// Read and verify row index.
			myRow = myDis.readInt();
			if (0 > myRow || myRow >= myHeight)
				{
				throw new PJGImageFileFormatException
					("Invalid PJG pixel data segment row index (= " +
					 myRow + ")");
				}

			// Read and verify column index.
			myCol = myDis.readInt();
			if (0 > myCol || myCol >= myWidth)
				{
				throw new PJGImageFileFormatException
					("Invalid PJG pixel data segment column index (= "
					 + myCol + ")");
				}

			// Read and verify row count.
			myRowLen = myDis.readInt();
			if (1 > myRowLen || myRow + myRowLen > myHeight)
				{
				throw new PJGImageFileFormatException
					("Invalid PJG pixel data segment row count (= " +
					 myRowLen + ")");
				}

			// Read and verify column count.
			myColLen = myDis.readInt();
			if (1 > myColLen || myCol + myColLen > myWidth)
				{
				throw new PJGImageFileFormatException
					("Invalid PJG pixel data segment column count (= " +
					 myColLen + ")");
				}
			}

		}

	}
