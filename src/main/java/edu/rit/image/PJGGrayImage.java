//******************************************************************************
//
// File:    PJGGrayImage.java
// Package: edu.rit.image
// Unit:    Class edu.rit.image.PJGGrayImage
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

import edu.rit.swing.Displayable;

import edu.rit.util.Range;

import java.awt.image.BufferedImage;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import java.util.Arrays;

/**
 * Class PJGGrayImage provides an 8-bit grayscale image that is read from or
 * written to a file in Parallel Java Graphics (PJG) format. For further
 * information about the PJG format, see class {@linkplain PJGImage}. The image
 * is layered on top of a byte matrix (type <TT>byte[][]</TT>). The height and
 * width of the image are equal to the number of rows and columns in the
 * underlying matrix.
 * <P>
 * To read and write the pixels of an image, use the <TT>getIntPixel()</TT>,
 * <TT>getPixel()</TT>, <TT>setIntPixel()</TT>, and <TT>setPixel()</TT> methods.
 * These methods represent a pixel value as an integer in the range 0 .. 255 or
 * as a floating point number in the range 0.0 .. 1.0. The image's
 * <I>interpretation</I> attribute specifies how to map between the numerical
 * pixel value and the shade of gray. If the interpretation is
 * <TT>ZERO_IS_WHITE</TT>, then a value of 0 or 0.0 represents white and 255 or
 * 1.0 represents black. If the interpretation is <TT>ZERO_IS_BLACK</TT>, then a
 * value of 0 or 0.0 represents black and 255 or 1.0 represents white. The
 * default interpretation is <TT>ZERO_IS_BLACK</TT>. To change the
 * interpretation, call the <TT>setInterpretation()</TT> method.
 * <P>
 * Changing the contents of the underlying matrix directly will also change the
 * image. In the underlying matrix, a value of 0 represents black and a value of
 * 255 represents white.
 * <P>
 * To write a PJGGrayImage object to a PJG image file, call the
 * <TT>prepareToWrite()</TT> method, specifying the output stream to write. The
 * <TT>prepareToWrite()</TT> method returns an instance of class {@linkplain
 * PJGImage.Writer}. Call the methods of the PJG image writer object to write
 * the pixel data, or sections of the pixel data, to the output stream. When
 * finished, close the PJG image writer.
 * <P>
 * When writing an image to a PJG file, class PJGGrayImage uses Huffman delta
 * encoded 8-bit grayscale pixel data segments (see class {@linkplain PJGImage})
 * to compress the pixel data. This technique will work well if the image has
 * many runs of consecutive identical or nearly identical pixel values.
 * <P>
 * To read a PJGGrayImage object from a PJG image file, call the
 * <TT>prepareToRead()</TT> method, specifying the input stream to read. The
 * <TT>prepareToRead()</TT> method returns an instance of class {@linkplain
 * PJGImage.Reader}. Call the methods of the PJG image reader object to read
 * the pixel data, or sections of the pixel data, from the input stream. When
 * finished, close the PJG image reader.
 * <P>
 * To get a BufferedImage object that uses the same underlying pixel data matrix
 * as the PJGGrayImage object, call the <TT>getBufferedImage()</TT> method. You
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
 * <I>Note:</I> Class PJGGrayImage is not multiple thread safe.
 * 
 * @author  Alan Kaminsky
 * @version 08-Apr-2008
 */
public class PJGGrayImage
	extends PJGImage
	{

// Exported enumerations.

	/**
	 * Enum Interpretation specifies how to interpret a numerical pixel value as
	 * a shade of gray in a grayscale image.
	 *
	 * @author  Alan Kaminsky
	 * @version 08-Nov-2007
	 */
	public static enum Interpretation
		{
		/**
		 * A pixel value of 0 or 0.0 is white, 255 or 1.0 is black.
		 */
		ZERO_IS_WHITE,

		/**
		 * A pixel value of 0 or 0.0 is black, 255 or 1.0 is white.
		 */
		ZERO_IS_BLACK
		}

	/**
	 * A pixel value of 0 or 0.0 is white, 255 or 1.0 is black.
	 */
	public static final Interpretation ZERO_IS_WHITE =
		Interpretation.ZERO_IS_WHITE;

	/**
	 * A pixel value of 0 or 0.0 is black, 255 or 1.0 is white.
	 */
	public static final Interpretation ZERO_IS_BLACK =
		Interpretation.ZERO_IS_BLACK;

// Hidden helper classes.

	/**
	 * Interface Transformation specifies the interface for an object that
	 * transforms a user-specified pixel value into the value stored in the
	 * underlying matrix.
	 *
	 * @author  Alan Kaminsky
	 * @version 08-Nov-2007
	 */
	static interface Transformation
		{
		/**
		 * Transform the given integer pixel value to a byte value to be stored
		 * in the underlying matrix.
		 *
		 * @param  value  Integer pixel value.
		 *
		 * @return  Byte value.
		 */
		public byte transformInt
			(int value);

		/**
		 * Transform the given float pixel value to a byte value to be stored in
		 * the underlying matrix.
		 *
		 * @param  value  Float pixel value.
		 *
		 * @return  Byte value.
		 */
		public byte transformFloat
			(float value);

		/**
		 * Transform the given byte value from the underlying matrix to an
		 * integer pixel value.
		 *
		 * @param  value  Byte value.
		 *
		 * @return  Integer pixel value.
		 */
		public int inverseTransformInt
			(byte value);

		/**
		 * Transform the given byte value from the underlying matrix to a float
		 * pixel value.
		 *
		 * @param  value  Byte value.
		 *
		 * @return  Float pixel value.
		 */
		public float inverseTransformFloat
			(byte value);
		}

	static final Transformation ZERO_IS_WHITE_TRANSFORMATION =
		new Transformation()
			{
			public byte transformInt
				(int value)
				{
				return (byte)
					(255 - Math.max (0, Math.min (value, 255)));
				}

			public byte transformFloat
				(float value)
				{
				return (byte)
					(255 - Math.max (0, Math.min ((int) (value*256.0f), 255)));
				}

			public int inverseTransformInt
				(byte value)
				{
				return 255 - (value & 0xFF);
				}

			public float inverseTransformFloat
				(byte value)
				{
				return (255 - (value & 0xFF)) / 256.0f;
				}
			};

	static final Transformation ZERO_IS_BLACK_TRANSFORMATION =
		new Transformation()
			{
			public byte transformInt
				(int value)
				{
				return (byte)
					(Math.max (0, Math.min (value, 255)));
				}

			public byte transformFloat
				(float value)
				{
				return (byte)
					(Math.max (0, Math.min ((int) (value*256.0f), 255)));
				}

			public int inverseTransformInt
				(byte value)
				{
				return value & 0xFF;
				}

			public float inverseTransformFloat
				(byte value)
				{
				return (value & 0xFF) / 256.0f;
				}
			};

// Hidden data members.

	byte[][] myMatrix;
	Transformation myTransformation = ZERO_IS_BLACK_TRANSFORMATION;

// Exported constructors.

	/**
	 * Construct a new PJG grayscale image. The image's height and width are
	 * uninitialized. Before accessing the image's pixels, specify the height
	 * and width by calling the <TT>setMatrix()</TT> method or by reading the
	 * image from an input stream. The image's interpretation is {@link
	 * #ZERO_IS_BLACK}.
	 */
	public PJGGrayImage()
		{
		super (IMAGE_TYPE_8_BIT_GRAY);
		}

	/**
	 * Construct a new PJG grayscale image with the given height, width, and
	 * underlying matrix. The image's interpretation is {@link #ZERO_IS_BLACK}.
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
	public PJGGrayImage
		(int theHeight,
		 int theWidth,
		 byte[][] theMatrix)
		{
		super (IMAGE_TYPE_8_BIT_GRAY);
		setMatrix (theHeight, theWidth, theMatrix);
		}

// Exported operations.

	/**
	 * Obtain this image's underlying matrix.
	 *
	 * @return  Underlying matrix, or null if none.
	 */
	public byte[][] getMatrix()
		{
		return myMatrix;
		}

	/**
	 * Set this image's height, width, and underlying matrix.
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
	public void setMatrix
		(int theHeight,
		 int theWidth,
		 byte[][] theMatrix)
		{
		setHeightAndWidth (theHeight, theWidth);
		if (theMatrix.length != myHeight)
			{
			throw new IllegalArgumentException
				("PJGGrayImage.setMatrix(): theMatrix.length (= " +
				 theMatrix.length + ") does not equal image height (= " +
				 myHeight + ")");
			}
		myMatrix = theMatrix;
		}

	/**
	 * Get this image's interpretation.
	 *
	 * @return  {@link #ZERO_IS_WHITE} or {@link #ZERO_IS_BLACK}.
	 */
	public Interpretation getInterpretation()
		{
		if (myTransformation == ZERO_IS_WHITE_TRANSFORMATION)
			{
			return ZERO_IS_WHITE;
			}
		else
			{
			return ZERO_IS_BLACK;
			}
		}

	/**
	 * Set this image's interpretation. If this method is not called, this
	 * image's default interpretation is {@link #ZERO_IS_BLACK}.
	 *
	 * @param  theInterpretation  {@link #ZERO_IS_WHITE} or
	 *                            {@link #ZERO_IS_BLACK}.
	 */
	public void setInterpretation
		(Interpretation theInterpretation)
		{
		switch (theInterpretation)
			{
			case ZERO_IS_WHITE:
				myTransformation = ZERO_IS_WHITE_TRANSFORMATION;
				break;
			case ZERO_IS_BLACK:
				myTransformation = ZERO_IS_BLACK_TRANSFORMATION;
				break;
			}
		}

	/**
	 * Obtain the pixel at the given row and column in this image.
	 *
	 * @param  r  Row index.
	 * @param  c  Column index.
	 *
	 * @return  Pixel[<I>r,c</I>] as an integer in the range 0 .. 255.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public int getIntPixel
		(int r,
		 int c)
		{
		return myTransformation.inverseTransformInt (myMatrix[r][c]);
		}

	/**
	 * Obtain the pixel at the given row and column in this image.
	 *
	 * @param  r  Row index.
	 * @param  c  Column index.
	 *
	 * @return  Pixel[<I>r,c</I>] as a floating point number in the range 0.0 ..
	 *          1.0.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public float getPixel
		(int r,
		 int c)
		{
		return myTransformation.inverseTransformFloat (myMatrix[r][c]);
		}

	/**
	 * Set the pixel at the given row and column in this image. If <TT>val</TT>
	 * is outside the range 0 .. 255, it is pinned to the appropriate boundary.
	 *
	 * @param  r    Row index.
	 * @param  c    Column index.
	 * @param  val  New pixel value as an integer in the range 0 .. 255.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setIntPixel
		(int r,
		 int c,
		 int val)
		{
		myMatrix[r][c] = myTransformation.transformInt (val);
		}

	/**
	 * Set the pixel at the given row and column in this image. If <TT>val</TT>
	 * is outside the range 0.0 .. 1.0, it is pinned to the appropriate
	 * boundary.
	 *
	 * @param  r    Row index.
	 * @param  c    Column index.
	 * @param  val  New pixel value as a floating point number in the range 0.0
	 *              .. 1.0.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>r</TT> or <TT>c</TT> is out of
	 *     bounds.
	 */
	public void setPixel
		(int r,
		 int c,
		 float val)
		{
		myMatrix[r][c] = myTransformation.transformFloat (val);
		}

	/**
	 * Set all pixels in this image to the given value. If <TT>val</TT> is
	 * outside the range 0 .. 255, it is pinned to the appropriate boundary.
	 *
	 * @param  val  New pixel value as an integer in the range 0 .. 255.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix has no storage
	 *     allocated.
	 */
	public void fill
		(int val)
		{
		byte byteval = myTransformation.transformInt (val);
		for (int r = 0; r < myHeight; ++ r)
			{
			Arrays.fill (myMatrix[r], byteval);
			}
		}

	/**
	 * Set all pixels in this image to the given value. If <TT>val</TT> is
	 * outside the range 0.0 .. 1.0, it is pinned to the appropriate boundary.
	 *
	 * @param  val  New pixel value as a floating point number in the range 0.0
	 *              .. 1.0.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if the underlying matrix has no storage
	 *     allocated.
	 */
	public void fill
		(float val)
		{
		byte byteval = myTransformation.transformFloat (val);
		for (int r = 0; r < myHeight; ++ r)
			{
			Arrays.fill (myMatrix[r], byteval);
			}
		}

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
		return new GrayWriter (theStream);
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
		return new GrayReader (theStream);
		}

	/**
	 * Obtain a BufferedImage whose pixel data comes from this image's
	 * underlying matrix.
	 *
	 * @return  BufferedImage.
	 */
	public BufferedImage getBufferedImage()
		{
		return new GrayBufferedImage (myHeight, myWidth, myMatrix);
		}

	/**
	 * Obtain a Displayable object with which to display this image in a Swing
	 * UI.
	 *
	 * @return  Displayable object.
	 */
	public Displayable getDisplayable()
		{
		return new GrayDisplayable (myHeight, myWidth, myMatrix);
		}

// Hidden operations.

	/**
	 * Write a pixel data segment to the given data output stream.
	 *
	 * @param  dos         Data output stream.
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
		 int row,
		 int col,
		 int rowlen,
		 int collen)
		throws IOException
		{
		// Write segment type, row index, column index, number of rows, number
		// of columns.
		dos.writeByte (SEGMENT_PIXEL_DATA_HDE_8_BIT_GRAY);
		dos.writeInt (row);
		dos.writeInt (col);
		dos.writeInt (rowlen);
		dos.writeInt (collen);

		// Prepare to write a bit stream.
		OutputBitStream obs = new OutputBitStream (dos);

		// Scan rows.
		int rowub = row + rowlen;
		int colub = col + collen;
		for (int r = row; r < rowub; ++ r)
			{
			byte[] matrix_r = myMatrix[r];

			// Scan columns and write Huffman encoded bit string.
			int prevPixel = 0;
			int pixel;
			int delta;
			for (int c = col; c < colub; ++ c)
				{
				// Get pixel and delta.
				pixel = matrix_r[c] & 0xFF;
				delta = pixel - prevPixel;
				prevPixel = pixel;

				// Encode bits and write them to output stream.
				if (delta < -42 || delta > +42)
					{
					obs.writeBits (0x0F00 | pixel, 12);
					}
				else if (delta < -10)
					{
					obs.writeBits (0x0380 | ((delta + 10) & 0x3F), 10);
					}
				else if (delta > +10)
					{
					obs.writeBits (0x0380 | ((delta - 11) & 0x3F), 10);
					}
				else if (delta < -2)
					{
					obs.writeBits (0x0060 | ((delta + 2) & 0x0F), 7);
					}
				else if (delta > +2)
					{
					obs.writeBits (0x0060 | ((delta - 3) & 0x0F), 7);
					}
				else if (delta < 0)
					{
					obs.writeBits (0x0008 | (delta & 0x03), 4);
					}
				else if (delta > 0)
					{
					obs.writeBits (0x0008 | ((delta - 1) & 0x03), 4);
					}
				else
					{
					obs.writeBits (0, 1);
					}
				}
			}

		// Dump any remaining bits to output stream.
		obs.flush();
		}

	/**
	 * Read a pixel data segment from the given data input stream. Assumes the
	 * segment type, row index, column index, row count, and column count have
	 * already been read.
	 *
	 * @param  dis      Data input stream.
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
		 int row,
		 int col,
		 int rowlen,
		 int collen)
		throws IOException
		{
		// Allocate storage for matrix if necessary.
		if (myMatrix == null) myMatrix = new byte [myHeight] [];

		// Prepare to read a bit stream.
		InputBitStream ibs = new InputBitStream (dis);

		// Scan rows.
		int rowub = row + rowlen;
		int colub = col + collen;
		for (int r = row; r < rowub; ++ r)
			{
			// Allocate storage for row if necessary.
			byte[] matrix_r = myMatrix[r];
			if (matrix_r == null)
				{
				matrix_r = new byte [colub];
				myMatrix[r] = matrix_r;
				}
			else if (matrix_r.length < colub)
				{
				byte[] new_matrix_r = new byte [colub];
				System.arraycopy
					(matrix_r, 0, new_matrix_r, 0, matrix_r.length);
				matrix_r = new_matrix_r;
				myMatrix[r] = new_matrix_r;
				}

			// Scan columns.
			int prevPixel = 0;
			int pixel;
			int delta;
			for (int c = col; c < colub; ++ c)
				{
				// Check high-order bits of Huffman code.
				if (ibs.peekBits (1) == 0)
					{
					// Delta = 0.
					ibs.skipBits (1);
					pixel = prevPixel;
					}

				else if (ibs.peekBits (2) == 0x02)
					{
					// Two-bit delta.
					ibs.skipBits (2);
					delta = ibs.readBits (2);
					if ((delta & 0x02) != 0)
						{
						delta = delta | 0xFFFFFFFC;
						}
					else
						{
						delta = delta + 1;
						}
					pixel = prevPixel + delta;
					}

				else if (ibs.peekBits (3) == 0x06)
					{
					// Four-bit delta.
					ibs.skipBits (3);
					delta = ibs.readBits (4);
					if ((delta & 0x08) != 0)
						{
						delta = (delta | 0xFFFFFFF0) - 2;
						}
					else
						{
						delta = delta + 3;
						}
					pixel = prevPixel + delta;
					}

				else if (ibs.peekBits (4) == 0x0E)
					{
					// Six-bit delta.
					ibs.skipBits (4);
					delta = ibs.readBits (6);
					if ((delta & 0x20) != 0)
						{
						delta = (delta | 0xFFFFFFC0) - 10;
						}
					else
						{
						delta = delta + 11;
						}
					pixel = prevPixel + delta;
					}

				else // (ibs.peekBits (4) == 0x0F)
					{
					// Eight-bit pixel.
					ibs.skipBits (4);
					pixel = ibs.readBits (8);
					}

				// Record pixel.
				matrix_r[c] = (byte) pixel;
				prevPixel = pixel;
				}
			}
		}

// Hidden helper classes.

	/**
	 * Class PJGGrayImage.GrayWriter provides an object with which to write a
	 * {@linkplain PJGGrayImage} to a file.
	 *
	 * @author  Alan Kaminsky
	 * @version 09-Nov-2007
	 */
	private class GrayWriter
		extends PJGImage.Writer
		{

	// Hidden constructors.

		/**
		 * Construct a new PJG gray image writer.
		 *
		 * @param  theStream  Output stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private GrayWriter
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
			writePixelData (myDos, 0, 0, myHeight, myWidth);
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
			writePixelData (myDos, rowlb, 0, rowlen, myWidth);
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
			writePixelData (myDos, 0, collb, myHeight, collen);
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
			writePixelData (myDos, rowlb, collb, rowlen, collen);
			}

		}

	/**
	 * Class PJGGrayImage.GrayReader provides an object with which to read a
	 * {@linkplain PJGGrayImage} from an input stream.
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
	 * <I>Note:</I> Class PJGGrayImage.GrayReader is not multiple thread safe.
	 *
	 * @author  Alan Kaminsky
	 * @version 02-Nov-2007
	 */
	private class GrayReader
		extends PJGImage.Reader
		{

	// Hidden data members.

		int myRow;
		int myCol;
		int myRowLen;
		int myColLen;

	// Hidden constructors.

		/**
		 * Construct a new PJG gray image reader.
		 *
		 * @param  theStream  Input stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		GrayReader
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
			readPixelData (myDis, myRow, myCol, myRowLen, myColLen);

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
			if (myNextSegmentType != SEGMENT_PIXEL_DATA_HDE_8_BIT_GRAY)
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
