//******************************************************************************
//
// File:    DoubleMatrixFile.java
// Package: edu.rit.io
// Unit:    Class edu.rit.io.DoubleMatrixFile
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

package edu.rit.io;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import edu.rit.util.Range;

/**
 * Class DoubleMatrixFile provides an object for reading or writing a double
 * matrix from or to a file. The matrix containing the data to read or write is
 * a separate object, specified as an argument to the constructor or the
 * <TT>setMatrix()</TT> method.
 * <P>
 * To write the matrix to a file, call the <TT>prepareToWrite()</TT> method,
 * specifying the output stream to write. The <TT>prepareToWrite()</TT> method
 * returns an instance of class {@linkplain DoubleMatrixFile.Writer}. Call the
 * methods of the writer object to write the matrix, or sections of the matrix,
 * to the output stream. When finished, close the writer.
 * <P>
 * To read the matrix from a file, call the <TT>prepareToRead()</TT> method,
 * specifying the input stream to read. The <TT>prepareToRead()</TT> method
 * returns an instance of class {@linkplain DoubleMatrixFile.Reader}. Call the
 * methods of the reader object to read the matrix, or sections of the matrix,
 * from the input stream. When finished, close the reader.
 * <P>
 * Class DoubleMatrixFile includes a main program to combine a group of double
 * matrix files into one double matrix file. You might use this main program
 * when the processes of a parallel program have computed slices of a matrix and
 * stored the slices in separate files, and you want to combine the slices
 * together into one file.
 * <P>
 * <B>Double Matrix File Format</B>
 * <P>
 * A double matrix file is a binary file containing the following items. Each
 * primitive item is written as though by java.io.DataOutput (<TT>int</TT>, four
 * bytes; <TT>double</TT>, eight bytes; most significant byte first).
 * <UL>
 * <LI>
 * Number of matrix rows, <I>R</I> (<TT>int</TT>). <I>R</I> &gt;= 0.
 * <LI>
 * Number of matrix columns, <I>C</I> (<TT>int</TT>). <I>C</I> &gt;= 0.
 * <LI>
 * Zero or more segments of matrix elements.
 * </UL>
 * <P>
 * Each segment contains the following items.
 * <UL>
 * <LI>
 * The segment's lower row index, <I>RL</I> (<TT>int</TT>). <I>RL</I> &gt;= 0.
 * <LI>
 * The segment's lower column index, <I>CL</I> (<TT>int</TT>). <I>CL</I> &gt;=
 * 0.
 * <LI>
 * Number of rows in the segment, <I>M</I> (<TT>int</TT>). <I>M</I> &gt;= 0.
 * <I>RL</I>+<I>M</I> &lt;= <I>R</I>.
 * <LI>
 * Number of columns in the segment, <I>N</I> (<TT>int</TT>). <I>N</I> &gt;= 0.
 * <I>CL</I>+<I>N</I> &lt;= <I>C</I>.
 * <LI>
 * The matrix elements in rows <I>RL</I>..<I>RL</I>+<I>M</I>-1 and columns
 * <I>CL</I>..<I>CL</I>+<I>N</I>-1 (<TT>double</TT>). Stored in row major order.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 06-Jan-2008
 */
public class DoubleMatrixFile
	{

// Hidden data members.

	private static final long BYTES_PER_ELEMENT = 8L;

	private int R = -1;
	private int C = -1;
	private double[][] myMatrix;

// Exported constructors.

	/**
	 * Construct a new double matrix file object. The matrix's number of rows
	 * and number of columns are uninitialized. Before using the matrix file,
	 * specify the number of rows and the number of columns by calling the
	 * <TT>setMatrix()</TT> method or by reading the matrix file from an input
	 * stream.
	 */
	public DoubleMatrixFile()
		{
		}

	/**
	 * Construct a new double matrix file object with the given number of rows,
	 * number of columns, and underlying matrix.
	 *
	 * @param  R          Number of rows.
	 * @param  C          Number of columns.
	 * @param  theMatrix  Underlying matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>R</TT> &lt; 0. Thrown if
	 *     <TT>C</TT> &lt; 0. Thrown if <TT>theMatrix.length</TT> does not
	 *     equal <TT>R</TT>.
	 */
	public DoubleMatrixFile
		(int R,
		 int C,
		 double[][] theMatrix)
		{
		setMatrix (R, C, theMatrix);
		}

// Exported operations.

	/**
	 * Returns the number of rows in this matrix file.
	 *
	 * @return  Number of rows <I>R</I>, or -1 if not initialized.
	 */
	public int getRowCount()
		{
		return R;
		}

	/**
	 * Returns the number of columns in this matrix file.
	 *
	 * @return  Number of columns, <I>C</I>, or -1 if not initialized.
	 */
	public int getColCount()
		{
		return C;
		}

	/**
	 * Obtain this matrix file's underlying matrix.
	 *
	 * @return  Underlying matrix, or null if not initialized.
	 */
	public double[][] getMatrix()
		{
		return myMatrix;
		}

	/**
	 * Set this matrix file's number of rows, number of columns, and underlying
	 * matrix.
	 *
	 * @param  R          Number of rows.
	 * @param  C          Number of columns.
	 * @param  theMatrix  Underlying matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>R</TT> &lt; 0. Thrown if
	 *     <TT>C</TT> &lt; 0. Thrown if <TT>theMatrix.length</TT> does not
	 *     equal <TT>R</TT>.
	 */
	public void setMatrix
		(int R,
		 int C,
		 double[][] theMatrix)
		{
		setRC (R, C);
		if (theMatrix.length != R)
			{
			throw new IllegalArgumentException
				("DoubleMatrixFile.setMatrix(): theMatrix.length (= " +
				 theMatrix.length + ") does not equal R (= " + R + ")");
			}
		myMatrix = theMatrix;
		}

	/**
	 * Prepare to write this matrix file to the given output stream. The number
	 * of rows and the number of columns in this matrix file are written to the
	 * output stream at this time. To write this matrix file's elements, call
	 * methods on the returned writer, then close the writer.
	 * <P>
	 * For improved performance, specify an output stream with buffering, such
	 * as an instance of class java.io.BufferedOutputStream.
	 *
	 * @param  theStream  Output stream.
	 *
	 * @return  Writer object with which to write this matrix file.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if this matrix file object is
	 *     uninitialized.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Writer prepareToWrite
		(OutputStream theStream)
		throws IOException
		{
		if (myMatrix == null)
			{
			throw new IllegalStateException
				("DoubleMatrixFile.prepareToWrite(): Not initialized");
			}
		return new Writer (theStream);
		}

	/**
	 * Prepare to read this matrix file from the given input stream. The number
	 * of rows and the number of columns are read from the input stream at this
	 * time. If this matrix file object is uninitialized, it becomes initialized
	 * with the given number of rows and columns, and storage is allocated for
	 * the underlying matrix's row references; call the <TT>getMatrix()</TT>
	 * method to obtain the underlying matrix. If this matrix file object is
	 * already initialized, the number of rows and columns read from the input
	 * stream must match those of this matrix file object, otherwise an
	 * exception is thrown. To read this matrix file's elements, call methods on
	 * the returned reader, then close the reader.
	 * <P>
	 * For improved performance, specify an input stream with buffering, such
	 * as an instance of class java.io.BufferedInputStream.
	 *
	 * @param  theStream  Input stream.
	 *
	 * @return  Reader object with which to read this matrix file.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 * @exception  InvalidMatrixFileException
	 *     (subclass of IOException) Thrown if the input stream's contents were
	 *     invalid.
	 */
	public Reader prepareToRead
		(InputStream theStream)
		throws IOException
		{
		return new Reader (theStream);
		}

	/**
	 * Main program to combine a group of double matrix files into one double
	 * matrix file. The program reads all the given input files into a matrix,
	 * then writes the entire matrix to the given output file as a single
	 * segment. Each input file must be for a matrix with the same number of
	 * rows and columns. There must be sufficient main memory to hold the entire
	 * matrix.
	 * <P>
	 * You might use this main program when the processes of a parallel program
	 * have computed slices of a matrix and stored the slices in separate files,
	 * and you want to combine the slices together into one file.
	 * <P>
	 * Usage: java edu.rit.io.DoubleMatrixFile <I>outfile</I> <I>infile1</I>
	 * [ <I>infile2</I> . . . ]
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Validate command line arguments.
		if (args.length < 2)
			{
			System.err.println ("Usage: java edu.rit.io.DoubleMatrixFile <outfile> <infile1> [<infile2> ...]");
			System.exit (1);
			}

		// Double matrix file object.
		DoubleMatrixFile dmf = new DoubleMatrixFile();

		// Read each input file.
		for (int i = 1; i < args.length; ++ i)
			{
			DoubleMatrixFile.Reader reader =
				dmf.prepareToRead
					(new BufferedInputStream
						(new FileInputStream (args[i])));
			reader.read();
			reader.close();
			}

		// Write output file.
		DoubleMatrixFile.Writer writer =
			dmf.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (args[0])));
		writer.write();
		writer.close();
		}

// Exported helper classes.

	/**
	 * Class DoubleMatrixFile.Writer provides an object with which to write a
	 * {@linkplain DoubleMatrixFile} to an output stream.
	 * <P>
	 * When a writer is created, the number of rows and number of columns are
	 * written to the output stream. Each time the <TT>write()</TT>,
	 * <TT>writeRowSlice()</TT>, <TT>writeColSlice()</TT>, or
	 * <TT>writePatch()</TT> method is called, one segment of matrix elements is
	 * written to the output stream. When finished, call the <TT>close()</TT>
	 * method.
	 * <P>
	 * <I>Note:</I> Class DoubleMatrixFile.Writer is not multiple thread safe.
	 *
	 * @author  Alan Kaminsky
	 * @version 06-Jan-2008
	 */
	public class Writer
		{

	// Hidden data members.

		private OutputStream myOs;
		private DataOutputStream myDos;

	// Hidden constructors.

		/**
		 * Construct a new double matrix file writer.
		 *
		 * @param  theStream  Output stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private Writer
			(OutputStream theStream)
			throws IOException
			{
			if (theStream == null)
				{
				throw new NullPointerException
					("DoubleMatrixFile.Writer(): theStream is null");
				}
			myOs = theStream;
			myDos = new DataOutputStream (theStream);
			myDos.writeInt (R); // Number of matrix rows
			myDos.writeInt (C); // Number of matrix columns
			}

	// Exported operations.

		/**
		 * Write all rows and columns of the matrix to the output stream.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void write()
			throws IOException
			{
			write (0, R-1, 0, C-1);
			}

		/**
		 * Write the given row slice of the matrix to the output stream.
		 * Elements in the given range of rows and in all columns are written.
		 * <P>
		 * <I>Note:</I> <TT>theRowRange</TT>'s stride must be 1.
		 *
		 * @param  theRowRange  Range of matrix rows.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theRowRange</TT>
		 *     is outside the range 0 .. <I>R</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writeRowSlice
			(Range theRowRange)
			throws IOException
			{
			if (theRowRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Writer.writeRowSlice(): theRowRange stride > 1");
				}
			int RL = theRowRange.lb();
			int RU = theRowRange.ub();
			if (0 > RL || RL + theRowRange.length() > R)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Writer.writeRowSlice(): theRowRange = " +
					 theRowRange + " out of bounds");
				}
			write (RL, RU, 0, C-1);
			}

		/**
		 * Write the given column slice of the matrix to the output stream.
		 * Elements in all rows and in the given range of columns are written.
		 * <P>
		 * <I>Note:</I> <TT>theColRange</TT>'s stride must be 1.
		 *
		 * @param  theColRange  Range of matrix columns.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theColRange</TT>
		 *     is outside the range 0 .. <I>C</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writeColSlice
			(Range theColRange)
			throws IOException
			{
			if (theColRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Writer.writeColSlice(): theColRange stride > 1");
				}
			int CL = theColRange.lb();
			int CU = theColRange.ub();
			if (0 > CL || CL + theColRange.length() > C)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Writer.writeColSlice(): theColRange = " +
					 theColRange + " out of bounds");
				}
			write (0, R-1, CL, CU);
			}

		/**
		 * Write the given patch of the matrix to the output stream. Elements in
		 * the given range of rows and in the given range of columns are
		 * written.
		 * <P>
		 * <I>Note:</I> <TT>theRowRange</TT>'s stride must be 1.
		 * <TT>theColRange</TT>'s stride must be 1.
		 *
		 * @param  theRowRange  Range of matrix rows.
		 * @param  theColRange  Range of matrix columns.
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
		 *     is outside the range 0 .. <I>R</I>-1. Thrown if any index in
		 *     <TT>theColRange</TT> is outside the range 0 .. <I>C</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writePatch
			(Range theRowRange,
			 Range theColRange)
			throws IOException
			{
			if (theRowRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Writer.writePatch(): theRowRange stride > 1");
				}
			if (theColRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Writer.writePatch(): theColRange stride > 1");
				}
			int RL = theRowRange.lb();
			int RU = theRowRange.ub();
			if (0 > RL || RL + theRowRange.length() > R)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Writer.writePatch(): theRowRange = " +
					 theRowRange + " out of bounds");
				}
			int CL = theColRange.lb();
			int CU = theColRange.ub();
			if (0 > CL || CL + theColRange.length() > C)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Writer.writePatch(): theColRange = " +
					 theColRange + " out of bounds");
				}
			write (RL, RU, CL, CU);
			}

		/**
		 * Close the output stream.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void close()
			throws IOException
			{
			myDos.close();
			}

	// Hidden operations.

		/**
		 * Write the given rows and columns of the matrix to the output stream.
		 *
		 * @param  RL  Segment lower row index.
		 * @param  RU  Segment upper row index.
		 * @param  CL  Segment lower column index.
		 * @param  CU  Segment upper column index.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private void write
			(int RL,
			 int RU,
			 int CL,
			 int CU)
			throws IOException
			{
			myDos.writeInt (RL); // Segment lower row index
			myDos.writeInt (CL); // Segment lower column index
			myDos.writeInt (RU-RL+1); // Number of rows in segment
			myDos.writeInt (CU-CL+1); // Number of columns in segment
			for (int r = RL; r <= RU; ++ r)
				{
				double[] myMatrix_r = myMatrix[r];
				for (int c = CL; c <= CU; ++ c)
					{
					myDos.writeDouble (myMatrix_r[c]);
					}
				}
			}
		}

	/**
	 * Class DoubleMatrixFile.Reader provides an object with which to read a
	 * {@linkplain DoubleMatrixFile} from an input stream.
	 * <P>
	 * When a reader is created, the number of rows and number of columns are
	 * read from the input stream. If the matrix file object is uninitialized,
	 * it becomes initialized with the given number of rows and columns, and
	 * storage is allocated for the underlying matrix's row references. If the
	 * matrix file object is already initialized, the number of rows and columns
	 * read from the input stream must match those of the matrix file object,
	 * otherwise an exception is thrown.
	 * <P>
	 * To read the matrix element segments one at a time, first call the
	 * <TT>getRowRange()</TT> and <TT>getColRange()</TT> methods to obtain the
	 * range of rows and columns in the next segment. At this point, allocate
	 * storage for the rows and columns in the underlying matrix if necessary.
	 * Then call the <TT>readSegment()</TT> method to read the actual matrix
	 * elements. Repeat these steps if there are additional segments.
	 * <P>
	 * To read all the matrix element segments (or all the remaining segments),
	 * call the <TT>read()</TT> method.
	 * <P>
	 * Methods are also provided to read the current segment, or all the
	 * remaining segments, through a <I>filter.</I> As the segment(s) are read,
	 * only those matrix elements within a given row range, a given column
	 * range, or a given row and column range are actually stored in the
	 * underlying matrix.
	 * <P>
	 * When finished, call the <TT>close()</TT> method.
	 * <P>
	 * <I>Note:</I> Class DoubleMatrixFile.Reader is not multiple thread safe.
	 *
	 * @author  Alan Kaminsky
	 * @version 07-Jan-2008
	 */
	public class Reader
		{

	// Hidden data members.

		private InputStream myIs;
		private DataInputStream myDis;
		private Range myRowRange;
		private Range myColRange;

	// Hidden constructors.

		/**
		 * Construct a new double matrix file reader.
		 *
		 * @param  theStream  Input stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		private Reader
			(InputStream theStream)
			throws IOException
			{
			if (theStream == null)
				{
				throw new NullPointerException
					("DoubleMatrixFile.Reader(): theStream is null");
				}
			myIs = theStream;
			myDis = new DataInputStream (theStream);

			int R = myDis.readInt(); // Number of matrix rows
			int C = myDis.readInt(); // Number of matrix columns
			if (myMatrix == null)
				{
				setRC (R, C);
				myMatrix = new double [R] [];
				}
			else if (DoubleMatrixFile.this.R != R)
				{
				throw new InvalidMatrixFileException
					("DoubleMatrixFile.Reader(): Number of rows from stream (" +
					 R + ") != number of rows in this matrix file (" +
					 DoubleMatrixFile.this.R + ")");
				}
			else if (DoubleMatrixFile.this.C != C)
				{
				throw new InvalidMatrixFileException
					("DoubleMatrixFile.Reader(): Number of columns from stream (" +
					 C + ") != number of columns in this matrix file (" +
					 DoubleMatrixFile.this.C + ")");
				}

			getNextSegment();
			}

	// Exported operations.

		/**
		 * Read all matrix element segments from the input stream. If some
		 * segments have already been read, the <TT>read()</TT> method reads all
		 * remaining segments. If there are no more segments, the
		 * <TT>read()</TT> method does nothing. If storage is not already
		 * allocated in the underlying matrix for the matrix elements, the
		 * <TT>read()</TT> method allocates the necessary storage.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void read()
			throws IOException
			{
			while (myRowRange != null) readSegment();
			}

		/**
		 * Read all matrix element segments from the input stream, storing only
		 * the matrix elements in the given row slice. If some segments have
		 * already been read, the <TT>readRowSlice()</TT> method reads all
		 * remaining segments. If there are no more segments, the
		 * <TT>readRowSlice()</TT> method does nothing. If storage is not
		 * already allocated in the underlying matrix for the matrix elements,
		 * the <TT>readRowSlice()</TT> method allocates the necessary storage.
		 *
		 * @param  theRowRange  Row range.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theRowRange</TT>
		 *     is outside the range 0 .. <I>R</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readRowSlice
			(Range theRowRange)
			throws IOException
			{
			while (myRowRange != null)
				{
				readSegmentRowSlice (theRowRange);
				}
			}

		/**
		 * Read all matrix element segments from the input stream, storing only
		 * the matrix elements in the given column slice. If some segments have
		 * already been read, the <TT>readColSlice()</TT> method reads all
		 * remaining segments. If there are no more segments, the
		 * <TT>readColSlice()</TT> method does nothing. If storage is not
		 * already allocated in the underlying matrix for the matrix elements,
		 * the <TT>readColSlice()</TT> method allocates the necessary storage.
		 *
		 * @param  theColRange  Column range.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theColRange</TT>
		 *     is outside the range 0 .. <I>C</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readColSlice
			(Range theColRange)
			throws IOException
			{
			while (myRowRange != null)
				{
				readSegmentColSlice (theColRange);
				}
			}

		/**
		 * Read all matrix element segments from the input stream, storing only
		 * the matrix elements in the given patch. If some segments have already
		 * been read, the <TT>readPatch()</TT> method reads all remaining
		 * segments. If there are no more segments, the <TT>readPatch()</TT>
		 * method does nothing. If storage is not already allocated in the
		 * underlying matrix for the matrix elements, the <TT>readPatch()</TT>
		 * method allocates the necessary storage.
		 *
		 * @param  theRowRange  Row range.
		 * @param  theColRange  Column range.
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
		 *     is outside the range 0 .. <I>R</I>-1. Thrown if any index in
		 *     <TT>theColRange</TT> is outside the range 0 .. <I>C</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readPatch
			(Range theRowRange,
			 Range theColRange)
			throws IOException
			{
			while (myRowRange != null)
				{
				readSegmentPatch (theRowRange, theColRange);
				}
			}

		/**
		 * Obtain the row range of the next matrix element segment in the input
		 * stream. If there are no more segments, null is returned.
		 *
		 * @return  Row range, or null.
		 */
		public Range getRowRange()
			{
			return myRowRange;
			}

		/**
		 * Obtain the column range of the next matrix element segment in the
		 * input stream. If there are no more segments, null is returned.
		 *
		 * @return  Column range, or null.
		 */
		public Range getColRange()
			{
			return myColRange;
			}

		/**
		 * Read the next matrix element segment from the input stream. If there
		 * are no more segments, the <TT>readSegment()</TT> method does nothing.
		 * If storage is not already allocated in the underlying matrix for the
		 * matrix elements, the <TT>readSegment()</TT> method allocates the
		 * necessary storage.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readSegment()
			throws IOException
			{
			readSegment (0, R-1, 0, C-1);
			}

		/**
		 * Read the next matrix element segment from the input stream, storing
		 * only the matrix elements in the given row slice. If there are no more
		 * segments, the <TT>readSegmentRowSlice()</TT> method does nothing. If
		 * storage is not already allocated in the underlying matrix for the
		 * matrix elements, the <TT>readSegmentRowSlice()</TT> method allocates
		 * the necessary storage.
		 *
		 * @param  theRowRange  Row range.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theRowRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theRowRange</TT>
		 *     is outside the range 0 .. <I>R</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readSegmentRowSlice
			(Range theRowRange)
			throws IOException
			{
			if (theRowRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Reader.readSegmentRowSlice(): theRowRange stride > 1");
				}
			int RL = theRowRange.lb();
			int RU = theRowRange.ub();
			if (0 > RL || RL + theRowRange.length() > R)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Reader.readSegmentRowSlice(): theRowRange = " +
					 theRowRange + " out of bounds");
				}
			readSegment (RL, RU, 0, C-1);
			}

		/**
		 * Read the next matrix element segment from the input stream, storing
		 * only the matrix elements in the given column slice. If there are no
		 * more segments, the <TT>readSegmentRowSlice()</TT> method does
		 * nothing. If storage is not already allocated in the underlying matrix
		 * for the matrix elements, the <TT>readSegmentRowSlice()</TT> method
		 * allocates the necessary storage.
		 *
		 * @param  theColRange  Column range.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT> is null.
		 * @exception  IllegalArgumentException
		 *     (unchecked exception) Thrown if <TT>theColRange</TT>'s stride is
		 *     greater than 1.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if any index in <TT>theColRange</TT>
		 *     is outside the range 0 .. <I>C</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readSegmentColSlice
			(Range theColRange)
			throws IOException
			{
			if (theColRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Reader.readSegmentColSlice(): theColRange stride > 1");
				}
			int CL = theColRange.lb();
			int CU = theColRange.ub();
			if (0 > CL || CL + theColRange.length() > C)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Reader.readSegmentColSlice(): theColRange = " +
					 theColRange + " out of bounds");
				}
			readSegment (0, R-1, CL, CU);
			}

		/**
		 * Read the next matrix element segment from the input stream, storing
		 * only the matrix elements in the given patch. If there are no more
		 * segments, the <TT>readSegmentPatch()</TT> method does nothing. If
		 * storage is not already allocated in the underlying matrix for the
		 * matrix elements, the <TT>readSegmentPatch()</TT> method allocates the
		 * necessary storage.
		 *
		 * @param  theRowRange  Row range.
		 * @param  theColRange  Column range.
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
		 *     is outside the range 0 .. <I>R</I>-1. Thrown if any index in
		 *     <TT>theColRange</TT> is outside the range 0 .. <I>C</I>-1.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		public void readSegmentPatch
			(Range theRowRange,
			 Range theColRange)
			throws IOException
			{
			if (theRowRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Reader.readSegmentPatch(): theRowRange stride > 1");
				}
			int RL = theRowRange.lb();
			int RU = theRowRange.ub();
			if (0 > RL || RL + theRowRange.length() > R)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Reader.readSegmentPatch(): theRowRange = " +
					 theRowRange + " out of bounds");
				}
			if (theColRange.stride() != 1)
				{
				throw new IllegalArgumentException
					("DoubleMatrixImage.Reader.readSegmentPatch(): theColRange stride > 1");
				}
			int CL = theColRange.lb();
			int CU = theColRange.ub();
			if (0 > CL || CL + theColRange.length() > C)
				{
				throw new IndexOutOfBoundsException
					("DoubleMatrixImage.Reader.readSegmentPatch(): theColRange = " +
					 theColRange + " out of bounds");
				}
			readSegment (RL, RU, CL, CU);
			}

		/**
		 * Close the input stream.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void close()
			throws IOException
			{
			myDis.close();
			}

	// Hidden operations.

		/**
		 * Read the bounds of the next segment from the input stream and store
		 * them in myRowRange and myColRange.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		private void getNextSegment()
			throws IOException
			{
			try
				{
				int RL = myDis.readInt(); // Segment lower row index
				int CL = myDis.readInt(); // Segment lower column index
				int M = myDis.readInt();  // Number of rows in segment
				int N = myDis.readInt();  // Number of columns in segment

				if (RL < 0)
					{
					throw new InvalidMatrixFileException
						("DoubleMatrixFile.Reader.getNextSegment(): Invalid segment lower row index (" +
						 RL + ")");
					}
				if (CL < 0)
					{
					throw new InvalidMatrixFileException
						("DoubleMatrixFile.Reader.getNextSegment(): Invalid segment lower column index (" +
						 CL + ")");
					}
				if (M < 0 || RL + M > R)
					{
					throw new InvalidMatrixFileException
						("DoubleMatrixFile.Reader.getNextSegment(): Invalid numer of rows in segment (" +
						 M + ")");
					}
				if (N < 0 || CL + N > C)
					{
					throw new InvalidMatrixFileException
						("DoubleMatrixFile.Reader.getNextSegment(): Invalid numer of columns in segment (" +
						 N + ")");
					}

				myRowRange = new Range (RL, RL+M-1);
				myColRange = new Range (CL, CL+N-1);
				}

			catch (EOFException exc)
				{
				myRowRange = null;
				myColRange = null;
				}
			}

		/**
		 * Read the next matrix element segment from the input stream, storing
		 * only the matrix elements within the given row and column index
		 * bounds. If there are no more segments, the <TT>readSegment()</TT>
		 * method does nothing. If storage is not already allocated in the
		 * underlying matrix for the matrix elements, the <TT>readSegment()</TT>
		 * method allocates the necessary storage.
		 *
		 * @param  RL  Lower row index.
		 * @param  RU  Upper row index.
		 * @param  CL  Lower column index.
		 * @param  CU  Upper column index.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 * @exception  InvalidMatrixFileException
		 *     (subclass of IOException) Thrown if the input stream's contents
		 *     were invalid.
		 */
		private void readSegment
			(int RL,
			 int RU,
			 int CL,
			 int CU)
			throws IOException
			{
			// Early return if no more segments.
			if (myRowRange == null) return;

			// Get row and column bounds of segment.
			int SRL = myRowRange.lb();
			int SRU = myRowRange.ub();
			int SCL = myColRange.lb();
			int SCU = myColRange.ub();

			// Number of bytes in an entire row.
			long rowBytes = BYTES_PER_ELEMENT * myColRange.length();

			// First row to read.
			int firstRow = Math.max (RL, SRL);

			// Last row to read.
			int lastRow = Math.min (RU, SRU);

			// Number of rows (bytes) to skip at the beginning of the segment.
			int preSkipRows = firstRow - SRL;
			long preSkipRowBytes = preSkipRows * rowBytes;

			// Number of rows (bytes) to skip at the end of the segment.
			int postSkipRows = SRU - lastRow;
			long postSkipRowBytes = postSkipRows * rowBytes;

			// First column to read.
			int firstCol = Math.max (CL, SCL);

			// Last column to read.
			int lastCol = Math.min (CU, SCU);

			// Number of columns (bytes) to skip at the beginning of the row.
			int preSkipCols = firstCol - SCL;
			long preSkipColBytes = preSkipCols * BYTES_PER_ELEMENT;

			// Number of columns (bytes) to skip at the end of the row.
			int postSkipCols = SCU - lastCol;
			long postSkipColBytes = postSkipCols * BYTES_PER_ELEMENT;

			// Skip rows at beginning of segment.
			skipFully (preSkipRowBytes);

			// Read rows.
			for (int r = firstRow; r <= lastRow; ++ r)
				{
				// Allocate storage for row if necessary.
				double[] myMatrix_r = myMatrix[r];
				if (myMatrix_r == null)
					{
					myMatrix_r = new double [C];
					myMatrix[r] = myMatrix_r;
					}

				// Skip columns at beginning of row.
				skipFully (preSkipColBytes);

				// Read columns.
				for (int c = firstCol; c <= lastCol; ++ c)
					{
					myMatrix_r[c] = myDis.readDouble();
					}

				// Skip columns at end of row.
				skipFully (postSkipColBytes);
				}

			// Skip rows at end of segment.
			skipFully (postSkipRowBytes);

			getNextSegment();
			}

		/**
		 * Skip over the given number of bytes in the input stream.
		 *
		 * @param  n  Number of bytes to skip.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private void skipFully
			(long n)
			throws IOException
			{
			while (n > 0L) n -= myDis.skip (n);
			}

		}

// Hidden operations.

	/**
	 * Set this matrix file's number of rows and number of columns.
	 *
	 * @param  R  Number of rows.
	 * @param  C  Number of columns.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>R</TT> &lt; 0. Thrown if
	 *     <TT>C</TT> &lt; 0.
	 */
	void setRC
		(int R,
		 int C)
		{
		if (R < 0)
			{
			throw new IllegalArgumentException
				("DoubleMatrixFile.setHeightAndWidth(): R = " + R + " illegal");
			}
		if (C < 0)
			{
			throw new IllegalArgumentException
				("DoubleMatrixFile.setHeightAndWidth(): C = " + C + " illegal");
			}
		this.R = R;
		this.C = C;
		}

	}
