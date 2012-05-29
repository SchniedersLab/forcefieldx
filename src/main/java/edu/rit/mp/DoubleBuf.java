//******************************************************************************
//
// File:    DoubleBuf.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.DoubleBuf
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

package edu.rit.mp;

import java.nio.ByteBuffer;

import edu.rit.mp.buf.DoubleArrayBuf;
import edu.rit.mp.buf.DoubleArrayBuf_1;
import edu.rit.mp.buf.DoubleItemBuf;
import edu.rit.mp.buf.DoubleMatrixBuf;
import edu.rit.mp.buf.DoubleMatrixBuf_1;
import edu.rit.mp.buf.EmptyDoubleBuf;
import edu.rit.mp.buf.SharedDoubleArrayBuf;
import edu.rit.mp.buf.SharedDoubleArrayBuf_1;
import edu.rit.mp.buf.SharedDoubleBuf;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class DoubleBuf is the abstract base class for a buffer of double items sent
 * or received using the Message Protocol (MP). In a message, a double item is
 * represented as eight bytes, most significant byte first.
 * <P>
 * A buffer may be used to send one or more messages at the same time in
 * multiple threads. If a buffer is being used to send a message or messages,
 * the buffer must not be used to receive a message at the same time.
 * <P>
 * A buffer may be used to receive one message at a time. If a buffer is being
 * used to receive a message, the buffer must not be used to receive another
 * message in a different thread, and the buffer must not be used to send a
 * message or messages.
 * <P>
 * A buffer is a conduit for retrieving and storing data in some underlying data
 * structure. If the underlying data structure is multiple thread safe, then one
 * thread can be retrieving or storing data via the buffer at the same time as
 * other threads are accessing the data structure. If the underlying data
 * structure is not multiple thread safe, then other threads must not access the
 * data structure while one thread is retrieving or storing data via the buffer.
 * <P>
 * To create a DoubleBuf, call one of the following static factory methods:
 * <UL>
 * <LI><TT>emptyBuffer()</TT>
 * <LI><TT>buffer()</TT>
 * <LI><TT>buffer (double)</TT>
 * <LI><TT>buffer (double[])</TT>
 * <LI><TT>sliceBuffer (double[], Range)</TT>
 * <LI><TT>sliceBuffers (double[], Range[])</TT>
 * <LI><TT>buffer (double[][])</TT>
 * <LI><TT>rowSliceBuffer (double[][], Range)</TT>
 * <LI><TT>rowSliceBuffers (double[][], Range[])</TT>
 * <LI><TT>colSliceBuffer (double[][], Range)</TT>
 * <LI><TT>colSliceBuffers (double[][], Range[])</TT>
 * <LI><TT>patchBuffer (double[][], Range, Range)</TT>
 * <LI><TT>patchBuffers (double[][], Range[], Range[])</TT>
 * <LI><TT>buffer (SharedDouble)</TT>
 * <LI><TT>buffer (SharedDoubleArray)</TT>
 * <LI><TT>sliceBuffer (SharedDoubleArray, Range)</TT>
 * <LI><TT>sliceBuffers (SharedDoubleArray, Range[])</TT>
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 03-May-2008
 */
public abstract class DoubleBuf
	extends Buf
	{

// Hidden constructors.

	/**
	 * Construct a new double buffer.
	 *
	 * @param  theLength     Number of items.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theLength</TT> &lt; 0.
	 */
	protected DoubleBuf
		(int theLength)
		{
		super (Constants.TYPE_DOUBLE, theLength);
		}

// Exported operations.

	/**
	 * Create an empty buffer. The buffer's length is 0. The buffer's item type
	 * is double.
	 *
	 * @return  Empty buffer.
	 */
	public static DoubleBuf emptyBuffer()
		{
		return new EmptyDoubleBuf();
		}

	/**
	 * Create a buffer for a double item. The item is stored in the
	 * <TT>item</TT> field of the buffer.
	 *
	 * @return  Buffer.
	 */
	public static DoubleItemBuf buffer()
		{
		return new DoubleItemBuf();
		}

	/**
	 * Create a buffer for a double item with the given initial value. The
	 * item is stored in the <TT>item</TT> field of the buffer.
	 *
	 * @param  item  Initial value of the <TT>item</TT> field.
	 *
	 * @return  Buffer.
	 */
	public static DoubleItemBuf buffer
		(double item)
		{
		return new DoubleItemBuf (item);
		}

	/**
	 * Create a buffer for the entire given double array. The returned buffer
	 * encompasses all the elements in <TT>theArray</TT>.
	 *
	 * @param  theArray  Array.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public static DoubleBuf buffer
		(double[] theArray)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("DoubleBuf.buffer(): theArray is null");
			}
		int nr = Arrays.length (theArray);
		return new DoubleArrayBuf_1 (theArray, new Range (0, nr-1));
		}

	/**
	 * Create a buffer for one slice of the given double array. The returned
	 * buffer encompasses <TT>theRange</TT> of elements in <TT>theArray</TT>.
	 * The range's stride may be 1 or greater than 1.
	 *
	 * @param  theArray  Array.
	 * @param  theRange  Range of elements to include.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null or
	 *     <TT>theRange</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> does not include
	 *     all the indexes in <TT>theRange</TT>.
	 */
	public static DoubleBuf sliceBuffer
		(double[] theArray,
		 Range theRange)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("DoubleBuf.sliceBuffer(): theArray is null");
			}
		int nr = Arrays.length (theArray);
		if (0 > theRange.lb() || theRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("DoubleBuf.sliceBuffer(): theArray index range = 0.." +
				 (nr-1) + ", theRange = " + theRange);
			}
		if (theRange.stride() == 1)
			{
			return new DoubleArrayBuf_1 (theArray, theRange);
			}
		else
			{
			return new DoubleArrayBuf (theArray, theRange);
			}
		}

	/**
	 * Create an array of buffers for multiple slices of the given double
	 * array. The returned buffer array has the same length as
	 * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
	 * encompasses the elements of <TT>theArray</TT> specified by
	 * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
	 *
	 * @param  theArray   Array.
	 * @param  theRanges  Array of ranges of elements to include.
	 *
	 * @return  Array of buffers.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null or
	 *     <TT>theRanges</TT> or any element thereof is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT>'s allocation does
	 *     not include any element of <TT>theRanges</TT>.
	 */
	public static DoubleBuf[] sliceBuffers
		(double[] theArray,
		 Range[] theRanges)
		{
		int n = theRanges.length;
		DoubleBuf[] result = new DoubleBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = sliceBuffer (theArray, theRanges[i]);
			}
		return result;
		}

	/**
	 * Create a buffer for the entire given double matrix. The returned
	 * buffer encompasses all the rows and all the columns in
	 * <TT>theMatrix</TT>.
	 *
	 * @param  theMatrix  Matrix.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null.
	 */
	public static DoubleBuf buffer
		(double[][] theMatrix)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("DoubleBuf.buffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		int nc = Arrays.colLength (theMatrix, 0);
		return new DoubleMatrixBuf_1
			(theMatrix, new Range (0, nr-1), new Range (0, nc-1));
		}

	/**
	 * Create a buffer for one row slice of the given double matrix. The
	 * returned buffer encompasses <TT>theRowRange</TT> of rows, and all the
	 * columns, in <TT>theMatrix</TT>. The range's stride may be 1 or greater
	 * than 1.
	 *
	 * @param  theMatrix    Matrix.
	 * @param  theRowRange  Range of rows to include.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null or
	 *     <TT>theRowRange</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT>'s allocation does
	 *     not include <TT>theRowRange</TT>.
	 */
	public static DoubleBuf rowSliceBuffer
		(double[][] theMatrix,
		 Range theRowRange)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("DoubleBuf.rowSliceBuffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		if (0 > theRowRange.lb() || theRowRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("DoubleBuf.rowSliceBuffer(): theMatrix row index range = 0.." +
				 (nr-1) + ", theRowRange = " + theRowRange);
			}
		int nc = Arrays.colLength (theMatrix, theRowRange.lb());
		if (theRowRange.stride() == 1)
			{
			return new DoubleMatrixBuf_1
				(theMatrix, theRowRange, new Range (0, nc-1));
			}
		else
			{
			return new DoubleMatrixBuf
				(theMatrix, theRowRange, new Range (0, nc-1));
			}
		}

	/**
	 * Create an array of buffers for multiple row slices of the given double
	 * matrix. The returned buffer array has the same length as
	 * <TT>theRowRanges</TT>. Each element [<I>i</I>] of the returned buffer
	 * array encompasses the rows of <TT>theMatrix</TT> specified by
	 * <TT>theRowRanges[i]</TT> and all the columns of <TT>theMatrix</TT>. Each
	 * range's stride may be 1 or greater than 1.
	 *
	 * @param  theMatrix     Matrix.
	 * @param  theRowRanges  Array of ranges of rows to include.
	 *
	 * @return  Array of buffers.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null or
	 *     <TT>theRowRanges</TT> or any element thereof is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT>'s allocation does
	 *     not include any element of <TT>theRowRanges</TT>.
	 */
	public static DoubleBuf[] rowSliceBuffers
		(double[][] theMatrix,
		 Range[] theRowRanges)
		{
		int n = theRowRanges.length;
		DoubleBuf[] result = new DoubleBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = rowSliceBuffer (theMatrix, theRowRanges[i]);
			}
		return result;
		}

	/**
	 * Create a buffer for one column slice of the given double matrix. The
	 * returned buffer encompasses all the rows, and <TT>theColRange</TT> of
	 * columns, in <TT>theMatrix</TT>. The range's stride may be 1 or greater
	 * than 1.
	 *
	 * @param  theMatrix    Matrix.
	 * @param  theColRange  Range of columns to include.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null or
	 *     <TT>theColRange</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT>'s allocation does
	 *     not include <TT>theColRange</TT>.
	 */
	public static DoubleBuf colSliceBuffer
		(double[][] theMatrix,
		 Range theColRange)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("DoubleBuf.colSliceBuffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		int nc = Arrays.colLength (theMatrix, 0);
		if (0 > theColRange.lb() || theColRange.ub() >= nc)
			{
			throw new IndexOutOfBoundsException
				("DoubleBuf.colSliceBuffer(): theMatrix column index range = 0.."
				 + (nc-1) + ", theColRange = " + theColRange);
			}
		if (theColRange.stride() == 1)
			{
			return new DoubleMatrixBuf_1
				(theMatrix, new Range (0, nr-1), theColRange);
			}
		else
			{
			return new DoubleMatrixBuf
				(theMatrix, new Range (0, nr-1), theColRange);
			}
		}

	/**
	 * Create an array of buffers for multiple column slices of the given
	 * double matrix. The returned buffer array has the same length as
	 * <TT>theColRanges</TT>. Each element [<I>i</I>] of the returned buffer
	 * array encompasses all the rows of <TT>theMatrix</TT> and the columns of
	 * <TT>theMatrix</TT> specified by <TT>theColRanges[i]</TT>. Each range's
	 * stride may be 1 or greater than 1.
	 *
	 * @param  theMatrix     Matrix.
	 * @param  theColRanges  Array of ranges of columns to include.
	 *
	 * @return  Array of buffers.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null or
	 *     <TT>theColRanges</TT> or any element thereof is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT>'s allocation does
	 *     not include any element of <TT>theColRanges</TT>.
	 */
	public static DoubleBuf[] colSliceBuffers
		(double[][] theMatrix,
		 Range[] theColRanges)
		{
		int n = theColRanges.length;
		DoubleBuf[] result = new DoubleBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = colSliceBuffer (theMatrix, theColRanges[i]);
			}
		return result;
		}

	/**
	 * Create a buffer for one patch of the given double matrix. The returned
	 * buffer encompasses <TT>theRowRange</TT> of rows, and <TT>theColRange</TT>
	 * of columns, in <TT>theMatrix</TT>. Each range's stride may be 1 or
	 * greater than 1.
	 *
	 * @param  theMatrix    Matrix.
	 * @param  theRowRange  Range of rows to include.
	 * @param  theColRange  Range of columns to include.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null,
	 *     <TT>theRowRange</TT> is null, or <TT>theColRange</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT>'s allocation does
	 *     not include <TT>theRowRange</TT> and <TT>theColRange</TT>.
	 */
	public static DoubleBuf patchBuffer
		(double[][] theMatrix,
		 Range theRowRange,
		 Range theColRange)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("DoubleBuf.patchBuffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		if (0 > theRowRange.lb() || theRowRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("DoubleBuf.patchBuffer(): theMatrix row index range = 0.."
				 + (nr-1) + ", theRowRange = " + theRowRange);
			}
		int nc = Arrays.colLength (theMatrix, theRowRange.lb());
		if (0 > theColRange.lb() || theColRange.ub() >= nc)
			{
			throw new IndexOutOfBoundsException
				("DoubleBuf.patchBuffer(): theMatrix column index range = 0.."
				 + (nc-1) + ", theColRange = " + theColRange);
			}
		if (theRowRange.stride() == 1 && theColRange.stride() == 1)
			{
			return new DoubleMatrixBuf_1 (theMatrix, theRowRange, theColRange);
			}
		else
			{
			return new DoubleMatrixBuf (theMatrix, theRowRange, theColRange);
			}
		}

	/**
	 * Create an array of buffers for multiple patches of the given double
	 * matrix. The length of the returned buffer array is equal to the length of
	 * <TT>theRowRanges</TT> times the length of <TT>theColRanges</TT>. Each
	 * element of the returned buffer array encompasses the rows given in one
	 * element of <TT>theRowRanges</TT> array, and the columns given in one
	 * element of <TT>theColRanges</TT> array, in all possible combinations, of
	 * <TT>theMatrix</TT>. Each range's stride may be 1 or greater than 1.
	 *
	 * @param  theMatrix     Matrix.
	 * @param  theRowRanges  Array of ranges of rows to include.
	 * @param  theColRanges  Array of ranges of columns to include.
	 *
	 * @return  Array of buffers.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null,
	 *     <TT>theRowRanges</TT> or any element thereof is null, or
	 *     <TT>theColRanges</TT> or any element thereof is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT>'s allocation does
	 *     not include any element of <TT>theRowRanges</TT> or
	 *     <TT>theColRanges</TT>.
	 */
	public static DoubleBuf[] patchBuffers
		(double[][] theMatrix,
		 Range[] theRowRanges,
		 Range[] theColRanges)
		{
		int m = theRowRanges.length;
		int n = theColRanges.length;
		DoubleBuf[] result = new DoubleBuf [m*n];
		int k = 0;
		for (int i = 0; i < m; ++ i)
			{
			Range rowrange = theRowRanges[i];
			for (int j = 0; j < n; ++ j)
				{
				result[k++] =
					patchBuffer (theMatrix, rowrange, theColRanges[j]);
				}
			}
		return result;
		}

	/**
	 * Create a buffer for a shared double item. The item is wrapped in an
	 * instance of class {@linkplain edu.rit.pj.reduction.SharedDouble
	 * SharedDouble}. Use the methods of the SharedDouble object to access
	 * the actual item.
	 *
	 * @param  item  SharedDouble object that wraps the item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>item</TT> is null.
	 */
	public static DoubleBuf buffer
		(SharedDouble item)
		{
		if (item == null)
			{
			throw new NullPointerException
				("DoubleBuf.buffer(): item is null");
			}
		return new SharedDoubleBuf (item);
		}

	/**
	 * Create a buffer for the entire given shared double array. The returned
	 * buffer encompasses all the elements in <TT>theArray</TT>.
	 *
	 * @param  theArray  Array.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public static DoubleBuf buffer
		(SharedDoubleArray theArray)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("DoubleBuf.buffer(): theArray is null");
			}
		int nr = theArray.length();
		return new SharedDoubleArrayBuf_1 (theArray, new Range (0, nr-1));
		}

	/**
	 * Create a buffer for one slice of the given shared double array. The
	 * returned buffer encompasses <TT>theRange</TT> of elements in
	 * <TT>theArray</TT>. The range's stride may be 1 or greater than 1.
	 *
	 * @param  theArray  Array.
	 * @param  theRange  Range of elements to include.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null or
	 *     <TT>theRange</TT> is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> does not include
	 *     all the indexes in <TT>theRange</TT>.
	 */
	public static DoubleBuf sliceBuffer
		(SharedDoubleArray theArray,
		 Range theRange)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("DoubleBuf.sliceBuffer(): theArray is null");
			}
		int nr = theArray.length();
		if (0 > theRange.lb() || theRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("DoubleBuf.sliceBuffer(): theArray index range = 0.."
				 + (nr-1) + ", theRange = " + theRange);
			}
		if (theRange.stride() == 1)
			{
			return new SharedDoubleArrayBuf_1 (theArray, theRange);
			}
		else
			{
			return new SharedDoubleArrayBuf (theArray, theRange);
			}
		}

	/**
	 * Create an array of buffers for multiple slices of the given shared
	 * double array. The returned buffer array has the same length as
	 * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
	 * encompasses the elements of <TT>theArray</TT> specified by
	 * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
	 *
	 * @param  theArray   Array.
	 * @param  theRanges  Array of ranges of elements to include.
	 *
	 * @return  Array of buffers.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null or
	 *     <TT>theRanges</TT> or any element thereof is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>theArray</TT>'s allocation does
	 *     not include any element of <TT>theRanges</TT>.
	 */
	public static DoubleBuf[] sliceBuffers
		(SharedDoubleArray theArray,
		 Range[] theRanges)
		{
		int n = theRanges.length;
		DoubleBuf[] result = new DoubleBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = sliceBuffer (theArray, theRanges[i]);
			}
		return result;
		}

	/**
	 * Obtain the given item from this buffer.
	 * <P>
	 * The <TT>get()</TT> method must not block the calling thread; if it does,
	 * all message I/O in MP will be blocked.
	 *
	 * @param  i  Item index in the range 0 .. <TT>length()</TT>-1.
	 *
	 * @return  Item at index <TT>i</TT>.
	 */
	public abstract double get
		(int i);

	/**
	 * Store the given item in this buffer.
	 * <P>
	 * The <TT>put()</TT> method must not block the calling thread; if it does,
	 * all message I/O in MP will be blocked.
	 *
	 * @param  i     Item index in the range 0 .. <TT>length()</TT>-1.
	 * @param  item  Item to be stored at index <TT>i</TT>.
	 */
	public abstract void put
		(int i,
		 double item);

	/**
	 * Copy items from the given buffer to this buffer. The number of items
	 * copied is this buffer's length or <TT>theSrc</TT>'s length, whichever is
	 * smaller. If <TT>theSrc</TT> is this buffer, the <TT>copy()</TT> method
	 * does nothing.
	 * <P>
	 * The default implementation of the <TT>copy()</TT> method calls the
	 * <TT>defaultCopy()</TT> method. A subclass can override the
	 * <TT>copy()</TT> method to use a more efficient algorithm.
	 *
	 * @param  theSrc  Source of items to copy into this buffer.
	 *
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if <TT>theSrc</TT>'s item data type is
	 *     not the same as this buffer's item data type.
	 */
	public void copy
		(Buf theSrc)
		{
		if (theSrc != this) defaultCopy ((DoubleBuf) theSrc, this);
		}

	/**
	 * Fill this buffer with the given item. The <TT>item</TT> is assigned to
	 * each element in this buffer.
	 * <P>
	 * The <TT>item</TT> must be an instance of class Double. If the
	 * <TT>item</TT> is null, 0.0 is assigned to each element in this buffer.
	 *
	 * @param  item  Item.
	 *
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if the <TT>item</TT>'s data type is not
	 *     the same as this buffer's item data type.
	 */
	public void fill
		(Object item)
		{
		double value = item == null ? 0.0 : ((Double) item).doubleValue();
		for (int i = 0; i < myLength; ++ i)
			{
			put (i, value);
			}
		}

	/**
	 * Create a temporary buffer with the same type of items and the same length
	 * as this buffer. The new buffer items are stored in a newly created array,
	 * separate from the storage for this buffer's items.
	 */
	public Buf getTemporaryBuf()
		{
		return buffer (new double [myLength]);
		}

// Hidden operations.

	/**
	 * Skip as many items as possible from the given byte buffer.
	 *
	 * @param  num     Number of items to skip.
	 * @param  buffer  Buffer.
	 *
	 * @return  Number of items actually skipped.
	 */
	int skipItems
		(int num,
		 ByteBuffer buffer)
		{
		int n = Math.min (num, buffer.remaining() / 8);
		buffer.position (buffer.position() + 8 * n);
		return n;
		}

	/**
	 * Copy items from the given source buffer to the given destination buffer.
	 * The number of items copied is <TT>theSrc</TT>'s length or
	 * <TT>theDst</TT>'s length, whichever is smaller. Each item is copied
	 * individually using the <TT>get()</TT> and <TT>put()</TT> methods. It is
	 * assumed that <TT>theSrc</TT> is not the same as <TT>theDst</TT>.
	 *
	 * @param  theSrc  Source of items to copy.
	 * @param  theDst  Destination of items to copy.
	 */
	protected static void defaultCopy
		(DoubleBuf theSrc,
		 DoubleBuf theDst)
		{
		int n = Math.min (theSrc.myLength, theDst.myLength);
		for (int i = 0; i < n; ++ i)
			{
			theDst.put (i, theSrc.get (i));
			}
		}

	}
