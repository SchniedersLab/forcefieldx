//******************************************************************************
//
// File:    ObjectBuf.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.ObjectBuf<T>
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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.nio.ByteBuffer;

import edu.rit.mp.buf.EmptyObjectBuf;
import edu.rit.mp.buf.ObjectArrayBuf;
import edu.rit.mp.buf.ObjectArrayBuf_1;
import edu.rit.mp.buf.ObjectItemBuf;
import edu.rit.mp.buf.ObjectMatrixBuf;
import edu.rit.mp.buf.ObjectMatrixBuf_1;
import edu.rit.mp.buf.SharedObjectArrayBuf;
import edu.rit.mp.buf.SharedObjectArrayBuf_1;
import edu.rit.mp.buf.SharedObjectBuf;
import edu.rit.pj.reduction.SharedObject;
import edu.rit.pj.reduction.SharedObjectArray;
import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class ObjectBuf is the abstract base class for a buffer of object items
 * sent or received using the Message Protocol (MP). In a message, an object
 * item is represented as a sequence of bytes using Java Object Serialization.
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
 * To create an ObjectBuf, call one of the following static factory methods:
 * <UL>
 * <LI><TT>emptyBuffer()</TT>
 * <LI><TT>buffer()</TT>
 * <LI><TT>buffer (T)</TT>
 * <LI><TT>buffer (T[])</TT>
 * <LI><TT>sliceBuffer (T[], Range)</TT>
 * <LI><TT>sliceBuffers (T[], Range[])</TT>
 * <LI><TT>objectBuffer (T[])</TT>
 * <LI><TT>buffer (T[][])</TT>
 * <LI><TT>rowSliceBuffer (T[][], Range)</TT>
 * <LI><TT>rowSliceBuffers (T[][], Range[])</TT>
 * <LI><TT>colSliceBuffer (T[][], Range)</TT>
 * <LI><TT>colSliceBuffers (T[][], Range[])</TT>
 * <LI><TT>patchBuffer (T[][], Range, Range)</TT>
 * <LI><TT>patchBuffers (T[][], Range[], Range[])</TT>
 * <LI><TT>objectBuffer (T[][])</TT>
 * <LI><TT>buffer (SharedObject&lt;T&gt;)</TT>
 * <LI><TT>buffer (SharedObjectArray&lt;T&gt;)</TT>
 * <LI><TT>sliceBuffer (SharedObjectArray&lt;T&gt;, Range)</TT>
 * <LI><TT>sliceBuffers (SharedObjectArray&lt;T&gt;, Range[])</TT>
 * </UL>
 * <P>
 * There are two ways to create a buffer for an array of objects (type
 * <TT>T[]</TT>):
 * <OL TYPE=1>
 * <LI>
 * With the <TT>buffer(T[])</TT>, <TT>sliceBuffer(T[],Range)</TT>, and
 * <TT>sliceBuffers(T[],Range[])</TT> methods. These methods create a buffer
 * that sends and receives the array elements as multiple separate objects of
 * type <TT>T</TT>. The receiver must allocate an array of the proper dimension
 * to receive the incoming objects and must create a buffer that receives the
 * array elements as separate objects.
 * <P><LI>
 * With the <TT>objectBuffer(T[])</TT> method. This method creates a buffer that
 * sends and receives the entire array as one object of type <TT>T[]</TT>. The
 * receiver must also create a buffer that receives the entire array as one
 * object; the buffer's <TT>item</TT> field is automatically set to an array of
 * the proper dimension.
 * </OL>
 * <P>
 * There are two ways to create a buffer for a matrix of objects (type
 * <TT>T[][]</TT>):
 * <OL TYPE=1>
 * <LI>
 * With the <TT>buffer(T[][])</TT>, <TT>rowSliceBuffer(T[][],Range)</TT>,
 * <TT>rowSliceBuffers(T[][],Range[])</TT>,
 * <TT>colSliceBuffer(T[][],Range)</TT>,
 * <TT>colSliceBuffers(T[][],Range[])</TT>, <TT>patchBuffer(T[][],Range)</TT>,
 * and <TT>patchBuffers(T[][],Range[])</TT> methods. These methods create a
 * buffer that sends and receives the matrix elements as multiple separate
 * objects of type <TT>T</TT>. The receiver must allocate a matrix of the proper
 * dimensions to receive the incoming objects and must create a buffer that
 * receives the matrix elements as separate objects.
 * <P><LI>
 * With the <TT>objectBuffer(T[][])</TT> method. This method creates a buffer
 * that sends and receives the entire matrix as one object of type
 * <TT>T[][]</TT>. The receiver must also create a buffer that receives the
 * matrix as one object; the buffer's <TT>item</TT> field is automatically set
 * to a matrix of the proper dimensions.
 * </OL>
 * <P>
 * <B><I>Important Note:</I></B> An ObjectBuf uses the protected field
 * <TT>mySerializedItems</TT> to store the serialized representation of the
 * objects in the buffer. If the buffer is used to receive a message, the
 * serialized representation of the received objects is cached in
 * <TT>mySerializedItems</TT>. If the buffer is used to send a message and
 * <TT>mySerializedItems</TT> is empty, the objects in the buffer are
 * serialized, the serialized representation is cached in
 * <TT>mySerializedItems</TT>, and the serialized representation is sent in the
 * message. If the buffer is used to send a message and
 * <TT>mySerializedItems</TT> is not empty, the objects in the buffer are
 * <I>not</I> serialized; rather, the cached serialized representation is sent.
 * This is done to avoid re-serializing the objects if the buffer is used to
 * send copies of a message to multiple destinations, or if the buffer is used
 * to receive and then immediately send a message. However, if the state of any
 * object in the buffer changes, the buffer's <TT>reset()</TT> method must be
 * called; this tells the buffer to discard the cached serialized representation
 * and re-serialize the objects in the buffer.
 *
 * @param  <T>  Data type of the objects in the buffer.
 *
 * @author  Alan Kaminsky
 * @version 03-Jul-2008
 */
public abstract class ObjectBuf<T>
	extends Buf
	{

// Hidden data members.

	/**
	 * Byte array containing this buffer's object items in serialized form. If
	 * null, the object items need to be serialized.
	 */
	protected byte[] mySerializedItems;

// Hidden constructors.

	/**
	 * Construct a new object buffer.
	 *
	 * @param  theLength     Number of items.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theLength</TT> &lt; 0.
	 */
	protected ObjectBuf
		(int theLength)
		{
		super (Constants.TYPE_OBJECT, theLength);
		}

// Exported operations.

	/**
	 * Create an empty buffer. The buffer's length is 0. The buffer's item type
	 * is Object.
	 *
	 * @return  Empty buffer.
	 */
	public static ObjectBuf<Object> emptyBuffer()
		{
		return new EmptyObjectBuf();
		}

	/**
	 * Create a buffer for an object item. The item is stored in the
	 * <TT>item</TT> field of the buffer.
	 *
	 * @param  <T>  Data type of the objects in the buffer.
	 *
	 * @return  Buffer.
	 */
	public static <T> ObjectItemBuf<T> buffer()
		{
		return new ObjectItemBuf<T>();
		}

	/**
	 * Create a buffer for an object item with the given initial value. The
	 * item is stored in the <TT>item</TT> field of the buffer.
	 *
	 * @param  <T>   Data type of the objects in the buffer.
	 * @param  item  Initial value of the <TT>item</TT> field.
	 *
	 * @return  Buffer.
	 */
	public static <T> ObjectItemBuf<T> buffer
		(T item)
		{
		return new ObjectItemBuf<T> (item);
		}

	/**
	 * Create a buffer for the entire given object array. The returned buffer
	 * encompasses all the elements in <TT>theArray</TT>. The array elements are
	 * sent and received as multiple separate objects of type <TT>T</TT>.
	 *
	 * @param  <T>       Data type of the objects in the buffer.
	 * @param  theArray  Array.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public static <T> ObjectBuf<T> buffer
		(T[] theArray)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("ObjectBuf.buffer(): theArray is null");
			}
		int nr = Arrays.length (theArray);
		return new ObjectArrayBuf_1<T> (theArray, new Range (0, nr-1));
		}

	/**
	 * Create a buffer for one slice of the given object array. The returned
	 * buffer encompasses <TT>theRange</TT> of elements in <TT>theArray</TT>.
	 * The range's stride may be 1 or greater than 1. The array elements are
	 * sent and received as multiple separate objects of type <TT>T</TT>.
	 *
	 * @param  <T>       Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T> sliceBuffer
		(T[] theArray,
		 Range theRange)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("ObjectBuf.sliceBuffer(): theArray is null");
			}
		int nr = Arrays.length (theArray);
		if (0 > theRange.lb() || theRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("ObjectBuf.sliceBuffer(): theArray index range = 0.." +
				 (nr-1) + ", theRange = " + theRange);
			}
		if (theRange.stride() == 1)
			{
			return new ObjectArrayBuf_1<T> (theArray, theRange);
			}
		else
			{
			return new ObjectArrayBuf<T> (theArray, theRange);
			}
		}

	/**
	 * Create an array of buffers for multiple slices of the given object
	 * array. The returned buffer array has the same length as
	 * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
	 * encompasses the elements of <TT>theArray</TT> specified by
	 * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
	 * The array elements are sent and received as multiple separate objects of
	 * type <TT>T</TT>.
	 *
	 * @param  <T>        Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T>[] sliceBuffers
		(T[] theArray,
		 Range[] theRanges)
		{
		int n = theRanges.length;
		ObjectBuf<T>[] result = (ObjectBuf<T>[]) new ObjectBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = sliceBuffer (theArray, theRanges[i]);
			}
		return result;
		}

	/**
	 * Create a buffer for the entire given object array. The returned buffer
	 * encompasses all the elements in <TT>theArray</TT>. The array is sent and
	 * received as a single object of type <TT>T[]</TT>.
	 *
	 * @param  <T>       Data type of the objects in the buffer.
	 * @param  theArray  Array. May be null.
	 *
	 * @return  Buffer.
	 */
	public static <T> ObjectItemBuf<T[]> objectBuffer
		(T[] theArray)
		{
		return new ObjectItemBuf<T[]> (theArray);
		}

	/**
	 * Create a buffer for the entire given object matrix. The returned
	 * buffer encompasses all the rows and all the columns in
	 * <TT>theMatrix</TT>. The matrix elements are sent and received as multiple
	 * separate objects of type <TT>T</TT>.
	 *
	 * @param  <T>        Data type of the objects in the buffer.
	 * @param  theMatrix  Matrix.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theMatrix</TT> is null.
	 */
	public static <T> ObjectBuf<T> buffer
		(T[][] theMatrix)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("ObjectBuf.buffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		int nc = Arrays.colLength (theMatrix, 0);
		return new ObjectMatrixBuf_1<T>
			(theMatrix, new Range (0, nr-1), new Range (0, nc-1));
		}

	/**
	 * Create a buffer for one row slice of the given object matrix. The
	 * returned buffer encompasses <TT>theRowRange</TT> of rows, and all the
	 * columns, in <TT>theMatrix</TT>. The range's stride may be 1 or greater
	 * than 1. The matrix elements are sent and received as multiple separate
	 * objects of type <TT>T</TT>.
	 *
	 * @param  <T>          Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T> rowSliceBuffer
		(T[][] theMatrix,
		 Range theRowRange)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("ObjectBuf.rowSliceBuffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		if (0 > theRowRange.lb() || theRowRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("ObjectBuf.rowSliceBuffer(): theMatrix row index range = 0.." +
				 (nr-1) + ", theRowRange = " + theRowRange);
			}
		int nc = Arrays.colLength (theMatrix, theRowRange.lb());
		if (theRowRange.stride() == 1)
			{
			return new ObjectMatrixBuf_1<T>
				(theMatrix, theRowRange, new Range (0, nc-1));
			}
		else
			{
			return new ObjectMatrixBuf<T>
				(theMatrix, theRowRange, new Range (0, nc-1));
			}
		}

	/**
	 * Create an array of buffers for multiple row slices of the given object
	 * matrix. The returned buffer array has the same length as
	 * <TT>theRowRanges</TT>. Each element [<I>i</I>] of the returned buffer
	 * array encompasses the rows of <TT>theMatrix</TT> specified by
	 * <TT>theRowRanges[i]</TT> and all the columns of <TT>theMatrix</TT>. Each
	 * range's stride may be 1 or greater than 1. The matrix elements are sent
	 * and received as multiple separate objects of type <TT>T</TT>.
	 *
	 * @param  <T>           Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T>[] rowSliceBuffers
		(T[][] theMatrix,
		 Range[] theRowRanges)
		{
		int n = theRowRanges.length;
		ObjectBuf<T>[] result = (ObjectBuf<T>[]) new ObjectBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = rowSliceBuffer (theMatrix, theRowRanges[i]);
			}
		return result;
		}

	/**
	 * Create a buffer for one column slice of the given object matrix. The
	 * returned buffer encompasses all the rows, and <TT>theColRange</TT> of
	 * columns, in <TT>theMatrix</TT>. The range's stride may be 1 or greater
	 * than 1. The matrix elements are sent and received as multiple separate
	 * objects of type <TT>T</TT>.
	 *
	 * @param  <T>          Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T> colSliceBuffer
		(T[][] theMatrix,
		 Range theColRange)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("ObjectBuf.colSliceBuffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		int nc = Arrays.colLength (theMatrix, 0);
		if (0 > theColRange.lb() || theColRange.ub() >= nc)
			{
			throw new IndexOutOfBoundsException
				("ObjectBuf.colSliceBuffer(): theMatrix column index range = 0.."
				 + (nc-1) + ", theColRange = " + theColRange);
			}
		if (theColRange.stride() == 1)
			{
			return new ObjectMatrixBuf_1<T>
				(theMatrix, new Range (0, nr-1), theColRange);
			}
		else
			{
			return new ObjectMatrixBuf<T>
				(theMatrix, new Range (0, nr-1), theColRange);
			}
		}

	/**
	 * Create an array of buffers for multiple column slices of the given
	 * object matrix. The returned buffer array has the same length as
	 * <TT>theColRanges</TT>. Each element [<I>i</I>] of the returned buffer
	 * array encompasses all the rows of <TT>theMatrix</TT> and the columns of
	 * <TT>theMatrix</TT> specified by <TT>theColRanges[i]</TT>. Each range's
	 * stride may be 1 or greater than 1. The matrix elements are sent and
	 * received as multiple separate objects of type <TT>T</TT>.
	 *
	 * @param  <T>           Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T>[] colSliceBuffers
		(T[][] theMatrix,
		 Range[] theColRanges)
		{
		int n = theColRanges.length;
		ObjectBuf<T>[] result = (ObjectBuf<T>[]) new ObjectBuf [n];
		for (int i = 0; i < n; ++ i)
			{
			result[i] = colSliceBuffer (theMatrix, theColRanges[i]);
			}
		return result;
		}

	/**
	 * Create a buffer for one patch of the given object matrix. The returned
	 * buffer encompasses <TT>theRowRange</TT> of rows, and <TT>theColRange</TT>
	 * of columns, in <TT>theMatrix</TT>. Each range's stride may be 1 or
	 * greater than 1. The matrix elements are sent and received as multiple
	 * separate objects of type <TT>T</TT>.
	 *
	 * @param  <T>          Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T> patchBuffer
		(T[][] theMatrix,
		 Range theRowRange,
		 Range theColRange)
		{
		if (theMatrix == null)
			{
			throw new NullPointerException
				("ObjectBuf.patchBuffer(): theMatrix is null");
			}
		int nr = Arrays.rowLength (theMatrix);
		if (0 > theRowRange.lb() || theRowRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("ObjectBuf.patchBuffer(): theMatrix row index range = 0.."
				 + (nr-1) + ", theRowRange = " + theRowRange);
			}
		int nc = Arrays.colLength (theMatrix, theRowRange.lb());
		if (0 > theColRange.lb() || theColRange.ub() >= nc)
			{
			throw new IndexOutOfBoundsException
				("ObjectBuf.patchBuffer(): theMatrix column index range = 0.."
				 + (nc-1) + ", theColRange = " + theColRange);
			}
		if (theRowRange.stride() == 1 && theColRange.stride() == 1)
			{
			return new ObjectMatrixBuf_1<T> (theMatrix, theRowRange, theColRange);
			}
		else
			{
			return new ObjectMatrixBuf<T> (theMatrix, theRowRange, theColRange);
			}
		}

	/**
	 * Create an array of buffers for multiple patches of the given object
	 * matrix. The length of the returned buffer array is equal to the length of
	 * <TT>theRowRanges</TT> times the length of <TT>theColRanges</TT>. Each
	 * element of the returned buffer array encompasses the rows given in one
	 * element of <TT>theRowRanges</TT> array, and the columns given in one
	 * element of <TT>theColRanges</TT> array, in all possible combinations, of
	 * <TT>theMatrix</TT>. Each range's stride may be 1 or greater than 1. The
	 * matrix elements are sent and received as multiple separate objects of
	 * type <TT>T</TT>.
	 *
	 * @param  <T>           Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T>[] patchBuffers
		(T[][] theMatrix,
		 Range[] theRowRanges,
		 Range[] theColRanges)
		{
		int m = theRowRanges.length;
		int n = theColRanges.length;
		ObjectBuf<T>[] result = (ObjectBuf<T>[]) new ObjectBuf [m*n];
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
	 * Create a buffer for the entire given object matrix. The returned buffer
	 * encompasses all the rows and all the columns in <TT>theMatrix</TT>. The
	 * matrix is sent and received as a single object of type <TT>T[][]</TT>.
	 *
	 * @param  <T>        Data type of the objects in the buffer.
	 * @param  theMatrix  Matrix. May be null.
	 *
	 * @return  Buffer.
	 */
	public static <T> ObjectItemBuf<T[][]> objectBuffer
		(T[][] theMatrix)
		{
		return new ObjectItemBuf<T[][]> (theMatrix);
		}

	/**
	 * Create a buffer for a shared object item. The item is wrapped in an
	 * instance of class {@linkplain edu.rit.pj.reduction.SharedObject
	 * SharedObject}. Use the methods of the SharedObject object to access
	 * the actual item.
	 *
	 * @param  <T>   Data type of the objects in the buffer.
	 * @param  item  SharedObject object that wraps the item.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>item</TT> is null.
	 */
	public static <T> ObjectBuf<T> buffer
		(SharedObject<T> item)
		{
		if (item == null)
			{
			throw new NullPointerException
				("ObjectBuf.buffer(): item is null");
			}
		return new SharedObjectBuf<T> (item);
		}

	/**
	 * Create a buffer for the entire given shared object array. The returned
	 * buffer encompasses all the elements in <TT>theArray</TT>.
	 *
	 * @param  <T>       Data type of the objects in the buffer.
	 * @param  theArray  Array.
	 *
	 * @return  Buffer.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theArray</TT> is null.
	 */
	public static <T> ObjectBuf<T> buffer
		(SharedObjectArray<T> theArray)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("ObjectBuf.buffer(): theArray is null");
			}
		int nr = theArray.length();
		return new SharedObjectArrayBuf_1<T> (theArray, new Range (0, nr-1));
		}

	/**
	 * Create a buffer for one slice of the given shared object array. The
	 * returned buffer encompasses <TT>theRange</TT> of elements in
	 * <TT>theArray</TT>. The range's stride may be 1 or greater than 1.
	 *
	 * @param  <T>       Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T> sliceBuffer
		(SharedObjectArray<T> theArray,
		 Range theRange)
		{
		if (theArray == null)
			{
			throw new NullPointerException
				("ObjectBuf.sliceBuffer(): theArray is null");
			}
		int nr = theArray.length();
		if (0 > theRange.lb() || theRange.ub() >= nr)
			{
			throw new IndexOutOfBoundsException
				("ObjectBuf.sliceBuffer(): theArray index range = 0.."
				 + (nr-1) + ", theRange = " + theRange);
			}
		if (theRange.stride() == 1)
			{
			return new SharedObjectArrayBuf_1<T> (theArray, theRange);
			}
		else
			{
			return new SharedObjectArrayBuf<T> (theArray, theRange);
			}
		}

	/**
	 * Create an array of buffers for multiple slices of the given shared
	 * object array. The returned buffer array has the same length as
	 * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
	 * encompasses the elements of <TT>theArray</TT> specified by
	 * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
	 *
	 * @param  <T>        Data type of the objects in the buffer.
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
	public static <T> ObjectBuf<T>[] sliceBuffers
		(SharedObjectArray<T> theArray,
		 Range[] theRanges)
		{
		int n = theRanges.length;
		ObjectBuf<T>[] result = (ObjectBuf<T>[]) new ObjectBuf [n];
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
	public abstract T get
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
		 T item);

	/**
	 * Copy items from the given buffer to this buffer. The number of items
	 * copied is this buffer's length or <TT>theSrc</TT>'s length, whichever is
	 * smaller. If <TT>theSrc</TT> is this buffer, the <TT>copy()</TT> method
	 * does nothing.
	 * <P>
	 * The default implementation of the <TT>copy()</TT> method calls the
	 * <TT>defaultCopy()</TT> method. A subclass can override the
	 * <TT>copy()</TT> method to use a more efficient algorithm.
	 * <P>
	 * The default implementation of the <TT>copy()</TT> method also calls the
	 * <TT>reset()</TT> method.
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
		if (theSrc != this)
			{
			defaultCopy ((ObjectBuf<T>) theSrc, this);
			reset();
			}
		}

	/**
	 * Fill this buffer with the given item. The <TT>item</TT> is assigned to
	 * each element in this buffer.
	 * <P>
	 * The <TT>item</TT> must be an instance of class T or a subclass thereof.
	 * The <TT>item</TT> may be null. Note that since <TT>item</TT> is
	 * <I>assigned</I> to every buffer element, every buffer element ends up
	 * referring to the same <TT>item</TT>.
	 * <P>
	 * The <TT>fill()</TT> method calls the <TT>reset()</TT> method.
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
		T value = (T) item;
		for (int i = 0; i < myLength; ++ i)
			{
			put (i, value);
			}
		reset();
		}

	/**
	 * Create a temporary buffer with the same type of items and the same length
	 * as this buffer. The new buffer items are stored in a newly created array,
	 * separate from the storage for this buffer's items.
	 */
	public Buf getTemporaryBuf()
		{
		return buffer ((T[]) new Object [myLength]);
		}

	/**
	 * Reset this buffer. Call <TT>reset()</TT> if the state of any object in
	 * this buffer changes.
	 */
	public void reset()
		{
		mySerializedItems = null;
		}

// Hidden operations.

	/**
	 * Called by the I/O thread before sending message items using this buffer.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void preSend()
		throws IOException
		{
		if (mySerializedItems == null)
			{
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			ObjectOutputStream oos = new ObjectOutputStream (baos);
			oos.writeInt (myLength);
			for (int i = 0; i < myLength; ++ i)
				{
				oos.writeObject (get (i));
				}
			oos.close();
			mySerializedItems = baos.toByteArray();
			myMessageLength = mySerializedItems.length;
			}
		}

	/**
	 * Send as many items as possible from this buffer to the given byte
	 * buffer.
	 * <P>
	 * The <TT>sendItems()</TT> method must not block the calling thread; if it
	 * does, all message I/O in MP will be blocked.
	 *
	 * @param  i       Index of first item to send, in the range 0 ..
	 *                 <TT>length</TT>-1.
	 * @param  buffer  Byte buffer.
	 *
	 * @return  Number of items sent.
	 */
	protected int sendItems
		(int i,
		 ByteBuffer buffer)
		{
		int len = Math.min (myMessageLength - i, buffer.remaining());
		buffer.put (mySerializedItems, i, len);
		return len;
		}

	/**
	 * Called by the I/O thread before receiving message items using this
	 * buffer.
	 *
	 * @param  theReadLength  Actual number of items in message.
	 */
	void preReceive
		(int theReadLength)
		{
		mySerializedItems = new byte [theReadLength];
		myMessageLength = theReadLength;
		}

	/**
	 * Receive as many items as possible from the given byte buffer to this
	 * buffer.
	 * <P>
	 * The <TT>receiveItems()</TT> method must not block the calling thread; if
	 * it does, all message I/O in MP will be blocked.
	 *
	 * @param  i       Index of first item to receive, in the range 0 ..
	 *                 <TT>length</TT>-1.
	 * @param  num     Maximum number of items to receive.
	 * @param  buffer  Byte buffer.
	 *
	 * @return  Number of items received.
	 */
	protected int receiveItems
		(int i,
		 int num,
		 ByteBuffer buffer)
		{
		int len = num;
		len = Math.min (len, myMessageLength - i);
		len = Math.min (len, buffer.remaining());
		buffer.get (mySerializedItems, i, len);
		return len;
		}

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
		int n = Math.min (num, buffer.remaining());
		buffer.position (buffer.position() + n);
		return n;
		}

	/**
	 * Called by the I/O thread after receiving message items using this
	 * buffer.
	 *
	 * @param  theStatus       Status object that will be returned for the
	 *                         message; its contents may be altered if
	 *                         necessary.
	 * @param  theClassLoader  Alternate class loader to be used when receiving
	 *                         objects, or null.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void postReceive
		(Status theStatus,
		 ClassLoader theClassLoader)
		throws IOException
		{
		try
			{
			byte[] savedSerializedItems = mySerializedItems;
			ByteArrayInputStream bais =
				new ByteArrayInputStream (mySerializedItems);
			ObjectInputStream ois =
				new MPObjectInputStream (bais, theClassLoader);
			int nmsg = ois.readInt();
			int n = Math.min (myLength, nmsg);
			for (int i = 0; i < n; ++ i)
				{
				put (i, (T) ois.readObject());
				}
			ois.close();
			theStatus.length = nmsg;
			mySerializedItems = savedSerializedItems;
			}
		catch (ClassNotFoundException exc)
			{
			IOException exc2 =
				new IOException ("ObjectBuf.postReceive(): Class not found");
			exc2.initCause (exc);
			throw exc2;
			}
		catch (ClassCastException exc)
			{
			IOException exc2 =
				new IOException ("ObjectBuf.postReceive(): Wrong type");
			exc2.initCause (exc);
			throw exc2;
			}
		}

	/**
	 * Copy items from the given source buffer to the given destination buffer.
	 * The number of items copied is <TT>theSrc</TT>'s length or
	 * <TT>theDst</TT>'s length, whichever is smaller. Each item is copied
	 * individually using the <TT>get()</TT> and <TT>put()</TT> methods. It is
	 * assumed that <TT>theSrc</TT> is not the same as <TT>theDst</TT>.
	 *
	 * @param  <T>     Data type of the objects in the buffer.
	 * @param  theSrc  Source of items to copy.
	 * @param  theDst  Destination of items to copy.
	 */
	protected static <T> void defaultCopy
		(ObjectBuf<T> theSrc,
		 ObjectBuf<T> theDst)
		{
		int n = Math.min (theSrc.myLength, theDst.myLength);
		for (int i = 0; i < n; ++ i)
			{
			theDst.put (i, theSrc.get (i));
			}
		}

	}
