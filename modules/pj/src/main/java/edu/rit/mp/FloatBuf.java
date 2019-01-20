//******************************************************************************
//
// File:    FloatBuf.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.FloatBuf
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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.mp;

import java.nio.ByteBuffer;

import edu.rit.mp.buf.EmptyFloatBuf;
import edu.rit.mp.buf.FloatArrayBuf;
import edu.rit.mp.buf.FloatArrayBuf_1;
import edu.rit.mp.buf.FloatItemBuf;
import edu.rit.mp.buf.FloatMatrixBuf;
import edu.rit.mp.buf.FloatMatrixBuf_1;
import edu.rit.mp.buf.SharedFloatArrayBuf;
import edu.rit.mp.buf.SharedFloatArrayBuf_1;
import edu.rit.mp.buf.SharedFloatBuf;
import edu.rit.pj.reduction.SharedFloat;
import edu.rit.pj.reduction.SharedFloatArray;
import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class FloatBuf is the abstract base class for a buffer of float items sent or
 * received using the Message Protocol (MP). In a message, a float item is
 * represented as four bytes, most significant byte first.
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
 * To create a FloatBuf, call one of the following static factory methods:
 * <UL>
 * <LI><code>emptyBuffer()</code>
 * <LI><code>buffer()</code>
 * <LI><code>buffer (float)</code>
 * <LI><code>buffer (float[])</code>
 * <LI><code>sliceBuffer (float[], Range)</code>
 * <LI><code>sliceBuffers (float[], Range[])</code>
 * <LI><code>buffer (float[][])</code>
 * <LI><code>rowSliceBuffer (float[][], Range)</code>
 * <LI><code>rowSliceBuffers (float[][], Range[])</code>
 * <LI><code>colSliceBuffer (float[][], Range)</code>
 * <LI><code>colSliceBuffers (float[][], Range[])</code>
 * <LI><code>patchBuffer (float[][], Range, Range)</code>
 * <LI><code>patchBuffers (float[][], Range[], Range[])</code>
 * <LI><code>buffer (SharedFloat)</code>
 * <LI><code>buffer (SharedFloatArray)</code>
 * <LI><code>sliceBuffer (SharedFloatArray, Range)</code>
 * <LI><code>sliceBuffers (SharedFloatArray, Range[])</code>
 * </UL>
 *
 * @author Alan Kaminsky
 * @version 03-May-2008
 */
public abstract class FloatBuf
        extends Buf {

// Hidden constructors.
    /**
     * Construct a new float buffer.
     *
     * @param theLength Number of items.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>theLength</code> &lt; 0.
     */
    protected FloatBuf(int theLength) {
        super(Constants.TYPE_FLOAT, theLength);
    }

// Exported operations.
    /**
     * Create an empty buffer. The buffer's length is 0. The buffer's item type
     * is float.
     *
     * @return Empty buffer.
     */
    public static FloatBuf emptyBuffer() {
        return new EmptyFloatBuf();
    }

    /**
     * Create a buffer for a float item. The item is stored in the
     * <code>item</code> field of the buffer.
     *
     * @return Buffer.
     */
    public static FloatItemBuf buffer() {
        return new FloatItemBuf();
    }

    /**
     * Create a buffer for a float item with the given initial value. The item
     * is stored in the <code>item</code> field of the buffer.
     *
     * @param item Initial value of the <code>item</code> field.
     * @return Buffer.
     */
    public static FloatItemBuf buffer(float item) {
        return new FloatItemBuf(item);
    }

    /**
     * Create a buffer for the entire given float array. The returned buffer
     * encompasses all the elements in <code>theArray</code>.
     *
     * @param theArray Array.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null.
     */
    public static FloatBuf buffer(float[] theArray) {
        if (theArray == null) {
            throw new NullPointerException("FloatBuf.buffer(): theArray is null");
        }
        int nr = Arrays.length(theArray);
        return new FloatArrayBuf_1(theArray, new Range(0, nr - 1));
    }

    /**
     * Create a buffer for one slice of the given float array. The returned
     * buffer encompasses <code>theRange</code> of elements in <code>theArray</code>.
     * The range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRange Range of elements to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null or
     * <code>theRange</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theArray</code> does not include all the indexes in <code>theRange</code>.
     */
    public static FloatBuf sliceBuffer(float[] theArray,
            Range theRange) {
        if (theArray == null) {
            throw new NullPointerException("FloatBuf.sliceBuffer(): theArray is null");
        }
        int nr = Arrays.length(theArray);
        if (0 > theRange.lb() || theRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("FloatBuf.sliceBuffer(): theArray index range = 0.."
                    + (nr - 1) + ", theRange = " + theRange);
        }
        if (theRange.stride() == 1) {
            return new FloatArrayBuf_1(theArray, theRange);
        } else {
            return new FloatArrayBuf(theArray, theRange);
        }
    }

    /**
     * Create an array of buffers for multiple slices of the given float array.
     * The returned buffer array has the same length as
     * <code>theRanges</code>. Each element [<I>i</I>] of the returned buffer array
     * encompasses the elements of <code>theArray</code> specified by
     * <code>theRanges[i]</code>. Each range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRanges Array of ranges of elements to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null or
     * <code>theRanges</code> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theArray</code>'s allocation does not include any element of
     * <code>theRanges</code>.
     */
    public static FloatBuf[] sliceBuffers(float[] theArray,
            Range[] theRanges) {
        int n = theRanges.length;
        FloatBuf[] result = new FloatBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = sliceBuffer(theArray, theRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for the entire given float matrix. The returned buffer
     * encompasses all the rows and all the columns in
     * <code>theMatrix</code>.
     *
     * @param theMatrix Matrix.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null.
     */
    public static FloatBuf buffer(float[][] theMatrix) {
        if (theMatrix == null) {
            throw new NullPointerException("FloatBuf.buffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        int nc = Arrays.colLength(theMatrix, 0);
        return new FloatMatrixBuf_1(theMatrix, new Range(0, nr - 1), new Range(0, nc - 1));
    }

    /**
     * Create a buffer for one row slice of the given float matrix. The returned
     * buffer encompasses <code>theRowRange</code> of rows, and all the columns, in
     * <code>theMatrix</code>. The range's stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRange Range of rows to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null or
     * <code>theRowRange</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theMatrix</code>'s allocation does not include <code>theRowRange</code>.
     */
    public static FloatBuf rowSliceBuffer(float[][] theMatrix,
            Range theRowRange) {
        if (theMatrix == null) {
            throw new NullPointerException("FloatBuf.rowSliceBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        if (0 > theRowRange.lb() || theRowRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("FloatBuf.rowSliceBuffer(): theMatrix row index range = 0.."
                    + (nr - 1) + ", theRowRange = " + theRowRange);
        }
        int nc = Arrays.colLength(theMatrix, theRowRange.lb());
        if (theRowRange.stride() == 1) {
            return new FloatMatrixBuf_1(theMatrix, theRowRange, new Range(0, nc - 1));
        } else {
            return new FloatMatrixBuf(theMatrix, theRowRange, new Range(0, nc - 1));
        }
    }

    /**
     * Create an array of buffers for multiple row slices of the given float
     * matrix. The returned buffer array has the same length as
     * <code>theRowRanges</code>. Each element [<I>i</I>] of the returned buffer
     * array encompasses the rows of <code>theMatrix</code> specified by
     * <code>theRowRanges[i]</code> and all the columns of <code>theMatrix</code>. Each
     * range's stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRanges Array of ranges of rows to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null or
     * <code>theRowRanges</code> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theMatrix</code>'s allocation does not include any element of
     * <code>theRowRanges</code>.
     */
    public static FloatBuf[] rowSliceBuffers(float[][] theMatrix,
            Range[] theRowRanges) {
        int n = theRowRanges.length;
        FloatBuf[] result = new FloatBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = rowSliceBuffer(theMatrix, theRowRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for one column slice of the given float matrix. The
     * returned buffer encompasses all the rows, and <code>theColRange</code> of
     * columns, in <code>theMatrix</code>. The range's stride may be 1 or greater
     * than 1.
     *
     * @param theMatrix Matrix.
     * @param theColRange Range of columns to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null or
     * <code>theColRange</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theMatrix</code>'s allocation does not include <code>theColRange</code>.
     */
    public static FloatBuf colSliceBuffer(float[][] theMatrix,
            Range theColRange) {
        if (theMatrix == null) {
            throw new NullPointerException("FloatBuf.colSliceBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        int nc = Arrays.colLength(theMatrix, 0);
        if (0 > theColRange.lb() || theColRange.ub() >= nc) {
            throw new IndexOutOfBoundsException("FloatBuf.colSliceBuffer(): theMatrix column index range = 0.."
                    + (nc - 1) + ", theColRange = " + theColRange);
        }
        if (theColRange.stride() == 1) {
            return new FloatMatrixBuf_1(theMatrix, new Range(0, nr - 1), theColRange);
        } else {
            return new FloatMatrixBuf(theMatrix, new Range(0, nr - 1), theColRange);
        }
    }

    /**
     * Create an array of buffers for multiple column slices of the given float
     * matrix. The returned buffer array has the same length as
     * <code>theColRanges</code>. Each element [<I>i</I>] of the returned buffer
     * array encompasses all the rows of <code>theMatrix</code> and the columns of
     * <code>theMatrix</code> specified by <code>theColRanges[i]</code>. Each range's
     * stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theColRanges Array of ranges of columns to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null or
     * <code>theColRanges</code> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theMatrix</code>'s allocation does not include any element of
     * <code>theColRanges</code>.
     */
    public static FloatBuf[] colSliceBuffers(float[][] theMatrix,
            Range[] theColRanges) {
        int n = theColRanges.length;
        FloatBuf[] result = new FloatBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = colSliceBuffer(theMatrix, theColRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for one patch of the given float matrix. The returned
     * buffer encompasses <code>theRowRange</code> of rows, and <code>theColRange</code>
     * of columns, in <code>theMatrix</code>. Each range's stride may be 1 or
     * greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRange Range of rows to include.
     * @param theColRange Range of columns to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null,
     * <code>theRowRange</code> is null, or <code>theColRange</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theMatrix</code>'s allocation does not include <code>theRowRange</code> and
     * <code>theColRange</code>.
     */
    public static FloatBuf patchBuffer(float[][] theMatrix,
            Range theRowRange,
            Range theColRange) {
        if (theMatrix == null) {
            throw new NullPointerException("FloatBuf.patchBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        if (0 > theRowRange.lb() || theRowRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("FloatBuf.patchBuffer(): theMatrix row index range = 0.."
                    + (nr - 1) + ", theRowRange = " + theRowRange);
        }
        int nc = Arrays.colLength(theMatrix, theRowRange.lb());
        if (0 > theColRange.lb() || theColRange.ub() >= nc) {
            throw new IndexOutOfBoundsException("FloatBuf.patchBuffer(): theMatrix column index range = 0.."
                    + (nc - 1) + ", theColRange = " + theColRange);
        }
        if (theRowRange.stride() == 1 && theColRange.stride() == 1) {
            return new FloatMatrixBuf_1(theMatrix, theRowRange, theColRange);
        } else {
            return new FloatMatrixBuf(theMatrix, theRowRange, theColRange);
        }
    }

    /**
     * Create an array of buffers for multiple patches of the given float
     * matrix. The length of the returned buffer array is equal to the length of
     * <code>theRowRanges</code> times the length of <code>theColRanges</code>. Each
     * element of the returned buffer array encompasses the rows given in one
     * element of <code>theRowRanges</code> array, and the columns given in one
     * element of <code>theColRanges</code> array, in all possible combinations, of
     * <code>theMatrix</code>. Each range's stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRanges Array of ranges of rows to include.
     * @param theColRanges Array of ranges of columns to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null,
     * <code>theRowRanges</code> or any element thereof is null, or
     * <code>theColRanges</code> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theMatrix</code>'s allocation does not include any element of
     * <code>theRowRanges</code> or
     * <code>theColRanges</code>.
     */
    public static FloatBuf[] patchBuffers(float[][] theMatrix,
            Range[] theRowRanges,
            Range[] theColRanges) {
        int m = theRowRanges.length;
        int n = theColRanges.length;
        FloatBuf[] result = new FloatBuf[m * n];
        int k = 0;
        for (int i = 0; i < m; ++i) {
            Range rowrange = theRowRanges[i];
            for (int j = 0; j < n; ++j) {
                result[k++]
                        = patchBuffer(theMatrix, rowrange, theColRanges[j]);
            }
        }
        return result;
    }

    /**
     * Create a buffer for a shared float item. The item is wrapped in an
     * instance of class {@linkplain edu.rit.pj.reduction.SharedFloat
     * SharedFloat}. Use the methods of the SharedFloat object to access the
     * actual item.
     *
     * @param item SharedFloat object that wraps the item.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>item</code> is null.
     * @return a {@link edu.rit.mp.FloatBuf} object.
     */
    public static FloatBuf buffer(SharedFloat item) {
        if (item == null) {
            throw new NullPointerException("FloatBuf.buffer(): item is null");
        }
        return new SharedFloatBuf(item);
    }

    /**
     * Create a buffer for the entire given shared float array. The returned
     * buffer encompasses all the elements in <code>theArray</code>.
     *
     * @param theArray Array.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null.
     */
    public static FloatBuf buffer(SharedFloatArray theArray) {
        if (theArray == null) {
            throw new NullPointerException("FloatBuf.buffer(): theArray is null");
        }
        int nr = theArray.length();
        return new SharedFloatArrayBuf_1(theArray, new Range(0, nr - 1));
    }

    /**
     * Create a buffer for one slice of the given shared float array. The
     * returned buffer encompasses <code>theRange</code> of elements in
     * <code>theArray</code>. The range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRange Range of elements to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null or
     * <code>theRange</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theArray</code> does not include all the indexes in <code>theRange</code>.
     */
    public static FloatBuf sliceBuffer(SharedFloatArray theArray,
            Range theRange) {
        if (theArray == null) {
            throw new NullPointerException("FloatBuf.sliceBuffer(): theArray is null");
        }
        int nr = theArray.length();
        if (0 > theRange.lb() || theRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("FloatBuf.sliceBuffer(): theArray index range = 0.."
                    + (nr - 1) + ", theRange = " + theRange);
        }
        if (theRange.stride() == 1) {
            return new SharedFloatArrayBuf_1(theArray, theRange);
        } else {
            return new SharedFloatArrayBuf(theArray, theRange);
        }
    }

    /**
     * Create an array of buffers for multiple slices of the given shared float
     * array. The returned buffer array has the same length as
     * <code>theRanges</code>. Each element [<I>i</I>] of the returned buffer array
     * encompasses the elements of <code>theArray</code> specified by
     * <code>theRanges[i]</code>. Each range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRanges Array of ranges of elements to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null or
     * <code>theRanges</code> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>theArray</code>'s allocation does not include any element of
     * <code>theRanges</code>.
     */
    public static FloatBuf[] sliceBuffers(SharedFloatArray theArray,
            Range[] theRanges) {
        int n = theRanges.length;
        FloatBuf[] result = new FloatBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = sliceBuffer(theArray, theRanges[i]);
        }
        return result;
    }

    /**
     * Obtain the given item from this buffer.
     * <P>
     * The <code>get()</code> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <code>length()</code>-1.
     * @return Item at index <code>i</code>.
     */
    public abstract float get(int i);

    /**
     * Store the given item in this buffer.
     * <P>
     * The <code>put()</code> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <code>length()</code>-1.
     * @param item Item to be stored at index <code>i</code>.
     */
    public abstract void put(int i,
            float item);

    /**
     * {@inheritDoc}
     *
     * Copy items from the given buffer to this buffer. The number of items
     * copied is this buffer's length or <code>theSrc</code>'s length, whichever is
     * smaller. If <code>theSrc</code> is this buffer, the <code>copy()</code> method
     * does nothing.
     * <P>
     * The default implementation of the <code>copy()</code> method calls the
     * <code>defaultCopy()</code> method. A subclass can override the
     * <code>copy()</code> method to use a more efficient algorithm.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <code>theSrc</code>'s item data type is not the same as this buffer's item
     * data type.
     */
    public void copy(Buf theSrc) {
        if (theSrc != this) {
            defaultCopy((FloatBuf) theSrc, this);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Fill this buffer with the given item. The <code>item</code> is assigned to
     * each element in this buffer.
     * <P>
     * The <code>item</code> must be an instance of class Float. If the
     * <code>item</code> is null, 0.0f is assigned to each element in this buffer.
     * @exception ClassCastException (unchecked exception) Thrown if the
     * <code>item</code>'s data type is not the same as this buffer's item data
     * type.
     */
    public void fill(Object item) {
        float value = item == null ? 0.0f : ((Float) item).floatValue();
        for (int i = 0; i < myLength; ++i) {
            put(i, value);
        }
    }

    /**
     * Create a temporary buffer with the same type of items and the same length
     * as this buffer. The new buffer items are stored in a newly created array,
     * separate from the storage for this buffer's items.
     *
     * @return a {@link edu.rit.mp.Buf} object.
     */
    public Buf getTemporaryBuf() {
        return buffer(new float[myLength]);
    }

// Hidden operations.
    /**
     * Skip as many items as possible from the given byte buffer.
     *
     * @param num Number of items to skip.
     * @param buffer Buffer.
     *
     * @return Number of items actually skipped.
     */
    int skipItems(int num,
            ByteBuffer buffer) {
        int n = Math.min(num, buffer.remaining() / 4);
        buffer.position(buffer.position() + 4 * n);
        return n;
    }

    /**
     * Copy items from the given source buffer to the given destination buffer.
     * The number of items copied is <code>theSrc</code>'s length or
     * <code>theDst</code>'s length, whichever is smaller. Each item is copied
     * individually using the <code>get()</code> and <code>put()</code> methods. It is
     * assumed that <code>theSrc</code> is not the same as <code>theDst</code>.
     *
     * @param theSrc Source of items to copy.
     * @param theDst Destination of items to copy.
     */
    protected static void defaultCopy(FloatBuf theSrc,
            FloatBuf theDst) {
        int n = Math.min(theSrc.myLength, theDst.myLength);
        for (int i = 0; i < n; ++i) {
            theDst.put(i, theSrc.get(i));
        }
    }

}
