//******************************************************************************
//
// File:    ShortBuf.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.ShortBuf
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

import edu.rit.mp.buf.EmptyShortBuf;
import edu.rit.mp.buf.SharedShortArrayBuf;
import edu.rit.mp.buf.SharedShortArrayBuf_1;
import edu.rit.mp.buf.SharedShortBuf;
import edu.rit.mp.buf.ShortArrayBuf;
import edu.rit.mp.buf.ShortArrayBuf_1;
import edu.rit.mp.buf.ShortItemBuf;
import edu.rit.mp.buf.ShortMatrixBuf;
import edu.rit.mp.buf.ShortMatrixBuf_1;
import edu.rit.pj.reduction.SharedShort;
import edu.rit.pj.reduction.SharedShortArray;
import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class ShortBuf is the abstract base class for a buffer of short items sent or
 * received using the Message Protocol (MP). In a message, a short item is
 * represented as two bytes, most significant byte first.
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
 * To create a ShortBuf, call one of the following static factory methods:
 * <UL>
 * <LI><TT>emptyBuffer()</TT>
 * <LI><TT>buffer()</TT>
 * <LI><TT>buffer (short)</TT>
 * <LI><TT>buffer (short[])</TT>
 * <LI><TT>sliceBuffer (short[], Range)</TT>
 * <LI><TT>sliceBuffers (short[], Range[])</TT>
 * <LI><TT>buffer (short[][])</TT>
 * <LI><TT>rowSliceBuffer (short[][], Range)</TT>
 * <LI><TT>rowSliceBuffers (short[][], Range[])</TT>
 * <LI><TT>colSliceBuffer (short[][], Range)</TT>
 * <LI><TT>colSliceBuffers (short[][], Range[])</TT>
 * <LI><TT>patchBuffer (short[][], Range, Range)</TT>
 * <LI><TT>patchBuffers (short[][], Range[], Range[])</TT>
 * <LI><TT>buffer (SharedShort)</TT>
 * <LI><TT>buffer (SharedShortArray)</TT>
 * <LI><TT>sliceBuffer (SharedShortArray, Range)</TT>
 * <LI><TT>sliceBuffers (SharedShortArray, Range[])</TT>
 * </UL>
 *
 * @author Alan Kaminsky
 * @version 03-May-2008
 */
public abstract class ShortBuf
        extends Buf {

// Hidden constructors.
    /**
     * Construct a new short buffer.
     *
     * @param theLength Number of items.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theLength</TT> &lt; 0.
     */
    protected ShortBuf(int theLength) {
        super(Constants.TYPE_SHORT, theLength);
    }

// Exported operations.
    /**
     * Create an empty buffer. The buffer's length is 0. The buffer's item type
     * is short.
     *
     * @return Empty buffer.
     */
    public static ShortBuf emptyBuffer() {
        return new EmptyShortBuf();
    }

    /**
     * Create a buffer for a short item. The item is stored in the
     * <TT>item</TT> field of the buffer.
     *
     * @return Buffer.
     */
    public static ShortItemBuf buffer() {
        return new ShortItemBuf();
    }

    /**
     * Create a buffer for a short item with the given initial value. The item
     * is stored in the <TT>item</TT> field of the buffer.
     *
     * @param item Initial value of the <TT>item</TT> field.
     * @return Buffer.
     */
    public static ShortItemBuf buffer(short item) {
        return new ShortItemBuf(item);
    }

    /**
     * Create a buffer for the entire given short array. The returned buffer
     * encompasses all the elements in <TT>theArray</TT>.
     *
     * @param theArray Array.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theArray</TT> is null.
     */
    public static ShortBuf buffer(short[] theArray) {
        if (theArray == null) {
            throw new NullPointerException("ShortBuf.buffer(): theArray is null");
        }
        int nr = Arrays.length(theArray);
        return new ShortArrayBuf_1(theArray, new Range(0, nr - 1));
    }

    /**
     * Create a buffer for one slice of the given short array. The returned
     * buffer encompasses <TT>theRange</TT> of elements in <TT>theArray</TT>.
     * The range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRange Range of elements to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theArray</TT> is null or
     * <TT>theRange</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theArray</TT> does not include all the indexes in <TT>theRange</TT>.
     */
    public static ShortBuf sliceBuffer(short[] theArray,
            Range theRange) {
        if (theArray == null) {
            throw new NullPointerException("ShortBuf.sliceBuffer(): theArray is null");
        }
        int nr = Arrays.length(theArray);
        if (0 > theRange.lb() || theRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("ShortBuf.sliceBuffer(): theArray index range = 0.."
                    + (nr - 1) + ", theRange = " + theRange);
        }
        if (theRange.stride() == 1) {
            return new ShortArrayBuf_1(theArray, theRange);
        } else {
            return new ShortArrayBuf(theArray, theRange);
        }
    }

    /**
     * Create an array of buffers for multiple slices of the given short array.
     * The returned buffer array has the same length as
     * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
     * encompasses the elements of <TT>theArray</TT> specified by
     * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRanges Array of ranges of elements to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theArray</TT> is null or
     * <TT>theRanges</TT> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theArray</TT>'s allocation does not include any element of
     * <TT>theRanges</TT>.
     */
    public static ShortBuf[] sliceBuffers(short[] theArray,
            Range[] theRanges) {
        int n = theRanges.length;
        ShortBuf[] result = new ShortBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = sliceBuffer(theArray, theRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for the entire given short matrix. The returned buffer
     * encompasses all the rows and all the columns in
     * <TT>theMatrix</TT>.
     *
     * @param theMatrix Matrix.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null.
     */
    public static ShortBuf buffer(short[][] theMatrix) {
        if (theMatrix == null) {
            throw new NullPointerException("ShortBuf.buffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        int nc = Arrays.colLength(theMatrix, 0);
        return new ShortMatrixBuf_1(theMatrix, new Range(0, nr - 1), new Range(0, nc - 1));
    }

    /**
     * Create a buffer for one row slice of the given short matrix. The returned
     * buffer encompasses <TT>theRowRange</TT> of rows, and all the columns, in
     * <TT>theMatrix</TT>. The range's stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRange Range of rows to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null or
     * <TT>theRowRange</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theMatrix</TT>'s allocation does not include <TT>theRowRange</TT>.
     */
    public static ShortBuf rowSliceBuffer(short[][] theMatrix,
            Range theRowRange) {
        if (theMatrix == null) {
            throw new NullPointerException("ShortBuf.rowSliceBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        if (0 > theRowRange.lb() || theRowRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("ShortBuf.rowSliceBuffer(): theMatrix row index range = 0.."
                    + (nr - 1) + ", theRowRange = " + theRowRange);
        }
        int nc = Arrays.colLength(theMatrix, theRowRange.lb());
        if (theRowRange.stride() == 1) {
            return new ShortMatrixBuf_1(theMatrix, theRowRange, new Range(0, nc - 1));
        } else {
            return new ShortMatrixBuf(theMatrix, theRowRange, new Range(0, nc - 1));
        }
    }

    /**
     * Create an array of buffers for multiple row slices of the given short
     * matrix. The returned buffer array has the same length as
     * <TT>theRowRanges</TT>. Each element [<I>i</I>] of the returned buffer
     * array encompasses the rows of <TT>theMatrix</TT> specified by
     * <TT>theRowRanges[i]</TT> and all the columns of <TT>theMatrix</TT>. Each
     * range's stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRanges Array of ranges of rows to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null or
     * <TT>theRowRanges</TT> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theMatrix</TT>'s allocation does not include any element of
     * <TT>theRowRanges</TT>.
     */
    public static ShortBuf[] rowSliceBuffers(short[][] theMatrix,
            Range[] theRowRanges) {
        int n = theRowRanges.length;
        ShortBuf[] result = new ShortBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = rowSliceBuffer(theMatrix, theRowRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for one column slice of the given short matrix. The
     * returned buffer encompasses all the rows, and <TT>theColRange</TT> of
     * columns, in <TT>theMatrix</TT>. The range's stride may be 1 or greater
     * than 1.
     *
     * @param theMatrix Matrix.
     * @param theColRange Range of columns to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null or
     * <TT>theColRange</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theMatrix</TT>'s allocation does not include <TT>theColRange</TT>.
     */
    public static ShortBuf colSliceBuffer(short[][] theMatrix,
            Range theColRange) {
        if (theMatrix == null) {
            throw new NullPointerException("ShortBuf.colSliceBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        int nc = Arrays.colLength(theMatrix, 0);
        if (0 > theColRange.lb() || theColRange.ub() >= nc) {
            throw new IndexOutOfBoundsException("ShortBuf.colSliceBuffer(): theMatrix column index range = 0.."
                    + (nc - 1) + ", theColRange = " + theColRange);
        }
        if (theColRange.stride() == 1) {
            return new ShortMatrixBuf_1(theMatrix, new Range(0, nr - 1), theColRange);
        } else {
            return new ShortMatrixBuf(theMatrix, new Range(0, nr - 1), theColRange);
        }
    }

    /**
     * Create an array of buffers for multiple column slices of the given short
     * matrix. The returned buffer array has the same length as
     * <TT>theColRanges</TT>. Each element [<I>i</I>] of the returned buffer
     * array encompasses all the rows of <TT>theMatrix</TT> and the columns of
     * <TT>theMatrix</TT> specified by <TT>theColRanges[i]</TT>. Each range's
     * stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theColRanges Array of ranges of columns to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null or
     * <TT>theColRanges</TT> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theMatrix</TT>'s allocation does not include any element of
     * <TT>theColRanges</TT>.
     */
    public static ShortBuf[] colSliceBuffers(short[][] theMatrix,
            Range[] theColRanges) {
        int n = theColRanges.length;
        ShortBuf[] result = new ShortBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = colSliceBuffer(theMatrix, theColRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for one patch of the given short matrix. The returned
     * buffer encompasses <TT>theRowRange</TT> of rows, and <TT>theColRange</TT>
     * of columns, in <TT>theMatrix</TT>. Each range's stride may be 1 or
     * greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRange Range of rows to include.
     * @param theColRange Range of columns to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null,
     * <TT>theRowRange</TT> is null, or <TT>theColRange</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theMatrix</TT>'s allocation does not include <TT>theRowRange</TT> and
     * <TT>theColRange</TT>.
     */
    public static ShortBuf patchBuffer(short[][] theMatrix,
            Range theRowRange,
            Range theColRange) {
        if (theMatrix == null) {
            throw new NullPointerException("ShortBuf.patchBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        if (0 > theRowRange.lb() || theRowRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("ShortBuf.patchBuffer(): theMatrix row index range = 0.."
                    + (nr - 1) + ", theRowRange = " + theRowRange);
        }
        int nc = Arrays.colLength(theMatrix, theRowRange.lb());
        if (0 > theColRange.lb() || theColRange.ub() >= nc) {
            throw new IndexOutOfBoundsException("ShortBuf.patchBuffer(): theMatrix column index range = 0.."
                    + (nc - 1) + ", theColRange = " + theColRange);
        }
        if (theRowRange.stride() == 1 && theColRange.stride() == 1) {
            return new ShortMatrixBuf_1(theMatrix, theRowRange, theColRange);
        } else {
            return new ShortMatrixBuf(theMatrix, theRowRange, theColRange);
        }
    }

    /**
     * Create an array of buffers for multiple patches of the given short
     * matrix. The length of the returned buffer array is equal to the length of
     * <TT>theRowRanges</TT> times the length of <TT>theColRanges</TT>. Each
     * element of the returned buffer array encompasses the rows given in one
     * element of <TT>theRowRanges</TT> array, and the columns given in one
     * element of <TT>theColRanges</TT> array, in all possible combinations, of
     * <TT>theMatrix</TT>. Each range's stride may be 1 or greater than 1.
     *
     * @param theMatrix Matrix.
     * @param theRowRanges Array of ranges of rows to include.
     * @param theColRanges Array of ranges of columns to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theMatrix</TT> is null,
     * <TT>theRowRanges</TT> or any element thereof is null, or
     * <TT>theColRanges</TT> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theMatrix</TT>'s allocation does not include any element of
     * <TT>theRowRanges</TT> or
     * <TT>theColRanges</TT>.
     */
    public static ShortBuf[] patchBuffers(short[][] theMatrix,
            Range[] theRowRanges,
            Range[] theColRanges) {
        int m = theRowRanges.length;
        int n = theColRanges.length;
        ShortBuf[] result = new ShortBuf[m * n];
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
     * Create a buffer for a shared short item. The item is wrapped in an
     * instance of class {@linkplain edu.rit.pj.reduction.SharedShort
     * SharedShort}. Use the methods of the SharedShort object to access the
     * actual item.
     *
     * @param item SharedShort object that wraps the item.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>item</TT> is null.
     * @return a {@link edu.rit.mp.ShortBuf} object.
     */
    public static ShortBuf buffer(SharedShort item) {
        if (item == null) {
            throw new NullPointerException("ShortBuf.buffer(): item is null");
        }
        return new SharedShortBuf(item);
    }

    /**
     * Create a buffer for the entire given shared short array. The returned
     * buffer encompasses all the elements in <TT>theArray</TT>.
     *
     * @param theArray Array.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theArray</TT> is null.
     */
    public static ShortBuf buffer(SharedShortArray theArray) {
        if (theArray == null) {
            throw new NullPointerException("ShortBuf.buffer(): theArray is null");
        }
        int nr = theArray.length();
        return new SharedShortArrayBuf_1(theArray, new Range(0, nr - 1));
    }

    /**
     * Create a buffer for one slice of the given shared short array. The
     * returned buffer encompasses <TT>theRange</TT> of elements in
     * <TT>theArray</TT>. The range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRange Range of elements to include.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theArray</TT> is null or
     * <TT>theRange</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theArray</TT> does not include all the indexes in <TT>theRange</TT>.
     */
    public static ShortBuf sliceBuffer(SharedShortArray theArray,
            Range theRange) {
        if (theArray == null) {
            throw new NullPointerException("ShortBuf.sliceBuffer(): theArray is null");
        }
        int nr = theArray.length();
        if (0 > theRange.lb() || theRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("ShortBuf.sliceBuffer(): theArray index range = 0.."
                    + (nr - 1) + ", theRange = " + theRange);
        }
        if (theRange.stride() == 1) {
            return new SharedShortArrayBuf_1(theArray, theRange);
        } else {
            return new SharedShortArrayBuf(theArray, theRange);
        }
    }

    /**
     * Create an array of buffers for multiple slices of the given shared short
     * array. The returned buffer array has the same length as
     * <TT>theRanges</TT>. Each element [<I>i</I>] of the returned buffer array
     * encompasses the elements of <TT>theArray</TT> specified by
     * <TT>theRanges[i]</TT>. Each range's stride may be 1 or greater than 1.
     *
     * @param theArray Array.
     * @param theRanges Array of ranges of elements to include.
     * @return Array of buffers.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theArray</TT> is null or
     * <TT>theRanges</TT> or any element thereof is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>theArray</TT>'s allocation does not include any element of
     * <TT>theRanges</TT>.
     */
    public static ShortBuf[] sliceBuffers(SharedShortArray theArray,
            Range[] theRanges) {
        int n = theRanges.length;
        ShortBuf[] result = new ShortBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = sliceBuffer(theArray, theRanges[i]);
        }
        return result;
    }

    /**
     * Obtain the given item from this buffer.
     * <P>
     * The <TT>get()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <TT>length()</TT>-1.
     * @return Item at index <TT>i</TT>.
     */
    public abstract short get(int i);

    /**
     * Store the given item in this buffer.
     * <P>
     * The <TT>put()</TT> method must not block the calling thread; if it does,
     * all message I/O in MP will be blocked.
     *
     * @param i Item index in the range 0 .. <TT>length()</TT>-1.
     * @param item Item to be stored at index <TT>i</TT>.
     * @param item Item to be stored at index <TT>i</TT>.
     */
    public abstract void put(int i,
            short item);

    /**
     * {@inheritDoc}
     *
     * Copy items from the given buffer to this buffer. The number of items
     * copied is this buffer's length or <TT>theSrc</TT>'s length, whichever is
     * smaller. If <TT>theSrc</TT> is this buffer, the <TT>copy()</TT> method
     * does nothing.
     * <P>
     * The default implementation of the <TT>copy()</TT> method calls the
     * <TT>defaultCopy()</TT> method. A subclass can override the
     * <TT>copy()</TT> method to use a more efficient algorithm.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>theSrc</TT>'s item data type is not the same as this buffer's item
     * data type.
     */
    public void copy(Buf theSrc) {
        if (theSrc != this) {
            defaultCopy((ShortBuf) theSrc, this);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Fill this buffer with the given item. The <TT>item</TT> is assigned to
     * each element in this buffer.
     * <P>
     * The <TT>item</TT> must be an instance of class Short. If the
     * <TT>item</TT> is null, 0 is assigned to each element in this buffer.
     * @exception ClassCastException (unchecked exception) Thrown if the
     * <TT>item</TT>'s data type is not the same as this buffer's item data
     * type.
     */
    public void fill(Object item) {
        short value = item == null ? (short) 0 : ((Short) item).shortValue();
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
        return buffer(new short[myLength]);
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
        int n = Math.min(num, buffer.remaining() / 2);
        buffer.position(buffer.position() + 2 * n);
        return n;
    }

    /**
     * Copy items from the given source buffer to the given destination buffer.
     * The number of items copied is <TT>theSrc</TT>'s length or
     * <TT>theDst</TT>'s length, whichever is smaller. Each item is copied
     * individually using the <TT>get()</TT> and <TT>put()</TT> methods. It is
     * assumed that <TT>theSrc</TT> is not the same as <TT>theDst</TT>.
     *
     * @param theSrc Source of items to copy.
     * @param theDst Destination of items to copy.
     */
    protected static void defaultCopy(ShortBuf theSrc,
            ShortBuf theDst) {
        int n = Math.min(theSrc.myLength, theDst.myLength);
        for (int i = 0; i < n; ++i) {
            theDst.put(i, theSrc.get(i));
        }
    }

}
