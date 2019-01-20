//******************************************************************************
//
// File:    BooleanBuf.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.BooleanBuf
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

import edu.rit.mp.buf.BooleanArrayBuf;
import edu.rit.mp.buf.BooleanArrayBuf_1;
import edu.rit.mp.buf.BooleanItemBuf;
import edu.rit.mp.buf.BooleanMatrixBuf;
import edu.rit.mp.buf.BooleanMatrixBuf_1;
import edu.rit.mp.buf.EmptyBooleanBuf;
import edu.rit.mp.buf.SharedBooleanArrayBuf;
import edu.rit.mp.buf.SharedBooleanArrayBuf_1;
import edu.rit.mp.buf.SharedBooleanBuf;
import edu.rit.pj.reduction.SharedBoolean;
import edu.rit.pj.reduction.SharedBooleanArray;
import edu.rit.util.Arrays;
import edu.rit.util.Range;

/**
 * Class BooleanBuf is the abstract base class for a buffer of Boolean items
 * sent or received using the Message Protocol (MP). In a message, a Boolean
 * item is represented as one byte, 0 for false, 1 for true.
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
 * To create a BooleanBuf, call one of the following static factory methods:
 * <UL>
 * <LI><code>emptyBuffer()</code>
 * <LI><code>buffer()</code>
 * <LI><code>buffer (boolean)</code>
 * <LI><code>buffer (boolean[])</code>
 * <LI><code>sliceBuffer (boolean[], Range)</code>
 * <LI><code>sliceBuffers (boolean[], Range[])</code>
 * <LI><code>buffer (boolean[][])</code>
 * <LI><code>rowSliceBuffer (boolean[][], Range)</code>
 * <LI><code>rowSliceBuffers (boolean[][], Range[])</code>
 * <LI><code>colSliceBuffer (boolean[][], Range)</code>
 * <LI><code>colSliceBuffers (boolean[][], Range[])</code>
 * <LI><code>patchBuffer (boolean[][], Range, Range)</code>
 * <LI><code>patchBuffers (boolean[][], Range[], Range[])</code>
 * <LI><code>buffer (SharedBoolean)</code>
 * <LI><code>buffer (SharedBooleanArray)</code>
 * <LI><code>sliceBuffer (SharedBooleanArray, Range)</code>
 * <LI><code>sliceBuffers (SharedBooleanArray, Range[])</code>
 * </UL>
 *
 * @author Alan Kaminsky
 * @version 03-May-2008
 */
public abstract class BooleanBuf
        extends Buf {

// Hidden constructors.
    /**
     * Construct a new Boolean buffer.
     *
     * @param theLength Number of items.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>theLength</code> &lt; 0.
     */
    protected BooleanBuf(int theLength) {
        super(Constants.TYPE_BOOLEAN, theLength);
    }

// Exported operations.
    /**
     * Create an empty buffer. The buffer's length is 0. The buffer's item type
     * is Boolean.
     *
     * @return Empty buffer.
     */
    public static BooleanBuf emptyBuffer() {
        return new EmptyBooleanBuf();
    }

    /**
     * Create a buffer for a Boolean item. The item is stored in the
     * <code>item</code> field of the buffer.
     *
     * @return Buffer.
     */
    public static BooleanItemBuf buffer() {
        return new BooleanItemBuf();
    }

    /**
     * Create a buffer for a Boolean item with the given initial value. The item
     * is stored in the <code>item</code> field of the buffer.
     *
     * @param item Initial value of the <code>item</code> field.
     * @return Buffer.
     */
    public static BooleanItemBuf buffer(boolean item) {
        return new BooleanItemBuf(item);
    }

    /**
     * Create a buffer for the entire given Boolean array. The returned buffer
     * encompasses all the elements in <code>theArray</code>.
     *
     * @param theArray Array.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null.
     */
    public static BooleanBuf buffer(boolean[] theArray) {
        if (theArray == null) {
            throw new NullPointerException("BooleanBuf.buffer(): theArray is null");
        }
        int nr = Arrays.length(theArray);
        return new BooleanArrayBuf_1(theArray, new Range(0, nr - 1));
    }

    /**
     * Create a buffer for one slice of the given Boolean array. The returned
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
    public static BooleanBuf sliceBuffer(boolean[] theArray,
            Range theRange) {
        if (theArray == null) {
            throw new NullPointerException("BooleanBuf.sliceBuffer(): theArray is null");
        }
        int nr = Arrays.length(theArray);
        if (0 > theRange.lb() || theRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("BooleanBuf.sliceBuffer(): theArray index range = 0.."
                    + (nr - 1) + ", theRange = " + theRange);
        }
        if (theRange.stride() == 1) {
            return new BooleanArrayBuf_1(theArray, theRange);
        } else {
            return new BooleanArrayBuf(theArray, theRange);
        }
    }

    /**
     * Create an array of buffers for multiple slices of the given Boolean
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
    public static BooleanBuf[] sliceBuffers(boolean[] theArray,
            Range[] theRanges) {
        int n = theRanges.length;
        BooleanBuf[] result = new BooleanBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = sliceBuffer(theArray, theRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for the entire given Boolean matrix. The returned buffer
     * encompasses all the rows and all the columns in <code>theMatrix</code>.
     *
     * @param theMatrix Matrix.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theMatrix</code> is null.
     */
    public static BooleanBuf buffer(boolean[][] theMatrix) {
        if (theMatrix == null) {
            throw new NullPointerException("BooleanBuf.buffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        int nc = Arrays.colLength(theMatrix, 0);
        return new BooleanMatrixBuf_1(theMatrix, new Range(0, nr - 1), new Range(0, nc - 1));
    }

    /**
     * Create a buffer for one row slice of the given Boolean matrix. The
     * returned buffer encompasses <code>theRowRange</code> of rows, and all the
     * columns, in <code>theMatrix</code>. The range's stride may be 1 or greater
     * than 1.
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
    public static BooleanBuf rowSliceBuffer(boolean[][] theMatrix,
            Range theRowRange) {
        if (theMatrix == null) {
            throw new NullPointerException("BooleanBuf.rowSliceBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        if (0 > theRowRange.lb() || theRowRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("BooleanBuf.rowSliceBuffer(): theMatrix row index range = 0.."
                    + (nr - 1) + ", theRowRange = " + theRowRange);
        }
        int nc = Arrays.colLength(theMatrix, theRowRange.lb());
        if (theRowRange.stride() == 1) {
            return new BooleanMatrixBuf_1(theMatrix, theRowRange, new Range(0, nc - 1));
        } else {
            return new BooleanMatrixBuf(theMatrix, theRowRange, new Range(0, nc - 1));
        }
    }

    /**
     * Create an array of buffers for multiple row slices of the given Boolean
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
    public static BooleanBuf[] rowSliceBuffers(boolean[][] theMatrix,
            Range[] theRowRanges) {
        int n = theRowRanges.length;
        BooleanBuf[] result = new BooleanBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = rowSliceBuffer(theMatrix, theRowRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for one column slice of the given Boolean matrix. The
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
    public static BooleanBuf colSliceBuffer(boolean[][] theMatrix,
            Range theColRange) {
        if (theMatrix == null) {
            throw new NullPointerException("BooleanBuf.colSliceBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        int nc = Arrays.colLength(theMatrix, 0);
        if (0 > theColRange.lb() || theColRange.ub() >= nc) {
            throw new IndexOutOfBoundsException("BooleanBuf.colSliceBuffer(): theMatrix column index range = 0.."
                    + (nc - 1) + ", theColRange = " + theColRange);
        }
        if (theColRange.stride() == 1) {
            return new BooleanMatrixBuf_1(theMatrix, new Range(0, nr - 1), theColRange);
        } else {
            return new BooleanMatrixBuf(theMatrix, new Range(0, nr - 1), theColRange);
        }
    }

    /**
     * Create an array of buffers for multiple column slices of the given
     * Boolean matrix. The returned buffer array has the same length as
     * <code>theColRanges</code>. Each element [<I>i</I>] of the returned buffer
     * array encompasses all the rows of <code>theMatrix</code> and the columns of
     * <code>theMatrix</code> specified by <code>theColRanges[i]</code>. Each range's
     * stride may be 1 or greater than 1. It is assumed that the rows and
     * columns of <code>theMatrix</code> are allocated and that each row of
     * <code>theMatrix</code> has the same number of columns.
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
    public static BooleanBuf[] colSliceBuffers(boolean[][] theMatrix,
            Range[] theColRanges) {
        int n = theColRanges.length;
        BooleanBuf[] result = new BooleanBuf[n];
        for (int i = 0; i < n; ++i) {
            result[i] = colSliceBuffer(theMatrix, theColRanges[i]);
        }
        return result;
    }

    /**
     * Create a buffer for one patch of the given Boolean matrix. The returned
     * buffer encompasses <code>theRowRange</code> of rows, and <code>theColRange</code>
     * of columns, in <code>theMatrix</code>. Each range's stride may be 1 or
     * greater than 1. It is assumed that the rows and columns of
     * <code>theMatrix</code> are allocated and that each row of <code>theMatrix</code>
     * has the same number of columns.
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
    public static BooleanBuf patchBuffer(boolean[][] theMatrix,
            Range theRowRange,
            Range theColRange) {
        if (theMatrix == null) {
            throw new NullPointerException("BooleanBuf.patchBuffer(): theMatrix is null");
        }
        int nr = Arrays.rowLength(theMatrix);
        if (0 > theRowRange.lb() || theRowRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("BooleanBuf.patchBuffer(): theMatrix row index range = 0.."
                    + (nr - 1) + ", theRowRange = " + theRowRange);
        }
        int nc = Arrays.colLength(theMatrix, theRowRange.lb());
        if (0 > theColRange.lb() || theColRange.ub() >= nc) {
            throw new IndexOutOfBoundsException("BooleanBuf.colSliceBuffer(): theMatrix column index range = 0.."
                    + (nc - 1) + ", theColRange = " + theColRange);
        }
        if (theRowRange.stride() == 1 && theColRange.stride() == 1) {
            return new BooleanMatrixBuf_1(theMatrix, theRowRange, theColRange);
        } else {
            return new BooleanMatrixBuf(theMatrix, theRowRange, theColRange);
        }
    }

    /**
     * Create an array of buffers for multiple patches of the given Boolean
     * matrix. The length of the returned buffer array is equal to the length of
     * <code>theRowRanges</code> times the length of <code>theColRanges</code>. Each
     * element of the returned buffer array encompasses the rows given in one
     * element of <code>theRowRanges</code> array, and the columns given in one
     * element of <code>theColRanges</code> array, in all possible combinations, of
     * <code>theMatrix</code>. Each range's stride may be 1 or greater than 1. It is
     * assumed that the rows and columns of <code>theMatrix</code> are allocated and
     * that each row of <code>theMatrix</code> has the same number of columns.
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
    public static BooleanBuf[] patchBuffers(boolean[][] theMatrix,
            Range[] theRowRanges,
            Range[] theColRanges) {
        int m = theRowRanges.length;
        int n = theColRanges.length;
        BooleanBuf[] result = new BooleanBuf[m * n];
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
     * Create a buffer for a shared Boolean item. The item is wrapped in an
     * instance of class {@linkplain edu.rit.pj.reduction.SharedBoolean
     * SharedBoolean}. Use the methods of the SharedBoolean object to access the
     * actual item.
     *
     * @param item SharedBoolean object that wraps the item.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>item</code> is null.
     * @return a {@link edu.rit.mp.BooleanBuf} object.
     */
    public static BooleanBuf buffer(SharedBoolean item) {
        if (item == null) {
            throw new NullPointerException("BooleanBuf.buffer(): item is null");
        }
        return new SharedBooleanBuf(item);
    }

    /**
     * Create a buffer for the entire given shared Boolean array. The returned
     * buffer encompasses all the elements in <code>theArray</code>.
     *
     * @param theArray Array.
     * @return Buffer.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null.
     */
    public static BooleanBuf buffer(SharedBooleanArray theArray) {
        if (theArray == null) {
            throw new NullPointerException("BooleanBuf.buffer(): theArray is null");
        }
        int nr = theArray.length();
        return new SharedBooleanArrayBuf_1(theArray, new Range(0, nr - 1));
    }

    /**
     * Create a buffer for one slice of the given shared Boolean array. The
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
    public static BooleanBuf sliceBuffer(SharedBooleanArray theArray,
            Range theRange) {
        if (theArray == null) {
            throw new NullPointerException("BooleanBuf.sliceBuffer(): theArray is null");
        }
        int nr = theArray.length();
        if (0 > theRange.lb() || theRange.ub() >= nr) {
            throw new IndexOutOfBoundsException("BooleanBuf.sliceBuffer(): theArray row index range = 0.."
                    + (nr - 1) + ", theRange = " + theRange);
        }
        if (theRange.stride() == 1) {
            return new SharedBooleanArrayBuf_1(theArray, theRange);
        } else {
            return new SharedBooleanArrayBuf(theArray, theRange);
        }
    }

    /**
     * Create an array of buffers for multiple slices of the given shared
     * Boolean array. The returned buffer array has the same length as
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
    public static BooleanBuf[] sliceBuffers(SharedBooleanArray theArray,
            Range[] theRanges) {
        int n = theRanges.length;
        BooleanBuf[] result = new BooleanBuf[n];
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
    public abstract boolean get(int i);

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
            boolean item);

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
            defaultCopy((BooleanBuf) theSrc, this);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Fill this buffer with the given item. The <code>item</code> is assigned to
     * each element in this buffer.
     * <P>
     * The <code>item</code> must be an instance of class Boolean. If the
     * <code>item</code> is null, false is assigned to each element in this buffer.
     * @exception ClassCastException (unchecked exception) Thrown if the
     * <code>item</code>'s data type is not the same as this buffer's item data
     * type.
     */
    public void fill(Object item) {
        boolean value = item == null ? false : ((Boolean) item).booleanValue();
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
        return buffer(new boolean[myLength]);
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
        int n = Math.min(num, buffer.remaining());
        buffer.position(buffer.position() + n);
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
    protected static void defaultCopy(BooleanBuf theSrc,
            BooleanBuf theDst) {
        int n = Math.min(theSrc.myLength, theDst.myLength);
        for (int i = 0; i < n; ++i) {
            theDst.put(i, theSrc.get(i));
        }
    }

}
