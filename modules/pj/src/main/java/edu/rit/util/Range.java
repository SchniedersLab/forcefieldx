//******************************************************************************
//
// File:    Range.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Range
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
package edu.rit.util;

import java.io.Externalizable;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Range provides a range of type <code>int</code>. A range object has the
 * following attributes: <B>lower bound</B> <I>L</I>, <B>upper bound</B>
 * <I>U</I>, <B>stride</B> <I>S</I>, and <B>length</B> <I>N</I>. A range object
 * represents the following set of integers: {<I>L</I>, <I>L</I>+<I>S</I>,
 * <I>L</I>+2*<I>S</I>, . . . , <I>L</I>+(<I>N</I>-1)*<I>S</I>}, where <I>U</I>
 * = <I>L</I>+(<I>N</I>-1)*<I>S</I>.
 * <P>
 * You construct a range object by specifying the lower bound, upper bound, and
 * stride. If the stride is omitted, the default stride is 1. The length is
 * determined automatically. If the lower bound is greater than the upper bound,
 * the range's length is 0 (an empty range).
 * <P>
 * You can use a range object to control a for loop like this:
 * <PRE>
 *     Range range = new Range (0, N-1);
 *     int lb = range.lb();
 *     int ub = range.ub();
 *     for (int i = lb; i &lt;= ub; ++ i)
 *         . . .
 * </PRE> Note that the range is from <code>lb()</code> to <code>ub()</code> inclusive,
 * so the appropriate test in the for loop is <code>i &lt;= ub</code>. Also note
 * that it usually reduces the running time to call <code>ub()</code> once, store
 * the result in a local variable, and use the local variable in the for loop
 * test, than to call <code>ub()</code> directly in the for loop test.
 * <P>
 * You can use a range object with a stride greater than 1 to control a for loop
 * like this:
 * <PRE>
 *     Range range = new Range (0, N-1, 2);
 *     int lb = range.lb();
 *     int ub = range.ub();
 *     int stride = range.stride();
 *     for (int i = lb; i &lt;= ub; i += stride)
 *         . . .
 * </PRE>
 *
 * @author Alan Kaminsky
 * @version 30-May-2007
 */
public class Range
        implements Externalizable {

// Hidden data members.
    private static final long serialVersionUID = -155434844308312282L;

    int lb;
    int stride;
    int length;
    int ub;

// Exported constructors.
    /**
     * Construct a new range object representing an empty range.
     */
    public Range() {
        this.lb = 0;
        this.stride = 1;
        this.length = 0;
        setUb();
    }

    /**
     * Construct a new range object with the given lower bound and upper bound.
     * The stride is 1. The range object represents the following set of
     * integers: {<I>L</I>, <I>L</I>+1, <I>L</I>+2, . . . , <I>U</I>}. The
     * range's length <I>N</I> is <I>U</I>-<I>L</I>+1.
     * <P>
     * <I>Note:</I> <I>L</I> &gt; <I>U</I> is allowed and stands for an empty
     * range.
     *
     * @param lb Lower bound <I>L</I>.
     * @param ub Upper bound <I>U</I>.
     */
    public Range(int lb,
            int ub) {
        this.lb = lb;
        this.stride = 1;
        this.length = Math.max(ub - lb + 1, 0);
        setUb();
    }

    /**
     * Construct a new range object with the given lower bound, upper bound, and
     * stride. The stride must be greater than or equal to 1. The range object
     * represents the following set of integers: {<I>L</I>, <I>L</I>+<I>S</I>,
     * <I>L</I>+2*<I>S</I>, . . . , <I>L</I>+(<I>N</I>-1)*<I>S</I>}, where the
     * range's length <I>N</I> is such that <I>L</I>+(<I>N</I>-1)*<I>S</I> is
     * the largest integer less than or equal to <I>U</I>. Note that the actual
     * upper bound of the range, <I>L</I>+(<I>N</I>-1)*<I>S</I>, may not be the
     * same as <I>U</I>.
     * <P>
     * <I>Note:</I> <I>L</I> &gt; <I>U</I> is allowed and stands for an empty
     * range.
     *
     * @param lb Lower bound <I>L</I>.
     * @param ub Upper bound <I>U</I>.
     * @param stride Stride <I>S</I> &gt;= 1.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>S</I> &lt; 1.
     */
    public Range(int lb,
            int ub,
            int stride) {
        if (stride < 1) {
            throw new IllegalArgumentException("Range(): stride = " + stride + " illegal");
        }
        this.lb = lb;
        this.stride = stride;
        this.length = Math.max((ub - lb + stride) / stride, 0);
        setUb();
    }

    /**
     * Construct a new range object that is a copy of the given range object.
     *
     * @param range Range object to copy.
     */
    public Range(Range range) {
        this.lb = range.lb;
        this.stride = range.stride;
        this.length = range.length;
        this.ub = range.ub;
    }

// Exported operations.
    /**
     * Returns this range's lower bound.
     *
     * @return Lower bound.
     */
    public int lb() {
        return lb;
    }

    /**
     * Returns this range's upper bound.
     *
     * @return Upper bound.
     */
    public int ub() {
        return ub;
    }

    /**
     * Returns this range's stride.
     *
     * @return Stride.
     */
    public int stride() {
        return stride;
    }

    /**
     * Returns this range's length.
     *
     * @return Length.
     */
    public int length() {
        return length;
    }

    /**
     * Determine if this range contains the given value. This range contains the
     * given value if <code>this.lb()</code> &lt;= <code>val</code> &lt;=
     * <code>this.ub()</code>. (The stride does not affect the outcome.)
     *
     * @param value Value to test.
     * @return True if this range contains the given <code>value</code>, false
     * otherwise.
     */
    public boolean contains(int value) {
        return this.lb <= value && value <= this.ub;
    }

    /**
     * Determine if this range contains the given range. This range contains the
     * given range if <code>this.lb()</code> &lt;= <code>range.lb()</code> and
     * <code>range.ub()</code> &lt;= <code>this.ub()</code>. (The strides do not affect
     * the outcome.)
     *
     * @param range Range to test.
     * @return True if this range contains the given <code>range</code>, false
     * otherwise.
     */
    public boolean contains(Range range) {
        return this.lb <= range.lb && range.ub <= this.ub;
    }

    /**
     * Partition this range and return one subrange. This range is split up into
     * subranges; the <code>size</code> argument specifies the number of subranges.
     * This range is divided as equally as possible among the subranges; the
     * lengths of the subranges differ by at most 1. The subranges are numbered
     * 0, 1, . . . <code>size-1</code>. This method returns the subrange whose
     * number is <code>rank</code>.
     * <P>
     * Note that if <code>size</code> is greater than the length of this range, the
     * returned subrange may be empty.
     *
     * @param size Number of subranges, <code>size</code> &gt;= 1.
     * @param rank Rank of the desired subrange, 0 &lt;= <code>rank</code> &lt;
     * <code>size</code>.
     * @return Subrange.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>size</code> or <code>rank</code> is out of bounds.
     */
    public Range subrange(int size,
            int rank) {
        // Verify preconditions.
        if (size < 1) {
            throw new IllegalArgumentException("Range.subrange(): size = " + size + " illegal");
        }
        if (0 > rank || rank >= size) {
            throw new IllegalArgumentException("Range.subrange(): rank = " + rank + " illegal");
        }

        // Split this range.
        Range result = new Range();
        int sublen = this.length / size;
        int subrem = this.length % size;
        if (rank < subrem) {
            ++sublen;
            result.lb = this.lb + (rank * sublen) * this.stride;
        } else {
            result.lb = this.lb + (subrem + rank * sublen) * this.stride;
        }
        result.stride = this.stride;
        result.length = sublen;
        result.setUb();
        return result;
    }

    /**
     * Partition this range and return all the subranges. This range is split up
     * into subranges; the <code>size</code> argument specifies the number of
     * subranges. This range is divided as equally as possible among the
     * subranges; the lengths of the subranges differ by at most 1. The
     * subranges are returned in an array with indexes 0, 1, . . .
     * <code>size-1</code>.
     * <P>
     * Note that if <code>size</code> is greater than the length of this range, some
     * of the returned subranges may be empty.
     *
     * @param size Number of subranges, size &gt;= 1.
     * @return Array of subranges.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>size</code> is out of bounds.
     */
    public Range[] subranges(int size) {
        // Verify preconditions.
        if (size < 1) {
            throw new IllegalArgumentException("Range.subranges(): size = " + size + " illegal");
        }

        // Allocate storage for subranges.
        Range[] result = new Range[size];

        // Compute subranges.
        int sublen = this.length / size;
        int subrem = this.length % size;
        int x = this.lb;
        ++sublen;
        for (int i = 0; i < subrem; ++i) {
            Range result_i = new Range();
            result_i.lb = x;
            x += sublen * this.stride;
            result_i.stride = this.stride;
            result_i.length = sublen;
            result_i.setUb();
            result[i] = result_i;
        }
        --sublen;
        for (int i = subrem; i < size; ++i) {
            Range result_i = new Range();
            result_i.lb = x;
            x += sublen * this.stride;
            result_i.stride = this.stride;
            result_i.length = sublen;
            result_i.setUb();
            result[i] = result_i;
        }

        return result;
    }

    /**
     * Slice off a chunk of this range and return the chunk. Considering this
     * range as a set of integers from the lower bound to the upper bound, the
     * first <code>N1</code> integers are sliced off and discarded, then the next
     * <code>N2</code> integers are sliced off to form a chunk, and the chunk is
     * returned. If after removing the first <code>N1</code> integers there are
     * fewer than <code>N2</code> integers left, a chunk consisting of all the
     * remaining integers is returned; this may be an empty chunk.
     *
     * @param N1 Number of integers to discard (must be &gt;= 0).
     * @param N2 Number of integers to include in the chunk (must be &gt;= 0).
     * @return Chunk.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>N1</code> or <code>N2</code> is out of bounds.
     */
    public Range chunk(int N1,
            int N2) {
        // Verify preconditions.
        if (N1 < 0) {
            throw new IllegalArgumentException("Range.chunk(): N1 = " + N1 + " illegal");
        }
        if (N2 < 0) {
            throw new IllegalArgumentException("Range.chunk(): N2 = " + N2 + " illegal");
        }

        Range result = new Range();
        result.lb = this.lb + N1 * this.stride;
        result.stride = this.stride;
        result.length = Math.min(N2, Math.max(0, this.length - N1));
        result.setUb();
        return result;
    }

    /**
     * {@inheritDoc}
     *
     * Determine if this range is equal to the given object. Two ranges are
     * equal if they both have the same lower bound, stride, and length.
     */
    public boolean equals(Object obj) {
        return obj instanceof Range
                && this.lb == ((Range) obj).lb
                && this.stride == ((Range) obj).stride
                && this.length == ((Range) obj).length;
    }

    /**
     * Returns a hash code for this range.
     *
     * @return a int.
     */
    public int hashCode() {
        return (((this.lb << 10) + this.stride) << 10) + this.length;
    }

    /**
     * Returns a string version of this range. If the stride is 1, the format is
     * <code>"<I>L</I>..<I>U</I>"</code>, where <I>L</I> is the lower bound and
     * <I>U</I> is the upper bound. If the stride is greater than 1, the format
     * is <code>"<I>L</I>..<I>U</I>;<I>S</I>"</code>, where <I>L</I> is the lower
     * bound, <I>U</I> is the upper bound, and <I>S</I> is the stride.
     *
     * @return a {@link java.lang.String} object.
     */
    public String toString() {
        StringBuilder b = new StringBuilder();
        b.append(lb);
        b.append("..");
        b.append(ub);
        if (stride > 1) {
            b.append(';');
            b.append(stride);
        }
        return b.toString();
    }

    /**
     * {@inheritDoc}
     *
     * Write this range to the given object output stream.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void writeExternal(ObjectOutput out)
            throws IOException {
        out.writeInt(lb);
        out.writeInt(stride);
        out.writeInt(length);
    }

    /**
     * {@inheritDoc}
     *
     * Read this range from the given object input stream.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void readExternal(ObjectInput in)
            throws IOException {
        lb = in.readInt();
        stride = in.readInt();
        length = in.readInt();
        setUb();
    }

// Hidden operations.
    /**
     * Set the upper bound of this range based on the lower bound, stride, and
     * length.
     */
    private void setUb() {
        ub = lb + (length - 1) * stride;
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length != 4)
//			{
//			System.err.println ("Usage: edu.rit.util.Range <lb> <ub> <stride> <size>");
//			System.exit (1);
//			}
//		int lb = Integer.parseInt (args[0]);
//		int ub = Integer.parseInt (args[1]);
//		int stride = Integer.parseInt (args[2]);
//		int size = Integer.parseInt (args[3]);
//		Range range = new Range (lb, ub, stride);
//		System.out.println
//			("Range = " + range + ", length = " + range.length());
//		for (int rank = 0; rank < size; ++ rank)
//			{
//			Range subrange = range.subrange (size, rank);
//			System.out.println
//				("Subrange rank " + rank + " = " + subrange +
//				 ", length = " + subrange.length());
//			}
//		Range[] subranges = range.subranges (size);
//		for (int rank = 0; rank < size; ++ rank)
//			{
//			System.out.println
//				("Subranges[" + rank + "] = " + subranges[rank] +
//				 ", length = " + subranges[rank].length());
//			}
//		}
}
