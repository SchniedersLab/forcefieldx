//******************************************************************************
//
// File:    ReduceArrays.java
// Package: edu.rit.pj.reduction
// Unit:    Class edu.rit.pj.reduction.ReduceArrays
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
package edu.rit.pj.reduction;

import edu.rit.util.Range;

/**
 * Class ReduceArrays provides static methods for reduction operations on arrays
 * and matrices of primitive types and object types.
 * <P>
 * <I>Note:</I> The operations in class ReduceArrays are not multiple thread
 * safe.
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2007
 */
public class ReduceArrays {

// Prevent construction.
    private ReduceArrays() {
    }

// Exported operations.
    /**
     * Combine a range of elements from one Boolean array with a range of
     * elements in another Boolean array. The number of elements combined is the
     * smaller of <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length.
     * Either or both of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may
     * be greater than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(boolean[] src,
            Range srcRange,
            boolean[] dst,
            Range dstRange,
            BooleanOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one byte array with a range of elements
     * in another byte array. The number of elements combined is the smaller of
     * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
     * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
     * than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(byte[] src,
            Range srcRange,
            byte[] dst,
            Range dstRange,
            ByteOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one character array with a range of
     * elements in another character array. The number of elements combined is
     * the smaller of <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length.
     * Either or both of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may
     * be greater than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(char[] src,
            Range srcRange,
            char[] dst,
            Range dstRange,
            CharacterOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one double array with a range of
     * elements in another double array. The number of elements combined is the
     * smaller of <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length.
     * Either or both of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may
     * be greater than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(double[] src,
            Range srcRange,
            double[] dst,
            Range dstRange,
            DoubleOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one float array with a range of elements
     * in another float array. The number of elements combined is the smaller of
     * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
     * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
     * than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(float[] src,
            Range srcRange,
            float[] dst,
            Range dstRange,
            FloatOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one integer array with a range of
     * elements in another integer array. The number of elements combined is the
     * smaller of <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length.
     * Either or both of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may
     * be greater than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(int[] src,
            Range srcRange,
            int[] dst,
            Range dstRange,
            IntegerOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one long array with a range of elements
     * in another long array. The number of elements combined is the smaller of
     * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
     * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
     * than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(long[] src,
            Range srcRange,
            long[] dst,
            Range dstRange,
            LongOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one short array with a range of elements
     * in another short array. The number of elements combined is the smaller of
     * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
     * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
     * than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static void reduce(short[] src,
            Range srcRange,
            short[] dst,
            Range dstRange,
            ShortOp op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one object array with a range of
     * elements in another object array. The number of elements combined is the
     * smaller of <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length.
     * Either or both of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may
     * be greater than 1.
     * <P>
     * For each destination array element <I>D</I> in the destination range and
     * each corresponding source array element <I>S</I> in the source range,
     * <I>D</I> is set to <I>D op S</I>.
     *
     * @param <DT> Destination array element data type.
     * @param src Source array.
     * @param srcRange Range of source elements.
     * @param srcRange Range of source elements.
     * @param dst Destination array.
     * @param dstRange Range of destination elements.
     * @param dstRange Range of destination elements.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRange</TT> is outside the bounds of the source array.
     * Thrown if any index in
     * <TT>dstRange</TT> is outside the bounds of the destination array.
     */
    public static <DT, ST extends DT> void reduce(ST[] src,
            Range srcRange,
            DT[] dst,
            Range dstRange,
            ObjectOp<DT> op) {
        int len = Math.min(srcRange.length(), dstRange.length());
        if (len == 0) {
            return;
        }
        int srcLower = srcRange.lb();
        int dstLower = dstRange.lb();
        int srcStride = srcRange.stride();
        int dstStride = dstRange.stride();
        int srcUpper = srcLower + (len - 1) * srcStride;
        int dstUpper = dstLower + (len - 1) * dstStride;
        if (0 > srcLower || srcUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src indexes = 0.." + (src.length - 1)
                    + ", srcRange = " + srcRange);
        }
        if (0 > dstLower || dstUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst indexes = 0.." + (dst.length - 1)
                    + ", dstRange = " + dstRange);
        }
        if (src != dst || srcLower > dstLower) {
            for (int i = srcLower, j = dstLower;
                    i <= srcUpper;
                    i += srcStride, j += dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        } else if (srcLower < dstLower) {
            for (int i = srcUpper, j = dstUpper;
                    i >= srcLower;
                    i -= srcStride, j -= dstStride) {
                dst[j] = op.op(dst[j], src[i]);
            }
        }
    }

    /**
     * Combine a range of elements from one Boolean matrix with a range of
     * elements in another Boolean matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(boolean[][] src,
            Range srcRowRange,
            Range srcColRange,
            boolean[][] dst,
            Range dstRowRange,
            Range dstColRange,
            BooleanOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one byte matrix with a range of elements
     * in another byte matrix. The number of rows combined is the smaller of
     * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
     * each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(byte[][] src,
            Range srcRowRange,
            Range srcColRange,
            byte[][] dst,
            Range dstRowRange,
            Range dstColRange,
            ByteOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one character matrix with a range of
     * elements in another character matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(char[][] src,
            Range srcRowRange,
            Range srcColRange,
            char[][] dst,
            Range dstRowRange,
            Range dstColRange,
            CharacterOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one double matrix with a range of
     * elements in another double matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(double[][] src,
            Range srcRowRange,
            Range srcColRange,
            double[][] dst,
            Range dstRowRange,
            Range dstColRange,
            DoubleOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one float matrix with a range of
     * elements in another float matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(float[][] src,
            Range srcRowRange,
            Range srcColRange,
            float[][] dst,
            Range dstRowRange,
            Range dstColRange,
            FloatOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one integer matrix with a range of
     * elements in another integer matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(int[][] src,
            Range srcRowRange,
            Range srcColRange,
            int[][] dst,
            Range dstRowRange,
            Range dstColRange,
            IntegerOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one long matrix with a range of elements
     * in another long matrix. The number of rows combined is the smaller of
     * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
     * each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(long[][] src,
            Range srcRowRange,
            Range srcColRange,
            long[][] dst,
            Range dstRowRange,
            Range dstColRange,
            LongOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one short matrix with a range of
     * elements in another short matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static void reduce(short[][] src,
            Range srcRowRange,
            Range srcColRange,
            short[][] dst,
            Range dstRowRange,
            Range dstColRange,
            ShortOp op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

    /**
     * Combine a range of elements from one object matrix with a range of
     * elements in another object matrix. The number of rows combined is the
     * smaller of <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s
     * length. Within each row, the number of columns combined is the smaller of
     * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
     * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
     * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
     * that the source matrix is fully allocated; each row in the source matrix
     * is the same length; the destination matrix is fully allocated; and each
     * row in the destination matrix is the same length.
     * <P>
     * For each destination matrix element <I>D</I> in the destination row and
     * column ranges and each corresponding source matrix element <I>S</I> in
     * the source row and column ranges, <I>D</I> is set to <I>D op S</I>.
     *
     * @param <DT> Destination matrix element data type.
     * @param src Source matrix.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param srcRowRange Range of source rows.
     * @param srcColRange Range of source columns.
     * @param dst Destination matrix.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param dstRowRange Range of destination rows.
     * @param dstColRange Range of destination columns.
     * @param op Binary operation.
     * @exception NullPointerException (unchecked exception) Thrown if any
     * argument is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if any
     * index in <TT>srcRowRange</TT> is outside the row bounds of the source
     * matrix. Thrown if any index in
     * <TT>srcColRange</TT> is outside the column bounds of the source matrix.
     * Thrown if any index in <TT>dstRowRange</TT> is outside the row bounds of
     * the destination matrix. Thrown if any index in
     * <TT>dstColRange</TT> is outside the column bounds of the destination
     * matrix.
     */
    public static <DT, ST extends DT> void reduce(ST[][] src,
            Range srcRowRange,
            Range srcColRange,
            DT[][] dst,
            Range dstRowRange,
            Range dstColRange,
            ObjectOp<DT> op) {
        int len = Math.min(srcRowRange.length(), dstRowRange.length());
        if (len == 0) {
            return;
        }
        int srcRowLower = srcRowRange.lb();
        int dstRowLower = dstRowRange.lb();
        int srcRowStride = srcRowRange.stride();
        int dstRowStride = dstRowRange.stride();
        int srcRowUpper = srcRowLower + (len - 1) * srcRowStride;
        int dstRowUpper = dstRowLower + (len - 1) * dstRowStride;
        if (0 > srcRowLower || srcRowUpper >= src.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): src row indexes = 0.."
                    + (src.length - 1) + ", srcRowRange = " + srcRowRange);
        }
        if (0 > dstRowLower || dstRowUpper >= dst.length) {
            throw new IndexOutOfBoundsException("ReduceArrays.reduce(): dst row indexes = 0.."
                    + (dst.length - 1) + ", dstRowRange = " + dstRowRange);
        }
        if (src != dst || srcRowLower > dstRowLower) {
            for (int i = srcRowLower, j = dstRowLower;
                    i <= srcRowUpper;
                    i += srcRowStride, j += dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        } else if (srcRowLower < dstRowLower) {
            for (int i = srcRowUpper, j = dstRowUpper;
                    i >= srcRowLower;
                    i -= srcRowStride, j -= dstRowStride) {
                reduce(src[i], srcColRange, dst[j], dstColRange, op);
            }
        }
    }

}
