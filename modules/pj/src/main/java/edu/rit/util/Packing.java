//******************************************************************************
//
// File:    Packing.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Packing
//
// This Java source file is copyright (C) 2012 by Alan Kaminsky. All rights
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

/**
 * Class Packing provides static methods for packing and unpacking arrays of
 * bytes into and out of integers, long integers, and arrays thereof.
 * <P>
 * <I>Note:</I> The operations in class Packing are not multiple thread safe.
 *
 * @author Alan Kaminsky
 * @version 26-Mar-2012
 */
public class Packing {

// Prevent construction.
    private Packing() {
    }

// Exported operations.
    /**
     * Pack bytes from the given array into an integer in little-endian order.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @return Elements <code>src[srcPos]</code> through <code>src[srcPos+3]</code>
     * packed into an integer in little-endian order.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds.
     */
    public static int packIntLittleEndian(byte[] src,
            int srcPos) {
        if (srcPos + 4 > src.length) {
            throw new IndexOutOfBoundsException();
        }
        int rv = 0;
        for (int i = 0; i <= 3; ++i) {
            rv |= (src[srcPos + i] & 0xFF) << (i * 8);
        }
        return rv;
    }

    /**
     * Pack bytes from the given array into an integer in big-endian order.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @return Elements <code>src[srcPos]</code> through <code>src[srcPos+3]</code>
     * packed into an integer in big-endian order.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds.
     */
    public static int packIntBigEndian(byte[] src,
            int srcPos) {
        if (srcPos + 4 > src.length) {
            throw new IndexOutOfBoundsException();
        }
        int rv = 0;
        for (int i = 0; i <= 3; ++i) {
            rv |= (src[srcPos + i] & 0xFF) << ((3 - i) * 8);
        }
        return rv;
    }

    /**
     * Pack bytes from the given array into the given array of integers in
     * little-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+4*len-1]</code> are packed into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+len-1]</code>.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @param dst Destination array of packed integers.
     * @param dstPos Index of first packed integer.
     * @param len Number of integers (not bytes!) to pack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds; in this case
     * <code>dst</code> is not altered.
     */
    public static void packIntLittleEndian(byte[] src,
            int srcPos,
            int[] dst,
            int dstPos,
            int len) {
        if (srcPos + 4 * len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            dst[dstPos + i] = packIntLittleEndian(src, srcPos + 4 * i);
        }
    }

    /**
     * Pack bytes from the given array into the given array of integers in
     * big-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+4*len-1]</code> are packed into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+len-1]</code>.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @param dst Destination array of packed integers.
     * @param dstPos Index of first packed integer.
     * @param len Number of integers (not bytes!) to pack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds; in this case
     * <code>dst</code> is not altered.
     */
    public static void packIntBigEndian(byte[] src,
            int srcPos,
            int[] dst,
            int dstPos,
            int len) {
        if (srcPos + 4 * len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            dst[dstPos + i] = packIntBigEndian(src, srcPos + 4 * i);
        }
    }

    /**
     * Pack bytes from the given array into a long integer in little-endian
     * order.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @return Elements <code>src[srcPos]</code> through <code>src[srcPos+7]</code>
     * packed into a long integer in little-endian order.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds.
     */
    public static long packLongLittleEndian(byte[] src,
            int srcPos) {
        if (srcPos + 8 > src.length) {
            throw new IndexOutOfBoundsException();
        }
        long rv = 0L;
        for (int i = 0; i <= 7; ++i) {
            rv |= (src[srcPos + i] & 0xFFL) << (i * 8);
        }
        return rv;
    }

    /**
     * Pack bytes from the given array into a long integer in big-endian order.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @return Elements <code>src[srcPos]</code> through <code>src[srcPos+7]</code>
     * packed into a long integer in big-endian order.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds.
     */
    public static long packLongBigEndian(byte[] src,
            int srcPos) {
        if (srcPos + 8 > src.length) {
            throw new IndexOutOfBoundsException();
        }
        long rv = 0L;
        for (int i = 0; i <= 7; ++i) {
            rv |= (src[srcPos + i] & 0xFFL) << ((7 - i) * 8);
        }
        return rv;
    }

    /**
     * Pack bytes from the given array into the given array of long integers in
     * little-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+8*len-1]</code> are packed into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+len-1]</code>.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @param dst Destination array of packed long integers.
     * @param dstPos Index of first packed long integer.
     * @param len Number of long integers (not bytes!) to pack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds; in this case
     * <code>dst</code> is not altered.
     */
    public static void packLongLittleEndian(byte[] src,
            int srcPos,
            long[] dst,
            int dstPos,
            int len) {
        if (srcPos + 8 * len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            dst[dstPos + i] = packLongLittleEndian(src, srcPos + 8 * i);
        }
    }

    /**
     * Pack bytes from the given array into the given array of long integers in
     * big-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+8*len-1]</code> are packed into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+len-1]</code>.
     *
     * @param src Source array of bytes to pack.
     * @param srcPos Index of first byte to pack.
     * @param dst Destination array of packed long integers.
     * @param dstPos Index of first packed long integer.
     * @param len Number of long integers (not bytes!) to pack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * packing would cause accessing array elements out of bounds; in this case
     * <code>dst</code> is not altered.
     */
    public static void packLongBigEndian(byte[] src,
            int srcPos,
            long[] dst,
            int dstPos,
            int len) {
        if (srcPos + 8 * len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            dst[dstPos + i] = packLongBigEndian(src, srcPos + 8 * i);
        }
    }

    /**
     * Unpack the given integer into the given array of bytes in little-endian
     * order. The integer is unpacked into elements <code>dst[dstPos]</code> through
     * <code>dst[dstPos+3]</code>.
     *
     * @param src Source integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackIntLittleEndian(int src,
            byte[] dst,
            int dstPos) {
        if (dstPos + 4 > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i <= 3; ++i) {
            dst[dstPos + i] = (byte) (src >> (i * 8));
        }
    }

    /**
     * Unpack the given integer into the given array of bytes in big-endian
     * order. The integer is unpacked into elements <code>dst[dstPos]</code> through
     * <code>dst[dstPos+3]</code>.
     *
     * @param src Source integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackIntBigEndian(int src,
            byte[] dst,
            int dstPos) {
        if (dstPos + 4 > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i <= 3; ++i) {
            dst[dstPos + i] = (byte) (src >> ((3 - i) * 8));
        }
    }

    /**
     * Unpack integers from the given array into the given array of bytes in
     * little-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+len-1]</code> are unpacked into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+4*len-1]</code>.
     *
     * @param src Source array of integers to unpack.
     * @param srcPos Index of first integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @param len Number of integers (not bytes!) to unpack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackIntLittleEndian(int[] src,
            int srcPos,
            byte[] dst,
            int dstPos,
            int len) {
        if (srcPos + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + 4 * len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            unpackIntLittleEndian(src[srcPos + i], dst, dstPos + 4 * i);
        }
    }

    /**
     * Unpack integers from the given array into the given array of bytes in
     * big-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+len-1]</code> are unpacked into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+4*len-1]</code>.
     *
     * @param src Source array of integers to unpack.
     * @param srcPos Index of first integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @param len Number of integers (not bytes!) to unpack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackIntBigEndian(int[] src,
            int srcPos,
            byte[] dst,
            int dstPos,
            int len) {
        if (srcPos + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + 4 * len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            unpackIntBigEndian(src[srcPos + i], dst, dstPos + 4 * i);
        }
    }

    /**
     * Unpack the given long integer into the given array of bytes in
     * little-endian order. The long integer is unpacked into elements
     * <code>dst[dstPos]</code> through <code>dst[dstPos+7]</code>.
     *
     * @param src Source long integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackLongLittleEndian(long src,
            byte[] dst,
            int dstPos) {
        if (dstPos + 8 > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i <= 7; ++i) {
            dst[dstPos + i] = (byte) (src >> (i * 8));
        }
    }

    /**
     * Unpack the given long integer into the given array of bytes in big-endian
     * order. The long integer is unpacked into elements <code>dst[dstPos]</code>
     * through <code>dst[dstPos+3]</code>.
     *
     * @param src Source long integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackLongBigEndian(long src,
            byte[] dst,
            int dstPos) {
        if (dstPos + 8 > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i <= 7; ++i) {
            dst[dstPos + i] = (byte) (src >> ((7 - i) * 8));
        }
    }

    /**
     * Unpack long integers from the given array into the given array of bytes
     * in little-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+len-1]</code> are unpacked into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+8*len-1]</code>.
     *
     * @param src Source array of long integers to unpack.
     * @param srcPos Index of first long integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @param len Number of integers (not bytes!) to unpack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackLongLittleEndian(long[] src,
            int srcPos,
            byte[] dst,
            int dstPos,
            int len) {
        if (srcPos + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + 8 * len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            unpackLongLittleEndian(src[srcPos + i], dst, dstPos + 8 * i);
        }
    }

    /**
     * Unpack long integers from the given array into the given array of bytes
     * in big-endian order. Elements <code>src[srcPos]</code> through
     * <code>src[srcPos+len-1]</code> are unpacked into <code>dst[dstPos]</code> through
     * <code>dst[dstPos+8*len-1]</code>.
     *
     * @param src Source array of long integers to unpack.
     * @param srcPos Index of first long integer to unpack.
     * @param dst Destination array to receive unpacked bytes.
     * @param dstPos Index of first unpacked byte.
     * @param len Number of integers (not bytes!) to unpack.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>src</code> is null. Thrown if
     * <code>dst</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * unpacking would cause accessing array elements out of bounds; in this
     * case <code>dst</code> is not altered.
     */
    public static void unpackLongBigEndian(long[] src,
            int srcPos,
            byte[] dst,
            int dstPos,
            int len) {
        if (srcPos + len > src.length) {
            throw new IndexOutOfBoundsException();
        }
        if (dstPos + 8 * len > dst.length) {
            throw new IndexOutOfBoundsException();
        }
        for (int i = 0; i < len; ++i) {
            unpackLongBigEndian(src[srcPos + i], dst, dstPos + 8 * i);
        }
    }

}
