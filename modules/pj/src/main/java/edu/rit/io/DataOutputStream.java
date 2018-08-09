//******************************************************************************
//
// File:    DataOutputStream.java
// Package: edu.rit.io
// Unit:    Class edu.rit.io.DataOutputStream
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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
package edu.rit.io;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Class DataOutputStream provides an output stream that writes primitive data
 * types and strings in binary form. It behaves similarly to class
 * java.io.DataOutputStream, except the methods for writing types byte, short,
 * char, int, long, and String are implemented differently. These methods write
 * an integer value using a variable number of bytes, as described below. This
 * can save space in the file if small integer values are written more
 * frequently than large integer values. The resulting byte stream can be read
 * using class {@linkplain DataInputStream}.
 * <P>
 * Note that class DataOutputStream does <I>not</I> implement interface
 * java.io.DataOutput, because the methods do not obey the contract specified in
 * that interface.
 *
 * @author Alan Kaminsky
 * @version 18-Dec-2009
 */
public class DataOutputStream
        extends FilterOutputStream {

// Exported constructors.
    /**
     * Construct a new data output stream.
     *
     * @param out Underlying output stream.
     */
    public DataOutputStream(OutputStream out) {
        super(out);
    }

// Exported operations.
    /**
     * Write the given Boolean value to this data output stream. One byte is
     * written, either 0 (if <TT>v</TT> is false) or 1 (if <TT>v</TT> is true).
     *
     * @param v Boolean value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeBoolean(boolean v)
            throws IOException {
        out.write(v ? 1 : 0);
    }

    /**
     * Write the given integer value to this data output stream. This method can
     * be used to write values of type byte, short, char, or int. From one to
     * five bytes are written, in big-endian order, as follows:
     * <UL>
     * <P>
     * <LI>
     * If &minus;64 &le; <TT>v</TT> &le; 63, then one byte is written,
     * containing 0 (1 bit) followed by <TT>v</TT> (7 bits).
     * <P>
     * <LI>
     * Else if &minus;8192 &le; <TT>v</TT> &le; 8191, then two bytes are
     * written, containing 10 (2 bits) followed by <TT>v</TT> (14 bits).
     * <P>
     * <LI>
     * Else if &minus;1048576 &le; <TT>v</TT> &le; 1048575, then three bytes are
     * written, containing 110 (3 bits) followed by <TT>v</TT> (21 bits).
     * <P>
     * <LI>
     * Else if &minus;134217728 &le; <TT>v</TT> &le; 134217727, then four bytes
     * are written, containing 1110 (4 bits) followed by <TT>v</TT> (28 bits).
     * <P>
     * <LI>
     * Else five bytes are written, containing 1111 (4 bits) followed by
     * <TT>v</TT> (sign-extended to 36 bits).
     * </UL>
     *
     * @param v Integer value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeInt(int v)
            throws IOException {
        if (-64 <= v && v <= 63) {
            out.write(v & 0x7F);
        } else if (-8192 <= v && v <= 8191) {
            out.write(((v >> 8) & 0x3F) | 0x80);
            out.write(v);
        } else if (-1048576 <= v && v <= 1048575) {
            out.write(((v >> 16) & 0x1F) | 0xC0);
            out.write(v >> 8);
            out.write(v);
        } else if (-134217728 <= v && v <= 134217727) {
            out.write(((v >> 24) & 0x0F) | 0xE0);
            out.write(v >> 16);
            out.write(v >> 8);
            out.write(v);
        } else {
            out.write(((v >> 32) & 0x0F) | 0xF0);
            out.write(v >> 24);
            out.write(v >> 16);
            out.write(v >> 8);
            out.write(v);
        }
    }

    /**
     * Write the given unsigned integer value to this data output stream. This
     * method can be used to write values of type byte, short, char, or int.
     * From one to five bytes are written, in big-endian order, as follows:
     * <UL>
     * <P>
     * <LI>
     * If 0 &le; <TT>v</TT> &le; 127, then one byte is written, containing 0 (1
     * bit) followed by <TT>v</TT> (7 bits).
     * <P>
     * <LI>
     * Else if 128 &le; <TT>v</TT> &le; 16383, then two bytes are written,
     * containing 10 (2 bits) followed by <TT>v</TT> (14 bits).
     * <P>
     * <LI>
     * Else if 16384 &le; <TT>v</TT> &le; 2097151, then three bytes are written,
     * containing 110 (3 bits) followed by <TT>v</TT> (21 bits).
     * <P>
     * <LI>
     * Else if 2097152 &le; <TT>v</TT> &le; 268435455, then four bytes are
     * written, containing 1110 (4 bits) followed by <TT>v</TT> (28 bits).
     * <P>
     * <LI>
     * Else five bytes are written, containing 1111 (4 bits) followed by
     * <TT>v</TT> (zero-extended to 36 bits).
     * </UL>
     *
     * @param v Integer value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeUnsignedInt(int v)
            throws IOException {
        if (0 <= v && v <= 127) {
            out.write(v);
        } else if (128 <= v && v <= 16383) {
            out.write((v >> 8) | 0x80);
            out.write(v);
        } else if (16384 <= v && v <= 2097151) {
            out.write((v >> 16) | 0xC0);
            out.write(v >> 8);
            out.write(v);
        } else if (2097152 <= v && v <= 268435455) {
            out.write((v >> 24) | 0xE0);
            out.write(v >> 16);
            out.write(v >> 8);
            out.write(v);
        } else {
            out.write(0xF0);
            out.write(v >> 24);
            out.write(v >> 16);
            out.write(v >> 8);
            out.write(v);
        }
    }

    /**
     * Write the given long value to this data output stream. From one to nine
     * bytes are written, in big-endian order, as follows:
     * <UL>
     * <P>
     * <LI>
     * If &minus;64 &le; <TT>v</TT> &le; 63, then one byte is written,
     * containing 0 (1 bit) followed by <TT>v</TT> (7 bits).
     * <P>
     * <LI>
     * Else if &minus;8192 &le; <TT>v</TT> &le; 8191, then two bytes are
     * written, containing 10 (2 bits) followed by <TT>v</TT> (14 bits).
     * <P>
     * <LI>
     * Else if &minus;1048576 &le; <TT>v</TT> &le; 1048575, then three bytes are
     * written, containing 110 (3 bits) followed by <TT>v</TT> (21 bits).
     * <P>
     * <LI>
     * Else if &minus;134217728 &le; <TT>v</TT> &le; 134217727, then four bytes
     * are written, containing 1110 (4 bits) followed by <TT>v</TT> (28 bits).
     * <P>
     * <LI>
     * Else if &minus;17179869184 &le; <TT>v</TT> &le; 17179869183, then five
     * bytes are written, containing 11110 (5 bits) followed by <TT>v</TT> (35
     * bits).
     * <P>
     * <LI>
     * Else if &minus;2199023255552 &le; <TT>v</TT> &le; 2199023255551, then six
     * bytes are written, containing 111110 (6 bits) followed by <TT>v</TT> (42
     * bits).
     * <P>
     * <LI>
     * Else if &minus;281474976710656 &le; <TT>v</TT> &le; 281474976710655, then
     * seven bytes are written, containing 1111110 (7 bits) followed by
     * <TT>v</TT> (49 bits).
     * <P>
     * <LI>
     * Else if &minus;36028797018963968 &le; <TT>v</TT> &le; 36028797018963967,
     * then eight bytes are written, containing 11111110 (8 bits) followed by
     * <TT>v</TT> (56 bits).
     * <P>
     * <LI>
     * Else nine bytes are written, containing 11111111 (8 bits) followed by
     * <TT>v</TT> (64 bits).
     * </UL>
     *
     * @param v Integer value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeLong(long v)
            throws IOException {
        if (-64L <= v && v <= 63L) {
            out.write((int) (v) & 0x7F);
        } else if (-8192L <= v && v <= 8191L) {
            out.write(((int) (v >> 8L) & 0x3F) | 0x80);
            out.write((int) (v));
        } else if (-1048576L <= v && v <= 1048575L) {
            out.write(((int) (v >> 16L) & 0x1F) | 0xC0);
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (-134217728L <= v && v <= 134217727L) {
            out.write(((int) (v >> 24L) & 0x0F) | 0xE0);
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (-17179869184L <= v && v <= 17179869183L) {
            out.write(((int) (v >> 32L) & 0x07) | 0xF0);
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (-2199023255552L <= v && v <= 2199023255551L) {
            out.write(((int) (v >> 40L) & 0x03) | 0xF8);
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (-281474976710656L <= v && v <= 281474976710656L) {
            out.write(((int) (v >> 48L) & 0x01) | 0xFC);
            out.write((int) (v >> 40L));
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (-36028797018963968L <= v && v <= 36028797018963968L) {
            out.write(0xFE);
            out.write((int) (v >> 48L));
            out.write((int) (v >> 40L));
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else {
            out.write(0xFF);
            out.write((int) (v >> 56L));
            out.write((int) (v >> 48L));
            out.write((int) (v >> 40L));
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        }
    }

    /**
     * Write the given unsigned long value to this data output stream. From one
     * to nine bytes are written, in big-endian order, as follows:
     * <UL>
     * <P>
     * <LI>
     * If 0 &le; <TT>v</TT> &le; 127, then one byte is written, containing 0 (1
     * bit) followed by <TT>v</TT> (7 bits).
     * <P>
     * <LI>
     * Else if 128 &le; <TT>v</TT> &le; 16383, then two bytes are written,
     * containing 10 (2 bits) followed by <TT>v</TT> (14 bits).
     * <P>
     * <LI>
     * Else if 16384 &le; <TT>v</TT> &le; 2097151, then three bytes are written,
     * containing 110 (3 bits) followed by <TT>v</TT> (21 bits).
     * <P>
     * <LI>
     * Else if 2097152 &le; <TT>v</TT> &le; 268435455, then four bytes are
     * written, containing 1110 (4 bits) followed by <TT>v</TT> (28 bits).
     * <P>
     * <LI>
     * Else if 268435456 &le; <TT>v</TT> &le; 34359738367, then five bytes are
     * written, containing 11110 (5 bits) followed by <TT>v</TT> (35 bits).
     * <P>
     * <LI>
     * Else if 34359738368 &le; <TT>v</TT> &le; 4398046511103, then six bytes
     * are written, containing 111110 (6 bits) followed by <TT>v</TT> (42 bits).
     * <P>
     * <LI>
     * Else if 4398046511104 &le; <TT>v</TT> &le; 562949953421311, then seven
     * bytes are written, containing 1111110 (7 bits) followed by <TT>v</TT> (49
     * bits).
     * <P>
     * <LI>
     * Else if 562949953421312 &le; <TT>v</TT> &le; 72057594037927935, then
     * eight bytes are written, containing 11111110 (8 bits) followed by
     * <TT>v</TT> (56 bits).
     * <P>
     * <LI>
     * Else nine bytes are written, containing 11111111 (8 bits) followed by
     * <TT>v</TT> (64 bits).
     * </UL>
     *
     * @param v Integer value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeUnsignedLong(long v)
            throws IOException {
        if (0L <= v && v <= 127L) {
            out.write((int) (v));
        } else if (128L <= v && v <= 16383L) {
            out.write((int) (v >> 8L) | 0x80);
            out.write((int) (v));
        } else if (16384L <= v && v <= 2097151L) {
            out.write((int) (v >> 16L) | 0xC0);
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (2097152L <= v && v <= 268435455L) {
            out.write((int) (v >> 24L) | 0xE0);
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (268435456L <= v && v <= 34359738367L) {
            out.write((int) (v >> 32L) | 0xF0);
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (34359738368L <= v && v <= 4398046511103L) {
            out.write((int) (v >> 40L) | 0xF8);
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (4398046511104L <= v && v <= 562949953421311L) {
            out.write((int) (v >> 48L) | 0xFC);
            out.write((int) (v >> 40L));
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else if (562949953421312L <= v && v <= 72057594037927935L) {
            out.write(0xFE);
            out.write((int) (v >> 48L));
            out.write((int) (v >> 40L));
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        } else {
            out.write(0xFF);
            out.write((int) (v >> 56L));
            out.write((int) (v >> 48L));
            out.write((int) (v >> 40L));
            out.write((int) (v >> 32L));
            out.write((int) (v >> 24L));
            out.write((int) (v >> 16L));
            out.write((int) (v >> 8L));
            out.write((int) (v));
        }
    }

    /**
     * Write the given float value to this data output stream. Four bytes are
     * written in big-endian order containing
     * <TT>Float.floatToRawIntBits(v)</TT>.
     *
     * @param v Float value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeFloat(float v)
            throws IOException {
        int x = Float.floatToRawIntBits(v);
        out.write(x >> 24);
        out.write(x >> 16);
        out.write(x >> 8);
        out.write(x);
    }

    /**
     * Write the given double value to this data output stream. Eight bytes are
     * written in big-endian order containing
     * <TT>Double.doubleToRawLongBits(v)</TT>.
     *
     * @param v Double value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeDouble(double v)
            throws IOException {
        long x = Double.doubleToRawLongBits(v);
        out.write((int) (x >> 56L));
        out.write((int) (x >> 48L));
        out.write((int) (x >> 40L));
        out.write((int) (x >> 32L));
        out.write((int) (x >> 24L));
        out.write((int) (x >> 16L));
        out.write((int) (x >> 8L));
        out.write((int) (x));
    }

    /**
     * Write the given string value to this data output stream. The length of
     * the string is written using <TT>writeUnsignedInt()</TT>, then each
     * character of the string is written using <TT>writeUnsignedInt()</TT>.
     *
     * @param v String value.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void writeString(String v)
            throws IOException {
        int n = v.length();
        writeUnsignedInt(n);
        for (int i = 0; i < n; ++i) {
            writeUnsignedInt(v.charAt(i));
        }
    }

}
