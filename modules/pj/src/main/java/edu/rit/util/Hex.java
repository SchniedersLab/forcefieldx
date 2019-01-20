//******************************************************************************
//
// File:    Hex.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Hex
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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
 * Class Hex provides static methods for converting between hexadecimal strings
 * and numbers of various kinds.
 *
 * @author Alan Kaminsky
 * @version 22-Apr-2006
 */
public class Hex {

// Hidden data members.
    private static final char[] int2hex = new char[]{'0', '1', '2', '3', '4', '5', '6', '7',
        '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

// Prevent construction.
    private Hex() {
    }

// Exported operations.
    /**
     * Convert the given <code>byte</code> value to a two-digit hexadecimal string.
     *
     * @param val Value.
     * @return Hexadecimal string.
     */
    public static String toString(byte val) {
        StringBuilder buf = new StringBuilder(2);
        buf.append(int2hex[(val >> 4) & 0xF]);
        buf.append(int2hex[(val) & 0xF]);
        return buf.toString();
    }

    /**
     * Convert the given <code>char</code> value to a four-digit hexadecimal string.
     *
     * @param val Value.
     * @return Hexadecimal string.
     */
    public static String toString(char val) {
        StringBuilder buf = new StringBuilder(4);
        buf.append(int2hex[(val >> 12) & 0xF]);
        buf.append(int2hex[(val >> 8) & 0xF]);
        buf.append(int2hex[(val >> 4) & 0xF]);
        buf.append(int2hex[(val) & 0xF]);
        return buf.toString();
    }

    /**
     * Convert the given <code>short</code> value to a four-digit hexadecimal
     * string.
     *
     * @param val Value.
     * @return Hexadecimal string.
     */
    public static String toString(short val) {
        StringBuilder buf = new StringBuilder(4);
        buf.append(int2hex[(val >> 12) & 0xF]);
        buf.append(int2hex[(val >> 8) & 0xF]);
        buf.append(int2hex[(val >> 4) & 0xF]);
        buf.append(int2hex[(val) & 0xF]);
        return buf.toString();
    }

    /**
     * Convert the given <code>int</code> value to an eight-digit hexadecimal
     * string.
     *
     * @param val Value.
     * @return Hexadecimal string.
     */
    public static String toString(int val) {
        StringBuilder buf = new StringBuilder(8);
        buf.append(int2hex[(val >> 28) & 0xF]);
        buf.append(int2hex[(val >> 24) & 0xF]);
        buf.append(int2hex[(val >> 20) & 0xF]);
        buf.append(int2hex[(val >> 16) & 0xF]);
        buf.append(int2hex[(val >> 12) & 0xF]);
        buf.append(int2hex[(val >> 8) & 0xF]);
        buf.append(int2hex[(val >> 4) & 0xF]);
        buf.append(int2hex[(val) & 0xF]);
        return buf.toString();
    }

    /**
     * Convert the given <code>long</code> value to a sixteen-digit hexadecimal
     * string.
     *
     * @param val Value.
     * @return Hexadecimal string.
     */
    public static String toString(long val) {
        StringBuilder buf = new StringBuilder(16);
        buf.append(int2hex[(int) (val >> 60) & 0xF]);
        buf.append(int2hex[(int) (val >> 56) & 0xF]);
        buf.append(int2hex[(int) (val >> 52) & 0xF]);
        buf.append(int2hex[(int) (val >> 48) & 0xF]);
        buf.append(int2hex[(int) (val >> 44) & 0xF]);
        buf.append(int2hex[(int) (val >> 40) & 0xF]);
        buf.append(int2hex[(int) (val >> 36) & 0xF]);
        buf.append(int2hex[(int) (val >> 32) & 0xF]);
        buf.append(int2hex[(int) (val >> 28) & 0xF]);
        buf.append(int2hex[(int) (val >> 24) & 0xF]);
        buf.append(int2hex[(int) (val >> 20) & 0xF]);
        buf.append(int2hex[(int) (val >> 16) & 0xF]);
        buf.append(int2hex[(int) (val >> 12) & 0xF]);
        buf.append(int2hex[(int) (val >> 8) & 0xF]);
        buf.append(int2hex[(int) (val >> 4) & 0xF]);
        buf.append(int2hex[(int) (val) & 0xF]);
        return buf.toString();
    }

    /**
     * Convert the given <code>byte</code> array to a hexadecimal string. Each byte
     * is converted to two hexadecimal digits.
     *
     * @param val Byte array.
     * @return Hexadecimal string.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>val</code> is null.
     */
    public static String toString(byte[] val) {
        return toString(val, 0, val.length);
    }

    /**
     * Convert a portion of the given <code>byte</code> array to a hexadecimal
     * string. Each byte is converted to two hexadecimal digits.
     *
     * @param val Byte array.
     * @param off Index of first byte to convert.
     * @param len Number of bytes to convert.
     * @return Hexadecimal string.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>val</code> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <code>off</code> &lt; 0, <code>len</code>
     * &lt; 0, or <code>off+len</code> &gt; <code>val.length</code>.
     */
    public static String toString(byte[] val,
            int off,
            int len) {
        if (off < 0 || len < 0 || off + len > val.length) {
            throw new IndexOutOfBoundsException();
        }
        StringBuilder buf = new StringBuilder(2 * len);
        while (len > 0) {
            buf.append(int2hex[(val[off] >> 4) & 0xF]);
            buf.append(int2hex[(val[off]) & 0xF]);
            ++off;
            --len;
        }
        return buf.toString();
    }

    /**
     * Convert the given hexadecimal string to a <code>byte</code> value. Characters
     * of the string from the highest index to the lowest index give the
     * hexadecimal digits of the value from least significant digit to most
     * significant digit. Any extra high-order digits in the string are omitted.
     *
     * @param str Hexadecimal string.
     * @return Value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static byte toByte(String str) {
        return (byte) toInt(str);
    }

    /**
     * Convert the given hexadecimal string to a <code>char</code> value. Characters
     * of the string from the highest index to the lowest index give the
     * hexadecimal digits of the value from least significant digit to most
     * significant digit. Any extra high-order digits in the string are omitted.
     *
     * @param str Hexadecimal string.
     * @return Value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static char toChar(String str) {
        return (char) toInt(str);
    }

    /**
     * Convert the given hexadecimal string to a <code>short</code> value.
     * Characters of the string from the highest index to the lowest index give
     * the hexadecimal digits of the value from least significant digit to most
     * significant digit. Any extra high-order digits in the string are omitted.
     *
     * @param str Hexadecimal string.
     * @return Value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static short toShort(String str) {
        return (short) toInt(str);
    }

    /**
     * Convert the given hexadecimal string to an <code>int</code> value. Characters
     * of the string from the highest index to the lowest index give the
     * hexadecimal digits of the value from least significant digit to most
     * significant digit. Any extra high-order digits in the string are omitted.
     *
     * @param str Hexadecimal string.
     * @return Value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static int toInt(String str) {
        int result = 0;
        int n = str.length();
        for (int i = 0; i < n; ++i) {
            result = (result << 4) + hex2int(str.charAt(i));
        }
        return result;
    }

    /**
     * Convert the given hexadecimal string to a <code>long</code> value. Characters
     * of the string from the highest index to the lowest index give the
     * hexadecimal digits of the value from least significant digit to most
     * significant digit. Any extra high-order digits in the string are omitted.
     *
     * @param str Hexadecimal string.
     * @return Value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static long toLong(String str) {
        long result = 0L;
        int n = str.length();
        for (int i = 0; i < n; ++i) {
            result = (result << 4) + hex2int(str.charAt(i));
        }
        return result;
    }

    /**
     * Convert the given hexadecimal string to a sequence of bytes stored in a
     * new <code>byte</code> array. Characters of the string from the highest index
     * to the lowest index give the hexadecimal digits of the value from least
     * significant digit to most significant digit. The value is stored in a
     * newly allocated byte array, with the least significant byte stored at the
     * highest index and the most significant byte stored at the lowest index.
     * The array's length is sufficient to hold the entire converted string. The
     * newly allocated byte array is returned.
     *
     * @param str Hexadecimal string.
     * @return Value (byte array).
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static byte[] toByteArray(String str) {
        int n = (str.length() + 1) / 2;
        byte[] val = new byte[n];
        toByteArray(str, val, 0, val.length);
        return val;
    }

    /**
     * Convert the given hexadecimal string to a sequence of bytes stored in the
     * given <code>byte</code> array. Characters of the string from the highest
     * index to the lowest index give the hexadecimal digits of the value from
     * least significant digit to most significant digit. The value is stored in
     * <code>val</code> from highest index to lowest index, with the least
     * significant byte stored at index <code>val.length-1</code>. If the converted
     * string requires more than <code>val.length</code> bytes, the extra digits at
     * the beginning of the string are not converted. If the converted string
     * requires fewer than <code>val.length</code> bytes, the extra bytes at the
     * beginning of <code>val</code> are set to 0.
     *
     * @param str Hexadecimal string.
     * @param val Byte array in which to store converted value.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null or <code>val</code>
     * is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> ..
     * <code>'f'</code>, or <code>'A'</code> .. <code>'F'</code>).
     */
    public static void toByteArray(String str,
            byte[] val) {
        toByteArray(str, val, 0, val.length);
    }

    /**
     * Convert the given hexadecimal string to a sequence of bytes stored in a
     * portion of the given <code>byte</code> array. Characters of the string from
     * the highest index to the lowest index give the hexadecimal digits of the
     * value from least significant digit to most significant digit. The value
     * is stored in <code>val</code> from highest index to lowest index, with the
     * least significant byte stored at index <code>off+len-1</code>. If the
     * converted string requires more than <code>len</code> bytes, the extra digits
     * at the beginning of the string are not converted. If the converted string
     * requires fewer than <code>len</code> bytes, the extra bytes starting at index
     * <code>off</code> of <code>val</code> are set to 0.
     *
     * @param str Hexadecimal string.
     * @param val Byte array in which to store converted value.
     * @param off Index of first byte to store.
     * @param len Number of bytes to store.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>str</code> is null or <code>val</code>
     * is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>off</code> &lt; 0, <code>len</code>
     * &lt; 0, or <code>off+len</code> &gt; <code>val.length</code>. Thrown if any
     * character of the string is not a hexadecimal digit (<code>'0'</code> ..
     * <code>'9'</code>, <code>'a'</code> .. <code>'f'</code>, or <code>'A'</code> ..
     * <code>'F'</code>).
     */
    public static void toByteArray(String str,
            byte[] val,
            int off,
            int len) {
        if (off < 0 || len < 0 || off + len > val.length) {
            throw new IndexOutOfBoundsException();
        }
        int stroff = str.length() - 1;
        int valoff = off + len - 1;
        int result;
        while (len > 0 && stroff >= 0) {
            result = hex2int(str.charAt(stroff));
            --stroff;
            if (stroff >= 0) {
                result += hex2int(str.charAt(stroff)) << 4;
            }
            --stroff;
            val[valoff] = (byte) result;
            --valoff;
            --len;
        }
        while (len > 0) {
            val[valoff] = (byte) 0;
            --valoff;
            --len;
        }
    }

// Hidden operations.
    /**
     * Convert the given hexadecimal digit to an <code>int</code> value.
     *
     * @param digit Hexadecimal digit.
     *
     * @return Value.
     *
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>digit</code> is not a hexadecimal digit (<code>'0'</code> .. <code>'9'</code>,
     * <code>'a'</code> .. <code>'f'</code>, or
     * <code>'A'</code> .. <code>'F'</code>).
     */
    private static int hex2int(char digit) {
        switch (digit) {
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                return digit - '0';
            case 'a':
            case 'b':
            case 'c':
            case 'd':
            case 'e':
            case 'f':
                return digit - 'a' + 10;
            case 'A':
            case 'B':
            case 'C':
            case 'D':
            case 'E':
            case 'F':
                return digit - 'A' + 10;
            default:
                throw new IllegalArgumentException("Not a hexadecimal digit: '" + digit + "'");
        }
    }

}
