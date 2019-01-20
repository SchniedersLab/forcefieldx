//******************************************************************************
//
// File:    DataInputStream.java
// Package: edu.rit.io
// Unit:    Class edu.rit.io.DataInputStream
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

import java.io.EOFException;
import java.io.FilterInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Class DataInputStream provides an input stream that reads primitive data
 * types and strings in binary form. It behaves similarly to class
 * java.io.DataInputStream, except the methods for reading types byte, short,
 * char, int, long, and String are implemented differently. These methods read
 * an integer value using a variable number of bytes, as described below. This
 * can save space in the file if small integer values are written more
 * frequently than large integer values. The byte stream being read must have
 * been written using class {@linkplain DataOutputStream}.
 * <P>
 * Note that class DataInputStream does <I>not</I> implement interface
 * java.io.DataInput, because the methods do not obey the contract specified in
 * that interface.
 *
 * @author Alan Kaminsky
 * @version 19-Dec-2009
 */
public class DataInputStream
        extends FilterInputStream {

// Exported constructors.
    /**
     * Construct a new data input stream.
     *
     * @param in Underlying input stream.
     */
    public DataInputStream(InputStream in) {
        super(in);
    }

// Exported operations.
    /**
     * Read a Boolean value from this data input stream. One byte is read; if it
     * is 0, false is returned; otherwise true is returned.
     *
     * @return Boolean value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public boolean readBoolean()
            throws IOException {
        int v = in.read();
        if (v == -1) {
            throw new EOFException();
        }
        return v != 0;
    }

    /**
     * Read a byte value from this data input stream. An int is read using
     * <code>readInt()</code>, the int is converted to a byte, and the byte is
     * returned.
     *
     * @return Byte value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public byte readByte()
            throws IOException {
        return (byte) readInt();
    }

    /**
     * Read a short value from this data input stream. An int is read using
     * <code>readInt()</code>, the int is converted to a short, and the short is
     * returned.
     *
     * @return Short value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public short readShort()
            throws IOException {
        return (short) readInt();
    }

    /**
     * Read a character value from this data input stream. An int is read using
     * <code>readInt()</code>, the int is converted to a char, and the char is
     * returned.
     *
     * @return Character value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public char readChar()
            throws IOException {
        return (char) readInt();
    }

    /**
     * Read an integer value from this data input stream. From one to five bytes
     * are read, in big-endian order, as follows:
     * <UL>
     *
     * <LI>
     * If the first byte's most significant bit is 0, then the least significant
     * 7 bits are sign-extended to an int and returned.
     *
     * <LI>
     * Else if the first byte's two most significant bits are 10, then one more
     * byte is read, and the least significant 14 bits are sign-extended to an
     * int and returned.
     *
     * <LI>
     * Else if the first byte's three most significant bits are 110, then two
     * more bytes are read, and the least significant 21 bits are sign-extended
     * to an int and returned.
     *
     * <LI>
     * Else if the first byte's four most significant bits are 1110, then three
     * more bytes are read, and the least significant 28 bits are sign-extended
     * to an int and returned.
     *
     * <LI>
     * Else four more bytes are read, and the least significant 32 bits are
     * returned.
     * </UL>
     *
     * @return Integer value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public int readInt()
            throws IOException {
        int v;
        int b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = b;
        if ((b & 0x80) == 0x00) {
            v <<= 25;
            v >>= 25;
            return v;
        } else if ((b & 0xC0) == 0x80) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            v <<= 18;
            v >>= 18;
            return v;
        } else if ((b & 0xE0) == 0xC0) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            v <<= 11;
            v >>= 11;
            return v;
        } else if ((b & 0xF0) == 0xE0) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            v <<= 4;
            v >>= 4;
            return v;
        } else {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            return v;
        }
    }

    /**
     * Read an unsigned byte value from this data input stream. An int is read
     * using <code>readUnsignedInt()</code>, the int is converted to a byte, and the
     * byte is returned.
     *
     * @return Byte value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public byte readUnsignedByte()
            throws IOException {
        return (byte) readUnsignedInt();
    }

    /**
     * Read an unsigned short value from this data input stream. An int is read
     * using <code>readUnsignedInt()</code>, the int is converted to a short, and
     * the short is returned.
     *
     * @return Short value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public short readUnsignedShort()
            throws IOException {
        return (short) readUnsignedInt();
    }

    /**
     * Read an unsigned character value from this data input stream. An int is
     * read using <code>readUnsignedInt()</code>, the int is converted to a char,
     * and the char is returned.
     *
     * @return Character value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public char readUnsignedChar()
            throws IOException {
        return (char) readUnsignedInt();
    }

    /**
     * Read an unsigned integer value from this data input stream. From one to
     * five bytes are read, in big-endian order, as follows:
     * <UL>
     *
     * <LI>
     * If the first byte's most significant bit is 0, then the least significant
     * 7 bits are zero-extended to an int and returned.
     *
     * <LI>
     * Else if the first byte's two most significant bits are 10, then one more
     * byte is read, and the least significant 14 bits are zero-extended to an
     * int and returned.
     *
     * <LI>
     * Else if the first byte's three most significant bits are 110, then two
     * more bytes are read, and the least significant 21 bits are zero-extended
     * to an int and returned.
     *
     * <LI>
     * Else if the first byte's four most significant bits are 1110, then three
     * more bytes are read, and the least significant 28 bits are zero-extended
     * to an int and returned.
     *
     * <LI>
     * Else four more bytes are read, and the least significant 32 bits are
     * returned.
     * </UL>
     *
     * @return Integer value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public int readUnsignedInt()
            throws IOException {
        int v;
        int b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = b;
        if ((b & 0x80) == 0x00) {
            return v;
        } else if ((b & 0xC0) == 0x80) {
            v &= 0x3F;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            return v;
        } else if ((b & 0xE0) == 0xC0) {
            v &= 0x1F;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            return v;
        } else if ((b & 0xF0) == 0xE0) {
            v &= 0x0F;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            return v;
        } else {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8) | b;
            return v;
        }
    }

    /**
     * Read a long value from this data input stream. From one to nine bytes are
     * read, in big-endian order, as follows:
     * <UL>
     *
     * <LI>
     * If the first byte's most significant bit is 0, then the least significant
     * 7 bits are sign-extended to a long and returned.
     *
     * <LI>
     * Else if the first byte's two most significant bits are 10, then one more
     * byte is read, and the least significant 14 bits are sign-extended to a
     * long and returned.
     *
     * <LI>
     * Else if the first byte's three most significant bits are 110, then two
     * more bytes are read, and the least significant 21 bits are sign-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's four most significant bits are 1110, then three
     * more bytes are read, and the least significant 28 bits are sign-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's five most significant bits are 11110, then four
     * more bytes are read, and the least significant 35 bits are sign-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's six most significant bits are 111110, then five
     * more bytes are read, and the least significant 42 bits are sign-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's seven most significant bits are 1111110, then
     * six more bytes are read, and the least significant 49 bits are
     * sign-extended to a long and returned.
     *
     * <LI>
     * Else if the first byte is 11111110, then seven more bytes are read, and
     * the least significant 56 bits are sign-extended to a long and returned.
     *
     * <LI>
     * Else eight more bytes are read, and the least significant 64 bits are
     * returned.
     * </UL>
     *
     * @return Long value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public long readLong()
            throws IOException {
        long v;
        int b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = b;
        if ((b & 0x80) == 0x00) {
            v <<= 57L;
            v >>= 57L;
            return v;
        } else if ((b & 0xC0) == 0x80) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 50L;
            v >>= 50L;
            return v;
        } else if ((b & 0xE0) == 0xC0) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 43L;
            v >>= 43L;
            return v;
        } else if ((b & 0xF0) == 0xE0) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 36L;
            v >>= 36L;
            return v;
        } else if ((b & 0xF8) == 0xF0) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 29L;
            v >>= 29L;
            return v;
        } else if ((b & 0xFC) == 0xF8) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 22L;
            v >>= 22L;
            return v;
        } else if ((b & 0xFE) == 0xFC) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 15L;
            v >>= 15L;
            return v;
        } else if (b == 0xFE) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            v <<= 8L;
            v >>= 8L;
            return v;
        } else {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        }
    }

    /**
     * Read an unsigned long value from this data input stream. From one to nine
     * bytes are read, in big-endian order, as follows:
     * <UL>
     *
     * <LI>
     * If the first byte's most significant bit is 0, then the least significant
     * 7 bits are zero-extended to a long and returned.
     *
     * <LI>
     * Else if the first byte's two most significant bits are 10, then one more
     * byte is read, and the least significant 14 bits are zero-extended to a
     * long and returned.
     *
     * <LI>
     * Else if the first byte's three most significant bits are 110, then two
     * more bytes are read, and the least significant 21 bits are zero-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's four most significant bits are 1110, then three
     * more bytes are read, and the least significant 28 bits are zero-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's five most significant bits are 11110, then four
     * more bytes are read, and the least significant 35 bits are zero-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's six most significant bits are 111110, then five
     * more bytes are read, and the least significant 42 bits are zero-extended
     * to a long and returned.
     *
     * <LI>
     * Else if the first byte's seven most significant bits are 1111110, then
     * six more bytes are read, and the least significant 49 bits are
     * zero-extended to a long and returned.
     *
     * <LI>
     * Else if the first byte is 11111110, then seven more bytes are read, and
     * the least significant 56 bits are zero-extended to a long and returned.
     *
     * <LI>
     * Else eight more bytes are read, and the least significant 64 bits are
     * returned.
     * </UL>
     *
     * @return Long value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public long readUnsignedLong()
            throws IOException {
        long v;
        int b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = b;
        if ((b & 0x80) == 0x00) {
            return v;
        } else if ((b & 0xC0) == 0x80) {
            v &= 0x3FL;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else if ((b & 0xE0) == 0xC0) {
            v &= 0x1FL;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else if ((b & 0xF0) == 0xE0) {
            v &= 0x0FL;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else if ((b & 0xF8) == 0xF0) {
            v &= 0x07L;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else if ((b & 0xFC) == 0xF8) {
            v &= 0x03L;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else if ((b & 0xFE) == 0xFC) {
            v &= 0x01L;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else if (b == 0xFE) {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        } else {
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            b = in.read();
            if (b == -1) {
                throw new EOFException();
            }
            v = (v << 8L) | b;
            return v;
        }
    }

    /**
     * Read a float value from this data input stream. Four bytes are read in
     * big-endian order yielding an int value <code>v</code>, then
     * <code>Float.intBitsToFloat(v)</code> is returned.
     *
     * @return Float value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public float readFloat()
            throws IOException {
        int v;
        int b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8) | b;
        return Float.intBitsToFloat(v);
    }

    /**
     * Read a double value from this data input stream. Eight bytes are read in
     * big-endian order yielding a long value <code>v</code>, then
     * <code>Double.longBitsToDouble(v)</code> is returned.
     *
     * @return Double value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public double readDouble()
            throws IOException {
        long v;
        int b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        b = in.read();
        if (b == -1) {
            throw new EOFException();
        }
        v = (v << 8L) | b;
        return Double.longBitsToDouble(v);
    }

    /**
     * Read a string value from this data input stream. The length of the string
     * is read using <code>readUnsignedInt()</code>, then each character of the
     * string is read using <code>readUnsignedInt()</code>.
     *
     * @return String value.
     * @exception EOFException Thrown if the end of the stream is encountered.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public String readString()
            throws IOException {
        int n = readUnsignedInt();
        char[] v = new char[n];
        for (int i = 0; i < n; ++i) {
            v[i] = (char) readUnsignedInt();
        }
        return new String(v);
    }

}
