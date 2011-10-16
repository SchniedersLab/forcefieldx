//******************************************************************************
//
// File:    AES256Cipher.java
// Package: edu.rit.crypto.blockcipher
// Unit:    Class edu.rit.crypto.blockcipher.AES256Cipher
//
// This Java source file is copyright (C) 2005 by Alan Kaminsky. All rights
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

package edu.rit.crypto.blockcipher;

/**
 * Class AES256Cipher provides a {@linkplain BlockCipher} using the Advanced
 * Encryption Standard (AES) with a 128-bit (16-byte) block length and a 256-bit
 * (32-byte) key length. The AES specification may be found at
 * <A HREF="http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf">http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf</A>.
 * <P>
 * <I>Note:</I> Class AES256Cipher is not multiple thread safe.
 *
 * @author  Alan Kaminsky
 * @version 06-Dec-2005
 */
public class AES256Cipher
	extends BlockCipher
	{

// Hidden tables.

	/**
	 * The number of round keys.
	 */
	private static final int ROUND_KEY_COUNT = 15;

	/**
	 * The fixed AES S-box. To apply the S-box to a byte i, write 
	 * <CODE>byteSub [i]</CODE>. 
	 */
	private static final short[] byteSub = new short[]
		{
		/* 0x00 */ (short) 0x63,
		/* 0x01 */ (short) 0x7c,
		/* 0x02 */ (short) 0x77,
		/* 0x03 */ (short) 0x7b,
		/* 0x04 */ (short) 0xf2,
		/* 0x05 */ (short) 0x6b,
		/* 0x06 */ (short) 0x6f,
		/* 0x07 */ (short) 0xc5,
		/* 0x08 */ (short) 0x30,
		/* 0x09 */ (short) 0x01,
		/* 0x0A */ (short) 0x67,
		/* 0x0B */ (short) 0x2b,
		/* 0x0C */ (short) 0xfe,
		/* 0x0D */ (short) 0xd7,
		/* 0x0E */ (short) 0xab,
		/* 0x0F */ (short) 0x76,
		/* 0x10 */ (short) 0xca,
		/* 0x11 */ (short) 0x82,
		/* 0x12 */ (short) 0xc9,
		/* 0x13 */ (short) 0x7d,
		/* 0x14 */ (short) 0xfa,
		/* 0x15 */ (short) 0x59,
		/* 0x16 */ (short) 0x47,
		/* 0x17 */ (short) 0xf0,
		/* 0x18 */ (short) 0xad,
		/* 0x19 */ (short) 0xd4,
		/* 0x1A */ (short) 0xa2,
		/* 0x1B */ (short) 0xaf,
		/* 0x1C */ (short) 0x9c,
		/* 0x1D */ (short) 0xa4,
		/* 0x1E */ (short) 0x72,
		/* 0x1F */ (short) 0xc0,
		/* 0x20 */ (short) 0xb7,
		/* 0x21 */ (short) 0xfd,
		/* 0x22 */ (short) 0x93,
		/* 0x23 */ (short) 0x26,
		/* 0x24 */ (short) 0x36,
		/* 0x25 */ (short) 0x3f,
		/* 0x26 */ (short) 0xf7,
		/* 0x27 */ (short) 0xcc,
		/* 0x28 */ (short) 0x34,
		/* 0x29 */ (short) 0xa5,
		/* 0x2A */ (short) 0xe5,
		/* 0x2B */ (short) 0xf1,
		/* 0x2C */ (short) 0x71,
		/* 0x2D */ (short) 0xd8,
		/* 0x2E */ (short) 0x31,
		/* 0x2F */ (short) 0x15,
		/* 0x30 */ (short) 0x04,
		/* 0x31 */ (short) 0xc7,
		/* 0x32 */ (short) 0x23,
		/* 0x33 */ (short) 0xc3,
		/* 0x34 */ (short) 0x18,
		/* 0x35 */ (short) 0x96,
		/* 0x36 */ (short) 0x05,
		/* 0x37 */ (short) 0x9a,
		/* 0x38 */ (short) 0x07,
		/* 0x39 */ (short) 0x12,
		/* 0x3A */ (short) 0x80,
		/* 0x3B */ (short) 0xe2,
		/* 0x3C */ (short) 0xeb,
		/* 0x3D */ (short) 0x27,
		/* 0x3E */ (short) 0xb2,
		/* 0x3F */ (short) 0x75,
		/* 0x40 */ (short) 0x09,
		/* 0x41 */ (short) 0x83,
		/* 0x42 */ (short) 0x2c,
		/* 0x43 */ (short) 0x1a,
		/* 0x44 */ (short) 0x1b,
		/* 0x45 */ (short) 0x6e,
		/* 0x46 */ (short) 0x5a,
		/* 0x47 */ (short) 0xa0,
		/* 0x48 */ (short) 0x52,
		/* 0x49 */ (short) 0x3b,
		/* 0x4A */ (short) 0xd6,
		/* 0x4B */ (short) 0xb3,
		/* 0x4C */ (short) 0x29,
		/* 0x4D */ (short) 0xe3,
		/* 0x4E */ (short) 0x2f,
		/* 0x4F */ (short) 0x84,
		/* 0x50 */ (short) 0x53,
		/* 0x51 */ (short) 0xd1,
		/* 0x52 */ (short) 0x00,
		/* 0x53 */ (short) 0xed,
		/* 0x54 */ (short) 0x20,
		/* 0x55 */ (short) 0xfc,
		/* 0x56 */ (short) 0xb1,
		/* 0x57 */ (short) 0x5b,
		/* 0x58 */ (short) 0x6a,
		/* 0x59 */ (short) 0xcb,
		/* 0x5A */ (short) 0xbe,
		/* 0x5B */ (short) 0x39,
		/* 0x5C */ (short) 0x4a,
		/* 0x5D */ (short) 0x4c,
		/* 0x5E */ (short) 0x58,
		/* 0x5F */ (short) 0xcf,
		/* 0x60 */ (short) 0xd0,
		/* 0x61 */ (short) 0xef,
		/* 0x62 */ (short) 0xaa,
		/* 0x63 */ (short) 0xfb,
		/* 0x64 */ (short) 0x43,
		/* 0x65 */ (short) 0x4d,
		/* 0x66 */ (short) 0x33,
		/* 0x67 */ (short) 0x85,
		/* 0x68 */ (short) 0x45,
		/* 0x69 */ (short) 0xf9,
		/* 0x6A */ (short) 0x02,
		/* 0x6B */ (short) 0x7f,
		/* 0x6C */ (short) 0x50,
		/* 0x6D */ (short) 0x3c,
		/* 0x6E */ (short) 0x9f,
		/* 0x6F */ (short) 0xa8,
		/* 0x70 */ (short) 0x51,
		/* 0x71 */ (short) 0xa3,
		/* 0x72 */ (short) 0x40,
		/* 0x73 */ (short) 0x8f,
		/* 0x74 */ (short) 0x92,
		/* 0x75 */ (short) 0x9d,
		/* 0x76 */ (short) 0x38,
		/* 0x77 */ (short) 0xf5,
		/* 0x78 */ (short) 0xbc,
		/* 0x79 */ (short) 0xb6,
		/* 0x7A */ (short) 0xda,
		/* 0x7B */ (short) 0x21,
		/* 0x7C */ (short) 0x10,
		/* 0x7D */ (short) 0xff,
		/* 0x7E */ (short) 0xf3,
		/* 0x7F */ (short) 0xd2,
		/* 0x80 */ (short) 0xcd,
		/* 0x81 */ (short) 0x0c,
		/* 0x82 */ (short) 0x13,
		/* 0x83 */ (short) 0xec,
		/* 0x84 */ (short) 0x5f,
		/* 0x85 */ (short) 0x97,
		/* 0x86 */ (short) 0x44,
		/* 0x87 */ (short) 0x17,
		/* 0x88 */ (short) 0xc4,
		/* 0x89 */ (short) 0xa7,
		/* 0x8A */ (short) 0x7e,
		/* 0x8B */ (short) 0x3d,
		/* 0x8C */ (short) 0x64,
		/* 0x8D */ (short) 0x5d,
		/* 0x8E */ (short) 0x19,
		/* 0x8F */ (short) 0x73,
		/* 0x90 */ (short) 0x60,
		/* 0x91 */ (short) 0x81,
		/* 0x92 */ (short) 0x4f,
		/* 0x93 */ (short) 0xdc,
		/* 0x94 */ (short) 0x22,
		/* 0x95 */ (short) 0x2a,
		/* 0x96 */ (short) 0x90,
		/* 0x97 */ (short) 0x88,
		/* 0x98 */ (short) 0x46,
		/* 0x99 */ (short) 0xee,
		/* 0x9A */ (short) 0xb8,
		/* 0x9B */ (short) 0x14,
		/* 0x9C */ (short) 0xde,
		/* 0x9D */ (short) 0x5e,
		/* 0x9E */ (short) 0x0b,
		/* 0x9F */ (short) 0xdb,
		/* 0xA0 */ (short) 0xe0,
		/* 0xA1 */ (short) 0x32,
		/* 0xA2 */ (short) 0x3a,
		/* 0xA3 */ (short) 0x0a,
		/* 0xA4 */ (short) 0x49,
		/* 0xA5 */ (short) 0x06,
		/* 0xA6 */ (short) 0x24,
		/* 0xA7 */ (short) 0x5c,
		/* 0xA8 */ (short) 0xc2,
		/* 0xA9 */ (short) 0xd3,
		/* 0xAA */ (short) 0xac,
		/* 0xAB */ (short) 0x62,
		/* 0xAC */ (short) 0x91,
		/* 0xAD */ (short) 0x95,
		/* 0xAE */ (short) 0xe4,
		/* 0xAF */ (short) 0x79,
		/* 0xB0 */ (short) 0xe7,
		/* 0xB1 */ (short) 0xc8,
		/* 0xB2 */ (short) 0x37,
		/* 0xB3 */ (short) 0x6d,
		/* 0xB4 */ (short) 0x8d,
		/* 0xB5 */ (short) 0xd5,
		/* 0xB6 */ (short) 0x4e,
		/* 0xB7 */ (short) 0xa9,
		/* 0xB8 */ (short) 0x6c,
		/* 0xB9 */ (short) 0x56,
		/* 0xBA */ (short) 0xf4,
		/* 0xBB */ (short) 0xea,
		/* 0xBC */ (short) 0x65,
		/* 0xBD */ (short) 0x7a,
		/* 0xBE */ (short) 0xae,
		/* 0xBF */ (short) 0x08,
		/* 0xC0 */ (short) 0xba,
		/* 0xC1 */ (short) 0x78,
		/* 0xC2 */ (short) 0x25,
		/* 0xC3 */ (short) 0x2e,
		/* 0xC4 */ (short) 0x1c,
		/* 0xC5 */ (short) 0xa6,
		/* 0xC6 */ (short) 0xb4,
		/* 0xC7 */ (short) 0xc6,
		/* 0xC8 */ (short) 0xe8,
		/* 0xC9 */ (short) 0xdd,
		/* 0xCA */ (short) 0x74,
		/* 0xCB */ (short) 0x1f,
		/* 0xCC */ (short) 0x4b,
		/* 0xCD */ (short) 0xbd,
		/* 0xCE */ (short) 0x8b,
		/* 0xCF */ (short) 0x8a,
		/* 0xD0 */ (short) 0x70,
		/* 0xD1 */ (short) 0x3e,
		/* 0xD2 */ (short) 0xb5,
		/* 0xD3 */ (short) 0x66,
		/* 0xD4 */ (short) 0x48,
		/* 0xD5 */ (short) 0x03,
		/* 0xD6 */ (short) 0xf6,
		/* 0xD7 */ (short) 0x0e,
		/* 0xD8 */ (short) 0x61,
		/* 0xD9 */ (short) 0x35,
		/* 0xDA */ (short) 0x57,
		/* 0xDB */ (short) 0xb9,
		/* 0xDC */ (short) 0x86,
		/* 0xDD */ (short) 0xc1,
		/* 0xDE */ (short) 0x1d,
		/* 0xDF */ (short) 0x9e,
		/* 0xE0 */ (short) 0xe1,
		/* 0xE1 */ (short) 0xf8,
		/* 0xE2 */ (short) 0x98,
		/* 0xE3 */ (short) 0x11,
		/* 0xE4 */ (short) 0x69,
		/* 0xE5 */ (short) 0xd9,
		/* 0xE6 */ (short) 0x8e,
		/* 0xE7 */ (short) 0x94,
		/* 0xE8 */ (short) 0x9b,
		/* 0xE9 */ (short) 0x1e,
		/* 0xEA */ (short) 0x87,
		/* 0xEB */ (short) 0xe9,
		/* 0xEC */ (short) 0xce,
		/* 0xED */ (short) 0x55,
		/* 0xEE */ (short) 0x28,
		/* 0xEF */ (short) 0xdf,
		/* 0xF0 */ (short) 0x8c,
		/* 0xF1 */ (short) 0xa1,
		/* 0xF2 */ (short) 0x89,
		/* 0xF3 */ (short) 0x0d,
		/* 0xF4 */ (short) 0xbf,
		/* 0xF5 */ (short) 0xe6,
		/* 0xF6 */ (short) 0x42,
		/* 0xF7 */ (short) 0x68,
		/* 0xF8 */ (short) 0x41,
		/* 0xF9 */ (short) 0x99,
		/* 0xFA */ (short) 0x2d,
		/* 0xFB */ (short) 0x0f,
		/* 0xFC */ (short) 0xb0,
		/* 0xFD */ (short) 0x54,
		/* 0xFE */ (short) 0xbb,
		/* 0xFF */ (short) 0x16,
		};

	/**
	 * Table for computing the xtime function. To apply the xtime function to a 
	 * short i, write <CODE>xtime [i]</CODE>. 
	 */
	private static final short[] xtime = new short[]
		{
		/* 0x00 */ (short) 0x00,
		/* 0x01 */ (short) 0x02,
		/* 0x02 */ (short) 0x04,
		/* 0x03 */ (short) 0x06,
		/* 0x04 */ (short) 0x08,
		/* 0x05 */ (short) 0x0A,
		/* 0x06 */ (short) 0x0C,
		/* 0x07 */ (short) 0x0E,
		/* 0x08 */ (short) 0x10,
		/* 0x09 */ (short) 0x12,
		/* 0x0A */ (short) 0x14,
		/* 0x0B */ (short) 0x16,
		/* 0x0C */ (short) 0x18,
		/* 0x0D */ (short) 0x1A,
		/* 0x0E */ (short) 0x1C,
		/* 0x0F */ (short) 0x1E,
		/* 0x10 */ (short) 0x20,
		/* 0x11 */ (short) 0x22,
		/* 0x12 */ (short) 0x24,
		/* 0x13 */ (short) 0x26,
		/* 0x14 */ (short) 0x28,
		/* 0x15 */ (short) 0x2A,
		/* 0x16 */ (short) 0x2C,
		/* 0x17 */ (short) 0x2E,
		/* 0x18 */ (short) 0x30,
		/* 0x19 */ (short) 0x32,
		/* 0x1A */ (short) 0x34,
		/* 0x1B */ (short) 0x36,
		/* 0x1C */ (short) 0x38,
		/* 0x1D */ (short) 0x3A,
		/* 0x1E */ (short) 0x3C,
		/* 0x1F */ (short) 0x3E,
		/* 0x20 */ (short) 0x40,
		/* 0x21 */ (short) 0x42,
		/* 0x22 */ (short) 0x44,
		/* 0x23 */ (short) 0x46,
		/* 0x24 */ (short) 0x48,
		/* 0x25 */ (short) 0x4A,
		/* 0x26 */ (short) 0x4C,
		/* 0x27 */ (short) 0x4E,
		/* 0x28 */ (short) 0x50,
		/* 0x29 */ (short) 0x52,
		/* 0x2A */ (short) 0x54,
		/* 0x2B */ (short) 0x56,
		/* 0x2C */ (short) 0x58,
		/* 0x2D */ (short) 0x5A,
		/* 0x2E */ (short) 0x5C,
		/* 0x2F */ (short) 0x5E,
		/* 0x30 */ (short) 0x60,
		/* 0x31 */ (short) 0x62,
		/* 0x32 */ (short) 0x64,
		/* 0x33 */ (short) 0x66,
		/* 0x34 */ (short) 0x68,
		/* 0x35 */ (short) 0x6A,
		/* 0x36 */ (short) 0x6C,
		/* 0x37 */ (short) 0x6E,
		/* 0x38 */ (short) 0x70,
		/* 0x39 */ (short) 0x72,
		/* 0x3A */ (short) 0x74,
		/* 0x3B */ (short) 0x76,
		/* 0x3C */ (short) 0x78,
		/* 0x3D */ (short) 0x7A,
		/* 0x3E */ (short) 0x7C,
		/* 0x3F */ (short) 0x7E,
		/* 0x40 */ (short) 0x80,
		/* 0x41 */ (short) 0x82,
		/* 0x42 */ (short) 0x84,
		/* 0x43 */ (short) 0x86,
		/* 0x44 */ (short) 0x88,
		/* 0x45 */ (short) 0x8A,
		/* 0x46 */ (short) 0x8C,
		/* 0x47 */ (short) 0x8E,
		/* 0x48 */ (short) 0x90,
		/* 0x49 */ (short) 0x92,
		/* 0x4A */ (short) 0x94,
		/* 0x4B */ (short) 0x96,
		/* 0x4C */ (short) 0x98,
		/* 0x4D */ (short) 0x9A,
		/* 0x4E */ (short) 0x9C,
		/* 0x4F */ (short) 0x9E,
		/* 0x50 */ (short) 0xA0,
		/* 0x51 */ (short) 0xA2,
		/* 0x52 */ (short) 0xA4,
		/* 0x53 */ (short) 0xA6,
		/* 0x54 */ (short) 0xA8,
		/* 0x55 */ (short) 0xAA,
		/* 0x56 */ (short) 0xAC,
		/* 0x57 */ (short) 0xAE,
		/* 0x58 */ (short) 0xB0,
		/* 0x59 */ (short) 0xB2,
		/* 0x5A */ (short) 0xB4,
		/* 0x5B */ (short) 0xB6,
		/* 0x5C */ (short) 0xB8,
		/* 0x5D */ (short) 0xBA,
		/* 0x5E */ (short) 0xBC,
		/* 0x5F */ (short) 0xBE,
		/* 0x60 */ (short) 0xC0,
		/* 0x61 */ (short) 0xC2,
		/* 0x62 */ (short) 0xC4,
		/* 0x63 */ (short) 0xC6,
		/* 0x64 */ (short) 0xC8,
		/* 0x65 */ (short) 0xCA,
		/* 0x66 */ (short) 0xCC,
		/* 0x67 */ (short) 0xCE,
		/* 0x68 */ (short) 0xD0,
		/* 0x69 */ (short) 0xD2,
		/* 0x6A */ (short) 0xD4,
		/* 0x6B */ (short) 0xD6,
		/* 0x6C */ (short) 0xD8,
		/* 0x6D */ (short) 0xDA,
		/* 0x6E */ (short) 0xDC,
		/* 0x6F */ (short) 0xDE,
		/* 0x70 */ (short) 0xE0,
		/* 0x71 */ (short) 0xE2,
		/* 0x72 */ (short) 0xE4,
		/* 0x73 */ (short) 0xE6,
		/* 0x74 */ (short) 0xE8,
		/* 0x75 */ (short) 0xEA,
		/* 0x76 */ (short) 0xEC,
		/* 0x77 */ (short) 0xEE,
		/* 0x78 */ (short) 0xF0,
		/* 0x79 */ (short) 0xF2,
		/* 0x7A */ (short) 0xF4,
		/* 0x7B */ (short) 0xF6,
		/* 0x7C */ (short) 0xF8,
		/* 0x7D */ (short) 0xFA,
		/* 0x7E */ (short) 0xFC,
		/* 0x7F */ (short) 0xFE,
		/* 0x80 */ (short) 0x1B,
		/* 0x81 */ (short) 0x19,
		/* 0x82 */ (short) 0x1F,
		/* 0x83 */ (short) 0x1D,
		/* 0x84 */ (short) 0x13,
		/* 0x85 */ (short) 0x11,
		/* 0x86 */ (short) 0x17,
		/* 0x87 */ (short) 0x15,
		/* 0x88 */ (short) 0x0B,
		/* 0x89 */ (short) 0x09,
		/* 0x8A */ (short) 0x0F,
		/* 0x8B */ (short) 0x0D,
		/* 0x8C */ (short) 0x03,
		/* 0x8D */ (short) 0x01,
		/* 0x8E */ (short) 0x07,
		/* 0x8F */ (short) 0x05,
		/* 0x90 */ (short) 0x3B,
		/* 0x91 */ (short) 0x39,
		/* 0x92 */ (short) 0x3F,
		/* 0x93 */ (short) 0x3D,
		/* 0x94 */ (short) 0x33,
		/* 0x95 */ (short) 0x31,
		/* 0x96 */ (short) 0x37,
		/* 0x97 */ (short) 0x35,
		/* 0x98 */ (short) 0x2B,
		/* 0x99 */ (short) 0x29,
		/* 0x9A */ (short) 0x2F,
		/* 0x9B */ (short) 0x2D,
		/* 0x9C */ (short) 0x23,
		/* 0x9D */ (short) 0x21,
		/* 0x9E */ (short) 0x27,
		/* 0x9F */ (short) 0x25,
		/* 0xA0 */ (short) 0x5B,
		/* 0xA1 */ (short) 0x59,
		/* 0xA2 */ (short) 0x5F,
		/* 0xA3 */ (short) 0x5D,
		/* 0xA4 */ (short) 0x53,
		/* 0xA5 */ (short) 0x51,
		/* 0xA6 */ (short) 0x57,
		/* 0xA7 */ (short) 0x55,
		/* 0xA8 */ (short) 0x4B,
		/* 0xA9 */ (short) 0x49,
		/* 0xAA */ (short) 0x4F,
		/* 0xAB */ (short) 0x4D,
		/* 0xAC */ (short) 0x43,
		/* 0xAD */ (short) 0x41,
		/* 0xAE */ (short) 0x47,
		/* 0xAF */ (short) 0x45,
		/* 0xB0 */ (short) 0x7B,
		/* 0xB1 */ (short) 0x79,
		/* 0xB2 */ (short) 0x7F,
		/* 0xB3 */ (short) 0x7D,
		/* 0xB4 */ (short) 0x73,
		/* 0xB5 */ (short) 0x71,
		/* 0xB6 */ (short) 0x77,
		/* 0xB7 */ (short) 0x75,
		/* 0xB8 */ (short) 0x6B,
		/* 0xB9 */ (short) 0x69,
		/* 0xBA */ (short) 0x6F,
		/* 0xBB */ (short) 0x6D,
		/* 0xBC */ (short) 0x63,
		/* 0xBD */ (short) 0x61,
		/* 0xBE */ (short) 0x67,
		/* 0xBF */ (short) 0x65,
		/* 0xC0 */ (short) 0x9B,
		/* 0xC1 */ (short) 0x99,
		/* 0xC2 */ (short) 0x9F,
		/* 0xC3 */ (short) 0x9D,
		/* 0xC4 */ (short) 0x93,
		/* 0xC5 */ (short) 0x91,
		/* 0xC6 */ (short) 0x97,
		/* 0xC7 */ (short) 0x95,
		/* 0xC8 */ (short) 0x8B,
		/* 0xC9 */ (short) 0x89,
		/* 0xCA */ (short) 0x8F,
		/* 0xCB */ (short) 0x8D,
		/* 0xCC */ (short) 0x83,
		/* 0xCD */ (short) 0x81,
		/* 0xCE */ (short) 0x87,
		/* 0xCF */ (short) 0x85,
		/* 0xD0 */ (short) 0xBB,
		/* 0xD1 */ (short) 0xB9,
		/* 0xD2 */ (short) 0xBF,
		/* 0xD3 */ (short) 0xBD,
		/* 0xD4 */ (short) 0xB3,
		/* 0xD5 */ (short) 0xB1,
		/* 0xD6 */ (short) 0xB7,
		/* 0xD7 */ (short) 0xB5,
		/* 0xD8 */ (short) 0xAB,
		/* 0xD9 */ (short) 0xA9,
		/* 0xDA */ (short) 0xAF,
		/* 0xDB */ (short) 0xAD,
		/* 0xDC */ (short) 0xA3,
		/* 0xDD */ (short) 0xA1,
		/* 0xDE */ (short) 0xA7,
		/* 0xDF */ (short) 0xA5,
		/* 0xE0 */ (short) 0xDB,
		/* 0xE1 */ (short) 0xD9,
		/* 0xE2 */ (short) 0xDF,
		/* 0xE3 */ (short) 0xDD,
		/* 0xE4 */ (short) 0xD3,
		/* 0xE5 */ (short) 0xD1,
		/* 0xE6 */ (short) 0xD7,
		/* 0xE7 */ (short) 0xD5,
		/* 0xE8 */ (short) 0xCB,
		/* 0xE9 */ (short) 0xC9,
		/* 0xEA */ (short) 0xCF,
		/* 0xEB */ (short) 0xCD,
		/* 0xEC */ (short) 0xC3,
		/* 0xED */ (short) 0xC1,
		/* 0xEE */ (short) 0xC7,
		/* 0xEF */ (short) 0xC5,
		/* 0xF0 */ (short) 0xFB,
		/* 0xF1 */ (short) 0xF9,
		/* 0xF2 */ (short) 0xFF,
		/* 0xF3 */ (short) 0xFD,
		/* 0xF4 */ (short) 0xF3,
		/* 0xF5 */ (short) 0xF1,
		/* 0xF6 */ (short) 0xF7,
		/* 0xF7 */ (short) 0xF5,
		/* 0xF8 */ (short) 0xEB,
		/* 0xF9 */ (short) 0xE9,
		/* 0xFA */ (short) 0xEF,
		/* 0xFB */ (short) 0xED,
		/* 0xFC */ (short) 0xE3,
		/* 0xFD */ (short) 0xE1,
		/* 0xFE */ (short) 0xE7,
		/* 0xFF */ (short) 0xE5,
		};

// Hidden helper classes.

	/**
	 * The round key: a 4x4 array of bytes (represented as ints). The bytes are
	 * declared as separate fields rather than an array to avoid the overhead of
	 * array indexing.
	 */
	private static class RoundKey
		{
		public int key00, key01, key02, key03;
		public int key10, key11, key12, key13;
		public int key20, key21, key22, key23;
		public int key30, key31, key32, key33;

		public void erase()
			{
			key00 = 0;
			key01 = 0;
			key02 = 0;
			key03 = 0;
			key10 = 0;
			key11 = 0;
			key12 = 0;
			key13 = 0;
			key20 = 0;
			key21 = 0;
			key22 = 0;
			key23 = 0;
			key30 = 0;
			key31 = 0;
			key32 = 0;
			key33 = 0;
			}
		}

// Hidden data members.

	/**
	 * The state: a 4x4 array of bytes (represented as ints). The bytes are
	 * declared as separate fields rather than an array to avoid the overhead of
	 * array indexing.
	 */
	private int state00, state01, state02, state03;
	private int state10, state11, state12, state13;
	private int state20, state21, state22, state23;
	private int state30, state31, state32, state33;

	/**
	 * The expanded key: an array of round keys derived from the encryption key
	 * using the key schedule.
	 */
	private RoundKey[] expandedKey;

// Exported constructors.

	/**
	 * Construct a new AES-256 cipher object with the given secret key. After
	 * the constructor returns, <TT>theKey</TT> is no longer needed and may be
	 * erased (set to all 0s). The length of <TT>theKey</TT> must be at least
	 * 32. Only the first 32 bytes are used.
	 *
	 * @param  theKey  Key (byte array).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theKey</TT> is null. 
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theKey.length</TT> &lt; 32.
	 */
	public AES256Cipher
		(byte[] theKey)
		{
		super (16, 32);
		setKey (theKey);
		}

// Exported operations.

	/**
	 * Set the key to be used for all subsequent encryptions and decryptions.
	 * After this method returns, <TT>theKey</TT> is no longer needed and may be
	 * erased (set to all 0s). The length of <TT>theKey</TT> must be at least
	 * 32. Only the first 32 bytes are used.
	 *
	 * @param  theKey  Key (byte array).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theKey</TT> is null. 
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theKey.length</TT> &lt; 32.
	 */
	public void setKey
		(byte[] theKey)
		{
		// Verify arguments.
		if (theKey.length < 32)
			{
			throw new IllegalArgumentException ("Key is not 32 bytes");
			}

		// Allocate storage for round keys if necessary.
		if (expandedKey == null)
			{
			expandedKey = new RoundKey [ROUND_KEY_COUNT];
			for (int i = 0; i < ROUND_KEY_COUNT; ++ i)
				{
				expandedKey[i] = new RoundKey();
				}
			}

		// Fill in first two round keys from theKey.
		RoundKey roundKey;

		roundKey = expandedKey[0];
		roundKey.key00 = theKey[ 0] & 0xFF;
		roundKey.key10 = theKey[ 1] & 0xFF;
		roundKey.key20 = theKey[ 2] & 0xFF;
		roundKey.key30 = theKey[ 3] & 0xFF;
		roundKey.key01 = theKey[ 4] & 0xFF;
		roundKey.key11 = theKey[ 5] & 0xFF;
		roundKey.key21 = theKey[ 6] & 0xFF;
		roundKey.key31 = theKey[ 7] & 0xFF;
		roundKey.key02 = theKey[ 8] & 0xFF;
		roundKey.key12 = theKey[ 9] & 0xFF;
		roundKey.key22 = theKey[10] & 0xFF;
		roundKey.key32 = theKey[11] & 0xFF;
		roundKey.key03 = theKey[12] & 0xFF;
		roundKey.key13 = theKey[13] & 0xFF;
		roundKey.key23 = theKey[14] & 0xFF;
		roundKey.key33 = theKey[15] & 0xFF;

		roundKey = expandedKey[1];
		roundKey.key00 = theKey[16] & 0xFF;
		roundKey.key10 = theKey[17] & 0xFF;
		roundKey.key20 = theKey[18] & 0xFF;
		roundKey.key30 = theKey[19] & 0xFF;
		roundKey.key01 = theKey[20] & 0xFF;
		roundKey.key11 = theKey[21] & 0xFF;
		roundKey.key21 = theKey[22] & 0xFF;
		roundKey.key31 = theKey[23] & 0xFF;
		roundKey.key02 = theKey[24] & 0xFF;
		roundKey.key12 = theKey[25] & 0xFF;
		roundKey.key22 = theKey[26] & 0xFF;
		roundKey.key32 = theKey[27] & 0xFF;
		roundKey.key03 = theKey[28] & 0xFF;
		roundKey.key13 = theKey[29] & 0xFF;
		roundKey.key23 = theKey[30] & 0xFF;
		roundKey.key33 = theKey[31] & 0xFF;

		// Fill in remaining round keys.
		int rcon = 1;
		for (int i = 2; i < 14; i += 2)
			{
			advanceRoundKey
				(expandedKey[i-2],
				 expandedKey[i-1],
				 expandedKey[i],
				 rcon);
			rcon = xtime [rcon];
			advanceRoundKey
				(expandedKey[i-1],
				 expandedKey[i],
				 expandedKey[i+1]);
			}
		advanceRoundKey
			(expandedKey[12],
			 expandedKey[13],
			 expandedKey[14],
			 rcon);
		}

	/**
	 * Encrypt the given plaintext block. On input, <TT>thePlaintext</TT>
	 * contains the plaintext. On output, the contents of <TT>theCiphertext</TT>
	 * have been replaced by the ciphertext. <TT>thePlaintext</TT> and
	 * <TT>theCiphertext</TT> may be the same block. The length of
	 * <TT>thePlaintext</TT> must be at least 16. Only the first 16 bytes are
	 * read. The length of <TT>theCiphertext</TT> must be at least 16. Only the
	 * first 16 bytes are written.
	 *
	 * @param  thePlaintext   Input plaintext block to be encrypted.
	 * @param  theCiphertext  Output ciphertext block.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlaintext</TT> is null or
	 *     <TT>theCiphertext</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>thePlaintext.length</TT> &lt; 16
	 *     or <TT>theCiphertext.length</TT> &lt; 16.
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if the key has been erased and not
	 *     re-set.
	 */
	public void encrypt
		(byte[] thePlaintext,
		 byte[] theCiphertext)
		{
		int i;

		// Verify preconditions.
		if (thePlaintext.length < 16)
			{
			throw new IllegalArgumentException ("thePlaintext is not 16 bytes");
			}
		if (theCiphertext.length < 16)
			{
			throw new IllegalArgumentException ("theCiphertext is not 16 bytes");
			}
		if (expandedKey == null)
			{
			throw new IllegalArgumentException ("Key is not set");
			}

		// Initialize the state.
		state00 = thePlaintext[ 0] & 0xFF;
		state10 = thePlaintext[ 1] & 0xFF;
		state20 = thePlaintext[ 2] & 0xFF;
		state30 = thePlaintext[ 3] & 0xFF;
		state01 = thePlaintext[ 4] & 0xFF;
		state11 = thePlaintext[ 5] & 0xFF;
		state21 = thePlaintext[ 6] & 0xFF;
		state31 = thePlaintext[ 7] & 0xFF;
		state02 = thePlaintext[ 8] & 0xFF;
		state12 = thePlaintext[ 9] & 0xFF;
		state22 = thePlaintext[10] & 0xFF;
		state32 = thePlaintext[11] & 0xFF;
		state03 = thePlaintext[12] & 0xFF;
		state13 = thePlaintext[13] & 0xFF;
		state23 = thePlaintext[14] & 0xFF;
		state33 = thePlaintext[15] & 0xFF;

		// Do the first 13 rounds.
		for (i = 0; i < 13; ++ i)
			{
			addRoundKeyByteSubShiftRow (expandedKey[i]);
			mixColumn();
			}

		// Do the 14th round.
		addRoundKeyByteSubShiftRow (expandedKey[13]);
		addRoundKey (expandedKey[14]);

		// Output encrypted message.
		theCiphertext[ 0] = (byte) state00;
		theCiphertext[ 1] = (byte) state10;
		theCiphertext[ 2] = (byte) state20;
		theCiphertext[ 3] = (byte) state30;
		theCiphertext[ 4] = (byte) state01;
		theCiphertext[ 5] = (byte) state11;
		theCiphertext[ 6] = (byte) state21;
		theCiphertext[ 7] = (byte) state31;
		theCiphertext[ 8] = (byte) state02;
		theCiphertext[ 9] = (byte) state12;
		theCiphertext[10] = (byte) state22;
		theCiphertext[11] = (byte) state32;
		theCiphertext[12] = (byte) state03;
		theCiphertext[13] = (byte) state13;
		theCiphertext[14] = (byte) state23;
		theCiphertext[15] = (byte) state33;
		}

	/**
	 * Erase this block cipher object's key material.
	 */
	public void erase()
		{
		state00 = 0;
		state01 = 0;
		state02 = 0;
		state03 = 0;
		state10 = 0;
		state11 = 0;
		state12 = 0;
		state13 = 0;
		state20 = 0;
		state21 = 0;
		state22 = 0;
		state23 = 0;
		state30 = 0;
		state31 = 0;
		state32 = 0;
		state33 = 0;

		if (expandedKey != null)
			{
			for (int i = 0; i < ROUND_KEY_COUNT; ++ i)
				{
				expandedKey[i].erase();
				expandedKey[i] = null;
				}
			expandedKey = null;
			}
		}

	/**
	 * Finalize this block cipher object. This block cipher object is erased if
	 * it is not already erased.
	 */
	protected void finalize()
		{
		erase();
		}

// Hidden operations.

	/**
	 * Compute the next round key from the previous round key and the round 
	 * constant. 
	 */
	private static void advanceRoundKey
		(RoundKey prevprev,
		 RoundKey prev,
		 RoundKey next,
		 int rcon)
		{
		// Column 0.
		next.key00 = prevprev.key00 ^ byteSub [prev.key13] ^ rcon;
		next.key10 = prevprev.key10 ^ byteSub [prev.key23];
		next.key20 = prevprev.key20 ^ byteSub [prev.key33];
		next.key30 = prevprev.key30 ^ byteSub [prev.key03];
		// Column 1.
		next.key01 = prevprev.key01 ^ next.key00;
		next.key11 = prevprev.key11 ^ next.key10;
		next.key21 = prevprev.key21 ^ next.key20;
		next.key31 = prevprev.key31 ^ next.key30;
		// Column 2.
		next.key02 = prevprev.key02 ^ next.key01;
		next.key12 = prevprev.key12 ^ next.key11;
		next.key22 = prevprev.key22 ^ next.key21;
		next.key32 = prevprev.key32 ^ next.key31;
		// Column 3.
		next.key03 = prevprev.key03 ^ next.key02;
		next.key13 = prevprev.key13 ^ next.key12;
		next.key23 = prevprev.key23 ^ next.key22;
		next.key33 = prevprev.key33 ^ next.key32;
		}

	/**
	 * Compute the next round key from the previous round key.
	 */
	private static void advanceRoundKey
		(RoundKey prevprev,
		 RoundKey prev,
		 RoundKey next)
		{
		// Column 0.
		next.key00 = prevprev.key00 ^ byteSub [prev.key03];
		next.key10 = prevprev.key10 ^ byteSub [prev.key13];
		next.key20 = prevprev.key20 ^ byteSub [prev.key23];
		next.key30 = prevprev.key30 ^ byteSub [prev.key33];
		// Column 1.
		next.key01 = prevprev.key01 ^ next.key00;
		next.key11 = prevprev.key11 ^ next.key10;
		next.key21 = prevprev.key21 ^ next.key20;
		next.key31 = prevprev.key31 ^ next.key30;
		// Column 2.
		next.key02 = prevprev.key02 ^ next.key01;
		next.key12 = prevprev.key12 ^ next.key11;
		next.key22 = prevprev.key22 ^ next.key21;
		next.key32 = prevprev.key32 ^ next.key31;
		// Column 3.
		next.key03 = prevprev.key03 ^ next.key02;
		next.key13 = prevprev.key13 ^ next.key12;
		next.key23 = prevprev.key23 ^ next.key22;
		next.key33 = prevprev.key33 ^ next.key32;
		}

	/**
	 * To each byte of the state, add the corresponding byte of the round key.
	 */
	private void addRoundKey
		(RoundKey roundKey)
		{
		state00 ^= roundKey.key00;
		state01 ^= roundKey.key01;
		state02 ^= roundKey.key02;
		state03 ^= roundKey.key03;
		state10 ^= roundKey.key10;
		state11 ^= roundKey.key11;
		state12 ^= roundKey.key12;
		state13 ^= roundKey.key13;
		state20 ^= roundKey.key20;
		state21 ^= roundKey.key21;
		state22 ^= roundKey.key22;
		state23 ^= roundKey.key23;
		state30 ^= roundKey.key30;
		state31 ^= roundKey.key31;
		state32 ^= roundKey.key32;
		state33 ^= roundKey.key33;
		}

	/**
	 * To each byte of the state, add the corresponding byte of the round key,
	 * apply the S-box, and shift the state row.
	 */
	private void addRoundKeyByteSubShiftRow
		(RoundKey roundKey)
		{
		int temp;
		// Row 0: No shift.
		state00 = byteSub [state00 ^ roundKey.key00];
		state01 = byteSub [state01 ^ roundKey.key01];
		state02 = byteSub [state02 ^ roundKey.key02];
		state03 = byteSub [state03 ^ roundKey.key03];
		// Row 1: Left circular shift 1 byte.
		temp = state10;
		state10 = byteSub [state11 ^ roundKey.key11];
		state11 = byteSub [state12 ^ roundKey.key12];
		state12 = byteSub [state13 ^ roundKey.key13];
		state13 = byteSub [temp    ^ roundKey.key10];
		// Row 2: Left circular shift 2 bytes.
		temp = state20;
		state20 = byteSub [state22 ^ roundKey.key22];
		state22 = byteSub [temp    ^ roundKey.key20];
		temp = state21;
		state21 = byteSub [state23 ^ roundKey.key23];
		state23 = byteSub [temp    ^ roundKey.key21];
		// Row 3: Left circular shift 3 bytes = right circular shift 1 byte.
		temp = state33;
		state33 = byteSub [state32 ^ roundKey.key32];
		state32 = byteSub [state31 ^ roundKey.key31];
		state31 = byteSub [state30 ^ roundKey.key30];
		state30 = byteSub [temp    ^ roundKey.key33];
		}

	/**
	 * Multiply each column of the state by a fixed polynomial,
	 * 0x03 x^3 + 0x01 x^2 + 0x01 x + 0x02.
	 */
	private void mixColumn()
		{
		int temp, sum;
		// Column 0.
		temp = state00;
		sum = state00 ^ state10 ^ state20 ^ state30;
		state00 ^= xtime [state00 ^ state10] ^ sum;
		state10 ^= xtime [state10 ^ state20] ^ sum;
		state20 ^= xtime [state20 ^ state30] ^ sum;
		state30 ^= xtime [state30 ^ temp   ] ^ sum;
		// Column 1.
		temp = state01;
		sum = state01 ^ state11 ^ state21 ^ state31;
		state01 ^= xtime [state01 ^ state11] ^ sum;
		state11 ^= xtime [state11 ^ state21] ^ sum;
		state21 ^= xtime [state21 ^ state31] ^ sum;
		state31 ^= xtime [state31 ^ temp   ] ^ sum;
		// Column 2.
		temp = state02;
		sum = state02 ^ state12 ^ state22 ^ state32;
		state02 ^= xtime [state02 ^ state12] ^ sum;
		state12 ^= xtime [state12 ^ state22] ^ sum;
		state22 ^= xtime [state22 ^ state32] ^ sum;
		state32 ^= xtime [state32 ^ temp   ] ^ sum;
		// Column 3.
		temp = state03;
		sum = state03 ^ state13 ^ state23 ^ state33;
		state03 ^= xtime [state03 ^ state13] ^ sum;
		state13 ^= xtime [state13 ^ state23] ^ sum;
		state23 ^= xtime [state23 ^ state33] ^ sum;
		state33 ^= xtime [state33 ^ temp   ] ^ sum;
		}

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		try
//			{
//			int i;
//
//			// Verify command line arguments.
//			if (args.length != 2)
//				{
//				System.err.println ("Usage: java edu.rit.crypto.blockcipher.AES256Cipher <key> <block>");
//				System.exit (1);
//				}
//			String keystring = args[0];
//			String blockstring = args[1];
//
//			// Convert key string to bytes.
//			byte[] key = new byte[32];
//			for (i = 0; i < 32; ++ i)
//				{
//				key[i] = (byte)
//					Integer.parseInt (keystring.substring (2*i, 2*i+2), 16);
//				}
//
//			// Convert block string to bytes.
//			byte[] block = new byte[16];
//			for (i = 0; i < 16; ++ i)
//				{
//				block[i] = (byte)
//					Integer.parseInt (blockstring.substring (2*i, 2*i+2), 16);
//				}
//
//			// Create AES cipher object.
//			AES256Cipher bc = new AES256Cipher (key);
//
//			// Dump key schedule.
//			System.out.println ("KEY SCHEDULE");
//			for (i = 0; i < ROUND_KEY_COUNT; ++ i)
//				{
//				RoundKey roundkey = bc.expandedKey[i];
//				System.out.print (4*i);
//				System.out.print ('\t');
//				System.out.print (hexdigit [(roundkey.key00 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key00     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key10 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key10     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key20 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key20     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key30 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key30     ) & 0xF]);
//				System.out.println();
//				System.out.print (4*i+1);
//				System.out.print ('\t');
//				System.out.print (hexdigit [(roundkey.key01 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key01     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key11 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key11     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key21 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key21     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key31 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key31     ) & 0xF]);
//				System.out.println();
//				System.out.print (4*i+2);
//				System.out.print ('\t');
//				System.out.print (hexdigit [(roundkey.key02 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key02     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key12 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key12     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key22 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key22     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key32 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key32     ) & 0xF]);
//				System.out.println();
//				System.out.print (4*i+3);
//				System.out.print ('\t');
//				System.out.print (hexdigit [(roundkey.key03 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key03     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key13 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key13     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key23 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key23     ) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key33 >> 4) & 0xF]);
//				System.out.print (hexdigit [(roundkey.key33     ) & 0xF]);
//				System.out.println();
//				}
//
//			// Dump the plaintext.
//			System.out.println();
//			System.out.println ("PLAINTEXT");
//			for (i = 0; i < 16; ++ i)
//				{
//				System.out.print (hexdigit [(block[i] >> 4) & 0xF]);
//				System.out.print (hexdigit [(block[i]     ) & 0xF]);
//				}
//			System.out.println();
//
//			// Encrypt the block.
//			bc.encrypt (block);
//
//			// Dump the ciphertext.
//			System.out.println();
//			System.out.println ("CIPHERTEXT");
//			for (i = 0; i < 16; ++ i)
//				{
//				System.out.print (hexdigit [(block[i] >> 4) & 0xF]);
//				System.out.print (hexdigit [(block[i]     ) & 0xF]);
//				}
//			System.out.println();
//			}
//
//		catch (Throwable exc)
//			{
//			System.err.println ("AES256Cipher: Uncaught exception");
//			exc.printStackTrace (System.err);
//			System.exit (1);
//			}
//		}
//
//	private static final char[] hexdigit = new char[]
//		{'0', '1', '2', '3', '4', '5', '6', '7',
//		 '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

	}
