//******************************************************************************
//
// File:    Encrypt.java
// Package: edu.rit.smp.keysearch
// Unit:    Class edu.rit.smp.keysearch.Encrypt
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package edu.rit.smp.keysearch;

import edu.rit.crypto.blockcipher.AES256Cipher;

import edu.rit.util.Hex;

/**
 * Class Encrypt prepares input for the AES partial key search programs. The
 * program's command line arguments are a message string to encrypt, the
 * encryption key, and <I>n</I>, the number of key bits to search for. The
 * program prints on the standard output the plaintext (128-bit hexadecimal
 * number), the ciphertext (128-bit hexadecimal number), the partial key with
 * the <I>n</I> least significant bits set to 0 (256-bit hexadecimal number),
 * and <I>n</I>, the number of key bits to search for.
 * <P>
 * Usage: java edu.rit.smp.keysearch.Encrypt <I>message</I> <I>key</I> <I>n</I>
 * <BR><I>message</I> = Message string to encrypt
 * <BR><I>key</I> = Encryption key (256-bit hexadecimal number)
 * <BR><I>n</I> = Number of key bits to search for
 *
 * @author  Alan Kaminsky
 * @version 04-Sep-2006
 */
public class Encrypt
	{

// Prevent construction.

	private Encrypt()
		{
		}

// Main program.

	/**
	 * Encrypt main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 3) usage();
		String message = args[0];
		byte[] key = Hex.toByteArray (args[1]);
		int n = Integer.parseInt (args[2]);

		// Set up plaintext block.
		byte[] msg = message.getBytes();
		byte[] block = new byte [16];
		System.arraycopy (msg, 0, block, 0, Math.min (msg.length, 16));
		System.out.println (Hex.toString (block));

		// Encrypt plaintext.
		AES256Cipher cipher = new AES256Cipher (key);
		cipher.encrypt (block);
		System.out.println (Hex.toString (block));

		// Wipe out n least significant bits of the key.
		int off = 31;
		int len = n;
		while (len >= 8)
			{
			key[off] = (byte) 0;
			-- off;
			len -= 8;
			}
		key[off] &= mask[len];
		System.out.println (Hex.toString (key));
		System.out.println (n);
		}

	private static final int[] mask = new int[]
		{0xff, 0xfe, 0xfc, 0xf8, 0xf0, 0xe0, 0xc0, 0x80};

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.smp.keysearch.Encrypt <message> <key> <n>");
		System.err.println ("<message> = Message string to encrypt");
		System.err.println ("<key> = Encryption key (256-bit hexadecimal number)");
		System.err.println ("<n> = Number of key bits to search for");
		System.exit (0);
		}

	}
