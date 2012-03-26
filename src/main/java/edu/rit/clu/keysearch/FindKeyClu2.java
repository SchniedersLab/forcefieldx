//******************************************************************************
//
// File:    FindKeyClu2.java
// Package: edu.rit.clu.keysearch
// Unit:    Class edu.rit.clu.keysearch.FindKeyClu2
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package edu.rit.clu.keysearch;

import edu.rit.crypto.blockcipher.AES256Cipher;

import edu.rit.pj.Comm;
import edu.rit.pj.WorkerLongForLoop;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;

import edu.rit.pj.reduction.BooleanOp;

import edu.rit.pj.replica.ReplicatedBoolean;

import edu.rit.util.Hex;

/**
 * Class FindKeyClu2 is a cluster parallel program that solves an AES partial
 * key search problem. The program's command line arguments are the plaintext
 * (128-bit hexadecimal number), the ciphertext (128-bit hexadecimal number),
 * the partial key with the <I>n</I> least significant bits set to 0 (256-bit
 * hexadecimal number), and <I>n</I>, the number of key bits to search for. The
 * ciphertext was created by encrypting the plaintext with the key; however, not
 * all bits of the key are provided. The problem is to find the complete key.
 * The program performs an exhaustive search, trying all possible values for the
 * missing key bits until it finds the key that reproduces the given ciphertext
 * from the given plaintext.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.keysearch.FindKeyClu2
 * <I>plaintext</I> <I>ciphertext</I> <I>partialkey</I> <I>n</I>
 * <BR><I>K</I> = Number of parallel processes
 * <BR><I>plaintext</I> = Plaintext (128-bit hexadecimal number)
 * <BR><I>ciphertext</I> = Ciphertext (128-bit hexadecimal number)
 * <BR><I>partialkey</I> = Partial key (256-bit hexadecimal number)
 * <BR><I>n</I> = Number of key bits to search for
 * <P>
 * Whereas class {@linkplain FindKeyClu} always tests all possible keys, class
 * FindKeyClu2 stops as soon as it finds the correct key.
 *
 * @author  Alan Kaminsky
 * @version 18-Jan-2012
 */
public class FindKeyClu2
	{

// Prevent construction.

	private FindKeyClu2()
		{
		}

// Shared variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static byte[] plaintext;
	static byte[] ciphertext;
	static byte[] partialkey;
	static int n;

	// Variables for doing trial encryptions.
	static long keylsbs;
	static long maxcounter;
	static byte[] foundkey;
	static byte[] trialkey;
	static byte[] trialciphertext;
	static AES256Cipher cipher;

	// For early loop exit.
	static ReplicatedBoolean found;

// Main program.

	/**
	 * AES partial key search main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize PJ middleware.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (args.length != 4) usage();
		plaintext = Hex.toByteArray (args[0]);
		ciphertext = Hex.toByteArray (args[1]);
		partialkey = Hex.toByteArray (args[2]);
		n = Integer.parseInt (args[3]);

		// Make sure n is not too small or too large.
		if (n < 0)
			{
			System.err.println ("n = " + n + " is too small");
			System.exit (1);
			}
		if (n > 63)
			{
			System.err.println ("n = " + n + " is too large");
			System.exit (1);
			}

		// Set up variables for doing trial encryptions.
		keylsbs =
			((partialkey[24] & 0xFFL) << 56) |
			((partialkey[25] & 0xFFL) << 48) |
			((partialkey[26] & 0xFFL) << 40) |
			((partialkey[27] & 0xFFL) << 32) |
			((partialkey[28] & 0xFFL) << 24) |
			((partialkey[29] & 0xFFL) << 16) |
			((partialkey[30] & 0xFFL) <<  8) |
			((partialkey[31] & 0xFFL)      );
		maxcounter = (1L << n) - 1L;
		trialkey = new byte [32];
		System.arraycopy (partialkey, 0, trialkey, 0, 32);
		trialciphertext = new byte [16];
		cipher = new AES256Cipher (trialkey);

		// Set up found flag replicated in all processes.
		found = new ReplicatedBoolean (BooleanOp.OR);

		// Try every possible combination of low-order key bits in parallel.
		new WorkerTeam().execute (new WorkerRegion()
			{
			public void run() throws Exception
				{
				execute (0L, maxcounter, new WorkerLongForLoop()
					{
					public void run (long first, long last) throws Exception
						{
						for (long counter = first; counter <= last; ++ counter)
							{
							// Fill in low-order key bits.
							long lsbs = keylsbs | counter;
							trialkey[24] = (byte) (lsbs >>> 56);
							trialkey[25] = (byte) (lsbs >>> 48);
							trialkey[26] = (byte) (lsbs >>> 40);
							trialkey[27] = (byte) (lsbs >>> 32);
							trialkey[28] = (byte) (lsbs >>> 24);
							trialkey[29] = (byte) (lsbs >>> 16);
							trialkey[30] = (byte) (lsbs >>>  8);
							trialkey[31] = (byte) (lsbs       );

							// Try the key.
							cipher.setKey (trialkey);
							cipher.encrypt (plaintext, trialciphertext);

							// If the result equals the ciphertext, we found the
							// key. Set the found flag in all processes.
							if (match (ciphertext, trialciphertext))
								{
								foundkey = new byte [32];
								System.arraycopy (trialkey, 0, foundkey, 0, 32);
								found.reduce (true);
								}

							// If key was found, exit loop.
							if (found.get()) break;
							}
						}
					});
				}
			});

		// If we found the key, print it.
		if (foundkey != null)
			{
			System.out.println (Hex.toString (foundkey));
			}

		// Stop timing.
		long t2 = System.currentTimeMillis();
		System.out.println ((t2-t1) + " msec (" + rank + ")");
		}

// Hidden operations.

	/**
	 * Returns true if the two byte arrays match.
	 */
	private static boolean match
		(byte[] a,
		 byte[] b)
		{
		boolean matchsofar = true;
		int n = a.length;
		for (int i = 0; i < n; ++ i)
			{
			matchsofar = matchsofar && a[i] == b[i];
			}
		return matchsofar;
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.keysearch.FindKeyClu2 <plaintext> <ciphertext> <partialkey> <n>");
		System.err.println ("<K> = Number of parallel processes");
		System.err.println ("<plaintext> = Plaintext (128-bit hexadecimal number)");
		System.err.println ("<ciphertext> = Ciphertext (128-bit hexadecimal number)");
		System.err.println ("<partialkey> = Partial key (256-bit hexadecimal number)");
		System.err.println ("<n> = Number of key bits to search for");
		System.exit (1);
		}

	}
