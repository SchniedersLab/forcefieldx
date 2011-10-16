//******************************************************************************
//
// File:    AesTestSeq.java
// Package: edu.rit.clu.monte
// Unit:    Class edu.rit.clu.monte.AesTestSeq
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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

package edu.rit.clu.monte;

import edu.rit.crypto.blockcipher.AES256Cipher;

import edu.rit.numeric.Statistics;

import edu.rit.pj.Comm;

import edu.rit.util.Hex;

import java.util.Arrays;

/**
 * Class AesTestSeq is a sequential program that tests the randomness of the AES
 * block cipher.
 * <P>
 * Usage: java edu.rit.clu.monte.AesTestSeq <I>key</I> <I>N</I>
 * <BR><I>key</I> = Block cipher key
 * <BR><I>N</I> = Number of blocks
 * <P>
 * The program:
 * <UL>
 * <LI>
 * Initializes the AES block cipher with the <I>key</I>.
 * <LI>
 * Generates <I>N</I> ciphertext blocks by encrypting the plaintext blocks 0, 1,
 * 2, and so on (counter mode).
 * <LI>
 * Takes the most significant halves of the 128-bit ciphertext blocks as a
 * series of 64-bit <TT>long</TT> values.
 * <LI>
 * Divides the <TT>long</TT> values by 2<SUP>64</SUP>, yielding a series of
 * supposedly-random <TT>double</TT> values in the range 0.0 to 1.0.
 * <LI>
 * Performs a Kolmogorov-Smirnov (K-S) test on the data to try to disprove the
 * null hypothesis that the data are drawn from a uniform distribution.
 * <LI>
 * Prints the K-S statistic and the p-value.
 * </UL>
 * <P>
 * The computation is performed sequentially in a single processor. The program
 * measures the computation's running time. This establishes a benchmark for
 * measuring the computation's running time on a parallel processor.
 *
 * @author  Alan Kaminsky
 * @version 03-May-2008
 */
public class AesTestSeq
	{

// Prevent construction.

	private AesTestSeq()
		{
		}

// Program shared variables.

	// Command line arguments.
	static byte[] key = new byte [32];
	static int N;

	// AES block cipher.
	static AES256Cipher cipher;

	// Plaintext and ciphertext blocks.
	static byte[] plaintext = new byte [16];
	static byte[] ciphertext = new byte [16];

	// Random data values.
	static double[] data;

	// 2^64.
	static double TWO_SUP_64;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long time = -System.currentTimeMillis();

		// Initialize middleware.
		Comm.init (args);

		// Validate command line arguments.
		if (args.length != 2) usage();
		Hex.toByteArray (args[0], key);
		N = Integer.parseInt (args[1]);

		// Set up AES block cipher.
		cipher = new AES256Cipher (key);

		// Allocate storage for random data values.
		data = new double [N];

		// Compute 2^64.
		TWO_SUP_64  = 2.0;        // 2^1
		TWO_SUP_64 *= TWO_SUP_64; // 2^2
		TWO_SUP_64 *= TWO_SUP_64; // 2^4
		TWO_SUP_64 *= TWO_SUP_64; // 2^8
		TWO_SUP_64 *= TWO_SUP_64; // 2^16
		TWO_SUP_64 *= TWO_SUP_64; // 2^32
		TWO_SUP_64 *= TWO_SUP_64; // 2^64

		// Generate N random data values.
		for (int i = 0; i < N; ++ i)
			{
			longToBytes (i, plaintext, 8);
			cipher.encrypt (plaintext, ciphertext);
			data[i] = bytesToDouble (ciphertext, 0);
			}

		// Compute the K-S statistic, D.
		Arrays.sort (data);
		double N_double = N;
		double D = 0.0;
		double F_lower = 0.0;
		double F_upper;
		double x;
		for (int i = 0; i < N; ++ i)
			{
			F_upper = (i+1) / N_double;
			x = data[i];
			D = Math.max (D, Math.abs (x - F_lower));
			D = Math.max (D, Math.abs (x - F_upper));
			F_lower = F_upper;
			}

		// Compute the p-value, P.
		double P = Statistics.ksPvalue (N, D);

		// Stop timing.
		time += System.currentTimeMillis();

		// Print results.
		System.out.println ("N = " + N);
		System.out.println ("D = " + D);
		System.out.println ("P = " + P);
		System.out.println (time + " msec");
		}

// Hidden operations.

	/**
	 * Convert the given <TT>long</TT> value to eight bytes stored starting at
	 * <TT>block[i]</TT>.
	 *
	 * @param  value  <TT>long</TT> value.
	 * @param  block  Plaintext block.
	 * @param  i      Starting index.
	 */
	private static void longToBytes
		(long value,
		 byte[] block,
		 int i)
		{
		for (int j = 7; j >= 0; -- j)
			{
			block[i+j] = (byte) (value & 0xFF);
			value >>>= 8;
			}
		}

	/**
	 * Convert the eight bytes starting at <TT>block[i]</TT> to a
	 * <TT>double</TT> value.
	 *
	 * @param  block  Ciphertext block.
	 * @param  i      Starting index.
	 *
	 * @return  <TT>double</TT> value.
	 */
	private static double bytesToDouble
		(byte[] block,
		 int i)
		{
		long result = 0L;
		for (int j = 0; j < 8; ++ j)
			{
			result = (result << 8) | (block[i+j] & 0xFF);
			}
		return result / TWO_SUP_64 + 0.5;
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.clu.monte.AesTestSeq <key> <N>");
		System.err.println ("<key> = Block cipher key");
		System.err.println ("<N> = Number of blocks");
		System.exit (1);
		}

	}
