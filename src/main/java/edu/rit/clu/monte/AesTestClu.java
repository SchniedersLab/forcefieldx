//******************************************************************************
//
// File:    AesTestClu.java
// Package: edu.rit.clu.monte
// Unit:    Class edu.rit.clu.monte.AesTestClu
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

import edu.rit.mp.DoubleBuf;
import edu.rit.mp.IntegerBuf;
import edu.rit.mp.LongBuf;

import edu.rit.mp.buf.DoubleItemBuf;
import edu.rit.mp.buf.LongItemBuf;

import edu.rit.numeric.Statistics;

import edu.rit.pj.Comm;

import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.LongOp;

import edu.rit.util.Hex;
import edu.rit.util.LongRange;
import edu.rit.util.Range;

import java.util.Arrays;

/**
 * Class AesTestClu is a cluster parallel program that tests the randomness of
 * the AES block cipher.
 * <P>
 * Usage: java -Dpj.np=<I>K</I> edu.rit.clu.monte.AesTestClu <I>key</I> <I>N</I>
 * <BR><I>K</I> = Number of parallel processes
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
 * The computation is performed in parallel in multiple processors. The program
 * measures the computation's running time.
 *
 * @author  Alan Kaminsky
 * @version 05-May-2008
 */
public class AesTestClu
	{

// Prevent construction.

	private AesTestClu()
		{
		}

// Program shared variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static byte[] key = new byte [32];
	static long N;

	// AES block cipher.
	static AES256Cipher cipher;

	// Plaintext and ciphertext blocks.
	static byte[] plaintext = new byte [16];
	static byte[] ciphertext = new byte [16];

	// Random data values, partitioned to be sent to all processes.
	static double[] sendData;

	// Number of data values sent from this process to each process, plus total.
	static int[] sendLength;
	static int sendN;

	// Number of data values received by this process from each process, plus
	// total.
	static int[] recvLength;
	static int recvN;

	// Index ranges in the sendData array from which to obtain data values sent
	// to each process.
	static Range[] sendRanges;

	// Index ranges in the data array in which to store data values received
	// from each process.
	static Range[] recvRanges;

	// Random data values received by this process.
	static double[] data;

	// Number of data values in lower-ranked processes.
	static long lowerN;

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
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Validate command line arguments.
		if (args.length != 2) usage();
		Hex.toByteArray (args[0], key);
		N = Long.parseLong (args[1]);

		// Set up AES block cipher.
		cipher = new AES256Cipher (key);

		// Compute 2^64.
		TWO_SUP_64  = 2.0;        // 2^1
		TWO_SUP_64 *= TWO_SUP_64; // 2^2
		TWO_SUP_64 *= TWO_SUP_64; // 2^4
		TWO_SUP_64 *= TWO_SUP_64; // 2^8
		TWO_SUP_64 *= TWO_SUP_64; // 2^16
		TWO_SUP_64 *= TWO_SUP_64; // 2^32
		TWO_SUP_64 *= TWO_SUP_64; // 2^64

		// Generate this process's subset of the N random data values.
		LongRange indexRange = new LongRange (0, N-1) .subrange (size, rank);
		long lb = indexRange.lb();
		long len = indexRange.length();
		sendData = new double [(int) len];
		for (long i = 0; i < len; ++ i)
			{
			longToBytes (lb+i, plaintext, 8);
			cipher.encrypt (plaintext, ciphertext);
			sendData[(int) i] = bytesToDouble (ciphertext, 0);
			}

		// If there's more than one process, do message passing.
		if (size > 1)
			{
			// Determine how many data values will be going to each process.
			Arrays.sort (sendData);
			sendLength = new int [size];
			int prevj = 0;
			int j = 0;
			for (int i = 0; i < size; ++ i)
				{
				double threshold = ((double) (i+1)) / ((double) size);
				while (j < len && sendData[j] < threshold) ++ j;
				sendLength[i] = j - prevj;
				prevj = j;
				}

			// Determine how many data values will be coming from each process.
			recvLength = new int [size];
			world.allToAll
				(IntegerBuf.sliceBuffers
					(sendLength, new Range (0, size-1) .subranges (size)),
				 IntegerBuf.sliceBuffers
					(recvLength, new Range (0, size-1) .subranges (size)));

			// Transfer data values.
			sendRanges = new Range [size];
			sendN = 0;
			recvRanges = new Range [size];
			recvN = 0;
			for (int i = 0; i < size; ++ i)
				{
				sendRanges[i] = new Range (sendN, sendN+sendLength[i]-1);
				sendN += sendLength[i];
				recvRanges[i] = new Range (recvN, recvN+recvLength[i]-1);
				recvN += recvLength[i];
				}
			data = new double [recvN];
			world.allToAll
				(DoubleBuf.sliceBuffers (sendData, sendRanges),
				 DoubleBuf.sliceBuffers (data, recvRanges));

			// Release storage for sent data values.
			sendData = null;

			// Determine how many data values ended up in lower-ranked
			// processes.
			LongItemBuf lowerNbuf = LongBuf.buffer (recvN);
			world.exclusiveScan (lowerNbuf, LongOp.SUM, 0L);
			lowerN = lowerNbuf.item;
			}

		// If there's only one process, don't bother with message passing.
		else
			{
			data = sendData;
			sendData = null;
			recvN = (int) len;
			lowerN = 0;
			}

		// Compute the K-S statistic, D, for this process's random data values.
		Arrays.sort (data);
		double N_double = N;
		double D = 0.0;
		double F_lower = lowerN / N_double;
		double F_upper;
		double x;
		for (int i = 0; i < recvN; ++ i)
			{
			F_upper = (lowerN+i+1) / N_double;
			x = data[i];
			D = Math.max (D, Math.abs (x - F_lower));
			D = Math.max (D, Math.abs (x - F_upper));
			F_lower = F_upper;
			}

		// Put the maximum of all processes' D values into process 0.
		DoubleItemBuf Dbuf = DoubleBuf.buffer (D);
		world.reduce (0, Dbuf, DoubleOp.MAXIMUM);
		D = Dbuf.item;

		// Compute the p-value, P.
		double P = Statistics.ksPvalue (N, D);

		// Stop timing.
		time += System.currentTimeMillis();

		// Print results in process 0.
		if (rank == 0)
			{
			System.out.println ("N = " + N);
			System.out.println ("D = " + D);
			System.out.println ("P = " + P);
			}
		System.out.println (time + " msec " + rank);
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
		System.err.println ("Usage: java -Dpj.np=<K> edu.rit.clu.monte.AesTestClu <key> <N>");
		System.err.println ("<K> = Number of parallel processes");
		System.err.println ("<key> = Block cipher key");
		System.err.println ("<N> = Number of blocks");
		System.exit (1);
		}

	}
