//******************************************************************************
//
// File:    Prime32File.java
// Package: edu.rit.hyb.prime
// Unit:    Class edu.rit.hyb.prime.Prime32File
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

package edu.rit.hyb.prime;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * Class Prime32File is a sequential program that stores a list of 32-bit prime
 * numbers in a file.
 * <P>
 * Usage: java edu.rit.hyb.prime.Prime32File <I>primefile</I>
 * <BR><I>primefile</I> = Output prime file
 * <P>
 * The program calculates all odd primes less than 2<SUP>32</SUP>. The program
 * stores a series of values in the output prime file. Each value is stored as
 * an unsigned byte. Each value is the difference between two consecutive odd
 * primes, divided by 2.
 * <P>
 * The odd primes are 3, 5, 7, 11, 13, 17, 19, 23, . . .
 * <P>
 * The differences between consecutive odd primes are 2, 2, 4, 2, 4, 2, 4, . . .
 * <P>
 * The bytes stored in the output prime file are the differences divided by 2,
 * namely 1, 1, 2, 1, 2, 1, 2, . . .
 *
 * @author  Alan Kaminsky
 * @version 03-Jun-2008
 */
public class Prime32File
	{

// Prevent construction.

	private Prime32File()
		{
		}

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length != 1) usage();
		File primefile = new File (args[0]);

		// Open output file.
		OutputStream out =
			new BufferedOutputStream (new FileOutputStream (primefile));

		// Find all odd primes less than 2^16.
		Sieve smallPrime = new Sieve (0, 65536);
		smallPrime.initialize();
		long p = 3;
		while (smallPrime.sieveOut (p))
			{
			do p += 2; while (! smallPrime.isPrime (p));
			}

		// Write differences of all odd primes less than 2^16 to output file.
		LongIterator iter = smallPrime.iterator();
		iter.next(); // First number returned is 3, omit.
		long prev_p = writePrimeDiffs (out, iter, 3);

		// Find all odd primes less than 2^32, in chunks of 65,536.
		Sieve largePrime = new Sieve (0, 65536);
		for (long lb = 65536; lb < 4294967296L; lb += 65536)
			{
			// Sieve multiples of small odd primes out of the chunk.
			largePrime.lb (lb);
			largePrime.sieveOut (smallPrime.iterator());

			// Write differences of odd primes in chunk to output file.
			prev_p = writePrimeDiffs (out, largePrime.iterator(), prev_p);
			}

		// Close output file.
		out.close();
		}

// Hidden operations.

	/**
	 * Write the differences between consecutive primes, divided by 2, to the
	 * given output stream.
	 *
	 * @param  out     Output stream.
	 * @param  iter    Iterator for primes to write.
	 * @param  prev_p  Previous prime written.
	 *
	 * @return  Last prime written.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private static long writePrimeDiffs
		(OutputStream out,
		 LongIterator iter,
		 long prev_p)
		throws IOException
		{
		long curr_p;
		while ((curr_p = iter.next()) != 0)
			{
			int diff_over_2 = (int) ((curr_p - prev_p) >> 1);
			out.write (diff_over_2);
			prev_p = curr_p;
			}
		return prev_p;
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.hyb.prime.Prime32File <primefile>");
		System.err.println ("<primefile> = Output prime file");
		System.exit (1);
		}

	}
