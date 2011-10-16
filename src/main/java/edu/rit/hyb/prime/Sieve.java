//******************************************************************************
//
// File:    Sieve.java
// Package: edu.rit.hyb.prime
// Unit:    Class edu.rit.hyb.prime.Sieve
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

import java.io.IOException;

/**
 * Class Sieve provides an object for finding prime numbers using the Sieve of
 * Eratosthenes. A sieve is an array of Boolean flags, of a certain length,
 * starting at a certain lower bound index. The flag at index <I>p</I> is true
 * if <I>p</I> is a prime number and is false otherwise.
 *
 * @author  Alan Kaminsky
 * @version 05-Jun-2008
 */
public class Sieve
	{

// Hidden data members.

	private boolean[] isPrime;
	private long lb;
	private int len;

	// Padding to avert cache interference.
	private long p0, p1, p2, p3, p4, p5, p6, p7;
	private long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.

	/**
	 * Construct a new sieve object.
	 *
	 * @param  lb   Lower bound index. Assumed to be a nonnegative even number.
	 * @param  len  Length. Assumed to be a nonnegative even number.
	 */
	public Sieve
		(long lb,
		 int len)
		{
		this.isPrime = new boolean [len+128]; // Padding
		this.lb = lb;
		this.len = len;
		}

// Exported operations.

	/**
	 * Get this sieve's lower bound index.
	 *
	 * @return  Lower bound index.
	 */
	public long lb()
		{
		return this.lb;
		}

	/**
	 * Set this sieve's lower bound index.
	 *
	 * @param  lb  Lower bound index. Assumed to be a nonnegative even number.
	 */
	public void lb
		(long lb)
		{
		this.lb = lb;
		}

	/**
	 * Get this sieve's length.
	 *
	 * @return  Length.
	 */
	public long length()
		{
		return this.len;
		}

	/**
	 * Set this sieve's length.
	 *
	 * @param  len  Length. Assumed to be a nonnegative even number.
	 */
	public void len
		(int len)
		{
		this.isPrime = new boolean [len+128]; // Padding
		this.len = len;
		}

	/**
	 * Initialize this sieve. Afterwards, all even-numbered flags are false and
	 * all odd-numbered flags are true.
	 */
	public void initialize()
		{
		for (int i = 1; i < len; i += 2) isPrime[i] = true;
		if (lb == 0) isPrime[1] = false; // 1 is not a prime
		}

	/**
	 * Sieve out the given prime. Afterwards, all flags corresponding to
	 * multiples of the given prime are false. It is assumed that <TT>p</TT> is
	 * an odd prime and that all multiples of smaller primes have previously
	 * been sieved out.
	 *
	 * @param  p  Prime to sieve out.
	 *
	 * @return  True if sieving with further primes is required, false
	 *          otherwise.
	 */
	public boolean sieveOut
		(long p)
		{
		// If p^2 is beyond the end of the array, report that further sieving is
		// not required.
		long psqr = p*p;
		if (psqr - lb >= len) return false;

		// Find the first odd multiple of p greater than or equal to lb.
		long m = (lb + p - 1) / p;
		if ((m & 1) == 0) ++ m;
		long mp = m*p;

		// Sieving begins at p^2 or mp, whichever is larger.
		mp = Math.max (mp, psqr);

		// Set all odd multiples of p to false and report that further sieving
		// is required.
		long two_p = 2*p;
		for (long i = mp - lb; i < len; i += two_p)
			{
			isPrime[(int) i] = false;
			}
		return true;
		}

	/**
	 * Sieve out all primes returned by the given iterator. It is assumed that
	 * the <TT>iterator</TT> returns a sequence of the odd primes in ascending
	 * order (3, 5, 7, 11, . . .). Sieving continues until no further sieving is
	 * required or until the end of the iterator's sequence of primes, whichever
	 * comes first.
	 *
	 * @param  iterator  Iterator for a sequence of odd primes.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void sieveOut
		(LongIterator iterator)
		throws IOException
		{
		long p;
		initialize();
		while ((p = iterator.next()) != 0 && sieveOut (p));
		}

	/**
	 * Determine if the given number is prime. It is assumed that all primes
	 * smaller than <TT>p</TT> have been sieved out.
	 *
	 * @param  p  Number to test.
	 *
	 * @return  True if <TT>p</TT> is prime, false otherwise.
	 */
	public boolean isPrime
		(long p)
		{
		return isPrime[(int)(p - lb)];
		}

	/**
	 * Obtain an iterator for the primes in this sieve. The iterator returns a
	 * sequence of the numbers whose flags are true in this sieve.
	 *
	 * @return  Iterator.
	 */
	public LongIterator iterator()
		{
		return new LongIterator()
			{
			private int i = 0;

			// Padding to avert cache interference.
			private long p0, p1, p2, p3, p4, p5, p6, p7;
			private long p8, p9, pa, pb, pc, pd, pe, pf;

			public long next()
				{
				do ++ i; while (i < len && ! isPrime[i]);
				return i < len ? lb + i : 0;
				}

			public void close()
				{
				}
			};
		}

	}
