//******************************************************************************
//
// File:    Prime32List.java
// Package: edu.rit.hyb.prime
// Unit:    Class edu.rit.hyb.prime.Prime32List
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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Class Prime32List encapsulates a list of 32-bit prime numbers. The list of
 * prime numbers is obtained by reading a sequence of values from a file. Each
 * value is an unsigned byte. Each value is the difference between two
 * consecutive odd primes, divided by 2. The file must be created by the
 * {@linkplain Prime32File} program.
 * <P>
 * Class Prime32List provides a method for creating an iterator. The numbers
 * returned by the iterator are the odd primes in ascending order: 3, 5, 7, 11,
 * and so on. The numbers are returned as type <TT>long</TT>. Once the iterator
 * has read all numbers from the prime file, the iterator returns 0 to signal
 * that there are no more numbers.
 *
 * @author  Alan Kaminsky
 * @version 05-Jun-2008
 */
public class Prime32List
	{

// Hidden data members.

	// Prime file.
	private File primeFile;

// Exported constructors.

	/**
	 * Construct a new 32-bit prime list. The prime numbers are contained in the
	 * given file.
	 *
	 * @param  primeFile  Prime file.
	 */
	public Prime32List
		(File primeFile)
		{
		this.primeFile = primeFile;
		}

// Exported operations.

	/**
	 * Obtain an iterator for the primes in this 32-bit prime list.
	 *
	 * @return  Iterator.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public LongIterator iterator()
		throws IOException
		{
		return new Prime32Iterator();
		}

// Hidden helper classes.

	/**
	 * Class Prime32Iterator encapsulates an iterator for the primes in this
	 * 32-bit prime list.
	 *
	 * @author  Alan Kaminsky
	 * @version 03-Jun-2008
	 */
	private class Prime32Iterator
		implements LongIterator
		{
		private long previousP = 2;
		private InputStream primeStream;

		// Padding to avert cache interference.
		private long p0, p1, p2, p3, p4, p5, p6, p7;
		private long p8, p9, pa, pb, pc, pd, pe, pf;

		/**
		 * Construct a new iterator.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public Prime32Iterator()
			throws IOException
			{
			primeStream =
				new BufferedInputStream (new FileInputStream (primeFile));
			}

		/**
		 * Returns the next number in the sequence.
		 *
		 * @return  Number, or 0 if there are no more numbers.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public long next()
			throws IOException
			{
			long currentP;
			if (previousP == 0)
				{
				currentP = 0;
				}
			else if (previousP == 2)
				{
				currentP = 3;
				}
			else
				{
				int diff = primeStream.read();
				if (diff == -1)
					{
					currentP = 0;
					}
				else
					{
					currentP = previousP + (diff << 1);
					}
				}
			previousP = currentP;
			return currentP;
			}

		/**
		 * Close this iterator. Call <TT>close()</TT> when done using this
		 * iterator to release resources.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void close()
			throws IOException
			{
			primeStream.close();
			}

		/**
		 * Finalize this iterator.
		 */
		protected void finalize()
			{
			try { close(); } catch (IOException exc) {}
			}

		}

// Unit test main program.

//	/**
//	 * Unit test main program. Prints all primes from the given prime file.
//	 * <P>
//	 * Usage: java edu.rit.hyb.prime.Prime32List <I>primefile</I>
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		if (args.length != 1)
//			{
//			System.err.println ("Usage: java edu.rit.hyb.prime.Prime32List <primefile>");
//			System.exit (1);
//			}
//		Prime32List list = new Prime32List (new File (args[0]));
//		LongIterator iter = list.iterator();
//		long p;
//		while ((p = iter.next()) != 0)
//			{
//			System.out.println (p);
//			}
//		iter.close();
//		}

	}
