//******************************************************************************
//
// File:    LongIterator.java
// Package: edu.rit.hyb.prime
// Unit:    Interface edu.rit.hyb.prime.LongIterator
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
 * Interface LongIterator specifies the interface for an iterator over a
 * sequence of positive numbers of type <TT>long</TT>.
 *
 * @author  Alan Kaminsky
 * @version 03-Jun-2008
 */
public interface LongIterator
	{

// Exported operations.

	/**
	 * Returns the next number in the sequence.
	 *
	 * @return  Number, or 0 if there are no more numbers.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public long next()
		throws IOException;

	/**
	 * Close this iterator. Call <TT>close()</TT> when done using this iterator,
	 * to release resources.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void close()
		throws IOException;

	}
