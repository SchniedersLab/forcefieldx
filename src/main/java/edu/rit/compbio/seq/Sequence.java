//******************************************************************************
//
// File:    Sequence.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.Sequence
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

package edu.rit.compbio.seq;

import java.util.Arrays;

/**
 * Class Sequence is the abstract base class for a biological sequence.
 * <P>
 * In a program, a sequence is represented as a byte array (type
 * <TT>byte[]</TT>). For a sequence of length <I>L</I>, the byte array contains
 * <I>L</I>+1 bytes. The byte at index 0 is unused and contains a value of -1.
 * The bytes at indexes 1 through <I>L</I> contain the elements of the sequence.
 * These are generally represented as small integers starting at 0, not as
 * characters.
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public abstract class Sequence
	{

// Hidden data members.

	String myDescription;
	byte[] mySequence;
	int myLength;

// Exported constructors.

	/**
	 * Construct a new sequence.
	 */
	public Sequence()
		{
		}

// Exported operations.

	/**
	 * Get this sequence's description.
	 *
	 * @return  Description string.
	 */
	public String description()
		{
		return myDescription;
		}

	/**
	 * Get this sequence's length <I>L</I>. This is the number of elements in
	 * this sequence.
	 *
	 * @return  Length <I>L</I>.
	 */
	public int length()
		{
		return myLength;
		}

	/**
	 * Get this sequence's elements. The return value is a byte array of length
	 * <I>L</I>+1. The byte at index 0 is unused and contains a value of -1. The
	 * bytes at indexes 1 through <I>L</I> contain the elements.
	 * <P>
	 * <I>Note:</I> Do not alter the contents of the returned byte array.
	 *
	 * @return  Array of elements.
	 */
	public byte[] sequence()
		{
		return mySequence;
		}

	/**
	 * Determine if this sequence is equal to the given object. Two sequences
	 * are equal if they have the same elements.
	 *
	 * @param  obj  Object to test.
	 *
	 * @return  True if <TT>obj</TT> is equal to this sequence, false otherwise.
	 */
	public boolean equals
		(Object obj)
		{
		return
			(obj instanceof Sequence) &&
			Arrays.equals (this.mySequence, ((Sequence) obj).mySequence);
		}

	/**
	 * Returns a hash code for this sequence.
	 *
	 * @return  Hash code.
	 */
	public int hashCode()
		{
		return Arrays.hashCode (mySequence);
		}

	/**
	 * Returns a character version of this sequence's element at the given
	 * index.
	 *
	 * @param  i  Index in the range 1 .. <I>L</I>.
	 *
	 * @return  Character corresponding to element <TT>i</TT>.
	 */
	public abstract char charAt
		(int i);

	/**
	 * Returns a string version of this sequence's elements.
	 *
	 * @return  String version.
	 */
	public String elementsToString()
		{
		char[] c = new char [myLength];
		for (int i = 0; i < myLength; ++ i)
			{
			c[i] = charAt(i+1);
			}
		return new String (c);
		}

	}
