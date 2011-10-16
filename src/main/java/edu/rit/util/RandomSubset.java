//******************************************************************************
//
// File:    RandomSubset.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.RandomSubset
//
// This Java source file is copyright (C) 2011 by Alan Kaminsky. All rights
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

package edu.rit.util;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

/**
 * Class RandomSubset provides an object that generates a random subset of a set
 * of integers. The original set consists of the integers from 0 to
 * <I>N</I>&minus;1 inclusive. The subset consists of integers chosen at random
 * without replacement from the original set; each member of the original set is
 * equally likely to be chosen. Class RandomSubset is an Iterator that visits
 * the elements of the subset; each <TT>next()</TT> method call returns another
 * element of the subset. The order in which elements are visited is not
 * specified.
 * <P>
 * Calling the <TT>remove(i)</TT> method one or more times removes the given
 * integers from the original set. Those integers will not be chosen for the
 * random subset.
 * <P>
 * Class RandomSubset is layered on top of a pseudorandom number generator
 * (PRNG), an instance of class {@linkplain edu.rit.util.Random Random}. Each
 * time the <TT>next()</TT> method is called, one random number is consumed from
 * the underlying PRNG.
 * <P>
 * Class RandomSubset's implementation is optimized for the case where the size
 * of the random subset (the number of <TT>next()</TT> method calls) is much
 * less than <I>N</I>, the size of the original set.
 *
 * @author  Alan Kaminsky
 * @version 28-Jun-2011
 */
public class RandomSubset
	implements Iterator<Integer>
	{

// Hidden data members.

	// The underlying PRNG.
	private Random prng;

	// The size of the original set.
	private int N;

	// The number of random subset elements returned so far.
	private int M;

	// A sparse array containing a permutation of the integers from 0 to N-1.
	// Implemented as a mapping from array index to array element. If an array
	// index is not in the map, the corresponding array element is the same as
	// the array index.
	private HashMap<Integer,Integer> permutation =
		new HashMap<Integer,Integer>();

// Hidden operations.

	/**
	 * Returns the element in the permutation array at index i.
	 */
	private int getElement
		(int i)
		{
		Integer element = permutation.get (i);
		return element == null ? i : element;
		}

	/**
	 * Sets the element in the permutation array at index i to the given value.
	 */
	private void setElement
		(int i,
		 int value)
		{
		if (value == i)
			{
			permutation.remove (i);
			}
		else
			{
			permutation.put (i, value);
			}
		}

	/**
	 * Swaps the elements in the permutation array at indexes i and j.
	 */
	private void swapElements
		(int i,
		 int j)
		{
		int tmp = getElement (i);
		setElement (i, getElement (j));
		setElement (j, tmp);
		}

	/**
	 * Returns the index in the permutation array at which the given value
	 * resides.
	 */
	private int indexOf
		(int value)
		{
		for (Map.Entry<Integer,Integer> entry : permutation.entrySet())
			{
			if (entry.getValue() == value) return entry.getKey();
			}
		return value;
		}

// Exported constructors.

	/**
	 * Construct a new random subset object for the original set consisting of
	 * the integers from 0 through <I>N</I>&minus;1 inclusive.
	 *
	 * @param  prng  Underlying PRNG.
	 * @param  N     Size of original set.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>prng</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>N</I> &lt; 0.
	 */
	public RandomSubset
		(Random prng,
		 int N)
		{
		if (prng == null)
			{
			throw new NullPointerException
				("RandomSubset(): prng is null");
			}
		if (N < 0)
			{
			throw new IllegalArgumentException
				("RandomSubset(): N = "+N+" illegal");
			}
		this.prng = prng;
		this.N = N;
		}

// Exported operations.

	/**
	 * Determine whether there are more integers in the random subset.
	 *
	 * @return  True if there are more integers in the random subset, false if
	 *          all the integers in the original set have been used up.
	 */
	public boolean hasNext()
		{
		return M < N;
		}

	/**
	 * Returns the next integer in the random subset.
	 *
	 * @return  Next integer in the random subset.
	 *
	 * @exception  NoSuchElementException
	 *     (unchecked exception) Thrown if all the integers in the original set
	 *     have been used up.
	 */
	public Integer next()
		{
		if (M >= N)
			{
			throw new NoSuchElementException
				("RandomSubset.next(): No further elements");
			}
		swapElements (M, M + prng.nextInt (N - M));
		++ M;
		return getElement (M - 1);
		}

	/**
	 * Unsupported operation.
	 *
	 * @exception  UnsupportedOperationException
	 *     (unchecked exception) Thrown always.
	 */
	public void remove()
		{
		throw new UnsupportedOperationException();
		}

	/**
	 * Remove the given integer from the original set.
	 * <P>
	 * If <TT>i</TT> has already been removed from the original set, either by a
	 * <TT>remove(i)</TT> method call, or by a <TT>next()</TT> method call that
	 * returned <TT>i</TT>, then this method does nothing.
	 *
	 * @param  i  Integer to remove.
	 *
	 * @return  This random subset object.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>i</TT> is not in the range 0
	 *     through <I>N</I>&minus;1 inclusive.
	 */
	public RandomSubset remove
		(int i)
		{
		if (0 > i || i >= N)
			{
			throw new IllegalArgumentException
				("RandomSubset.remove(): i = "+i+" illegal");
			}
		int j = indexOf (i);
		if (j >= M)
			{
			swapElements (M, j);
			++ M;
			}
		return this;
		}

// Unit test main program.

//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.util.RandomSubset <I>seed</I> <I>N</I> <I>M</I> [
//	 * <I>i</I> ... ]
//	 * <BR><I>seed</I> = Random seed
//	 * <BR><I>N</I> = Size of original set
//	 * <BR><I>M</I> = Size of random subset
//	 * <BR><I>i</I> = Integer(s) to remove from original set
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length < 3) usage();
//		long seed = Long.parseLong (args[0]);
//		int N = Integer.parseInt (args[1]);
//		int M = Integer.parseInt (args[2]);
//		RandomSubset rs = new RandomSubset (Random.getInstance (seed), N);
//		for (int j = 3; j < args.length; ++ j)
//			{
//			rs.remove (Integer.parseInt (args[j]));
//			}
//		for (int j = 0; j < M; ++ j)
//			{
//			System.out.printf ("%d  ", rs.next());
//			}
//		System.out.println();
//		}
//
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.util.RandomSubset <seed> <N> <M> [<i> ...]");
//		System.err.println ("<seed> = Random seed");
//		System.err.println ("<N> = Size of original set");
//		System.err.println ("<M> = Size of random subset");
//		System.err.println ("<i> = Integer(s) to remove from original set");
//		System.exit (1);
//		}

	}
