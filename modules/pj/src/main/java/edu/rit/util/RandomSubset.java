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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
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
 * the elements of the original set in a random order; each <code>next()</code>
 * method call returns another element of the subset.
 * <P>
 * Calling the <code>remove(i)</code> method one or more times removes the given
 * integers from the original set. Those integers will not be chosen for the
 * random subset.
 * <P>
 * Class RandomSubset is layered on top of a pseudorandom number generator
 * (PRNG), an instance of class {@linkplain edu.rit.util.Random Random}. Each
 * time the <code>next()</code> method is called, one random number is consumed from
 * the underlying PRNG.
 * <P>
 * Class RandomSubset has two different implementations with different storage
 * requirements. The implementation is specified with a constructor argument.
 * <UL>
 *
 * <LI>
 * The <I>sparse</I> implementation requires less storage when the size of the
 * random subset (the number of <code>next()</code> method calls) is a small
 * fraction of the size of the original set (<I>N</I>).
 *
 * <LI>
 * The <I>dense</I> implementation requires less storage when the size of the
 * random subset is a large fraction of the size of the original set.
 * </UL>
 *
 * @author Alan Kaminsky
 * @version 22-Nov-2011
 */
public class RandomSubset
        implements Iterator<Integer> {

// Hidden data members.
    // The underlying PRNG.
    private Random prng;

    // The size of the original set.
    private int N;

    // The number of random subset elements returned so far.
    private int M;

    // Helper class.
    private Helper helper;

// Hidden helper classes.
    // Helper abstract base class.
    private abstract class Helper {

        // Returns the element in the permutation array at index i.
        public abstract int getElement(int i);

        // Sets the element in the permutation array at index i to the given
        // value.
        public abstract void setElement(int i,
                int value);

        // Swaps the elements in the permutation array at indexes i and j.
        public abstract void swapElements(int i,
                int j);

        // Returns the index in the permutation array at which the given value
        // resides.
        public abstract int indexOf(int value);

        // Restarts the iteration.
        public abstract void restart();
    }

    // Sparse implementation helper class.
    private class SparseHelper
            extends Helper {
        // A sparse array containing a permutation of the integers from 0 to
        // N-1. Implemented as a mapping from array index to array element. If
        // an array index is not in the map, the corresponding array element is
        // the same as the array index.

        private HashMap<Integer, Integer> permutation
                = new HashMap<Integer, Integer>();

        // Returns the element in the permutation array at index i.
        public int getElement(int i) {
            Integer element = permutation.get(i);
            return element == null ? i : element;
        }

        // Sets the element in the permutation array at index i to the given
        // value.
        public void setElement(int i,
                int value) {
            if (value == i) {
                permutation.remove(i);
            } else {
                permutation.put(i, value);
            }
        }

        // Swaps the elements in the permutation array at indexes i and j.
        public void swapElements(int i,
                int j) {
            int tmp = getElement(i);
            setElement(i, getElement(j));
            setElement(j, tmp);
        }

        // Returns the index in the permutation array at which the given value
        // resides.
        public int indexOf(int value) {
            for (Map.Entry<Integer, Integer> entry : permutation.entrySet()) {
                if (entry.getValue() == value) {
                    return entry.getKey();
                }
            }
            return value;
        }

        // Restarts the iteration.
        public void restart() {
            permutation.clear();
        }
    }

    // Dense implementation helper class.
    private class DenseHelper
            extends Helper {
        // A dense array containing a permutation of the integers from 0 to
        // N-1.

        private int[] permutation = new int[N];

        // Construct a new dense helper.
        public DenseHelper() {
            restart();
        }

        // Returns the element in the permutation array at index i.
        public int getElement(int i) {
            return permutation[i];
        }

        // Sets the element in the permutation array at index i to the given
        // value.
        public void setElement(int i,
                int value) {
            permutation[i] = value;
        }

        // Swaps the elements in the permutation array at indexes i and j.
        public void swapElements(int i,
                int j) {
            int tmp = permutation[i];
            permutation[i] = permutation[j];
            permutation[j] = tmp;
        }

        // Returns the index in the permutation array at which the given value
        // resides.
        public int indexOf(int value) {
            int i = 0;
            while (permutation[i] != value) {
                ++i;
            }
            return i;
        }

        // Restarts the iteration.
        public void restart() {
            for (int i = 0; i < N; ++i) {
                permutation[i] = i;
            }
        }
    }

// Exported constructors.
    /**
     * Construct a new random subset object for the original set consisting of
     * the integers from 0 through <I>N</I>&minus;1 inclusive. The sparse
     * implementation is used.
     *
     * @param prng Underlying PRNG.
     * @param N Size of original set.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>prng</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>N</I> &lt; 0.
     */
    public RandomSubset(Random prng,
            int N) {
        this(prng, N, false);
    }

    /**
     * Construct a new random subset object for the original set consisting of
     * the integers from 0 through <I>N</I>&minus;1 inclusive. If <code>dense</code>
     * is false, the sparse implementation is used. If <code>dense</code> is true,
     * the dense implementation is used.
     *
     * @param prng Underlying PRNG.
     * @param N Size of original set.
     * @param dense False to use the sparse implementation, true to use the
     * dense implementation.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>prng</code> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>N</I> &lt; 0.
     */
    public RandomSubset(Random prng,
            int N,
            boolean dense) {
        if (prng == null) {
            throw new NullPointerException("RandomSubset(): prng is null");
        }
        if (N < 0) {
            throw new IllegalArgumentException("RandomSubset(): N = " + N + " illegal");
        }
        this.prng = prng;
        this.N = N;
        this.helper = dense ? new DenseHelper() : new SparseHelper();
    }

// Exported operations.
    /**
     * Determine whether there are more integers in the random subset.
     *
     * @return True if there are more integers in the random subset, false if
     * all the integers in the original set have been used up.
     */
    public boolean hasNext() {
        return M < N;
    }

    /**
     * Returns the next integer in the random subset.
     *
     * @return Next integer in the random subset.
     * @exception NoSuchElementException (unchecked exception) Thrown if all the
     * integers in the original set have been used up.
     */
    public Integer next() {
        if (M >= N) {
            throw new NoSuchElementException("RandomSubset.next(): No further elements");
        }
        helper.swapElements(M, M + prng.nextInt(N - M));
        ++M;
        return helper.getElement(M - 1);
    }

    /**
     * Unsupported operation.
     *
     * @exception UnsupportedOperationException (unchecked exception) Thrown
     * always.
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Remove the given integer from the original set.
     * <P>
     * If <code>i</code> has already been removed from the original set, either by a
     * <code>remove(i)</code> method call, or by a <code>next()</code> method call that
     * returned <code>i</code>, then this method does nothing.
     *
     * @param i Integer to remove.
     * @return This random subset object.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>i</code> is not in the range 0 through <I>N</I>&minus;1 inclusive.
     */
    public RandomSubset remove(int i) {
        if (0 > i || i >= N) {
            throw new IllegalArgumentException("RandomSubset.remove(): i = " + i + " illegal");
        }
        int j = helper.indexOf(i);
        if (j >= M) {
            helper.swapElements(M, j);
            ++M;
        }
        return this;
    }

    /**
     * Restart this random subset's iteration. The original set is reset to
     * contain all the integers from 0 through <I>N</I>&minus;1 inclusive. If
     * integers had been removed from the original set by calling the
     * <code>remove(i)</code> method, those integers must be removed again by
     * calling the <code>remove(i)</code> method again. The iteration is restarted,
     * and calling the <code>next()</code> method will yield a different random
     * subset.
     * <P>
     * The <code>restart()</code> method lets you generate multiple random subsets
     * from the same original set using the same RandomSubset object. This
     * avoids the multiple storage allocations that would take place when
     * creating multiple RandomSubset objects. This in turn can reduce the
     * running time, especially when <I>N</I> is large.
     */
    public void restart() {
        M = 0;
        helper.restart();
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.util.RandomSubset <I>impl</I> <I>seed</I> <I>N</I>
//	 * <I>M</I> [ <I>i</I> ... ]
//	 * <BR><I>impl</I> = "sparse" or "dense"
//	 * <BR><I>seed</I> = Random seed
//	 * <BR><I>N</I> = Size of original set
//	 * <BR><I>M</I> = Size of random subset
//	 * <BR><I>i</I> = Integer(s) to remove from original set
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length < 3) usage();
//		boolean dense = args[0].equals ("dense");
//		long seed = Long.parseLong (args[1]);
//		int N = Integer.parseInt (args[2]);
//		int M = Integer.parseInt (args[3]);
//		RandomSubset rs =
//			new RandomSubset (Random.getInstance (seed), N, dense);
//		for (int r = 0; r < 4; ++ r)
//			{
//			rs.restart();
//			for (int j = 3; j < args.length; ++ j)
//				{
//				rs.remove (Integer.parseInt (args[j]));
//				}
//			for (int j = 0; j < M; ++ j)
//				{
//				System.out.printf ("%d  ", rs.next());
//				}
//			System.out.println();
//			}
//		}
//
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.util.RandomSubset <impl> <seed> <N> <M> [<i> ...]");
//		System.err.println ("<impl> = \"sparse\" or \"dense\"");
//		System.err.println ("<seed> = Random seed");
//		System.err.println ("<N> = Size of original set");
//		System.err.println ("<M> = Size of random subset");
//		System.err.println ("<i> = Integer(s) to remove from original set");
//		System.exit (1);
//		}
}
