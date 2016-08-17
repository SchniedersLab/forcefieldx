//******************************************************************************
//
// File:    RandomSample.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.RandomSample
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Class RandomSample provides objects that generate random samples from
 * discrete sets.
 *
 * @author Alan Kaminsky
 * @version 15-Feb-2010
 */
public class RandomSample {

// Prevent construction.
    private RandomSample() {
    }

// Exported operations.
    /**
     * Create an iterator that generates a random sample of <TT>int</TT>s
     * without replacement. The set to be sampled consists of the integers from
     * 0 through <I>N</I>&minus;1 inclusive. The sample consists of <I>n</I>
     * items. Each item in the set is equally likely to occur in the sample. The
     * iterator returns the items in the sample in ascending order. The iterator
     * uses the given <TT>prng</TT> to sample items at random. For each item
     * returned, the iterator consumes one random number from the <TT>prng</TT>.
     * <P>
     * As a special case, if <I>n</I> &ge; <I>N</I>, then the iterator returns
     * all items in the set in ascending order, and the iterator does not use
     * the <TT>prng</TT>.
     * <P>
     * The iterator uses the algorithm in A. Bissell, Ordered random selection
     * without replacement, <I>Applied Statistics,</I> 35(1):73-75, 1986.
     *
     * @param prng Pseudorandom number generator.
     * @param n Number of items in the sample.
     * @param N Number of items in the set.
     * @return Iterator over the sample.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>n</I> &lt; 0 or <I>N</I> &lt; 0.
     */
    public static Iterator<Integer> withoutReplacement(final Random prng,
            final int n,
            final int N) {
        if (n < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): n < 0 illegal");
        } else if (N < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): N < 0 illegal");
        } else if (n < N) {
            return new Iterator<Integer>() {
                private int i = 0;
                private int M = N;
                private int r = N - n;

                public boolean hasNext() {
                    return i < n;
                }

                public Integer next() {
                    if (i >= n) {
                        throw new NoSuchElementException();
                    }
                    ++i;
                    double x = prng.nextDouble();
                    double p = 1.0;
                    for (;;) {
                        p = p * r / M;
                        if (p <= x) {
                            int result = N - M;
                            --M;
                            return result;
                        } else {
                            --M;
                            --r;
                        }
                    }
                }

                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };
        } else {
            return new Iterator<Integer>() {
                private int i = 0;

                public boolean hasNext() {
                    return i < N;
                }

                public Integer next() {
                    if (i >= N) {
                        throw new NoSuchElementException();
                    }
                    return i++;
                }

                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };
        }
    }

    /**
     * Generate a random sample of <TT>int</TT>s without replacement. The set to
     * be sampled consists of the integers from 0 through <I>N</I>&minus;1
     * inclusive. The sample consists of <I>n</I> items. Each item in the set is
     * equally likely to occur in the sample. The items in the sample are stored
     * in ascending order in the given <TT>buf</TT> starting at index 0. The
     * iterator uses the given <TT>prng</TT> to sample items at random. For each
     * item sampled, the iterator consumes one random number from the
     * <TT>prng</TT>.
     * <P>
     * As a special case, if <I>n</I> &ge; <I>N</I>, then all items in the set
     * are stored in ascending order in the given <TT>buf</TT> starting at index
     * 0, and the iterator does not use the <TT>prng</TT>.
     * <P>
     * This method uses the algorithm in A. Bissell, Ordered random selection
     * without replacement, <I>Applied Statistics,</I> 35(1):73-75, 1986.
     *
     * @param prng Pseudorandom number generator.
     * @param n Number of items in the sample.
     * @param N Number of items in the set.
     * @param buf Array in which to store the sampled items.
     * @return Number of sampled items actually stored in <TT>buf</TT>.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>n</I> &lt; 0 or <I>N</I> &lt; 0.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>buf.length</TT> &lt; <TT>n</TT>.
     */
    public static int withoutReplacement(Random prng,
            int n,
            int N,
            int[] buf) {
        if (n < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): n < 0 illegal");
        } else if (N < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): N < 0 illegal");
        } else if (n < N) {
            int M = N;
            int r = N - n;
            for (int i = 0; i < n; ++i) {
                double x = prng.nextDouble();
                double p = 1.0;
                for (;;) {
                    p = p * r / M;
                    if (p <= x) {
                        buf[i] = N - M;
                        --M;
                        break;
                    } else {
                        --M;
                        --r;
                    }
                }
            }
            return n;
        } else {
            for (int i = 0; i < N; ++i) {
                buf[i] = i;
            }
            return N;
        }
    }

    /**
     * Create an iterator that generates a random sample of <TT>long</TT>s
     * without replacement. The set to be sampled consists of the integers from
     * 0 through <I>N</I>&minus;1 inclusive. The sample consists of <I>n</I>
     * items. Each item in the set is equally likely to occur in the sample. The
     * iterator returns the items in the sample in ascending order. The iterator
     * uses the given <TT>prng</TT> to sample items at random. For each item
     * returned, the iterator consumes one random number from the <TT>prng</TT>.
     * <P>
     * As a special case, if <I>n</I> &ge; <I>N</I>, then the iterator returns
     * all items in the set in ascending order, and the iterator does not use
     * the <TT>prng</TT>.
     * <P>
     * The iterator uses the algorithm in A. Bissell, Ordered random selection
     * without replacement, <I>Applied Statistics,</I> 35(1):73-75, 1986.
     *
     * @param prng Pseudorandom number generator.
     * @param n Number of items in the sample.
     * @param N Number of items in the set.
     * @return Iterator over the sample.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>n</I> &lt; 0 or <I>N</I> &lt; 0.
     */
    public static Iterator<Long> withoutReplacement(final Random prng,
            final int n,
            final long N) {
        if (n < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): n < 0 illegal");
        } else if (N < 0L) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): N < 0 illegal");
        } else if (n < N) {
            return new Iterator<Long>() {
                private int i = 0;
                private long M = N;
                private long r = N - n;

                public boolean hasNext() {
                    return i < n;
                }

                public Long next() {
                    if (i >= n) {
                        throw new NoSuchElementException();
                    }
                    ++i;
                    double x = prng.nextDouble();
                    double p = 1.0;
                    for (;;) {
                        p = p * r / M;
                        if (p <= x) {
                            long result = N - M;
                            --M;
                            return result;
                        } else {
                            --M;
                            --r;
                        }
                    }
                }

                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };
        } else {
            return new Iterator<Long>() {
                private long i = 0L;

                public boolean hasNext() {
                    return i < N;
                }

                public Long next() {
                    if (i >= N) {
                        throw new NoSuchElementException();
                    }
                    return i++;
                }

                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };
        }
    }

    /**
     * Generate a random sample of <TT>long</TT>s without replacement. The set
     * to be sampled consists of the integers from 0 through <I>N</I>&minus;1
     * inclusive. The sample consists of <I>n</I> items. Each item in the set is
     * equally likely to occur in the sample. The items in the sample are stored
     * in ascending order in the given <TT>buf</TT> starting at index 0. The
     * iterator uses the given <TT>prng</TT> to sample items at random. For each
     * item sampled, the iterator consumes one random number from the
     * <TT>prng</TT>.
     * <P>
     * As a special case, if <I>n</I> &ge; <I>N</I>, then all items in the set
     * are stored in ascending order in the given <TT>buf</TT> starting at index
     * 0, and the iterator does not use the <TT>prng</TT>.
     * <P>
     * This method uses the algorithm in A. Bissell, Ordered random selection
     * without replacement, <I>Applied Statistics,</I> 35(1):73-75, 1986.
     *
     * @param prng Pseudorandom number generator.
     * @param n Number of items in the sample.
     * @param N Number of items in the set.
     * @param buf Array in which to store the sampled items.
     * @return Number of sampled items actually stored in <TT>buf</TT>.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>n</I> &lt; 0 or <I>N</I> &lt; 0.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null.
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>buf.length</TT> &lt; <TT>n</TT>.
     */
    public static int withoutReplacement(Random prng,
            int n,
            long N,
            long[] buf) {
        if (n < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): n < 0 illegal");
        } else if (N < 0) {
            throw new IllegalArgumentException("RandomSample.withoutReplacement(): N < 0 illegal");
        } else if (n < N) {
            long M = N;
            long r = N - n;
            for (int i = 0; i < n; ++i) {
                double x = prng.nextDouble();
                double p = 1.0;
                for (;;) {
                    p = p * r / M;
                    if (p <= x) {
                        buf[i] = N - M;
                        --M;
                        break;
                    } else {
                        --M;
                        --r;
                    }
                }
            }
            return n;
        } else {
            for (int i = 0; i < N; ++i) {
                buf[i] = i;
            }
            return (int) N;
        }
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length != 3)
//			{
//			System.err.println ("Usage: java edu.rit.util.RandomSample <n> <N> <seed>");
//			System.exit (1);
//			}
//		int n = Integer.parseInt (args[0]);
//		int N = Integer.parseInt (args[1]);
//		long seed = Long.parseLong (args[2]);
//		Random prng = Random.getInstance (seed);
//		Iterator<Integer> iter = withoutReplacement (prng, n, N);
//		while (iter.hasNext()) System.out.printf (" %d", iter.next());
//		System.out.println();
//		}
}
