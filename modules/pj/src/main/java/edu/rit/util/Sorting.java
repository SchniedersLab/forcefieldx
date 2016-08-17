//******************************************************************************
//
// File:    Sorting.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Sorting
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

/**
 * Class Sorting provides static methods for sorting arrays of primitive types
 * and object types.
 * <P>
 * <I>Note:</I> The operations in class Sorting are not multiple thread safe.
 *
 * @author Alan Kaminsky
 * @version 02-Nov-2011
 */
public class Sorting {

// Prevent construction.
    private Sorting() {
    }

// Exported helper classes.
    /**
     * Class Sorting.Byte is the base class for a helper object used to sort an
     * array of type <TT>byte[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Byte {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(byte[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(byte[] x,
                int a,
                int b) {
            byte t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Character is the base class for a helper object used to
     * sort an array of type <TT>char[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Character {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(char[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(char[] x,
                int a,
                int b) {
            char t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Short is the base class for a helper object used to sort an
     * array of type <TT>short[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Short {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(short[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(short[] x,
                int a,
                int b) {
            short t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Integer is the base class for a helper object used to sort
     * an array of type <TT>int[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Integer {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(int[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(int[] x,
                int a,
                int b) {
            int t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Long is the base class for a helper object used to sort an
     * array of type <TT>long[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Long {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(long[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(long[] x,
                int a,
                int b) {
            long t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Float is the base class for a helper object used to sort an
     * array of type <TT>float[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Float {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(float[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(float[] x,
                int a,
                int b) {
            float t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Double is the base class for a helper object used to sort
     * an array of type <TT>double[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static class Double {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         * <P>
         * The default implementation returns true if <TT>x[a] &lt; x[b]</TT>,
         * which sorts the array into ascending order. A subclass can override
         * this method to obtain a different ordering criterion; for example,
         * descending order.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public boolean comesBefore(double[] x,
                int a,
                int b) {
            return x[a] < x[b];
        }

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(double[] x,
                int a,
                int b) {
            double t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    /**
     * Class Sorting.Object is the abstract base class for a helper object used
     * to sort an array of objects of type <TT>T[]</TT>.
     *
     * @param <T> Data type of the array elements.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2010
     */
    public static abstract class Object<T> {

        /**
         * Compare two elements in the given array. This determines the order of
         * the elements in the sorted array.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being compared.
         * @param b Index of second array element being compared.
         *
         * @return True if <TT>x[a]</TT> comes before <TT>x[b]</TT> in the
         * desired ordering, false otherwise.
         */
        public abstract boolean comesBefore(T[] x,
                int a,
                int b);

        /**
         * Swap two elements in the given array.
         * <P>
         * The default implementation swaps <TT>x[a]</TT> with <TT>x[b]</TT>. A
         * subclass can override this method to do something different; for
         * example, to swap the elements of other arrays in addition to
         * <TT>x</TT>.
         *
         * @param x Array being sorted.
         * @param a Index of first array element being swapped.
         * @param b Index of second array element being swapped.
         */
        public void swap(T[] x,
                int a,
                int b) {
            T t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

// Exported operations.
    /**
     * Sort the given array of type <TT>byte[]</TT>. The given helper object is
     * used to determine the desired ordering of the array elements and to swap
     * the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>) heapsort
     * algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static byte[] sort(byte[] x,
            Sorting.Byte helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(byte[] x,
            int c, // 1-based index
            Sorting.Byte helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(byte[] x,
            int n, // 1-based index
            Sorting.Byte helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given array of type <TT>char[]</TT>. The given helper object is
     * used to determine the desired ordering of the array elements and to swap
     * the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>) heapsort
     * algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static char[] sort(char[] x,
            Sorting.Character helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(char[] x,
            int c, // 1-based index
            Sorting.Character helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(char[] x,
            int n, // 1-based index
            Sorting.Character helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given array of type <TT>short[]</TT>. The given helper object is
     * used to determine the desired ordering of the array elements and to swap
     * the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>) heapsort
     * algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static short[] sort(short[] x,
            Sorting.Short helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(short[] x,
            int c, // 1-based index
            Sorting.Short helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(short[] x,
            int n, // 1-based index
            Sorting.Short helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given array of type <TT>int[]</TT>. The given helper object is
     * used to determine the desired ordering of the array elements and to swap
     * the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>) heapsort
     * algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static int[] sort(int[] x,
            Sorting.Integer helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(int[] x,
            int c, // 1-based index
            Sorting.Integer helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(int[] x,
            int n, // 1-based index
            Sorting.Integer helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given array of type <TT>long[]</TT>. The given helper object is
     * used to determine the desired ordering of the array elements and to swap
     * the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>) heapsort
     * algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static long[] sort(long[] x,
            Sorting.Long helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(long[] x,
            int c, // 1-based index
            Sorting.Long helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(long[] x,
            int n, // 1-based index
            Sorting.Long helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given array of type <TT>float[]</TT>. The given helper object is
     * used to determine the desired ordering of the array elements and to swap
     * the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>) heapsort
     * algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static float[] sort(float[] x,
            Sorting.Float helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(float[] x,
            int c, // 1-based index
            Sorting.Float helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(float[] x,
            int n, // 1-based index
            Sorting.Float helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given array of type <TT>double[]</TT>. The given helper object
     * is used to determine the desired ordering of the array elements and to
     * swap the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>)
     * heapsort algorithm is used.
     *
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static double[] sort(double[] x,
            Sorting.Double helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static void siftUp(double[] x,
            int c, // 1-based index
            Sorting.Double helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static void siftDown(double[] x,
            int n, // 1-based index
            Sorting.Double helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

    /**
     * Sort the given object array of type <TT>T[]</TT>. The given helper object
     * is used to determine the desired ordering of the array elements and to
     * swap the array elements. An <I>O</I>(<I>n</I>&nbsp;log&nbsp;<I>n</I>)
     * heapsort algorithm is used.
     *
     * @param <T> Data type of the array elements.
     * @param x Array to be sorted.
     * @param helper Helper object.
     * @return The array that was sorted (<TT>x</TT>).
     */
    public static <T> T[] sort(T[] x,
            Sorting.Object<T> helper) {
        int n = x.length;
        for (int i = 2; i <= n; ++i) {
            siftUp(x, i, helper);
        }
        for (int i = n; i >= 2; --i) {
            helper.swap(x, 0, i - 1);
            siftDown(x, i - 1, helper);
        }
        return x;
    }

    private static <T> void siftUp(T[] x,
            int c, // 1-based index
            Sorting.Object<T> helper) {
        int p = c >> 1; // 1-based index
        while (p >= 1) {
            if (helper.comesBefore(x, p - 1, c - 1)) {
                helper.swap(x, p - 1, c - 1);
            } else {
                break;
            }
            c = p;
            p = c >> 1;
        }
    }

    private static <T> void siftDown(T[] x,
            int n, // 1-based index
            Sorting.Object<T> helper) {
        int p = 1; // 1-based index
        int ca = 2; // 1-based index
        int cb = 3; // 1-based index
        while (ca <= n) {
            if (cb <= n && helper.comesBefore(x, ca - 1, cb - 1)) {
                if (helper.comesBefore(x, p - 1, cb - 1)) {
                    helper.swap(x, p - 1, cb - 1);
                    p = cb;
                } else {
                    break;
                }
            } else {
                if (helper.comesBefore(x, p - 1, ca - 1)) {
                    helper.swap(x, p - 1, ca - 1);
                    p = ca;
                } else {
                    break;
                }
            }
            ca = p << 1;
            cb = ca + 1;
        }
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.util.Sorting <I>n</I> <I>seed</I>
//	 * <BR><I>n</I> = Array length
//	 * <BR><I>seed</I> = Random seed
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length != 2) usage();
//		int n = java.lang.Integer.parseInt (args[0]);
//		long seed = java.lang.Long.parseLong (args[1]);
//		byte[] x = new byte [n];
//		Random prng = Random.getInstance (seed);
//		for (int i = 0; i < n; ++ i) x[i] = prng.nextByte();
//		for (int i = 0; i < n; ++ i) System.out.printf ("%d  ", x[i]);
//		System.out.println();
//		Sorting.sort (x, new Sorting.Byte());
//		for (int i = 0; i < n; ++ i) System.out.printf ("%d  ", x[i]);
//		System.out.println();
//		}
//
//	/**
//	 * Print a usage message and exit.
//	 */
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.util.Sorting <n> <seed>");
//		System.err.println ("<n> = Array length");
//		System.err.println ("<seed> = Random seed");
//		System.exit (1);
//		}
}
