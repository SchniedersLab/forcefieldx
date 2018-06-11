//******************************************************************************
//
// File:    Searching.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Searching
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

import java.util.Comparator;

/**
 * Class Searching provides static methods for searching arrays of primitive
 * types and object types.
 * <P>
 * <I>Note:</I> The operations in class Searching are not multiple thread safe.
 *
 * @author Alan Kaminsky
 * @version 22-Nov-2011
 */
@SuppressWarnings("unchecked")
public class Searching {

// Prevent construction.
    private Searching() {
    }

// Exported helper classes.
    /**
     * Class Searching.Byte is the base class for a helper object used to search
     * an array of type <TT>byte[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Byte {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(byte a,
                byte b) {
            return a - b;
        }
    }

    /**
     * Class Searching.Character is the base class for a helper object used to
     * search an array of type <TT>char[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Character {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(char a,
                char b) {
            return a - b;
        }
    }

    /**
     * Class Searching.Short is the base class for a helper object used to
     * search an array of type <TT>short[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Short {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(short a,
                short b) {
            return a - b;
        }
    }

    /**
     * Class Searching.Integer is the base class for a helper object used to
     * search an array of type <TT>int[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Integer {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(int a,
                int b) {
            return a - b;
        }
    }

    /**
     * Class Searching.Long is the base class for a helper object used to search
     * an array of type <TT>long[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Long {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(long a,
                long b) {
            long d = a - b;
            return d < 0L ? -1 : d > 0L ? 1 : 0;
        }
    }

    /**
     * Class Searching.Float is the base class for a helper object used to
     * search an array of type <TT>float[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Float {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(float a,
                float b) {
            float d = a - b;
            return d < 0.0f ? -1 : d > 0.0f ? 1 : 0;
        }
    }

    /**
     * Class Searching.Double is the base class for a helper object used to
     * search an array of type <TT>double[]</TT>.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static class Double {

        /**
         * Compare two elements according to the desired ordering criterion.
         * <P>
         * The default implementation compares <TT>a</TT> and <TT>b</TT> using
         * ascending order. A subclass can override this method to obtain a
         * different ordering criterion; for example, descending order.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public int compare(double a,
                double b) {
            double d = a - b;
            return d < 0.0 ? -1 : d > 0.0 ? 1 : 0;
        }
    }

    /**
     * Class Searching.Object is the base class for a helper object used to
     * search an array of type <TT>T[]</TT>.
     *
     * @param <T> Array element data type.
     *
     * @author Alan Kaminsky
     * @version 22-Nov-2011
     */
    public static abstract class Object<T> {

        /**
         * Compare two elements according to the desired ordering criterion.
         *
         * @param a First element being compared.
         * @param b Second element being compared.
         *
         * @return A number less than, equal to, or greater than 0 if
         * <TT>a</TT> comes before, is the same as, or comes after
         * <TT>b</TT>, respectively.
         */
        public abstract int compare(T a,
                T b);
    }

// Exported operations.
    /**
     * Search the given unordered array of type <TT>byte[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(byte[] x,
            byte a,
            Searching.Byte helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>byte[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(byte[] x,
            byte a,
            Searching.Byte helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>char[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(char[] x,
            char a,
            Searching.Character helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>char[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(char[] x,
            char a,
            Searching.Character helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>short[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(short[] x,
            short a,
            Searching.Short helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>short[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(short[] x,
            short a,
            Searching.Short helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>int[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(int[] x,
            int a,
            Searching.Integer helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>int[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(int[] x,
            int a,
            Searching.Integer helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>long[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(long[] x,
            long a,
            Searching.Long helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>long[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(long[] x,
            long a,
            Searching.Long helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>float[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(float[] x,
            float a,
            Searching.Float helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>float[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(float[] x,
            float a,
            Searching.Float helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>double[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchUnsorted(double[] x,
            double a,
            Searching.Double helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>double[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static int searchSorted(double[] x,
            double a,
            Searching.Double helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>T[]</TT> for the given
     * element. The given helper object is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param <T> Array element data type.
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static <T> int searchUnsorted(T[] x,
            T a,
            Searching.Object<T> helper) {
        for (int i = 0; i < x.length; ++i) {
            if (helper.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>T[]</TT> for the given
     * element. The given helper object is used to compare elements for order
     * and equality. It is assumed that the array is sorted in the order
     * determined by the helper object; otherwise, the <TT>searchSorted()</TT>
     * method's behavior is not specified. An <I>O</I>(log <I>n</I>) binary
     * search algorithm is used.
     *
     * @param <T> Array element data type.
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param helper Helper object.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static <T> int searchSorted(T[] x,
            T a,
            Searching.Object<T> helper) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = helper.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = helper.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = helper.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>T[]</TT> for the given
     * element. The given comparator is used to compare elements for equality
     * only. An <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param <T> Array element data type.
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param comp Comparator.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static <T> int searchUnsorted(T[] x,
            T a,
            Comparator<T> comp) {
        for (int i = 0; i < x.length; ++i) {
            if (comp.compare(x[i], a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>T[]</TT> for the given
     * element. The given comparator is used to compare elements for order and
     * equality. It is assumed that the array is sorted in the order determined
     * by the comparator; otherwise, the <TT>searchSorted()</TT> method's
     * behavior is not specified. An <I>O</I>(log <I>n</I>) binary search
     * algorithm is used.
     *
     * @param <T> Array element data type.
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @param comp Comparator.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static <T> int searchSorted(T[] x,
            T a,
            Comparator<T> comp) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = comp.compare(x[lo], a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = comp.compare(x[hi], a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = comp.compare(x[mid], a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

    /**
     * Search the given unordered array of type <TT>T[]</TT> for the given
     * element. The array element's natural ordering (<TT>compareTo()</TT>
     * method) is used to compare elements for equality only. An
     * <I>O</I>(<I>n</I>) linear search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static <T extends Comparable> int searchUnsorted(T[] x,
            T a) {
        for (int i = 0; i < x.length; ++i) {
            if (x[i].compareTo(a) == 0) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Search the given ordered array of type <TT>T[]</TT> for the given
     * element. The array element's natural ordering (<TT>compareTo()</TT>
     * method) is used to compare elements for order and equality. It is assumed
     * that the array is sorted in the order determined by the natural ordering;
     * otherwise, the <TT>searchSorted()</TT> method's behavior is not
     * specified. An <I>O</I>(log <I>n</I>) binary search algorithm is used.
     *
     * @param x Array to be searched.
     * @param a Element to be searched for.
     * @return If an element the same as <TT>a</TT> exists in <TT>x</TT>, then
     * the index of that element is returned. Otherwise, &minus;1 is returned.
     */
    public static <T extends Comparable> int searchSorted(T[] x,
            T a) {
        // Establish loop invariant.
        if (x.length == 0) {
            return -1;
        }

        int lo = 0;
        int locomp = x[lo].compareTo(a);
        if (locomp == 0) {
            return lo;
        } else if (locomp > 0) {
            return -1;
        }

        int hi = x.length - 1;
        int hicomp = x[hi].compareTo(a);
        if (hicomp == 0) {
            return hi;
        } else if (hicomp < 0) {
            return -1;
        }

        // Loop invariant: x[lo] comes before a; x[hi] comes after a.
        while (hi - lo > 1) {
            int mid = (hi + lo) / 2;
            int midcomp = x[mid].compareTo(a);
            if (midcomp == 0) {
                return mid;
            } else if (midcomp < 0) {
                lo = mid;
                locomp = midcomp;
            } else {
                hi = mid;
                hicomp = midcomp;
            }
        }

        return locomp == 0 ? lo : hicomp == 0 ? hi : -1;
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.util.Searching <I>a</I> <I>x_elements</I>
//	 * <BR><I>a</I> = Element to be searched for
//	 * <BR><I>x_elements</I> = Array elements to be searched
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length < 1) usage();
//		int a = java.lang.Integer.parseInt (args[0]);
//		int n = args.length - 1;
//		int[] x = new int [n];
//		for (int i = 0; i < n; ++ i)
//			x[i] = java.lang.Integer.parseInt (args[i+1]);
//		Searching.Integer helper = new Searching.Integer();
//		System.out.printf ("searchUnsorted() returns %d%n",
//			Searching.searchUnsorted (x, a, helper));
//		System.out.printf ("searchSorted() returns %d%n",
//			Searching.searchSorted (x, a, helper));
//		}
//
//	/**
//	 * Print a usage message and exit.
//	 */
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.util.Searching <a> <x_elements>");
//		System.err.println ("<a> = Element to be searched for");
//		System.err.println ("<x_elements> = Array elements to be searched");
//		System.exit (1);
//		}
}
