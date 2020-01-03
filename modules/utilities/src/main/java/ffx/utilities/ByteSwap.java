//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.utilities;

/**
 * <p>
 * ByteSwap class.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class ByteSwap {

    /**
     * Byte swap a single short value.
     *
     * @param value Value to byte swap.
     * @return Byte swapped representation.
     */
    public static short swap(short value) {
        int b1 = value & 0xff;
        int b2 = (value >> 8) & 0xff;

        return (short) (b1 << 8 | b2);
    }

    /**
     * Byte swap a single int value.
     *
     * @param value Value to byte swap.
     * @return Byte swapped representation.
     */
    public static int swap(int value) {
        int b1 = (value) & 0xff;
        int b2 = (value >> 8) & 0xff;
        int b3 = (value >> 16) & 0xff;
        int b4 = (value >> 24) & 0xff;

        return b1 << 24 | b2 << 16 | b3 << 8 | b4;
    }

    /**
     * Byte swap a single long value.
     *
     * @param value Value to byte swap.
     * @return Byte swapped representation.
     */
    public static long swap(long value) {
        long b1 = (value) & 0xff;
        long b2 = (value >> 8) & 0xff;
        long b3 = (value >> 16) & 0xff;
        long b4 = (value >> 24) & 0xff;
        long b5 = (value >> 32) & 0xff;
        long b6 = (value >> 40) & 0xff;
        long b7 = (value >> 48) & 0xff;
        long b8 = (value >> 56) & 0xff;

        return b1 << 56 | b2 << 48 | b3 << 40 | b4 << 32
                | b5 << 24 | b6 << 16 | b7 << 8 | b8;
    }

    /**
     * Byte swap a single float value.
     *
     * @param value Value to byte swap.
     * @return Byte swapped representation.
     */
    public static float swap(float value) {
        int intValue = Float.floatToIntBits(value);
        intValue = swap(intValue);

        return Float.intBitsToFloat(intValue);
    }

    /**
     * Byte swap a single double value.
     *
     * @param value Value to byte swap.
     * @return Byte swapped representation.
     */
    public static double swap(double value) {
        long longValue = Double.doubleToLongBits(value);
        longValue = swap(longValue);

        return Double.longBitsToDouble(longValue);
    }

    /**
     * Byte swap an array of shorts. The result of the swapping is put back into
     * the specified array.
     *
     * @param array Array of values to swap
     */
    public static void swap(short[] array) {
        int n = array.length;
        for (int i = 0; i < n; i++) {
            array[i] = swap(array[i]);
        }
    }

    /**
     * Byte swap an array of ints. The result of the swapping is put back into
     * the specified array.
     *
     * @param array Array of values to swap
     */
    public static void swap(int[] array) {
        int n = array.length;
        for (int i = 0; i < n; i++) {
            array[i] = swap(array[i]);
        }
    }

    /**
     * Byte swap an array of longs. The result of the swapping is put back into
     * the specified array.
     *
     * @param array Array of values to swap
     */
    public static void swap(long[] array) {
        int n = array.length;
        for (int i = 0; i < n; i++) {
            array[i] = swap(array[i]);
        }
    }

    /**
     * Byte swap an array of floats. The result of the swapping is put back into
     * the specified array.
     *
     * @param array Array of values to swap
     */
    public static void swap(float[] array) {
        int n = array.length;
        for (int i = 0; i < n; i++) {
            array[i] = swap(array[i]);
        }
    }

    /**
     * Byte swap an array of doubles. The result of the swapping is put back
     * into the specified array.
     *
     * @param array Array of values to swap
     */
    public static void swap(double[] array) {
        int n = array.length;
        for (int i = 0; i < n; i++) {
            array[i] = swap(array[i]);
        }
    }
}
