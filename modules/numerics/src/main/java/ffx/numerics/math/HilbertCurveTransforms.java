// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
// ******************************************************************************
package ffx.numerics.math;

import org.apache.commons.math3.util.FastMath;

import java.math.BigInteger;
import java.util.Arrays;

/**
 * HilbertCurveTransforms is a class that provides static methods for converting
 * between Hilbert indices and coordinates. This is used in the torsion scan of
 * crystals. This implementation is based on the one in OpenMM, which in turn
 * is based on Rice Universities implementation. GPT4.0 was used to convert
 * the c++ code to Java. The copyright is in the ForceFieldX code base under
 * Licenses/openmm-license/hilbert-curve-license.txt.
 */

public class HilbertCurveTransforms {
    private static int adjust_rotation(int rotation, int nDims, int bits) {
        long nd1Ones = (ones(nDims) >> 1);
        bits &= (int) (-bits & nd1Ones);
        while (bits != 0) {
            bits >>= 1;
            ++rotation;
        }
        if (rotation >= nDims) {
            rotation -= nDims;
        }
        return rotation;
    }

    private static long ones(int k) {
        return ((2L << (k - 1)) - 1);
    }

    private static long rotateLeft(long arg, int nRots, int nDims) {
        return (((arg << nRots) | (arg >> (nDims - nRots))) & ones(nDims));
    }

    private static long rotateRight(long arg, int nRots, int nDims) {
        return ((arg >> nRots) | (arg << (nDims - nRots))) & ones(nDims);
    }

    private static long bitTranspose(int nDims, int nBits, long inCoords) {
        int nDims1 = nDims - 1;
        int inB = nBits;
        int utB;
        long inFieldEnds = 1;
        long inMask = ones(inB);
        long coords = 0;

        while ((utB = inB / 2) != 0) {
            int shiftAmt = nDims1 * utB;
            long utFieldEnds = inFieldEnds | (inFieldEnds << (shiftAmt + utB));
            long utMask = (utFieldEnds << utB) - utFieldEnds;
            long utCoords = 0;
            int d;
            if (inB % 2 == 1) {
                long inFieldStarts = inFieldEnds << (inB - 1);
                int oddShift = 2 * shiftAmt;
                for (d = 0; d < nDims; ++d) {
                    long in = inCoords & inMask;
                    inCoords >>= inB;
                    coords |= (in & inFieldStarts) << oddShift++;
                    in &= ~inFieldStarts;
                    in = (in | (in << shiftAmt)) & utMask;
                    utCoords |= in << (d * utB);
                }
            } else {
                for (d = 0; d < nDims; ++d) {
                    long in = inCoords & inMask;
                    inCoords >>= inB;
                    in = (in | (in << shiftAmt)) & utMask;
                    utCoords |= in << (d * utB);
                }
            }
            inCoords = utCoords;
            inB = utB;
            inFieldEnds = utFieldEnds;
            inMask = utMask;
        }
        coords |= inCoords;
        return coords;
    }

    public static long[] hilbertIndexToCoordinates(int nBonds, int nBitsPerDim, long index) {
        long[] coord = new long[nBonds];

        if (nBonds > 1) {
            long coords;
            int nbOnes = (int) ones(nBitsPerDim);
            int d;

            if (nBitsPerDim > 1) {
                int nDimsBits = nBonds * nBitsPerDim;
                int ndOnes = (int) ones(nBonds);
                int nd1Ones = ndOnes >> 1;
                int b = nDimsBits;
                int rotation = 0;
                int flipBit = 0;
                long nthbits = ones(nDimsBits) / ndOnes;
                index ^= (index ^ nthbits) >> 1;
                coords = 0;

                do {
                    int bits = (int) ((index >> (b -= nBonds)) & ndOnes);
                    coords <<= nBonds;
                    coords |= rotateLeft(bits, rotation, nBonds) ^ flipBit;
                    flipBit = 1 << rotation;
                    rotation = adjust_rotation(rotation, nBonds, bits);
                } while (b > 0);

                for (b = nBonds; b < nDimsBits; b *= 2) {
                    coords ^= coords >> b;
                }
                coords = bitTranspose(nBitsPerDim, nBonds, coords);
            } else {
                coords = index ^ (index >> 1);
            }

            for (d = 0; d < nBonds; ++d) {
                coord[d] = coords & nbOnes;
                coords >>= nBitsPerDim;
            }
        } else {
            coord[0] = index;
        }

        return coord;
    }

    public static long coordinatesToHilbertIndex(int nDims, int nBits, long[] coord) {
        long index;
        if (nDims > 1) {
            int nDimsBits = nDims * nBits;
            long coords = 0;
            for (int d = nDims; d-- > 0; ) {
                coords <<= nBits;
                coords |= coord[d];
            }

            if (nBits > 1) {
                long ndOnes = ones(nDims);
                long nd1Ones = ndOnes >> 1;
                int b = nDimsBits;
                int rotation = 0;
                long flipBit = 0;
                long nthbits = ones(nDimsBits) / ndOnes;
                coords = bitTranspose(nDims, nBits, coords);
                coords ^= coords >> nDims;
                index = 0;

                do {
                    long bits = (coords >> (b -= nDims)) & ndOnes;
                    bits = rotateRight(flipBit ^ bits, rotation, nDims);
                    index <<= nDims;
                    index |= bits;
                    flipBit = 1L << rotation;
                    rotation = adjust_rotation(rotation, nDims, (int) bits);
                } while (b > 0);

                index ^= nthbits >> 1;
            } else {
                index = coords;
            }

            for (int d = 1; d < nDimsBits; d *= 2) {
                index ^= index >> d;
            }

        } else {
            index = coord[0];
        }

        return index;
    }


    public static void main(String[] args) {
        int nBonds = 2; // Dimensions of the space
        int nTorsions = 4; // Bits per dimension
        int nBits = (int) Math.ceil(FastMath.log(2, nTorsions));
        // Calculate the maximum index of number of configurations using BigInteger
        //BigInteger maxIndex = BigInteger.valueOf(nTorsions).pow(nBonds).subtract(BigInteger.ONE);


        // Max index allowing nTorsions >= nBonds
        BigInteger maxIndex = BigInteger.valueOf(2).pow(nBonds * nBits).subtract(BigInteger.ONE);
        BigInteger numConfigs = BigInteger.valueOf(nTorsions).pow(nBonds);

        System.out.println("Maximum index: " + maxIndex);
        System.out.println("Number of configurations: " + numConfigs);
        // Iterate over all indices
        BigInteger index = BigInteger.ZERO;
        int counter = 0;
        while (index.longValue() <= maxIndex.longValue()) {
            long[] coordinates = hilbertIndexToCoordinates(nBonds, nBits, index.longValue());
            long convertedIndex = coordinatesToHilbertIndex(nBonds, nBits, coordinates);
            boolean valid = true;

            for (long coord : coordinates) {
                if (coord > nTorsions-1) {
                    valid = false;
                    break;
                }
            }

            // Print out the coordinates
            System.out.println("Hilbert index: " + counter +  "; Coordinates: " + Arrays.toString(coordinates));
            System.out.println("Converted index: " + convertedIndex);
            counter++;
            index = index.add(BigInteger.ONE);
        }
        System.out.println("Number of valid configurations: " + counter);
    }
}
