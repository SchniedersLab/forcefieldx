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
package ffx.numerics.fft;

/**
 * Compute the 3D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs.
 * <p>
 * The location of the input point [i, j, k] within the input array must be:
 * <br>
 * double real = input[x*nextX + y*nextY + z*nextZ]
 * <br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + 1]
 * <br>
 * where
 * <br> int nextX = 2
 * <br> int nextY = 2*nX
 * <br> int nextZ = 2*nX*nY<br>
 * <p>
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex3D {

    private final int nX, nY, nZ;
    private final int nZ2;
    private final int nextX, nextY, nextZ;
    private final double[] work;
    private final Complex fftX, fftY, fftZ;
    private final double[] recip;

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX X-dimension.
     * @param nY Y-dimension.
     * @param nZ Z-dimension.
     */
    public Complex3D(int nX, int nY, int nZ) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        nZ2 = 2 * nZ;
        nextX = 2;
        nextY = 2 * nX;
        nextZ = nextY * nY;

        work = new double[nZ2];
        recip = new double[nX * nY * nZ];

        fftX = new Complex(nX);
        fftY = new Complex(nY);
        fftZ = new Complex(nZ);
    }

    /**
     * Compute the 3D FFT.
     *
     * @param input The input array must be of size 2 * nX * nY * nZ.
     */
    public void fft(final double[] input) {
        for (int z = 0; z < nZ; z++) {
            for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.fft(input, offset, nextX);
            }
            for (int offset = z * nextZ, x = 0; x < nX; x++, offset += nextX) {
                fftY.fft(input, offset, nextY);
            }
        }
        for (int x = 0, offset = 0; x < nX; x++) {
            for (int y = 0; y < nY; y++, offset += nextX) {
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    work[i] = input[z];
                    work[i + 1] = input[z + 1];
                }
                fftZ.fft(work, 0, 2);
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    input[z] = work[i];
                    input[z + 1] = work[i + 1];
                }
            }
        }
    }

    /**
     * Compute the inverse 3D FFT.
     *
     * @param input The input array must be of size 2 * nX * nY * nZ.
     */
    public void ifft(final double[] input) {
        for (int offset = 0, x = 0; x < nX; x++) {
            for (int y = 0; y < nY; y++, offset += nextX) {
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    work[i] = input[z];
                    work[i + 1] = input[z + 1];
                }
                fftZ.ifft(work, 0, 2);
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    input[z] = work[i];
                    input[z + 1] = work[i + 1];
                }
            }
        }
        for (int z = 0; z < nZ; z++) {
            for (int offset = z * nextZ, x = 0; x < nX; x++, offset += nextX) {
                fftY.ifft(input, offset, nextY);
            }
            for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.ifft(input, offset, nextX);
            }
        }
    }

    /**
     * <p>
     * Setter for the field <code>recip</code>.</p>
     *
     * @param recip an array of double.
     */
    public void setRecip(double[] recip) {
        // Reorder the reciprocal space data into the order it is needed by the convolution routine.
        for (int offset = 0, index = 0, y = 0; y < nY; y++) {
            for (int x = 0; x < nX; x++, offset += 1) {
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nX * nY) {
                    this.recip[index++] = recip[z];
                }
            }
        }
    }

    /**
     * <p>
     * convolution</p>
     *
     * @param input an array of double.
     */
    public void convolution(final double[] input) {
        for (int z = 0; z < nZ; z++) {
            for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.fft(input, offset, nextX);
            }
            for (int offset = z * nextZ, x = 0; x < nX; x++, offset += nextX) {
                fftY.fft(input, offset, nextY);
            }
        }
        for (int offset = 0, index = 0, y = 0; y < nY; y++) {
            for (int x = 0; x < nX; x++, offset += nextX) {
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    work[i] = input[z];
                    work[i + 1] = input[z + 1];
                }
                fftZ.fft(work, 0, 2);
                for (int i = 0; i < nZ2; i += 2) {
                    double r = recip[index++];
                    work[i] *= r;
                    work[i + 1] *= r;
                }
                fftZ.ifft(work, 0, 2);
                for (int i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    input[z] = work[i];
                    input[z + 1] = work[i + 1];
                }
            }
        }
        for (int z = 0; z < nZ; z++) {
            for (int offset = z * nextZ, x = 0; x < nX; x++, offset += nextX) {
                fftY.ifft(input, offset, nextY);
            }
            for (int offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.ifft(input, offset, nextX);
            }
        }
    }

    /**
     * <p>
     * iComplex3D</p>
     *
     * @param i  a int.
     * @param j  a int.
     * @param k  a int.
     * @param nX a int.
     * @param nY a int.
     * @return a int.
     */
    public static int iComplex3D(int i, int j, int k, int nX, int nY) {
        return 2 * (i + nX * (j + nY * k));
    }
}
