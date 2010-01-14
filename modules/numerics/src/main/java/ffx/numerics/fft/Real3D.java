/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.numerics.fft;

/**
 * Compute the 3D FFT of real, double precision input of arbitrary dimensions.<p>
 * 
 * @author Michal J. Schnieders
 */
public class Real3D {

    private final int n, nextX, nextY, nextZ;
    private final int nX, nY, nZ;
    private final int nX1, nZ2;
    private final double work[];
    private final double recip[];
    private final Real fftX;
    private final Complex fftY, fftZ;

    /**
     * Initialize the 3D FFT for complex 3D matrix.
     *
     * @param nX
     *            X-dimension.
     * @param nY
     *            Y-dimension.
     * @param nZ
     *            Z-dimension.
     */
    public Real3D(int nX, int nY, int nZ) {
        this.n = nX;
        this.nX = nX / 2;
        this.nY = nY;
        this.nZ = nZ;
        nX1 = this.nX + 1;
        nZ2 = 2 * nZ;
        nextX = 2;
        nextY = n + 2;
        nextZ = nextY * nY;
        work = new double[nZ2];
        recip = new double[nX1 * nY * nZ];
        fftX = new Real(n);
        fftY = new Complex(nY);
        fftZ = new Complex(nZ);
    }

    /**
     * Compute the 3D FFT.
     *
     * @param input
     *            The input array must be of size (nX + 2) * nY * nZ.
     */
    public void fft(final double input[]) {
        int i, x, y, z, offset, stride;
        for (z = 0; z < nZ; z++) {
            for (offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.fft(input, offset);
            }
            for (offset = z * nextZ, stride = nextY, x = 0; x < nX1; x++, offset += nextX) {
                fftY.fft(input, offset, stride);
            }
        }
        for (x = 0; x < nX1; x++) {
            for (offset = x * 2, y = 0; y < nY; y++, offset += nextY) {
                for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    work[i] = input[z];
                    work[i + 1] = input[z + 1];
                }
                fftZ.fft(work, 0, 2);
                for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    input[z] = work[i];
                    input[z + 1] = work[i + 1];
                }
            }
        }
    }

    /**
     * Compute the inverese 3D FFT.
     *
     * @param input
     *            The input array must be of size (nX + 2) * nY * nZ.
     */
    public void ifft(final double input[]) {
        int i, x, y, z, stride, offset;
        for (x = 0; x < nX1; x++) {
            for (offset = x * 2, y = 0; y < nY; y++, offset += nextY) {
                for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    work[i] = input[z];
                    work[i + 1] = input[z + 1];
                }
                fftZ.ifft(work, 0, 2);
                for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    input[z] = work[i];
                    input[z + 1] = work[i + 1];
                }
            }
        }
        for (z = 0; z < nZ; z++) {
            for (stride = nextY, offset = z * nextZ, x = 0; x < nX1; x++, offset += nextX) {
                fftY.ifft(input, offset, stride);
            }
            for (offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.ifft(input, offset);
            }
        }
    }

    public void setRecip(double recip[]) {
        int offset, y, x, z, i;

        /**
         * Reorder the reciprocal space data into the order it is needed
         * by the convolution routine.
         */
        int index = 0;
        for (index = 0, offset = 0, y = 0; y < nY; y++) {
            for (x = 0; x < nX1; x++, offset += 1) {
                for (i = 0, z = offset; i < nZ; i++, z += nX1 * nY) {
                    this.recip[index++] = recip[z];
                }
            }
        }
    }

    public void convolution(final double input[]) {
        int i, j, x, y, z, index, offset, stride;
        for (z = 0; z < nZ; z++) {
            for (offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.fft(input, offset);
            }
            for (offset = z * nextZ, stride = nextY, x = 0; x < nX1; x++, offset += nextX) {
                fftY.fft(input, offset, stride);
            }
        }
        for (index = 0, offset = 0, y = 0; y < nY; y++) {
            for (x = 0; x < nX1; x++, offset += nextX) {
                for (i = 0, j = 0, z = offset; i < nZ2; i += 2, j++, z += nextZ) {
                    work[i] = input[z];
                    work[i + 1] = input[z + 1];
                }
                fftZ.fft(work, 0, 2);
                for (i = 0; i < nZ2; i += 2) {
                    double r = recip[index++];
                    work[i] *= r;
                    work[i + 1] *= r;
                }
                fftZ.ifft(work, 0, 2);
                for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
                    input[z] = work[i];
                    input[z + 1] = work[i + 1];
                }
            }
        }
        for (z = 0; z < nZ; z++) {
            for (stride = nextY, offset = z * nextZ, x = 0; x < nX1; x++, offset += nextX) {
                fftY.ifft(input, offset, stride);
            }
            for (offset = z * nextZ, y = 0; y < nY; y++, offset += nextY) {
                fftX.ifft(input, offset);
            }
        }
    }
}
