/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Biophysics Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
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
 * Compute the 3D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs.
 * <p>
 * The location of the input point [i, j, k] within the input array must be:<br>
 * <br>
 * double real = input[x*nextX + y*nextY + z*nextZ]<br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + 1]<br>
 * <br>
 * where<br>
 * int nextX = 2<br>
 * int nextY = 2*nX<br>
 * int nextZ = 2*nX*nY<br>
 * <p>
 * @author Michal J. Schnieders
 */
public class Complex3D {

    private final int nX, nY, nZ;
    private final int nZ2;
    private final int nextX, nextY, nextZ;
    private final double work[];
    private final Complex fftX, fftY, fftZ;

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
    public Complex3D(int nX, int nY, int nZ) {
        this.nX = nX;
        this.nY = nY;
        this.nZ = nZ;
        nZ2 = 2 * nZ;
        nextX = 2;
        nextY = 2 * nX;
        nextZ = nextY * nY;
        work = new double[nZ2];
        fftX = new Complex(nX);
        fftY = new Complex(nY);
        fftZ = new Complex(nZ);
    }

    /**
     * Compute the 3D FFT.
     *
     * @param input
     *            The input array must be of size 2 * nX * nY * nZ.
     */
    public void fft(final double input[]) {
        int x, y, z, i, stride, offset;

        for (z = 0; z < nZ; z++) {
            for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                fftX.fft(input, offset, stride);
            }
            for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                fftY.fft(input, offset, stride);
            }
        }
        for (x = 0, offset = 0; x < nX; x++) {
            for (y = 0; y < nY; y++, offset += nextX) {
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
     * Compute the inverse 3D FFT.
     *
     * @param input
     *            The input array must be of size 2 * nX * nY * nZ.
     */
    public void ifft(final double input[]) {
        int x, y, z, i, stride, offset;
        for (offset = 0, x = 0; x < nX; x++) {
            for (y = 0; y < nY; y++, offset += nextX) {
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
            for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                fftY.ifft(input, offset, stride);
            }
            for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                fftX.ifft(input, offset, stride);
            }
        }
    }

    public void convolution(final double input[], final double recip[]) {
        int x, y, z, i, j, index, stride, offset;
        for (z = 0; z < nZ; z++) {
            for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                fftX.fft(input, offset, stride);
            }
            for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                fftY.fft(input, offset, stride);
            }
        }
        index = 0;
        for (offset = 0, y = 0; y < nY; y++) {
            for (x = 0; x < nX; x++, offset += nextX) {
                for (i = 0, z = offset; i < nZ2; i += 2, z += nextZ) {
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
            for (offset = z * nextZ, stride = nextY, x = 0; x < nX; x++, offset += nextX) {
                fftY.ifft(input, offset, stride);
            }
            for (offset = z * nextZ, stride = nextX, y = 0; y < nY; y++, offset += nextY) {
                fftX.ifft(input, offset, stride);
            }
        }
    }
}
