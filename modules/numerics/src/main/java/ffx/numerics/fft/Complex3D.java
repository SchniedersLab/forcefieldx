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
package ffx.numerics.fft;

/**
 * Compute the 3D FFT of complex, double precision input of arbitrary dimensions via 1D Mixed Radix
 * FFTs.
 *
 * <p>
 * For interleaved data, the location of the input point [x, y, z] within the input array must be:
 * <br>
 * double real = input[x*nextX + y*nextY + z*nextZ] <br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + 1] <br>
 * where <br>
 * int nextX = 2 <br>
 * int nextY = 2*nX <br>
 * int nextZ = 2*nX*nY <br>
 *
 * <p>
 * For blocked data along x, the location of the input point [x, y, z] within the input array must be:
 * <br>
 * double real = input[x*nextX + y*nextY + z*nextZ] <br>
 * double imag = input[x*nextX + y*nextY + z*nextZ + nX] <br>
 * where <br>
 * int nextX = 1 <br>
 * int nextY = 2*nX <br>
 * int nextZ = 2*nX*nY <br>
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex3D {

  private final int nX, nY, nZ;
  private final int im;

  private final int nextX, nextY, nextZ;
  private final double[] work;
  private final Complex fftX, fftY, fftZ;
  private final double[] recip;

  private final int internalImZ;
  private final int internalNextZ;

  /**
   * Initialize the 3D FFT for complex 3D matrix using interleaved data layout.
   *
   * @param nX X-dimension.
   * @param nY Y-dimension.
   * @param nZ Z-dimension.
   */
  public Complex3D(int nX, int nY, int nZ) {
    this(nX, nY, nZ, DataLayout3D.INTERLEAVED);
  }

  /**
   * Initialize the 3D FFT for complex 3D matrix.
   *
   * @param nX     X-dimension.
   * @param nY     Y-dimension.
   * @param nZ     Z-dimension.
   * @param layout Data layout.
   */
  public Complex3D(int nX, int nY, int nZ, DataLayout3D layout) {
    this.nX = nX;
    this.nY = nY;
    this.nZ = nZ;

    DataLayout1D dataLayout1D;
    switch (layout) {
      default:
      case INTERLEAVED:
        // Interleaved data layout.
        im = 1;
        nextX = 2;
        nextY = 2 * nX;
        nextZ = 2 * nX * nY;
        // Internal 1D FFTs will be performed in interleaved format.
        dataLayout1D = DataLayout1D.INTERLEAVED;
        // Transforms along the Z-axis will be repacked into 1D interleaved format.
        internalImZ = 1;
        internalNextZ = 2;
        break;
      case BLOCKED_X:
        // Blocking is along the X-axis.
        im = nX;
        nextX = 1;
        nextY = 2 * nX;
        nextZ = 2 * nX * nY;
        // Internal 1D FFTs will be performed in blocked format.
        dataLayout1D = DataLayout1D.BLOCKED;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalNextZ = 1;
        internalImZ = nZ;
        break;
      case BLOCKED_XY:
        // Blocking is based on 2D XY-planes.
        im = nX * nY;
        nextX = 1;
        nextY = nX;
        nextZ = 2 * nY * nX;
        // Internal 1D FFTs will be performed in blocked format.
        dataLayout1D = DataLayout1D.BLOCKED;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalNextZ = 1;
        internalImZ = nZ;
        break;
      case BLOCKED_XYZ:
        // Blocking is based on 3D XYZ-volume with all real values followed by all imaginary.
        im = nX * nY * nZ;
        nextX = 1;
        nextY = nX;
        nextZ = nY * nX;
        // Internal 1D FFTs will be performed in blocked format.
        dataLayout1D = DataLayout1D.BLOCKED;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalNextZ = 1;
        internalImZ = nZ;
        break;
    }

    // Allocate memory for the reciprocal space data to be repacked into the order needed by the convolution routine.
    recip = new double[nX * nY * nZ];
    fftX = new Complex(nX, dataLayout1D, im);
    fftY = new Complex(nY, dataLayout1D, im);
    // Z-data is always repacked into a local work array.
    work = new double[2 * nZ];
    fftZ = new Complex(nZ, dataLayout1D, internalImZ);
  }

  /**
   * Determine the index of the complex number in the 1D array from the X, Y and Z indices.
   *
   * @param i  the index along the X-axis.
   * @param j  the index along the Y-axis.
   * @param k  the index along the Z-axis.
   * @param nX the number of points along the X-axis.
   * @param nY the number of points along the Y-axis.
   * @return the index of the complex number in the 1D array.
   */
  public static int interleavedIndex(int i, int j, int k, int nX, int nY) {
    return 2 * (i + nX * (j + nY * k));
  }

  /**
   * Determine the index of the complex number in the blocked 1D array from the X, Y and Z indices.
   *
   * @param i  the index along the X-axis.
   * @param j  the index along the Y-axis.
   * @param k  the index along the Z-axis.
   * @param nX the number of points along the X-axis.
   * @param nY the number of points along the Y-axis.
   * @return the index of the complex number in the 1D array using blocked data layout.
   */
  public static int index3D(int i, int j, int k, int nX, int nY, DataLayout3D layout) {
    return switch (layout) {
      case INTERLEAVED -> interleavedIndex(i, j, k, nX, nY);
      case BLOCKED_X -> i + 2 * (nX * (j + nY * k));
      case BLOCKED_XY -> i + nX * (j + 2 * nY * k);
      case BLOCKED_XYZ -> i + nX * (j + nY * k);
    };
  }

  /**
   * Perform a convolution.
   *
   * @param input The input array.
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
        for (int i = 0, z = offset; i < nZ; i++, z += nextZ) {
          int w = i * nextX;
          work[w] = input[z];
          work[w + internalImZ] = input[z + im];
        }
        fftZ.fft(work, 0, internalNextZ);
        for (int i = 0; i < nZ; i++) {
          double r = recip[index++];
          int w = i * nextX;
          work[w] *= r;
          work[w + internalImZ] *= r;
        }
        fftZ.ifft(work, 0, internalNextZ);
        for (int i = 0, z = offset; i < nZ; i++, z += nextZ) {
          int w = i * nextX;
          input[z] = work[w];
          input[z + im] = work[w + internalImZ];
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
        for (int i = 0, z = offset; i < nZ; i++, z += nextZ) {
          int w = i * nextX;
          work[w] = input[z];
          work[w + internalImZ] = input[z + im];
        }
        fftZ.fft(work, 0, im);
        for (int i = 0, z = offset; i < nZ; i++, z += nextZ) {
          int w = i * nextX;
          input[z] = work[w];
          input[z + im] = work[w + internalImZ];
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
        for (int i = 0, z = offset; i < nZ; i++, z += nextZ) {
          int w = i * nextX;
          work[w] = input[z];
          work[w + internalImZ] = input[z + im];
        }
        fftZ.ifft(work, 0, internalNextZ);
        for (int i = 0, z = offset; i < nZ; i++, z += nextZ) {
          int w = i * nextX;
          input[z] = work[w];
          input[z + im] = work[w + internalImZ];
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
   * Setter for the field <code>recip</code>.
   *
   * @param recip an array of double.
   */
  public void setRecip(double[] recip) {
    // Reorder the reciprocal space data into the order it is needed by the convolution routine.
    int index = 0;
    for (int offset = 0, y = 0; y < nY; y++) {
      for (int x = 0; x < nX; x++, offset++) {
        for (int i = 0, z = offset; i < nZ; i++, z += nX * nY) {
          this.recip[index++] = recip[z];
        }
      }
    }
  }
}
