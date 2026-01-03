// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
 * double imag = input[x*nextX + y*nextY + z*nextZ + im] <br>
 * where <br>
 * int nextX = 1 <br>
 * int nextY = 2*nX <br>
 * int nextZ = 2*nX*nY <br>
 * int im = nX for BLOCKED_X<br>
 * int im = nX*nY for BLOCKED_XY<br>
 * int im = nX*nY*nZ for BLOCKED_XYZ<br>
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex3D {

  /**
   * The number of points along the X-axis.
   */
  private final int nX;
  /**
   * The number of points along the Y-axis.
   */
  private final int nY;
  /**
   * The number of points along the Z-axis.
   */
  private final int nZ;
  /**
   * The offset from each real point to the corresponding imaginary point.
   */
  private final int im;
  /**
   * The offset from each real point to the next real point (either 1 or 2).
   */
  private final int ii;
  /**
   * The offset from each real point to the next real point along the X-dimension (either 1 or 2).
   */
  private final int nextX;
  /**
   * The offset from each real point to the next real point along the Y-dimension (either nX or 2*nX).
   */
  private final int nextY;
  /**
   * The offset from each real point to the next real point along the Z-dimension (either nX*nY or 2*nX*nY).
   */
  private final int nextZ;
  /**
   * The 2D FFT for the XY-plane.
   */
  private final Complex2D fftXY;
  /**
   * The 1D FFT for the Z-axis.
   */
  private final Complex fftZN;
  /**
   * The work array holding a 2D YZ-plane for transforms along the Z-axis.
   */
  private final double[] work;
  /**
   * The offset from each real point to its corresponding imaginary point for the 2D YZ-plane.
   */
  private final int internalImZ;
  /**
   * The offset from each real point to the next real point for the 2D YZ-plane (1 or 2).
   */
  private final int internalNextZ;
  /**
   * The reciprocal space data.
   */
  private final double[] recip;

  private boolean useSIMD;
  private boolean packFFTs;

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

    DataLayout2D dataLayoutXY;
    DataLayout1D dataLayoutZ;
    switch (layout) {
      default:
      case INTERLEAVED:
        // Interleaved data layout.
        im = 1;
        ii = 2;
        nextX = 2;
        nextY = 2 * nX;
        nextZ = 2 * nX * nY;
        dataLayoutXY = DataLayout2D.INTERLEAVED;
        dataLayoutZ = DataLayout1D.INTERLEAVED;
        // Transforms along the Z-axis will be repacked into 1D interleaved format.
        internalNextZ = 2;
        internalImZ = 1;
        break;
      case BLOCKED_X:
        // Blocking is along the X-axis.
        im = nX;
        ii = 1;
        nextX = 1;
        nextY = 2 * nX;
        nextZ = 2 * nX * nY;
        dataLayoutXY = DataLayout2D.BLOCKED_X;
        dataLayoutZ = DataLayout1D.BLOCKED;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalNextZ = 1;
        internalImZ = nY * nZ;
        break;
      case BLOCKED_XY:
        // Blocking is based on 2D XY-planes.
        im = nX * nY;
        ii = 1;
        nextX = 1;
        nextY = nX;
        nextZ = 2 * nY * nX;
        dataLayoutXY = DataLayout2D.BLOCKED_XY;
        dataLayoutZ = DataLayout1D.BLOCKED;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalNextZ = 1;
        internalImZ = nY * nZ;
        break;
      case BLOCKED_XYZ:
        // Blocking is based on 3D XYZ-volume with all real values followed by all imaginary.
        im = nX * nY * nZ;
        ii = 1;
        nextX = 1;
        nextY = nX;
        nextZ = nY * nX;
        dataLayoutXY = DataLayout2D.BLOCKED_XY;
        dataLayoutZ = DataLayout1D.BLOCKED;
        // Transforms along the Z-axis will be repacked into 1D blocked format.
        internalNextZ = 1;
        internalImZ = nY * nZ;
        break;
    }

    // Use SIMD by default.
    useSIMD = true;
    String simd = System.getProperty("fft.simd", Boolean.toString(useSIMD));
    try {
      useSIMD = Boolean.parseBoolean(simd);
    } catch (Exception e) {
      useSIMD = false;
    }

    packFFTs = true;
    String pack = System.getProperty("fft.pack", Boolean.toString(packFFTs));
    try {
      packFFTs = Boolean.parseBoolean(pack);
    } catch (Exception e) {
      packFFTs = false;
    }

    // Allocate memory for the reciprocal space data to be repacked into the order needed by the convolution routine.
    recip = new double[nX * nY * nZ];
    fftXY = new Complex2D(nX, nY, dataLayoutXY, im);
    fftXY.setPackFFTs(packFFTs);
    fftXY.setUseSIMD(useSIMD);
    fftZN = new Complex(nZ, dataLayoutZ, internalImZ, nY);
    fftZN.setUseSIMD(useSIMD);

    // Transforms along Z-data are repacked into a local work array.
    work = new double[2 * nZ * nY];
  }

  /**
   * Set the 2D transform to use SIMD instructions.
   *
   * @param useSIMD True to use SIMD instructions.
   */
  public void setUseSIMD(boolean useSIMD) {
    this.useSIMD = useSIMD;
    fftXY.setUseSIMD(useSIMD);
    fftZN.setUseSIMD(useSIMD);
  }

  /**
   * Set the 2D transform to pack FFTs.
   *
   * @param packFFTs True to pack FFTs.
   */
  public void setPackFFTs(boolean packFFTs) {
    this.packFFTs = packFFTs;
    fftXY.setPackFFTs(packFFTs);
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
   * @param i      the index along the X-axis.
   * @param j      the index along the Y-axis.
   * @param k      the index along the Z-axis.
   * @param nX     the number of points along the X-axis.
   * @param nY     the number of points along the Y-axis.
   * @param layout the data layout.
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
   * Compute the 3D FFT.
   *
   * @param input The input array must be of size 2 * nX * nY * nZ.
   */
  public void fft(final double[] input) {
    for (int z = 0; z < nZ; z++) {
      fftXY.fft(input, z * nextZ);
    }
    for (int x = 0, offset = 0; x < nX; x++) {
      selectYZPlane(offset, input);
      fftZN.fft(work, 0, ii);
      replaceYZPlane(offset, input);
      offset += nextX;
    }
  }

  /**
   * Compute the inverse 3D FFT.
   *
   * @param input The input array must be of size 2 * nX * nY * nZ.
   */
  public void ifft(final double[] input) {
    for (int offset = 0, x = 0; x < nX; x++) {
      selectYZPlane(offset, input);
      fftZN.ifft(work, 0, ii);
      replaceYZPlane(offset, input);
      offset += nextX;
    }
    for (int z = 0; z < nZ; z++) {
      fftXY.ifft(input, z * nextZ);
    }
  }

  /**
   * Perform a convolution.
   *
   * @param input The input array.
   */
  public void convolution(final double[] input) {
    for (int z = 0; z < nZ; z++) {
      fftXY.fft(input, z * nextZ);
    }
    for (int x = 0, offset = 0; x < nX; x++) {
      selectYZPlane(offset, input);
      fftZN.fft(work, 0, internalNextZ);
      recipConv(x, work);
      fftZN.ifft(work, 0, internalNextZ);
      replaceYZPlane(offset, input);
      offset += nextX;
    }
    for (int z = 0; z < nZ; z++) {
      fftXY.ifft(input, z * nextZ);
    }
  }

  /**
   * Setter for the field <code>recip</code>.
   *
   * @param recip an array of double.
   */
  public void setRecip(double[] recip) {
    // Reorder the reciprocal space data for convolution.
    // Input
    // trNextY = ii
    // trNextZ = nY*ii
    // real[y, z] = work[y*trNextY + z*trNextZ]
    // imag[y, z] = work[y*trNextY + z*trNextZ + internalImZ]
    int recipNextY = nX;
    int recipNextZ = nY * nX;
    int index = 0;
    for (int x = 0; x < nX; x++) {
      int dx = x;
      for (int z = 0; z < nZ; z++) {
        int dz = dx + z * recipNextZ;
        for (int y = 0; y < nY; y++) {
          int conv = y * recipNextY + dz;
          this.recip[index] = recip[conv];
          index++;
        }
      }
    }
  }

  /**
   * Select a YZ-plane at fixed X into a contiguous block of memory.
   *
   * @param offset The offset into the input array.
   * @param input  The input array.
   */
  private void selectYZPlane(int offset, double[] input) {
    // Input
    // real[x, y, z] = input[offset + y*nextY + z*nextZ]
    // imag[x, y, z] = input[offset + y*nextY + z*nextZ + im]
    // Collect a Y-Z plane at fixed X.
    // Output
    // trNextY = ii
    // trNextZ = nY*ii
    // real[y, z] = work[y*trNextY + z*trNextZ]
    // imag[y, z] = work[y*trNextY + z*trNextZ + internalImZ]
    int index = 0;
    for (int z = 0; z < nZ; z++) {
      int dz = offset + z * nextZ;
      for (int y = 0; y < nY; y++) {
        double real = input[y * nextY + dz];
        double imag = input[y * nextY + dz + im];
        work[index] = real;
        work[index + internalImZ] = imag;
        index += ii;
      }
    }
  }

  /**
   * Replace the Y-Z plane at fixed X back into the input array.
   *
   * @param offset The offset into the output array.
   * @param output The input array.
   */
  private void replaceYZPlane(int offset, double[] output) {
    // Input
    // trNextY = ii
    // trNextZ = nY*ii
    // real[y, z] = work[y*trNextY + z*trNextZ]
    // imag[y, z] = work[y*trNextY + z*trNextZ + internalImZ]
    // Output
    // real[x, y, z] = input[offset + y*nextY + z*nextZ]
    // imag[x, y, z] = input[offset + y*nextY + z*nextZ + im]
    int index = 0;
    for (int z = 0; z < nZ; z++) {
      int dzOut = offset + z * nextZ;
      for (int y = 0; y < nY; y++) {
        int dyOut = y * nextY;
        double real = work[index];
        double imag = work[index + internalImZ];
        index += ii;
        output[dyOut + dzOut] = real;
        output[dyOut + dzOut + im] = imag;
      }
    }
  }

  /**
   * Perform a multiplication by the reciprocal space data.
   *
   * @param x    The X-value for this Y-Z plane.
   * @param work The input array.
   */
  private void recipConv(int x, double[] work) {
    int index = 0;
    int rindex = x * (nY * nZ);
    for (int z = 0; z < nZ; z++) {
      for (int y = 0; y < nY; y++) {
        double r = recip[rindex++];
        work[index] *= r;
        work[index + internalImZ] *= r;
        index += ii;
      }
    }
  }

}
