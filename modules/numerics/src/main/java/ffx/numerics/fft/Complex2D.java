// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
 * Compute the 2D FFT of complex, double precision input of arbitrary dimensions
 * via 1D Mixed Radix FFTs.
 * <p>
 * For interleaved data, the location of the input point [x, y] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY]<br>
 * double imag = input[x*nextX + y*nextY + im]<br>
 * where <br>
 * nextX = 2
 * nextY = 2*nX
 * im = 1
 * <p>
 * For blocked data along x, the location of the input point [x, y] within the input array must be: <br>
 * double real = input[x*nextX + y*nextY]<br>
 * double imag = input[x*nextX + y*nextY + im]<br>
 * where for BLOCKED_X <br>
 * nextX = 1<br>
 * nextY = nX<br>
 * im = nX<br>
 * and for BLOCKED_XY <br>
 * nextX = 1<br>
 * nextY = nX*nY<br>
 * im = nX*nY<br>
 *
 * @author Michal J. Schnieders
 * @see Complex
 * @since 1.0
 */
public class Complex2D {

  /**
   * The 2D data layout.
   */
  private final DataLayout2D layout;
  /**
   * The external imaginary offset.
   */
  private final int externalIm;
  /**
   * The X-dimension.
   */
  private final int nX;
  /**
   * The Y-dimension.
   */
  private final int nY;
  /**
   * The next real value along the X dimension.
   */
  private final int nextX;
  /**
   * The next real value along the Y dimension.
   */
  private final int nextY;
  /**
   * Compute FFTs along X one at a time.
   */
  private final Complex fftX;
  /**
   * Compute FFTs along Y one at a time.
   */
  private final Complex fftY;
  /**
   * If true, pack FFTs along X (or Y) into a contiguous array to compute all FFTs along X (or Y) at once.
   */
  private boolean packFFTs;
  /**
   * If true, use SIMD instructions.
   */
  private boolean useSIMD;
  /**
   * Compute nY FFTs along the X dimension all at once.
   */
  private final Complex packedFFTX;
  /**
   * Compute nX FFTs along the Y dimension all at once.
   */
  private final Complex packedFFTY;
  /**
   * Working array for packed FFTs.
   */
  private final double[] packedData;
  /**
   * The offset between real values in the packed data.
   */
  private final int ii;
  /**
   * The offset between any real value and its corresponding imaginary value for the packed data.
   */
  private final int im;
  /**
   * The offset between real values along the X-dimension in the transposed packed data.
   */
  private final int trNextX;
  /**
   * The offset between real values along the Y-dimension in the transposed packed data.
   */
  private final int trNextY;

  /**
   * Create a new 2D Complex FFT for interleaved data.
   *
   * @param nX The number of points in the X dimension.
   * @param nY The number of points in the Y dimension.
   */
  public Complex2D(int nX, int nY) {
    this(nX, nY, DataLayout2D.INTERLEAVED, 1);
  }

  /**
   * Create a new 2D Complex FFT.
   *
   * @param nX       The number of points in the X dimension.
   * @param nY       The number of points in the Y dimension.
   * @param layout   The data layout.
   * @param imOffset The offset between real and imaginary values.
   */
  public Complex2D(int nX, int nY, DataLayout2D layout, int imOffset) {
    this.nX = nX;
    this.nY = nY;
    this.externalIm = imOffset;
    this.layout = layout;
    if (layout == DataLayout2D.INTERLEAVED) {
      im = 1;
      ii = 2;
      nextX = 2;
      nextY = 2 * nX;
      trNextY = 2;
      trNextX = 2 * nY;
    } else if (layout == DataLayout2D.BLOCKED_X) {
      im = nX;
      ii = 1;
      nextX = 1;
      nextY = 2 * nX;
      trNextY = 1;
      trNextX = 2 * nY;
    } else if (layout == DataLayout2D.BLOCKED_XY) {
      im = nX * nY;
      ii = 1;
      nextX = 1;
      nextY = nX;
      trNextY = 1;
      trNextX = nY;
    } else {
      throw new IllegalArgumentException("Unsupported data layout: " + layout);
    }

    if (this.externalIm != im) {
      throw new IllegalArgumentException("Unsupported im offset: " + imOffset);
    }

    // Do not use SIMD by default for now.
    useSIMD = false;
    String simd = System.getProperty("fft.useSIMD", Boolean.toString(useSIMD));
    try {
      useSIMD = Boolean.parseBoolean(simd);
    } catch (Exception e) {
      useSIMD = false;
    }

    packFFTs = false;
    String pack = System.getProperty("fft.packFFTs", Boolean.toString(packFFTs));
    try {
      packFFTs = Boolean.parseBoolean(pack);
    } catch (Exception e) {
      packFFTs = false;
    }

    DataLayout1D layout1D;
    if (layout == DataLayout2D.INTERLEAVED) {
      layout1D = DataLayout1D.INTERLEAVED;
    } else {
      layout1D = DataLayout1D.BLOCKED;
    }

    // Create 1D FFTs that will be used for computing 1 FFT at a time.
    fftX = new Complex(nX, layout1D, externalIm);
    fftX.setUseSIMD(useSIMD);
    fftY = new Complex(nY, layout1D, externalIm);
    fftY.setUseSIMD(useSIMD);

    // Create 1D FFTs that will be used for computing all FFTs at once.
    packedFFTY = new Complex(nY, layout1D, externalIm, nX);
    packedFFTY.setUseSIMD(useSIMD);
    packedFFTX = new Complex(nX, layout1D, im, nY);
    packedFFTX.setUseSIMD(useSIMD);
    packedData = new double[2 * nX * nY];
  }

  /**
   * Get the Data Layout.
   *
   * @return The layout.
   */
  public DataLayout2D getLayout() {
    return layout;
  }

  /**
   * Set the 2D transform to use SIMD instructions.
   *
   * @param useSIMD True to use SIMD instructions.
   */
  public void setUseSIMD(boolean useSIMD) {
    this.useSIMD = useSIMD;
    fftX.setUseSIMD(useSIMD);
    fftY.setUseSIMD(useSIMD);
    packedFFTX.setUseSIMD(useSIMD);
    packedFFTY.setUseSIMD(useSIMD);
  }

  /**
   * Set the 2D transform to pack FFTs into a contiguous array to compute all FFTs at once.
   *
   * @param packFFTs True to pack FFTs.
   */
  public void setPackFFTs(boolean packFFTs) {
    this.packFFTs = packFFTs;
  }

  /**
   * Compute the 2D FFT.
   *
   * @param input The input array must be of size 2 * nX * nY.
   * @param index The offset into the input array of the first element.
   */
  public void fft(final double[] input, int index) {
    if (!packFFTs) {
      for (int offset = index, y = 0; y < nY; y++, offset += nextY) {
        fftX.fft(input, offset, nextX);
      }
      for (int offset = index, x = 0; x < nX; x++, offset += nextX) {
        fftY.fft(input, offset, nextY);
      }
    } else {
      // For FFTs along Y, the input data is already contiguous.
      packedFFTY.fft(input, index, ii);
      transpose(input, index);
      packedFFTX.fft(packedData, 0, ii);
      untranspose(input, index);
    }
  }

  /**
   * Pack the input array for Fourier transforms along the X dimension.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + im]
   * Output order:
   * real(x,y) = packedData[y*trNextY + x*trNextX]
   * imag(x,y) = packedData[y*trNextY + x*trXextX + im]
   *
   * @param input  The input data.
   * @param offset The offset into the input data.
   */
  private void transpose(final double[] input, int offset) {
    int index = 0;
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * nextX;
      // Inner loop over the Y dimension (the number of FFTs).
      for (int y = 0; y < nY; y++) {
        double real = input[dx + y * nextY];
        double imag = input[dx + y * nextY + externalIm];
        // Contiguous storage into the packed array.
        packedData[index] = real;
        packedData[index + im] = imag;
        index += ii;
      }
    }
  }

  /**
   * Unpack the output array after Fourier transforms.
   * <p>
   * Input order:
   * real_xy = packedData[y*trNextY + x*trNextX]
   * imag_xy = packedData[y*trNextY + x*trXextX + im]
   * Output order:
   * real_xy = output[offset + x*nextX + y*nextY]
   * imag_xy = output[offset + x*nextX + y*nextY + im]
   *
   * @param output The output data.
   * @param offset The offset into the output data.
   */
  private void untranspose(final double[] output, int offset) {
    int index = offset;
    // Outer loop over the Y dimension.
    for (int y = 0; y < nY; y++) {
      int dy = y * trNextY;
      // Inner loop over the X dimension.
      for (int x = 0; x < nX; x++) {
        double real = packedData[dy + x * trNextX];
        double imag = packedData[dy + x * trNextX + im];
        // Contiguous storage into the output array.
        output[index] = real;
        output[index + externalIm] = imag;
        index += ii;
      }
    }
  }

  /**
   * Compute the 2D IFFT.
   *
   * @param input The input array must be of size 2 * nX * nY.
   * @param index The offset into the input array of the first element.
   */
  public void ifft(final double[] input, int index) {
    if (!packFFTs) {
      for (int offset = index, y = 0; y < nY; y++, offset += nextY) {
        fftX.ifft(input, offset, nextX);
      }
      for (int offset = index, x = 0; x < nX; x++, offset += nextX) {
        fftY.ifft(input, offset, nextY);
      }
    } else {
      // For FFTs along Y, the input data is already contiguous.
      packedFFTY.ifft(input, index, ii);
      transpose(input, index);
      packedFFTX.ifft(packedData, 0, ii);
      untranspose(input, index);
    }
  }
}
