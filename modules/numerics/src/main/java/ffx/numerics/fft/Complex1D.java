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

import java.util.Random;

/**
 * This algorithm factors a size n FFT into nX * nY,
 * computes nY inner FFTs of size nX and nX inner FFTs of size nY,
 * then combines the results to get the final answer.
 *
 * This is incomplete.
 */
public class Complex1D {

  /**
   * Overall FFT length.
   */
  private final int n;
  /**
   * The overall FFT length n = nX * nY.
   */
  private final int nX;
  /**
   * Compute FFTs of length nX.
   */
  private final Complex fftX;
  /**
   * The next real value along the X dimension.
   */
  private final int nextX;
  /**
   * The overall FFT length n = nX * nY.
   */
  private final int nY;
  /**
   * Compute FFTs of length height.
   */
  private final Complex fftY;
  /**
   * The next real value along the Y dimension.
   */
  private final int nextY;
  /**
   * The input data layout as interleaved or blocked.
   */
  private final DataLayout1D dataLayout;
  /**
   * Offset to the imaginary part of the input
   * This is 1 for interleaved data.
   * For blocked data, a typical offset is n in 1-dimension (or nX*nY + n in 2D).
   */
  private final int externalIm;
  /**
   * Internal data format.
   */
  private final DataLayout1D internalDataLayout;
  /**
   * Internal offset to the next real value (2 for interleaved, 1 for blocked).
   */
  private final int ii;
  /**
   * Internally use an interleaved data format.
   */
  private final int internalIm;
  /**
   * Internal buffer to store rearranged data.
   */
  private final double[] buffer;
  /**
   * The offset between real values along the X-dimension in the transposed packed data.
   */
  private final int trNextX;
  /**
   * The offset between real values along the Y-dimension in the transposed packed data.
   */
  private final int trNextY;

  /**
   * Construct a Complex instance for interleaved data of length n. Factorization of n is designed to use special
   * methods for small factors, and a general routine for large odd prime factors. Scratch memory is
   * created of length 2*n, which is reused each time a transform is computed.
   *
   * @param n Number of complex numbers (n .GT. 1).
   */
  public Complex1D(int n) {
    this(n, DataLayout1D.INTERLEAVED, 1);
  }

  /**
   * Construct a Complex instance for data of length n.
   * The offset to each imaginary part relative to the real part is given by im.
   * Factorization of n is designed to use special methods for small factors.
   * Scratch memory is created of length 2*n, which is reused each time a transform is computed.
   *
   * @param n          Number of complex numbers (n .GT. 1).
   * @param dataLayout Data layout (interleaved or blocked).
   * @param imOffset   Offset to the imaginary part of each complex number relative to its real part.
   */
  public Complex1D(int n, DataLayout1D dataLayout, int imOffset) {
    this.n = n;
    this.dataLayout = dataLayout;
    this.externalIm = imOffset;

    // Determine nX * nY = n.
    int[] factors = Complex.factor(n);

    // Only the width transform will be used.
    if (factors.length == 1) {
      this.nX = n;
      this.nY = 1;
      fftX = new Complex(nX, dataLayout, imOffset);
      // The variables below are not used in this case.
      fftY = null;
      internalDataLayout = null;
      buffer = null;
      internalIm = -1;
      ii = -1;
      nextX = -1;
      nextY = -1;
      trNextX = -1;
      trNextY = 1;
    } else {
      int n1 = 1;
      int n2 = 1;
      for (int i = 0; i < factors.length; i++) {
        // Even factors contribute to n1.
        // Odd factors contribute to n2.
        if (i % 2 == 0) {
          n1 *= factors[i];
        } else {
          n2 *= factors[i];
        }
      }
      nX = n1;
      nY = n2;
      buffer = new double[n * 2];
      internalDataLayout = dataLayout;
      if (dataLayout == DataLayout1D.INTERLEAVED) {
        nextX = 2;
        nextY = 2 * nX;
        trNextY = 2;
        trNextX = 2 * nY;
        ii = 2;
        internalIm = 1;
      } else {
        nextX = 1;
        nextY = nX;
        trNextY = 1;
        trNextX = nY;
        ii = 1;
        internalIm = nX * nY;
      }
      fftY = new Complex(nY, dataLayout, externalIm, nX);
      fftX = new Complex(nX, internalDataLayout, internalIm, nY);
    }
  }

  /**
   * Compute the Fast Fourier Transform of data leaving the result in data. The array data must
   * contain the data points in the following locations:
   *
   * <PRE>
   * Re(d[i]) = data[offset + stride*i]
   * Im(d[i]) = data[offset + stride*i + im]
   * </PRE>
   * <p>
   * where im is 1 for interleaved data or a constant set when the class was constructed.
   *
   * @param data   an array of double.
   * @param offset the offset to the beginning of the data.
   * @param stride the stride between data points.
   */
  public void fft(double[] data, int offset, int stride) {
    // STEP 1: We need to pack the data if the stride is greater than 2. Skip for now.
    // To Do.

    // STEP 2: Perform nX FFTs of size nY.
    fftY.fft(data, offset, stride);

    // STEP 3: Apply twiddle factors
    // To Do.

    // STEP 4: Transpose
    // transpose(data, offset);

    // STEP 5: Perform nY FFTs of size nX.
    fftX.fft(buffer, 0, 2);

    // STEP 6: Un-Transpose
    // unTranspose(data, offset);
  }

  /**
   * Compute the (un-normalized) inverse FFT of data, leaving it in place. The frequency domain data
   * must be in wrap-around order, and be stored in the following locations:
   *
   * <PRE>
   * Re(D[i]) = data[offset + stride*i]
   * Im(D[i]) = data[offset + stride*i + im]
   * </PRE>
   *
   * @param data   an array of double.
   * @param offset the offset to the beginning of the data.
   * @param stride the stride between data points.
   */
  public void ifft(double[] data, int offset, int stride) {
    // TO DO
  }

  /**
   * Transpose the input array for Fourier transforms of length width.
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
        buffer[index] = real;
        buffer[index + internalIm] = imag;
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
  private void unTranspose(final double[] output, int offset) {
    int index = offset;
    // Outer loop over the Y dimension.
    for (int y = 0; y < nY; y++) {
      int dy = y * trNextY;
      // Inner loop over the X dimension.
      for (int x = 0; x < nX; x++) {
        double real = buffer[dy + x * trNextX];
        double imag = buffer[dy + x * trNextX + internalIm];
        // Contiguous storage into the output array.
        output[index] = real;
        output[index + externalIm] = imag;
        index += ii;
      }
    }
  }

  /**
   * Test the Complex FFT.
   *
   * @param args an array of {@link java.lang.String} objects.
   * @throws java.lang.Exception if any.
   * @since 1.0
   */
  public static void main(String[] args) throws Exception {
    int dimNotFinal = 128;
    int reps = 5;
    try {
      dimNotFinal = Integer.parseInt(args[0]);
      if (dimNotFinal < 1) {
        dimNotFinal = 128;
      }
      reps = Integer.parseInt(args[1]);
      if (reps < 1) {
        reps = 5;
      }
    } catch (Exception e) {
      //
    }
    final int dim = dimNotFinal;
    System.out.printf("Initializing a 1D array of length %d.\n"
        + "The best timing out of %d repetitions will be used.%n", dim, reps);
    Complex1D complex = new Complex1D(dim);
    final double[] data = new double[dim * 2];
    Random random = new Random(1);
    for (int i = 0; i < dim; i++) {
      data[2 * i] = random.nextDouble();
    }
    double toSeconds = 0.000000001;
    long seqTime = Long.MAX_VALUE;
    for (int i = 0; i < reps; i++) {
      System.out.printf("Iteration %d%n", i + 1);
      long time = System.nanoTime();
      complex.fft(data, 0, 2);
      complex.ifft(data, 0, 2);
      time = (System.nanoTime() - time);
      System.out.printf("Sequential: %12.9f%n", toSeconds * time);
      if (time < seqTime) {
        seqTime = time;
      }
    }
    System.out.printf("Best Sequential Time:  %12.9f%n", toSeconds * seqTime);
  }

}
