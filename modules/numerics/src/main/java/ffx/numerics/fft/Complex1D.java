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
import java.util.logging.Logger;

/**
 * This algorithm factors a size n FFT into nX * nY,
 * computes nY inner FFTs of size nX and nX inner FFTs of size nY,
 * then combines the results to get the final answer.
 *
 * This is incomplete.
 */
public class Complex1D {

  private static final Logger logger = Logger.getLogger(Complex1D.class.getName());

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
   * Cached real twiddle factors for indices x*nY + y.
   */
  private final double[] twiddleRe;
  /**
   * Cached imaginary twiddle factors for indices x*nY + y.
   */
  private final double[] twiddleIm;
  /**
   * Use SIMD operators.
   */
  private boolean useSIMD;

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
      twiddleRe = null;
      twiddleIm = null;
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
      twiddleRe = new double[n];
      twiddleIm = new double[n];
      precomputeTwiddleFactors();
    }

    // Use SIMD by default.
    useSIMD = true;
    String simd = System.getProperty("fft.simd", Boolean.toString(useSIMD));
    try {
      useSIMD = Boolean.parseBoolean(simd);
    } catch (Exception e) {
      logger.info(" Invalid value for fft.simd: " + simd);
      useSIMD = false;
    }
    fftX.setUseSIMD(useSIMD);
    if (fftY != null) {
      fftY.setUseSIMD(useSIMD);
    }
  }

  /**
   * Configure use of SIMD operators.
   *
   * @param useSIMD True to use SIMD operators.
   */
  public void setUseSIMD(boolean useSIMD) {
    fftX.setUseSIMD(useSIMD);
    if (fftY != null) {
      fftY.setUseSIMD(useSIMD);
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
    if (nY == 1) {
      fftX.fft(data, offset, stride);
      return;
    }

    // STEP 1: We need to pack the data if the stride is greater than the native layout. Skip for now.
    if (stride != ii) {
      throw new UnsupportedOperationException(
          "Complex1D.fft currently requires stride == " + ii + " for this data layout.");
    }

    // STEP 2: Perform nX FFTs of size nY.
    fftY.fft(data, offset, stride);

    // STEP 3: Apply twiddle factors while transposing.
    transpose(data, offset);

    // STEP 4: Perform nY FFTs of size nX.
    fftX.fft(buffer, 0, ii);

    // STEP 5: Un-Transpose
    unTranspose(data, offset);
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
    if (nY == 1) {
      fftX.ifft(data, offset, stride);
      return;
    }

    // STEP 1: We need to pack the data if the stride is greater than the native layout. Skip for now.
    if (stride != ii) {
      throw new UnsupportedOperationException(
          "Complex1D.ifft currently requires stride == " + ii + " for this data layout.");
    }

    // STEP 2: Pack the frequency-domain data for nY inverse FFTs of size nX.
    transposeFrequencyDomain(data, offset);

    // STEP 3: Perform nY inverse FFTs of size nX.
    fftX.ifft(buffer, 0, ii);

    // STEP 4: Scatter back to the intermediate layout while applying inverse twiddle factors.
    unTransposeIntermediate(data, offset);

    // STEP 5: Perform nX inverse FFTs of size nY.
    fftY.ifft(data, offset, stride);
  }

  /**
   * Precompute the Cooley-Tukey twiddle factors between the inner FFT stages.
   */
  private void precomputeTwiddleFactors() {
    final double twoPiOverN = -2.0 * Math.PI / n;
    for (int x = 0; x < nX; x++) {
      int offset = x * nY;
      twiddleRe[offset] = 1.0;
      twiddleIm[offset] = 0.0;
      for (int y = 1; y < nY; y++) {
        double theta = twoPiOverN * x * y;
        twiddleRe[offset + y] = Math.cos(theta);
        twiddleIm[offset + y] = Math.sin(theta);
      }
    }
  }

  /**
   * Pack the twiddled, transposed data into the point-major layout used by batched Complex FFTs.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*nextX + y*nextY]
   * imag(x,y) = input[offset + x*nextX + y*nextY + im]
   * Output order:
   * real(x,y) = buffer[x*trNextX + y*trNextY]
   * imag(x,y) = buffer[x*trNextX + y*trNextY + internalIm]
   *
   * @param input  The input data.
   * @param offset The offset into the input data.
   */
  private void transpose(final double[] input, int offset) {
    int index = 0;
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * nextX;
      int twiddleOffset = x * nY;
      // Inner loop over the Y dimension (the number of FFTs).
      for (int y = 0; y < nY; y++) {
        double real = input[dx + y * nextY];
        double imag = input[dx + y * nextY + externalIm];
        if (y == 0) {
          buffer[index] = real;
          buffer[index + internalIm] = imag;
        } else {
          double cosine = twiddleRe[twiddleOffset + y];
          double sine = twiddleIm[twiddleOffset + y];
          buffer[index] = real * cosine - imag * sine;
          buffer[index + internalIm] = real * sine + imag * cosine;
        }
        index += ii;
      }
    }
  }

  /**
   * Pack contiguous frequency-domain input into the point-major layout used by batched inverse FFTs
   * of length nX.
   * <p>
   * Input order:
   * real(x,y) = input[offset + x*trNextX + y*trNextY]
   * imag(x,y) = input[offset + x*trNextX + y*trNextY + im]
   * Output order:
   * real(x,y) = buffer[x*trNextX + y*trNextY]
   * imag(x,y) = buffer[x*trNextX + y*trNextY + internalIm]
   *
   * @param input  The input data.
   * @param offset The offset into the input data.
   */
  private void transposeFrequencyDomain(final double[] input, int offset) {
    int index = 0;
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * trNextX;
      // Inner loop over the Y dimension.
      for (int y = 0; y < nY; y++) {
        int indexXY = dx + y * trNextY;
        buffer[index] = input[indexXY];
        buffer[index + internalIm] = input[indexXY + externalIm];
        index += ii;
      }
    }
  }

  /**
   * Unpack the point-major FFT results back into contiguous output order.
   * <p>
   * Input order:
   * real_xy = buffer[x*trNextX + y*trNextY]
   * imag_xy = buffer[x*trNextX + y*trNextY + internalIm]
   * Output order:
   * real_xy = output[offset + (x*nY + y)*ii]
   * imag_xy = output[offset + (x*nY + y)*ii + im]
   *
   * @param output The output data.
   * @param offset The offset into the output data.
   */
  private void unTranspose(final double[] output, int offset) {
    int outputIndex = offset;
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = x * trNextX;
      // Inner loop over the Y dimension.
      for (int y = 0; y < nY; y++) {
        int indexXY = dx + y * trNextY;
        output[outputIndex] = buffer[indexXY];
        output[outputIndex + externalIm] = buffer[indexXY + internalIm];
        outputIndex += ii;
      }
    }
  }

  /**
   * Unpack point-major inverse FFT results back into the intermediate x + nX*y layout used by the
   * final inverse FFTs of length nY, while applying the inverse twiddle factors.
   * <p>
   * Input order:
   * real_xy = buffer[x*trNextX + y*trNextY]
   * imag_xy = buffer[x*trNextX + y*trNextY + internalIm]
   * Output order:
   * real_xy = output[offset + x*nextX + y*nextY]
   * imag_xy = output[offset + x*nextX + y*nextY + im]
   *
   * @param output The output data.
   * @param offset The offset into the output data.
   */
  private void unTransposeIntermediate(final double[] output, int offset) {
    int index = 0;
    // Outer loop over the X dimension.
    for (int x = 0; x < nX; x++) {
      int dx = offset + x * nextX;
      int twiddleOffset = x * nY;
      // Inner loop over the Y dimension.
      for (int y = 0; y < nY; y++) {
        int indexXY = dx + y * nextY;
        double real = buffer[index];
        double imag = buffer[index + internalIm];
        if (y == 0) {
          output[indexXY] = real;
          output[indexXY + externalIm] = imag;
        } else {
          double cosine = twiddleRe[twiddleOffset + y];
          double sine = twiddleIm[twiddleOffset + y];
          output[indexXY] = real * cosine + imag * sine;
          output[indexXY + externalIm] = imag * cosine - real * sine;
        }
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
