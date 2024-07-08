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

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

import java.util.Random;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.fma;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Compute the FFT of complex, double precision data of arbitrary length n. This class uses a mixed
 * radix method and has special methods to handle factors [2, 3, 4, 5, 6, 7] and a general method for
 * larger prime factors.
 *
 * @author Michal J. Schnieders<br> Derived from: <br> Bruce R. Miller (bruce.miller@nist.gov) <br>
 * Contribution of the National Institute of Standards and Technology, not subject to copyright.
 * <br>
 * Derived from:<br> GSL (Gnu Scientific Library) FFT Code by Brian Gough
 * (bjg@network-theory.co.uk)
 * @see <ul>
 * <li><a href="http://dx.doi.org/10.1016/0021-9991(83)90013-X" target="_blank"> Clive
 * Temperton. Self-sorting mixed-radix fast fourier transforms. Journal of Computational
 * Physics, 52(1):1-23, 1983. </a>
 * <li><a href="http://www.jstor.org/stable/2003354" target="_blank"> J. W. Cooley and J. W.
 * Tukey, Mathematics of Computation 19 (90), 297 (1965) </a>
 * <li><a href="http://en.wikipedia.org/wiki/Fast_Fourier_transform" target="_blank">FFT at
 * Wikipedia </a>
 * </ul>
 * @since 1.0
 */
public class Complex {

  private static final Logger logger = Logger.getLogger(Complex.class.getName());
  // TINKER v. 5.0 factors to achieve exact numerical agreement.
  // private static final int[] availableFactors = {5, 4, 3, 2};
  // private static final int firstUnavailablePrime = 7;
  private static final int[] availableFactors = {7, 6, 5, 4, 3, 2};
  private static final int firstUnavailablePrime = 11;
  private static final double oneThird = 1.0 / 3.0;
  private static final double sqrt3_2 = sqrt(3.0) / 2.0;
  private static final double sqrt5_4 = sqrt(5.0) / 4.0;
  private static final double sinPI_5 = sin(PI / 5.0);
  private static final double sin2PI_5 = sin(2.0 * PI / 5.0);
  private static final double sin2PI_7 = sin(2.0 * PI / 7.0);
  private static final double sin4PI_7 = sin(4.0 * PI / 7.0);
  private static final double sin6PI_7 = sin(6.0 * PI / 7.0);
  private static final double cos2PI_7 = cos(2.0 * PI / 7.0);
  private static final double cos4PI_7 = cos(4.0 * PI / 7.0);
  private static final double cos6PI_7 = cos(6.0 * PI / 7.0);

  /**
   * The preferred SIMD vector size.
   */
  private static final VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_PREFERRED;
  // private static final VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_512;
  /**
   * Vector used to change the sign of the imaginary members of the vector via multiplication.
   */
  private static final DoubleVector negateIm;
  /**
   * Vector used to change the sign of the real members of the vector via multiplication.
   */
  private static final DoubleVector negateRe;
  /**
   * Shuffle used to swap real and imaginary members of the vector.
   */
  private static final VectorShuffle<Double> shuffleReIm;
  /**
   * The number of contiguous elements that will be read from the input data array.
   */
  private static final int SPECIES_LENGTH = DOUBLE_SPECIES.length();
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  private static final int LOOP_INCREMENT = SPECIES_LENGTH / 2;

  static {
    // Assume that 512 is the largest vector size.
    assert (SPECIES_LENGTH <= 8);
    double[] negateReal = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
    int[] shuffleMask = {1, 0, 3, 2, 5, 4, 7, 6};
    negateRe = DoubleVector.fromArray(DOUBLE_SPECIES, negateReal, 0);
    negateIm = negateRe.mul(-1.0);
    shuffleReIm = VectorShuffle.fromArray(DOUBLE_SPECIES, shuffleMask, 0);
  }

  /**
   * Number of complex numbers in the transform.
   */
  private final int n;
  /**
   * Factorization of n.
   */
  private final int[] factors;
  /**
   * Twiddle factors.
   */
  private final double[][][] twiddle;
  /**
   * Packing of non-contiguous data.
   */
  private final double[] packedData;
  /**
   * Scratch space for the transform.
   */
  private final double[] scratch;
  /**
   * References to the input and output data arrays.
   */
  private final PassData[] passData;
  /**
   * Constants for each pass.
   */
  private final PassConstants[] passConstants;
  /**
   * Sign of negative -1 is for forward transform. Sign of 1 is for inverse transform.
   */
  private int sign = -1;

  private boolean useSIMD;

  /**
   * Construct a Complex instance for data of length n. Factorization of n is designed to use special
   * methods for small factors, and a general routine for large odd prime factors. Scratch memory is
   * created of length 2*n, which is reused each time a transform is computed.
   *
   * @param n Number of complex numbers (n .GT. 1).
   */
  public Complex(int n) {
    assert (n > 1);
    this.n = n;
    factors = factor();
    twiddle = wavetable();
    packedData = new double[2 * n];
    scratch = new double[2 * n];
    passData = new PassData[2];

    // Compute the constants for each pass.
    passConstants = new PassConstants[factors.length];
    int product = 1;
    for (int i = 0; i < factors.length; i++) {
      final int factor = factors[i];
      product *= factor;
      // Create pass constants only for factors [2, 3, 4, 5, 6, 7].
      if (factor > 7) {
        passConstants[i] = null;
        continue;
      }
      final int outerLoopLimit = n / product;
      final int innerLoopLimit = product / factor;
      final int nextInput = n / factor;
      final int di = 2 * nextInput;
      final int dj = 2 * innerLoopLimit;
      passConstants[i] = new PassConstants(factor, outerLoopLimit, innerLoopLimit, nextInput, di, dj, twiddle[i]);
    }

    // Use SIMD by default only for AVX-512.
    useSIMD = SPECIES_LENGTH == 8;
    String simd = System.getProperty("fft.useSIMD", Boolean.toString(useSIMD));
    try {
      useSIMD = Boolean.parseBoolean(simd);
    } catch (Exception e) {
      logger.info(" Invalid value for fft.useSIMD: " + simd);
      useSIMD = false;
    }
  }

  /**
   * Configure use of SIMD operators.
   *
   * @param useSIMD True to use SIMD operators.
   */
  public void setUseSIMD(boolean useSIMD) {
    this.useSIMD = useSIMD;
  }

  /**
   * Check if a dimension is a preferred dimension.
   *
   * @param dim the dimension to check.
   * @return true if the dimension is a preferred dimension.
   */
  public static boolean preferredDimension(int dim) {
    if (dim < 2) {
      return false;
    }

    // Apply preferred factors.
    for (int factor : availableFactors) {
      while ((dim % factor) == 0) {
        dim /= factor;
      }
    }
    return dim <= 1;
  }

  /**
   * Getter for the field <code>factors</code>.
   *
   * @return an array of int.
   */
  public int[] getFactors() {
    return factors;
  }

  /**
   * Compute the Fast Fourier Transform of data leaving the result in data. The array data must
   * contain the data points in the following locations:
   *
   * <PRE>
   * Re(d[i]) = data[offset + stride*i]
   * Im(d[i]) = data[offset + stride*i+1]
   * </PRE>
   *
   * @param data   an array of double.
   * @param offset the offset to the beginning of the data.
   * @param stride the stride between data points.
   */
  public void fft(double[] data, int offset, int stride) {
    transformInternal(data, offset, stride, -1);
  }

  /**
   * Compute the (un-normalized) inverse FFT of data, leaving it in place. The frequency domain data
   * must be in wrap-around order, and be stored in the following locations:
   *
   * <PRE>
   * Re(D[i]) = data[offset + stride*i]
   * Im(D[i]) = data[offset + stride*i+1]
   * </PRE>
   *
   * @param data   an array of double.
   * @param offset the offset to the beginning of the data.
   * @param stride the stride between data points.
   */
  public void ifft(double[] data, int offset, int stride) {
    transformInternal(data, offset, stride, +1);
  }

  /**
   * Compute the normalized inverse FFT of data, leaving it in place. The frequency domain data must
   * be stored in the following locations:
   *
   * <PRE>
   * Re(D[i]) = data[offset + stride*i]
   * Im(D[i]) = data[offset + stride*i+1]
   * </PRE>
   *
   * @param data   an array of double.
   * @param offset the offset to the beginning of the data.
   * @param stride the stride between data points.
   */
  public void inverse(double[] data, int offset, int stride) {
    ifft(data, offset, stride);

    // Normalize inverse FFT with 1/n.
    double norm = normalization();
    for (int i = 0; i < n; i++) {
      final int index = offset + stride * i;
      data[index] *= norm;
      data[index + 1] *= norm;
    }
  }

  /**
   * Factor the data length into preferred factors (those with special methods), falling back to odd
   * primes that the general routine must handle.
   *
   * @return integer factors
   */
  private int[] factor() {
    if (n < 2) {
      return null;
    }
    Vector<Integer> v = new Vector<>();
    int nTest = n;

    // Use the preferred factors first
    for (int factor : availableFactors) {
      while ((nTest % factor) == 0) {
        nTest /= factor;
        v.add(factor);
      }
    }

    // Unavailable odd prime factors.
    int factor = firstUnavailablePrime;
    while (nTest > 1) {
      while ((nTest % factor) != 0) {
        factor += 2;
      }
      nTest /= factor;
      v.add(factor);
    }
    int product = 1;
    int nf = v.size();
    int[] ret = new int[nf];
    for (int i = 0; i < nf; i++) {
      ret[i] = v.get(i);
      product *= ret[i];
    }

    // Report a failed factorization.
    if (product != n) {
      StringBuilder sb = new StringBuilder(" FFT factorization failed for " + n + "\n");
      for (int i = 0; i < nf; i++) {
        sb.append(" ");
        sb.append(ret[i]);
      }
      sb.append("\n");
      sb.append(" Factor product = ");
      sb.append(product);
      sb.append("\n");
      logger.severe(sb.toString());
      System.exit(-1);
    } else {
      if (logger.isLoggable(Level.FINEST)) {
        StringBuilder sb = new StringBuilder(" FFT factorization for " + n + " = ");
        for (int i = 0; i < nf - 1; i++) {
          sb.append(ret[i]);
          sb.append(" * ");
        }
        sb.append(ret[nf - 1]);
        logger.finest(sb.toString());
      }
    }
    return ret;
  }

  /**
   * References to the input and output data arrays.
   *
   * @param in        Input data for the current pass.
   * @param inOffset  Offset into the input data.
   * @param out       Output data for the current pass.
   * @param outOffset Offset into output array.
   */
  private record PassData(double[] in, int inOffset, double[] out, int outOffset) {
    // Empty.
  }

  /**
   * Constant factors needed for each pass.
   *
   * @param factor         The factor.
   * @param outerLoopLimit The outer loop limit (n / product).
   * @param innerLoopLimit The inner loop limit (product / factor).
   * @param nextInput      The next input (n / factor).
   * @param di             Twice the next input to account for complex numbers.
   * @param dj             Twice the inner loop limit to account for complex numbers.
   * @param twiddles       The twiddle factors for this pass.
   */
  private record PassConstants(int factor, int outerLoopLimit, int innerLoopLimit, int nextInput,
                               int di, int dj, double[][] twiddles) {
    // Empty.
  }

  /**
   * Compute the Fast Fourier Transform of data leaving the result in data.
   *
   * @param data   data an array of double.
   * @param offset the offset to the beginning of the data.
   * @param stride the stride between data points.
   * @param sign   the sign to apply.
   */
  private void transformInternal(final double[] data,
                                 final int offset, final int stride, final int sign) {

    boolean packed = false;
    if (stride != 2) {
      // Pack non-contiguous (stride > 2) data into a contiguous array.
      packed = true;
      for (int i = 0, i2 = 0, index = offset; i < n; i++, i2 += 2, index += stride) {
        packedData[i2] = data[index];
        packedData[i2 + 1] = data[index + 1];
      }
      passData[0] = new PassData(packedData, 0, scratch, 0);
      passData[1] = new PassData(scratch, 0, packedData, 0);
    } else {
      passData[0] = new PassData(data, offset, scratch, 0);
      passData[1] = new PassData(scratch, 0, data, offset);
    }

    this.sign = sign;
    int product = 1;

    final int nfactors = factors.length;
    for (int i = 0; i < nfactors; i++) {
      final int pass = i % 2;
      final int factor = factors[i];
      product *= factor;
      if (useSIMD) {
        switch (factor) {
          case 2 -> pass2SIMD(passConstants[i], passData[pass]);
          case 3 -> pass3SIMD(passConstants[i], passData[pass]);
          case 4 -> pass4SIMD(passConstants[i], passData[pass]);
          case 5 -> pass5SIMD(passConstants[i], passData[pass]);
          case 6 -> pass6SIMD(passConstants[i], passData[pass]);
          case 7 -> pass7SIMD(passConstants[i], passData[pass]);
          default -> passOdd(factor, product, passData[pass], twiddle[i]);
        }
      } else {
        switch (factor) {
          case 2 -> pass2(passConstants[i], passData[pass]);
          case 3 -> pass3(passConstants[i], passData[pass]);
          case 4 -> pass4(passConstants[i], passData[pass]);
          case 5 -> pass5(passConstants[i], passData[pass]);
          case 6 -> pass6(passConstants[i], passData[pass]);
          case 7 -> pass7(passConstants[i], passData[pass]);
          default -> passOdd(factor, product, passData[pass], twiddle[i]);
        }
      }
    }

    // If the number of factors is odd, the final result is in the scratch array.
    if (nfactors % 2 == 1) {
      // Copy the scratch array to the data array.
      if (stride == 2) {
        arraycopy(scratch, 0, data, offset, 2 * n);
      } else {
        for (int i = 0, i2 = 0, index = offset; i < n; i++, i2 += 2, index += stride) {
          data[index] = scratch[i2];
          data[index + 1] = scratch[i2 + 1];
        }
      }
      // If the number of factors is even, the data may need to be unpacked.
    } else if (packed) {
      for (int i = 0, i2 = 0, index = offset; i < n; i++, i2 += 2, index += stride) {
        data[index] = packedData[i2];
        data[index + 1] = packedData[i2 + 1];
      }
    }
  }

  /**
   * Return the normalization factor. Multiply the elements of the back-transformed data to get the
   * normalized inverse.
   *
   * @return a double.
   */
  private double normalization() {
    return 1.0 / n;
  }

  /**
   * Handle factors of 2.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass2(PassConstants passConstants, PassData passData) {
    final int innerLoopLimit = passConstants.innerLoopLimit;
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      final double w_r = twids[0];
      final double w_i = -sign * twids[1];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        final double z0_r = data[i];
        final double z0_i = data[i + 1];
        final int idi = i + di;
        final double z1_r = data[idi];
        final double z1_i = data[idi + 1];
        ret[j] = z0_r + z1_r;
        ret[j + 1] = z0_i + z1_i;
        final double x_r = z0_r - z1_r;
        final double x_i = z0_i - z1_i;
        final int jdj = j + dj;
        ret[jdj] = fma(w_r, x_r, -w_i * x_i);
        ret[jdj + 1] = fma(w_r, x_i, w_i * x_r);
      }
    }
  }

  /**
   * Handle factors of 2.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass2SIMD(PassConstants passConstants, PassData passData) {
    final int innerLoopLimit = passConstants.innerLoopLimit;
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      // System.out.printf("Scalar %d product=%d innerLoopLimit=%d increment=%d%n",
      // factor, product, innerLoopLimit, LOOP_INCREMENT);
      pass2(passConstants, passData);
      return;
    }
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(negateIm);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT, i += SPECIES_LENGTH, j += SPECIES_LENGTH) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(shuffleReIm)).intoArray(ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 3.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass3(PassConstants passConstants, PassData passData) {
    final int factor = 3;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final double[][] twiddles = passConstants.twiddles;
    final double tau = sign * sqrt3_2;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final int jstep = (factor - 1) * dj;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      final double w1_r = twids[0];
      final double w1_i = -sign * twids[1];
      final double w2_r = twids[2];
      final double w2_i = -sign * twids[3];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        final double z0_r = data[i];
        final double z0_i = data[i + 1];
        int idi = i + di;
        final double z1_r = data[idi];
        final double z1_i = data[idi + 1];
        idi += di;
        final double z2_r = data[idi];
        final double z2_i = data[idi + 1];
        final double t1_r = z1_r + z2_r;
        final double t1_i = z1_i + z2_i;
        final double t2_r = fma(-0.5, t1_r, z0_r);
        final double t2_i = fma(-0.5, t1_i, z0_i);
        final double t3_r = tau * (z1_r - z2_r);
        final double t3_i = tau * (z1_i - z2_i);
        ret[j] = z0_r + t1_r;
        ret[j + 1] = z0_i + t1_i;
        double x_r = t2_r - t3_i;
        double x_i = t2_i + t3_r;
        int jdj = j + dj;
        ret[jdj] = fma(w1_r, x_r, -w1_i * x_i);
        ret[jdj + 1] = fma(w1_r, x_i, w1_i * x_r);
        x_r = t2_r + t3_i;
        x_i = t2_i - t3_r;
        jdj += dj;
        ret[jdj] = fma(w2_r, x_r, -w2_i * x_i);
        ret[jdj + 1] = fma(w2_r, x_i, w2_i * x_r);
      }
    }
  }

  /**
   * Handle factors of 3.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass3SIMD(PassConstants passConstants, PassData passData) {
    final int factor = 3;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      pass3(passConstants, passData);
      return;
    }
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final double tau = sign * sqrt3_2;
    final int di2 = 2 * di;
    final int dj2 = 2 * dj;
    final int jstep = (factor - 1) * dj;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(negateIm),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(negateIm);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT, i += SPECIES_LENGTH, j += SPECIES_LENGTH) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2);
        DoubleVector
            t1 = z1.add(z2),
            t2 = t1.mul(-0.5).add(z0),
            t3 = z1.sub(z2).mul(tau).rearrange(shuffleReIm);
        z0.add(t1).intoArray(ret, j);
        DoubleVector x = t2.add(t3.mul(negateRe));
        w1r.fma(x, x.mul(w1i).rearrange(shuffleReIm)).intoArray(ret, j + dj);
        x = t2.add(t3.mul(negateIm));
        w2r.fma(x, x.mul(w2i).rearrange(shuffleReIm)).intoArray(ret, j + dj2);
      }
    }
  }

  /**
   * Handle factors of 4.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass4(PassConstants passConstants, PassData passData) {
    final int factor = 4;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final int jstep = (factor - 1) * dj;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      final double w1_r = twids[0];
      final double w1_i = -sign * twids[1];
      final double w2_r = twids[2];
      final double w2_i = -sign * twids[3];
      final double w3_r = twids[4];
      final double w3_i = -sign * twids[5];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        final double z0_r = data[i];
        final double z0_i = data[i + 1];
        int idi = i + di;
        final double z1_r = data[idi];
        final double z1_i = data[idi + 1];
        idi += di;
        final double z2_r = data[idi];
        final double z2_i = data[idi + 1];
        idi += di;
        final double z3_r = data[idi];
        final double z3_i = data[idi + 1];
        final double t1_r = z0_r + z2_r;
        final double t1_i = z0_i + z2_i;
        final double t2_r = z1_r + z3_r;
        final double t2_i = z1_i + z3_i;
        final double t3_r = z0_r - z2_r;
        final double t3_i = z0_i - z2_i;
        final double t4_r = sign * (z1_r - z3_r);
        final double t4_i = sign * (z1_i - z3_i);
        ret[j] = t1_r + t2_r;
        ret[j + 1] = t1_i + t2_i;
        double x_r = t3_r - t4_i;
        double x_i = t3_i + t4_r;
        int jdj = j + dj;
        ret[jdj] = fma(w1_r, x_r, -w1_i * x_i);
        ret[jdj + 1] = fma(w1_r, x_i, w1_i * x_r);
        x_r = t1_r - t2_r;
        x_i = t1_i - t2_i;
        jdj += dj;
        ret[jdj] = fma(w2_r, x_r, -w2_i * x_i);
        ret[jdj + 1] = fma(w2_r, x_i, w2_i * x_r);
        x_r = t3_r + t4_i;
        x_i = t3_i - t4_r;
        jdj += dj;
        ret[jdj] = fma(w3_r, x_r, -w3_i * x_i);
        ret[jdj + 1] = fma(w3_r, x_i, w3_i * x_r);
      }
    }
  }

  /**
   * Handle factors of 4.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass4SIMD(PassConstants passConstants, PassData passData) {
    final int factor = 4;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      // System.out.printf("Scalar %d product=%d innerLoopLimit=%d increment=%d%n",
      // factor, product, innerLoopLimit, LOOP_INCREMENT);
      pass4(passConstants, passData);
      return;
    }
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final int di2 = 2 * di;
    final int di3 = 3 * di;
    final int dj2 = 2 * dj;
    final int dj3 = 3 * dj;
    final int jstep = (factor - 1) * dj;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(negateIm),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(negateIm),
          w3r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(negateIm);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT, i += SPECIES_LENGTH, j += SPECIES_LENGTH) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3);
        DoubleVector
            t1 = z0.add(z2),
            t2 = z1.add(z3),
            t3 = z0.sub(z2),
            t4 = z1.sub(z3).mul(sign).rearrange(shuffleReIm);
        t1.add(t2).intoArray(ret, j);
        DoubleVector x = t3.add(t4.mul(negateRe));
        w1r.fma(x, x.mul(w1i).rearrange(shuffleReIm)).intoArray(ret, j + dj);
        x = t1.sub(t2);
        w2r.fma(x, x.mul(w2i).rearrange(shuffleReIm)).intoArray(ret, j + dj2);
        x = t3.add(t4.mul(negateIm));
        w3r.fma(x, x.mul(w3i).rearrange(shuffleReIm)).intoArray(ret, j + dj3);
      }
    }
  }

  /**
   * Handle factors of 5.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass5(PassConstants passConstants, PassData passData) {
    final int factor = 5;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final double tau = sqrt5_4;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    final int jstep = (factor - 1) * dj;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      final double w1r = twids[0];
      final double w1i = -sign * twids[1];
      final double w2r = twids[2];
      final double w2i = -sign * twids[3];
      final double w3r = twids[4];
      final double w3i = -sign * twids[5];
      final double w4r = twids[6];
      final double w4i = -sign * twids[7];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        final double z0r = data[i];
        final double z0i = data[i + 1];
        int idi = i + di;
        final double z1r = data[idi];
        final double z1i = data[idi + 1];
        idi += di;
        final double z2r = data[idi];
        final double z2i = data[idi + 1];
        idi += di;
        final double z3r = data[idi];
        final double z3i = data[idi + 1];
        idi += di;
        final double z4r = data[idi];
        final double z4i = data[idi + 1];
        final double t1r = z1r + z4r;
        final double t1i = z1i + z4i;
        final double t2r = z2r + z3r;
        final double t2i = z2i + z3i;
        final double t3r = z1r - z4r;
        final double t3i = z1i - z4i;
        final double t4r = z2r - z3r;
        final double t4i = z2i - z3i;
        final double t5r = t1r + t2r;
        final double t5i = t1i + t2i;
        final double t6r = tau * (t1r - t2r);
        final double t6i = tau * (t1i - t2i);
        final double t7r = fma(-0.25, t5r, z0r);
        final double t7i = fma(-0.25, t5i, z0i);
        final double t8r = t7r + t6r;
        final double t8i = t7i + t6i;
        final double t9r = t7r - t6r;
        final double t9i = t7i - t6i;
        final double t10r = fma(sin2PI_5s, t3r, sinPI_5s * t4r);
        final double t10i = fma(sin2PI_5s, t3i, sinPI_5s * t4i);
        final double t11r = fma(-sin2PI_5s, t4r, sinPI_5s * t3r);
        final double t11i = fma(-sin2PI_5s, t4i, sinPI_5s * t3i);
        ret[j] = z0r + t5r;
        ret[j + 1] = z0i + t5i;
        double xr = t8r - t10i;
        double xi = t8i + t10r;
        int jdj = j + dj;
        ret[jdj] = fma(w1r, xr, -w1i * xi);
        ret[jdj + 1] = fma(w1r, xi, w1i * xr);
        xr = t9r - t11i;
        xi = t9i + t11r;
        jdj += dj;
        ret[jdj] = fma(w2r, xr, -w2i * xi);
        ret[jdj + 1] = fma(w2r, xi, w2i * xr);
        xr = t9r + t11i;
        xi = t9i - t11r;
        jdj += dj;
        ret[jdj] = fma(w3r, xr, -w3i * xi);
        ret[jdj + 1] = fma(w3r, xi, w3i * xr);
        xr = t8r + t10i;
        xi = t8i - t10r;
        jdj += dj;
        ret[jdj] = fma(w4r, xr, -w4i * xi);
        ret[jdj + 1] = fma(w4r, xi, w4i * xr);
      }
    }
  }

  /**
   * Handle factors of 5.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass5SIMD(PassConstants passConstants, PassData passData) {
    final int factor = 5;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      pass5(passConstants, passData);
      return;
    }
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final int di2 = 2 * di;
    final int di3 = 3 * di;
    final int di4 = 4 * di;
    final int dj2 = 2 * dj;
    final int dj3 = 3 * dj;
    final int dj4 = 4 * dj;
    final double tau = sqrt5_4;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    final int jstep = (factor - 1) * dj;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(negateIm),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(negateIm),
          w3r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(negateIm),
          w4r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[7]).mul(negateIm);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT, i += SPECIES_LENGTH, j += SPECIES_LENGTH) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4);
        DoubleVector
            t1 = z1.add(z4),
            t2 = z2.add(z3),
            t3 = z1.sub(z4),
            t4 = z2.sub(z3),
            t5 = t1.add(t2),
            t6 = t1.sub(t2).mul(tau),
            t7 = t5.mul(-0.25).add(z0),
            t8 = t7.add(t6),
            t9 = t7.sub(t6),
            t10 = t3.mul(sin2PI_5s).add(t4.mul(sinPI_5s)).rearrange(shuffleReIm),
            t11 = t4.mul(-sin2PI_5s).add(t3.mul(sinPI_5s)).rearrange(shuffleReIm);
        z0.add(t5).intoArray(ret, j);
        DoubleVector x = t8.add(t10.mul(negateRe));
        w1r.mul(x).add(w1i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj);
        x = t9.add(t11.mul(negateRe));
        w2r.mul(x).add(w2i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj2);
        x = t9.add(t11.mul(negateIm));
        w3r.mul(x).add(w3i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj3);
        x = t8.add(t10.mul(negateIm));
        w4r.mul(x).add(w4i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj4);
      }
    }
  }

  /**
   * Handle factors of 6.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass6(PassConstants passConstants, PassData passData) {
    final int factor = 6;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final int jstep = (factor - 1) * dj;
    final double tau = sign * sqrt3_2;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      final double w1r = twids[0];
      final double w1i = -sign * twids[1];
      final double w2r = twids[2];
      final double w2i = -sign * twids[3];
      final double w3r = twids[4];
      final double w3i = -sign * twids[5];
      final double w4r = twids[6];
      final double w4i = -sign * twids[7];
      final double w5r = twids[8];
      final double w5i = -sign * twids[9];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        final double z0r = data[i];
        final double z0i = data[i + 1];
        int idi = i + di;
        final double z1r = data[idi];
        final double z1i = data[idi + 1];
        idi += di;
        final double z2r = data[idi];
        final double z2i = data[idi + 1];
        idi += di;
        final double z3r = data[idi];
        final double z3i = data[idi + 1];
        idi += di;
        final double z4r = data[idi];
        final double z4i = data[idi + 1];
        idi += di;
        final double z5r = data[idi];
        final double z5i = data[idi + 1];
        final double ta1r = z2r + z4r;
        final double ta1i = z2i + z4i;
        final double ta2r = fma(-0.5, ta1r, z0r);
        final double ta2i = fma(-0.5, ta1i, z0i);
        final double ta3r = tau * (z2r - z4r);
        final double ta3i = tau * (z2i - z4i);
        final double a0r = z0r + ta1r;
        final double a0i = z0i + ta1i;
        final double a1r = ta2r - ta3i;
        final double a1i = ta2i + ta3r;
        final double a2r = ta2r + ta3i;
        final double a2i = ta2i - ta3r;
        final double tb1r = z5r + z1r;
        final double tb1i = z5i + z1i;
        final double tb2r = fma(-0.5, tb1r, z3r);
        final double tb2i = fma(-0.5, tb1i, z3i);
        final double tb3r = tau * (z5r - z1r);
        final double tb3i = tau * (z5i - z1i);
        final double b0r = z3r + tb1r;
        final double b0i = z3i + tb1i;
        final double b1r = tb2r - tb3i;
        final double b1i = tb2i + tb3r;
        final double b2r = tb2r + tb3i;
        final double b2i = tb2i - tb3r;
        ret[j] = a0r + b0r;
        ret[j + 1] = a0i + b0i;
        double xr = a1r - b1r;
        double xi = a1i - b1i;
        int jdj = j + dj;
        ret[jdj] = fma(w1r, xr, -w1i * xi);
        ret[jdj + 1] = fma(w1r, xi, w1i * xr);
        xr = a2r + b2r;
        xi = a2i + b2i;
        jdj += dj;
        ret[jdj] = fma(w2r, xr, -w2i * xi);
        ret[jdj + 1] = fma(w2r, xi, w2i * xr);
        xr = a0r - b0r;
        xi = a0i - b0i;
        jdj += dj;
        ret[jdj] = fma(w3r, xr, -w3i * xi);
        ret[jdj + 1] = fma(w3r, xi, w3i * xr);
        xr = a1r + b1r;
        xi = a1i + b1i;
        jdj += dj;
        ret[jdj] = fma(w4r, xr, -w4i * xi);
        ret[jdj + 1] = fma(w4r, xi, w4i * xr);
        xr = a2r - b2r;
        xi = a2i - b2i;
        jdj += dj;
        ret[jdj] = fma(w5r, xr, -w5i * xi);
        ret[jdj + 1] = fma(w5r, xi, w5i * xr);
      }
    }
  }

  /**
   * Handle factors of 6.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass6SIMD(PassConstants passConstants, PassData passData) {
    final int factor = 6;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      pass6(passConstants, passData);
      return;
    }
    final int outerLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final int di2 = 2 * di;
    final int di3 = 3 * di;
    final int di4 = 4 * di;
    final int di5 = 5 * di;
    final int dj2 = 2 * dj;
    final int dj3 = 3 * dj;
    final int dj4 = 4 * dj;
    final int dj5 = 5 * dj;
    final double tau = sign * sqrt3_2;
    final int jstep = (factor - 1) * dj;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(negateIm),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(negateIm),
          w3r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(negateIm),
          w4r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[7]).mul(negateIm),
          w5r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[8]),
          w5i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[9]).mul(negateIm);
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4),
            z5 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di5);
        DoubleVector
            ta1 = z2.add(z4),
            ta2 = ta1.mul(-0.5).add(z0),
            ta3 = z2.sub(z4).mul(tau).rearrange(shuffleReIm),
            a0 = z0.add(ta1),
            a1 = ta2.add(ta3.mul(negateRe)),
            a2 = ta2.add(ta3.mul(negateIm)),
            tb1 = z5.add(z1),
            tb2 = tb1.mul(-0.5).add(z3),
            tb3 = z5.sub(z1).mul(tau).rearrange(shuffleReIm),
            b0 = z3.add(tb1),
            b1 = tb2.add(tb3.mul(negateRe)),
            b2 = tb2.add(tb3.mul(negateIm));
        a0.add(b0).intoArray(ret, j);
        DoubleVector x = a1.sub(b1);
        w1r.mul(x).add(w1i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj);
        x = a2.add(b2);
        w2r.mul(x).add(w2i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj2);
        x = a0.sub(b0);
        w3r.mul(x).add(w3i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj3);
        x = a1.add(b1);
        w4r.mul(x).add(w4i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj4);
        x = a2.sub(b2);
        w5r.mul(x).add(w5i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj5);
      }
    }
  }

  /**
   * Handle factors of 7.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass7(PassConstants passConstants, PassData passData) {
    final int factor = 7;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    final int outLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final double c1 = cos2PI_7;
    final double c2 = cos4PI_7;
    final double c3 = cos6PI_7;
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v1 = (c1 + c2 + c3) * oneThird - 1.0;
    final double v2 = (2.0 * c1 - c2 - c3) * oneThird;
    final double v3 = (c1 - 2.0 * c2 + c3) * oneThird;
    final double v4 = (c1 + c2 - 2.0 * c3) * oneThird;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    final int jstep = (factor - 1) * dj;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    for (int k = 0; k < outLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      final double w1r = twids[0];
      final double w1i = -sign * twids[1];
      final double w2r = twids[2];
      final double w2i = -sign * twids[3];
      final double w3r = twids[4];
      final double w3i = -sign * twids[5];
      final double w4r = twids[6];
      final double w4i = -sign * twids[7];
      final double w5r = twids[8];
      final double w5i = -sign * twids[9];
      final double w6r = twids[10];
      final double w6i = -sign * twids[11];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        final double z0r = data[i];
        final double z0i = data[i + 1];
        int idi = i + di;
        final double z1r = data[idi];
        final double z1i = data[idi + 1];
        idi += di;
        final double z2r = data[idi];
        final double z2i = data[idi + 1];
        idi += di;
        final double z3r = data[idi];
        final double z3i = data[idi + 1];
        idi += di;
        final double z4r = data[idi];
        final double z4i = data[idi + 1];
        idi += di;
        final double z5r = data[idi];
        final double z5i = data[idi + 1];
        idi += di;
        final double z6r = data[idi];
        final double z6i = data[idi + 1];
        final double t0r = z1r + z6r;
        final double t0i = z1i + z6i;
        final double t1r = z1r - z6r;
        final double t1i = z1i - z6i;
        final double t2r = z2r + z5r;
        final double t2i = z2i + z5i;
        final double t3r = z2r - z5r;
        final double t3i = z2i - z5i;
        final double t4r = z4r + z3r;
        final double t4i = z4i + z3i;
        final double t5r = z4r - z3r;
        final double t5i = z4i - z3i;
        final double t6r = t2r + t0r;
        final double t6i = t2i + t0i;
        final double t7r = t5r + t3r;
        final double t7i = t5i + t3i;
        final double b0r = z0r + t6r + t4r;
        final double b0i = z0i + t6i + t4i;
        final double b1r = v1 * (t6r + t4r);
        final double b1i = v1 * (t6i + t4i);
        final double b2r = v2 * (t0r - t4r);
        final double b2i = v2 * (t0i - t4i);
        final double b3r = v3 * (t4r - t2r);
        final double b3i = v3 * (t4i - t2i);
        final double b4r = v4 * (t2r - t0r);
        final double b4i = v4 * (t2i - t0i);
        final double b5r = v5 * (t7r + t1r);
        final double b5i = v5 * (t7i + t1i);
        final double b6r = v6 * (t1r - t5r);
        final double b6i = v6 * (t1i - t5i);
        final double b7r = v7 * (t5r - t3r);
        final double b7i = v7 * (t5i - t3i);
        final double b8r = v8 * (t3r - t1r);
        final double b8i = v8 * (t3i - t1i);
        final double u0r = b0r + b1r;
        final double u0i = b0i + b1i;
        final double u1r = b2r + b3r;
        final double u1i = b2i + b3i;
        final double u2r = b4r - b3r;
        final double u2i = b4i - b3i;
        final double u3r = -b2r - b4r;
        final double u3i = -b2i - b4i;
        final double u4r = b6r + b7r;
        final double u4i = b6i + b7i;
        final double u5r = b8r - b7r;
        final double u5i = b8i - b7i;
        final double u6r = -b8r - b6r;
        final double u6i = -b8i - b6i;
        final double u7r = u0r + u1r;
        final double u7i = u0i + u1i;
        final double u8r = u0r + u2r;
        final double u8i = u0i + u2i;
        final double u9r = u0r + u3r;
        final double u9i = u0i + u3i;
        final double u10r = u4r + b5r;
        final double u10i = u4i + b5i;
        final double u11r = u5r + b5r;
        final double u11i = u5i + b5i;
        final double u12r = u6r + b5r;
        final double u12i = u6i + b5i;
        ret[j] = b0r;
        ret[j + 1] = b0i;
        double xr = u7r + u10i;
        double xi = u7i - u10r;
        int jdj = j + dj;
        ret[jdj] = fma(w1r, xr, -w1i * xi);
        ret[jdj + 1] = fma(w1r, xi, w1i * xr);
        xr = u9r + u12i;
        xi = u9i - u12r;
        jdj += dj;
        ret[jdj] = fma(w2r, xr, -w2i * xi);
        ret[jdj + 1] = fma(w2r, xi, w2i * xr);
        xr = u8r - u11i;
        xi = u8i + u11r;
        jdj += dj;
        ret[jdj] = fma(w3r, xr, -w3i * xi);
        ret[jdj + 1] = fma(w3r, xi, w3i * xr);
        xr = u8r + u11i;
        xi = u8i - u11r;
        jdj += dj;
        ret[jdj] = fma(w4r, xr, -w4i * xi);
        ret[jdj + 1] = fma(w4r, xi, w4i * xr);
        xr = u9r - u12i;
        xi = u9i + u12r;
        jdj += dj;
        ret[jdj] = fma(w5r, xr, -w5i * xi);
        ret[jdj + 1] = fma(w5r, xi, w5i * xr);
        xr = u7r - u10i;
        xi = u7i + u10r;
        jdj += dj;
        ret[jdj] = fma(w6r, xr, -w6i * xi);
        ret[jdj + 1] = fma(w6r, xi, w6i * xr);
      }
    }
  }

  /**
   * Handle factors of 7.
   *
   * @param passConstants the constants.
   * @param passData      the data.
   */
  private void pass7SIMD(PassConstants passConstants, PassData passData) {
    final int factor = 7;
    final int innerLoopLimit = passConstants.innerLoopLimit;
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      pass7(passConstants, passData);
      return;
    }
    final int outLoopLimit = passConstants.outerLoopLimit;
    final int di = passConstants.di;
    final int dj = passConstants.dj;
    final double[][] twiddles = passConstants.twiddles;
    final double c1 = cos2PI_7;
    final double c2 = cos4PI_7;
    final double c3 = cos6PI_7;
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v1 = (c1 + c2 + c3) * oneThird - 1.0;
    final double v2 = (2.0 * c1 - c2 - c3) * oneThird;
    final double v3 = (c1 - 2.0 * c2 + c3) * oneThird;
    final double v4 = (c1 + c2 - 2.0 * c3) * oneThird;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    final int di2 = 2 * di;
    final int di3 = 3 * di;
    final int di4 = 4 * di;
    final int di5 = 5 * di;
    final int di6 = 6 * di;
    final int dj2 = 2 * dj;
    final int dj3 = 3 * dj;
    final int dj4 = 4 * dj;
    final int dj5 = 5 * dj;
    final int dj6 = 6 * dj;
    final int jstep = (factor - 1) * dj;
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int i = passData.inOffset;
    int j = passData.outOffset;
    for (int k = 0; k < outLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(negateIm),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(negateIm),
          w3r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(negateIm),
          w4r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[7]).mul(negateIm),
          w5r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[8]),
          w5i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[9]).mul(negateIm),
          w6r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[10]),
          w6i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[11]).mul(negateIm);
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4),
            z5 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di5),
            z6 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di6);
        DoubleVector
            t0 = z1.add(z6),
            t1 = z1.sub(z6),
            t2 = z2.add(z5),
            t3 = z2.sub(z5),
            t4 = z4.add(z3),
            t5 = z4.sub(z3),
            t6 = t2.add(t0),
            t7 = t5.add(t3);
        DoubleVector
            b0 = z0.add(t6).add(t4),
            b1 = t6.add(t4).mul(v1),
            b2 = t0.sub(t4).mul(v2),
            b3 = t4.sub(t2).mul(v3),
            b4 = t2.sub(t0).mul(v4),
            b5 = t7.add(t1).mul(v5),
            b6 = t1.sub(t5).mul(v6),
            b7 = t5.sub(t3).mul(v7),
            b8 = t3.sub(t1).mul(v8);
        DoubleVector
            u0 = b0.add(b1),
            u1 = b2.add(b3),
            u2 = b4.sub(b3),
            u3 = b2.add(b4).neg(),
            u4 = b6.add(b7),
            u5 = b8.sub(b7),
            u6 = b8.add(b6).neg(),
            u7 = u0.add(u1),
            u8 = u0.add(u2),
            u9 = u0.add(u3),
            u10 = u4.add(b5).rearrange(shuffleReIm),
            u11 = u5.add(b5).rearrange(shuffleReIm),
            u12 = u6.add(b5).rearrange(shuffleReIm);
        b0.intoArray(ret, j);
        DoubleVector x = u7.add(u10.mul(negateIm));
        w1r.mul(x).add(w1i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj);
        x = u9.add(u12.mul(negateIm));
        w2r.mul(x).add(w2i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj2);
        x = u8.add(u11.mul(negateRe));
        w3r.mul(x).add(w3i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj3);
        x = u8.add(u11.mul(negateIm));
        w4r.mul(x).add(w4i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj4);
        x = u9.add(u12.mul(negateRe));
        w5r.mul(x).add(w5i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj5);
        x = u7.add(u10.mul(negateRe));
        w6r.mul(x).add(w6i.mul(x).rearrange(shuffleReIm)).intoArray(ret, j + dj6);
      }
    }
  }

  /**
   * Note that passOdd is only intended for odd factors (and fails for even factors).
   *
   * @param factor   Factor to apply.
   * @param product  Product to apply.
   * @param passData the data.
   * @param twiddles the twiddle factors.
   */
  private void passOdd(int factor, int product, PassData passData, double[][] twiddles) {
    final double[] data = passData.in;
    final int dataOffset = passData.inOffset;
    final double[] ret = passData.out;
    final int retOffset = passData.outOffset;
    final int m = n / factor;
    final int q = n / product;
    final int p_1 = product / factor;
    final int jump = (factor - 1) * p_1;
    for (int i = 0; i < m; i++) {
      ret[retOffset + 2 * i] = data[dataOffset + 2 * i];
      ret[retOffset + 2 * i + 1] = data[dataOffset + 2 * i + 1];
    }
    for (int e = 1; e < (factor - 1) / 2 + 1; e++) {
      for (int i = 0; i < m; i++) {
        int idx = i + e * m;
        int idxc = i + (factor - e) * m;
        ret[retOffset + 2 * idx] = data[dataOffset + 2 * idx] + data[dataOffset + 2 * idxc];
        ret[retOffset + 2 * idx + 1] =
            data[dataOffset + 2 * idx + 1] + data[dataOffset + 2 * idxc + 1];
        ret[retOffset + 2 * idxc] = data[dataOffset + 2 * idx] - data[dataOffset + 2 * idxc];
        ret[retOffset + 2 * idxc + 1] =
            data[dataOffset + 2 * idx + 1] - data[dataOffset + 2 * idxc + 1];
      }
    }
    for (int i = 0; i < m; i++) {
      data[dataOffset + 2 * i] = ret[retOffset + 2 * i];
      data[dataOffset + 2 * i + 1] = ret[retOffset + 2 * i + 1];
    }
    for (int e1 = 1; e1 < (factor - 1) / 2 + 1; e1++) {
      for (int i = 0; i < m; i++) {
        int i1 = retOffset + 2 * (i + e1 * m);
        data[dataOffset + 2 * i] += ret[i1];
        data[dataOffset + 2 * i + 1] += ret[i1 + 1];
      }
    }
    double[] twiddl = twiddles[q];
    for (int e = 1; e < (factor - 1) / 2 + 1; e++) {
      int idx = e;
      double wr, wi;
      int em = e * m;
      int ecm = (factor - e) * m;
      for (int i = 0; i < m; i++) {
        data[dataOffset + 2 * (i + em)] = ret[retOffset + 2 * i];
        data[dataOffset + 2 * (i + em) + 1] = ret[retOffset + 2 * i + 1];
        data[dataOffset + 2 * (i + ecm)] = ret[retOffset + 2 * i];
        data[dataOffset + 2 * (i + ecm) + 1] = ret[retOffset + 2 * i + 1];
      }
      for (int e1 = 1; e1 < (factor - 1) / 2 + 1; e1++) {
        if (idx == 0) {
          wr = 1;
          wi = 0;
        } else {
          wr = twiddl[2 * (idx - 1)];
          wi = -sign * twiddl[2 * (idx - 1) + 1];
        }
        for (int i = 0; i < m; i++) {
          int i1 = retOffset + 2 * (i + e1 * m);
          int i2 = retOffset + 2 * (i + (factor - e1) * m);
          double ap = wr * ret[i1];
          double am = wi * ret[i2 + 1];
          double bp = wr * ret[i1 + 1];
          double bm = wi * ret[i2];
          data[dataOffset + 2 * (i + em)] += (ap - am);
          data[dataOffset + 2 * (i + em) + 1] += (bp + bm);
          data[dataOffset + 2 * (i + ecm)] += (ap + am);
          data[dataOffset + 2 * (i + ecm) + 1] += (bp - bm);
        }
        idx += e;
        idx %= factor;
      }
    }
    /* k = 0 */
    for (int k1 = 0; k1 < p_1; k1++) {
      ret[retOffset + 2 * k1] = data[dataOffset + 2 * k1];
      ret[retOffset + 2 * k1 + 1] = data[dataOffset + 2 * k1 + 1];
    }
    for (int e1 = 1; e1 < factor; e1++) {
      for (int k1 = 0; k1 < p_1; k1++) {
        int i = retOffset + 2 * (k1 + e1 * p_1);
        int i1 = dataOffset + 2 * (k1 + e1 * m);
        ret[i] = data[i1];
        ret[i + 1] = data[i1 + 1];
      }
    }
    int i = p_1;
    int j = product;
    for (int k = 1; k < q; k++) {
      for (int k1 = 0; k1 < p_1; k1++) {
        ret[retOffset + 2 * j] = data[dataOffset + 2 * i];
        ret[retOffset + 2 * j + 1] = data[dataOffset + 2 * i + 1];
        i++;
        j++;
      }
      j += jump;
    }
    i = p_1;
    j = product;
    for (int k = 1; k < q; k++) {
      twiddl = twiddles[k];
      for (int k1 = 0; k1 < p_1; k1++) {
        for (int e1 = 1; e1 < factor; e1++) {
          int i1 = dataOffset + 2 * (i + e1 * m);
          double xr = data[i1];
          double xi = data[i1 + 1];
          double wr = twiddl[2 * (e1 - 1)];
          double wi = -sign * twiddl[2 * (e1 - 1) + 1];
          int i2 = retOffset + 2 * (j + e1 * p_1);
          ret[i2] = fma(wr, xr, -wi * xi);
          ret[i2 + 1] = fma(wr, xi, wi * xr);
        }
        i++;
        j++;
      }
      j += jump;
    }
  }

  /**
   * Compute twiddle factors. These are trigonometric constants that depend on the factoring of n.
   *
   * @return twiddle factors.
   */
  private double[][][] wavetable() {
    if (n < 2) {
      return null;
    }
    final double d_theta = -2.0 * PI / n;
    final double[][][] ret = new double[factors.length][][];
    int product = 1;
    for (int i = 0; i < factors.length; i++) {
      int factor = factors[i];
      int product_1 = product;
      product *= factor;
      final int q = n / product;
      ret[i] = new double[q + 1][2 * (factor - 1)];
      final double[][] twid = ret[i];
      for (int j = 0; j < factor - 1; j++) {
        twid[0][2 * j] = 1.0;
        twid[0][2 * j + 1] = 0.0;
      }
      for (int k = 1; k <= q; k++) {
        int m = 0;
        for (int j = 0; j < factor - 1; j++) {
          m += k * product_1;
          m %= n;
          final double theta = d_theta * m;
          twid[k][2 * j] = cos(theta);
          twid[k][2 * j + 1] = sin(theta);
        }
      }
    }
    return ret;
  }

  /**
   * Static DFT method used to test the FFT.
   *
   * @param in  input array.
   * @param out output array.
   */
  public static void dft(double[] in, double[] out) {
    int n = in.length / 2;
    for (int k = 0; k < n; k++) { // For each output element
      double sumReal = 0;
      double simImag = 0;
      for (int t = 0; t < n; t++) { // For each input element
        double angle = (2 * PI * t * k) / n;
        int re = 2 * t;
        int im = 2 * t + 1;
        sumReal = fma(in[re], cos(angle), sumReal);
        sumReal = fma(in[im], sin(angle), sumReal);
        simImag = fma(-in[re], sin(angle), simImag);
        simImag = fma(in[im], cos(angle), simImag);
      }
      int re = 2 * k;
      int im = 2 * k + 1;
      out[re] = sumReal;
      out[im] = simImag;
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
        dimNotFinal = 100;
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
    Complex complex = new Complex(dim);
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
