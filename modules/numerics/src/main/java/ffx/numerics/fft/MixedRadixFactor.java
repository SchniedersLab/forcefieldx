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

import static java.lang.Math.fma;

public abstract class MixedRadixFactor {

  private static final double[] negateReal = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  private static final int[] shuffleMask = {1, 0, 3, 2, 5, 4, 7, 6};

  protected static final VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_PREFERRED;
  /**
   * Vector used to change the sign of the imaginary members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_IM;
  /**
   * Vector used to change the sign of the real members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_RE;
  /**
   * Shuffle used to swap real and imaginary members of the vector.
   */
  protected static final VectorShuffle<Double> SHUFFLE_RE_IM;
  /**
   * The number of contiguous elements that will be read from the input data array.
   */
  protected static final int LENGTH = DOUBLE_SPECIES.length();
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int LOOP = LENGTH / 2;
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int BLOCK_LOOP = LENGTH;

  /**
   * Vector used to change the sign of the imaginary members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_IM_128;
  /**
   * Vector used to change the sign of the real members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_RE_128;
  /**
   * Shuffle used to swap real and imaginary members of the vector.
   */
  protected static final VectorShuffle<Double> SHUFFLE_RE_IM_128;
  /**
   * The number of contiguous elements that will be read from the input data array.
   */
  protected static final int LENGTH_128 = DoubleVector.SPECIES_128.length();
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int LOOP_128 = LENGTH_128 / 2;
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int BLOCK_LOOP_128 = LENGTH_128;

  /**
   * Vector used to change the sign of the imaginary members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_IM_256;
  /**
   * Vector used to change the sign of the real members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_RE_256;
  /**
   * Shuffle used to swap real and imaginary members of the vector.
   */
  protected static final VectorShuffle<Double> SHUFFLE_RE_IM_256;
  /**
   * The number of contiguous elements that will be read from the input data array.
   */
  protected static final int LENGTH_256 = DoubleVector.SPECIES_256.length();
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int LOOP_256 = LENGTH_256 / 2;
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int BLOCK_LOOP_256 = LENGTH_256;

  /**
   * Vector used to change the sign of the imaginary members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_IM_512;
  /**
   * Vector used to change the sign of the real members of the vector via multiplication.
   */
  protected static final DoubleVector NEGATE_RE_512;
  /**
   * Shuffle used to swap real and imaginary members of the vector.
   */
  protected static final VectorShuffle<Double> SHUFFLE_RE_IM_512;
  /**
   * The number of contiguous elements that will be read from the input data array.
   */
  protected static final int LENGTH_512 = DoubleVector.SPECIES_512.length();
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int LOOP_512 = LENGTH_512 / 2;
  /**
   * The number of complex elements that will be processed in each inner loop iteration.
   * The number of elements to process in the inner loop must be evenly divisible by this loop increment.
   */
  protected static final int BLOCK_LOOP_512 = LENGTH_512;

  static {
    // Assume that 512 is the largest vector size.
    if (LENGTH > 8) {
      throw new IllegalStateException("Unsupported SIMD vector size: " + LENGTH);
    }

    NEGATE_RE_128 = DoubleVector.fromArray(DoubleVector.SPECIES_128, negateReal, 0);
    NEGATE_IM_128 = NEGATE_RE_128.mul(-1.0);
    SHUFFLE_RE_IM_128 = VectorShuffle.fromArray(DoubleVector.SPECIES_128, shuffleMask, 0);

    NEGATE_RE_256 = DoubleVector.fromArray(DoubleVector.SPECIES_256, negateReal, 0);
    NEGATE_IM_256 = NEGATE_RE_256.mul(-1.0);
    SHUFFLE_RE_IM_256 = VectorShuffle.fromArray(DoubleVector.SPECIES_256, shuffleMask, 0);

    NEGATE_RE_512 = DoubleVector.fromArray(DoubleVector.SPECIES_512, negateReal, 0);
    NEGATE_IM_512 = NEGATE_RE_512.mul(-1.0);
    SHUFFLE_RE_IM_512 = VectorShuffle.fromArray(DoubleVector.SPECIES_512, shuffleMask, 0);

    switch (LENGTH) {
      case 2:
        NEGATE_RE = NEGATE_RE_128;
        NEGATE_IM = NEGATE_IM_128;
        SHUFFLE_RE_IM = SHUFFLE_RE_IM_128;
        break;
      case 4:
        NEGATE_RE = NEGATE_RE_256;
        NEGATE_IM = NEGATE_IM_256;
        SHUFFLE_RE_IM = SHUFFLE_RE_IM_256;
        break;
      case 8:
        NEGATE_RE = NEGATE_RE_512;
        NEGATE_IM = NEGATE_IM_512;
        SHUFFLE_RE_IM = SHUFFLE_RE_IM_512;
        break;
      default:
        throw new IllegalStateException("Unsupported SIMD DoubleVector size: " + LENGTH);
    }
  }

  /**
   * The size of the input.
   */
  protected final int n;
  /**
   * The mixed radix factor.
   */
  protected final int factor;
  /**
   * The product of all factors applied so far.
   */
  protected final int product;
  /**
   * The outer loop limit (n / product).
   */
  protected final int outerLoopLimit;
  /**
   * The inner loop limit (product / factor).
   */
  protected final int innerLoopLimit;
  /**
   * The next input (n / factor).
   * This is the separation between the input data for each pass.
   */
  protected final int nextInput;
  /**
   * Equal to 2 * nextInput for interleaved complex data.
   * Equal to nextInput for separate real and imaginary arrays.
   */
  protected final int di;
  /**
   * Equal to 2 * innerLoopLimit for interleaved complex data.
   * Equal to innerLoopLimit for separate real and imaginary arrays.
   */
  protected final int dj;
  /**
   * The twiddle factors for this pass.
   */
  protected final double[][] twiddles;
  /**
   * The offset for the imaginary part of the input.
   * For interleaved complex data, this is 1.
   * For separate real and imaginary arrays, this is n (the size of the input).
   */
  protected final int im;
  /**
   * The increment for input data within the inner loop.
   * This is equal to 2 for interleaved complex data.
   * This is equal to 1 for separate real and imaginary arrays.
   */
  protected final int ii;
  /**
   * Increment for the inner loop.
   */
  protected final int jstep;

  protected double[] data;
  protected double[] ret;
  protected int sign;
  protected int i;
  protected int j;

  private boolean useSIMD = false;
  private int minSIMDLoopLength = 1;

  public MixedRadixFactor(PassConstants passConstants) {
    n = passConstants.n();
    factor = passConstants.factor();
    product = passConstants.product();
    im = passConstants.im();
    twiddles = passConstants.twiddles();
    outerLoopLimit = n / product;
    innerLoopLimit = product / factor;
    nextInput = n / factor;
    if (im == 1) {
      ii = 2;
      // For interleaved complex data, the di and dj offsets are doubled.
      di = 2 * nextInput;
      dj = 2 * innerLoopLimit;
    } else {
      ii = 1;
      // For separate real and imaginary arrays, the di and dj offsets
      // are the same as the next input and inner loop limit.
      di = nextInput;
      dj = innerLoopLimit;
    }
    jstep = (factor - 1) * dj;
  }

  /**
   * Apply the mixed radix factor.
   * SIMD operations will be used if enabled.
   *
   * @param passData the pass data.
   */
  protected void pass(PassData passData) {
    data = passData.in();
    ret = passData.out();
    sign = passData.sign();
    i = passData.inOffset();
    j = passData.outOffset();
    initRadixSpecificConstants();
    if (useSIMD && innerLoopLimit >= minSIMDLoopLength) {
      passSIMD();
    } else {
      passScalar();
    }
  }

  /**
   * Initialize radix factor specific constants.
   */
  protected void initRadixSpecificConstants() {
    // Override this method to initialize any constants specific to the radix factor.
  }

  /**
   * Apply the mixed radix factor using scalar operations.
   */
  protected abstract void passScalar();

  /**
   * Apply the mixed radix factor using SIMD operations.
   */
  protected abstract void passSIMD();

  /**
   * Minimum SIMD inner loop length.
   * For interleaved data, the default minimum SIMD loop length is 1 using AVX-128 (1 Real and 1 Imaginary per load).
   * For blocked data, the default minimum SIMD loop length is 1 using AVX-128 (2 Real or 2 Imaginary per load).
   * Setting this value above 1 reverts to scalar operations for short inner loop lengths.
   *
   * @param minSIMDLoopLength the minimum SIMD inner loop length.
   */
  public void setMinSIMDLoopLength(int minSIMDLoopLength) {
    this.minSIMDLoopLength = minSIMDLoopLength;
  }

  /**
   * Use SIMD instructions.
   *
   * @param useSIMD true to use SIMD instructions, false to use scalar instructions.
   */
  public void setUseSIMD(boolean useSIMD) {
    this.useSIMD = useSIMD;
  }

  /**
   * Multiply two complex numbers [x_r, x_i] and [w_r, w_i] and store the result.
   *
   * @param x_r the real part of the complex number.
   * @param x_i the imaginary part of the complex number.
   * @param w_r the real part of the twiddle factor.
   * @param w_i the imaginary part of the twiddle factor.
   * @param ret the array to store the result.
   * @param re  the real part index in the result array.
   * @param im  the imaginary part index in the result array.
   */
  protected static void multiplyAndStore(double x_r, double x_i, double w_r, double w_i, double[] ret, int re, int im) {
    ret[re] = fma(w_r, x_r, -w_i * x_i);
    ret[im] = fma(w_r, x_i, w_i * x_r);
  }

  /*
   * Multiply two complex numbers [x_r, x_i] and [w_r, w_i] and store the result.
   *
   * @param x_r the real part of the complex number.
   * @param x_i the imaginary part of the complex number.
   * @param w_r the real part of the twiddle factor.
   * @param w_i the imaginary part of the twiddle factor.
   * @param ret the array to store the result.
   * @param re the real part index in the result array.
   * @param im the imaginary part index in the result array.
   */
  protected static void multiplyAndStoreBlocked(DoubleVector x_r, DoubleVector x_i, DoubleVector w_r, DoubleVector w_i, double[] ret, int re, int im) {
    w_r.mul(x_r).add(w_i.neg().mul(x_i)).intoArray(ret, re);
    w_r.mul(x_i).add(w_i.mul(x_r)).intoArray(ret, im);
  }

  /*
   * Multiply a complex number x and with a complex twiddle w using interleaved data and the preferred SIMD vector size.
   *
   * @param x a vector of complex numbers [x0_r, x0_i].
   * @param wr real part of the twiddle factor in all lanes.
   * @param wi imag part of the twiddle factor in all lanes (the odd lanes should be negated by the caller).
   * @param ret the array to store the result.
   * @param index the index in the result array.
   */
  protected static void multiplyAndStoreInterleaved(DoubleVector x, DoubleVector w_r, DoubleVector w_i, double[] ret, int index) {
    x.mul(w_r).add(x.mul(w_i).rearrange(SHUFFLE_RE_IM)).intoArray(ret, index);
  }

  /*
   * Multiply a complex number x and with a complex twiddle w.
   *
   * @param x a vector of complex numbers [x0_r, x0_i].
   * @param wr real part of the twiddle factor in all lanes.
   * @param wi imag part of the twiddle factor in all lanes (the odd lanes should be negated by the caller).
   * @param ret the array to store the result.
   * @param index the index in the result array.
   */
  protected static void multiplyAndStore128(DoubleVector x, DoubleVector w_r, DoubleVector w_i, double[] ret, int index) {
    x.mul(w_r).add(x.mul(w_i).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, index);
  }

  /*
   * Multiply a complex number x and with a complex twiddle w.
   *
   * @param x a vector of complex numbers [x0_r, x0_i, x1_r, x1_i].
   * @param wr real part of the twiddle factor in all lanes.
   * @param wi imag part of the twiddle factor in all lanes (the odd lanes should be negated by the caller).
   * @param ret the array to store the result.
   * @param index the index in the result array.
   */
  protected static void multiplyAndStore256(DoubleVector x, DoubleVector w_r, DoubleVector w_i, double[] ret, int index) {
    x.mul(w_r).add(x.mul(w_i).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, index);
  }

  /*
   * Multiply a complex number x and with a complex twiddle w.
   *
   * @param x a vector of complex numbers [x0_r, x0_i, x1_r, x1_i, x2_r, x2_i, x3_r, x3_i].
   * @param wr real part of the twiddle factor in all lanes.
   * @param wi imag part of the twiddle factor in all lanes (the odd lanes should be negated by the caller).
   * @param ret the array to store the result.
   * @param index the index in the result array.
   */
  protected static void multiplyAndStore512(DoubleVector x, DoubleVector w_r, DoubleVector w_i, double[] ret, int index) {
    x.mul(w_r).add(x.mul(w_i).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, index);
  }
}
