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

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorShuffle;
import jdk.incubator.vector.VectorSpecies;

import static java.lang.Math.fma;

/**
 * Mixed radix factor is extended by the pass classes to apply the mixed radix factor.
 */
public abstract class MixedRadixFactor {

  private static final double[] negateReal = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  private static final int[] shuffleMask = {1, 0, 3, 2, 5, 4, 7, 6};

  /**
   * The preferred vector species for double precision.
   */
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
   * The number of complex elements that will be processed in each inner loop iteration for
   * interleaved data and the preferred SIMD species. The number of elements to process in the inner loop
   * must be evenly divisible by this loop increment.
   */
  protected static final int INTERLEAVED_LOOP = LENGTH / 2;
  /**
   * The number of complex elements that will be processed in each inner loop iteration for block data
   * and the preferred SIMD species. The number of elements to process in the inner loop must be
   * evenly divisible by this loop increment.
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
   * 1 complex element will be processed in each inner loop iteration for
   * interleaved data and a 128-bit SIMD width.
   */
  protected static final int INTERLEAVED_LOOP_128 = LENGTH_128 / 2;
  /**
   * 2 complex elements will be processed in each inner loop iteration for block data
   * and a 128-bit SIMD width. The number of elements to process in the inner loop must be
   * evenly divisible by 2.
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
   * 2 complex elements will be processed in each inner loop iteration for
   * interleaved data and the 256-bit SIMD width. The number of elements to process in the inner loop
   * must be evenly divisible by 2.
   */
  protected static final int INTERLEAVED_LOOP_256 = LENGTH_256 / 2;
  /**
   * 4 complex elements will be processed in each inner loop iteration for block data
   * and a 256-bit SIMD width. The number of elements to process in the inner loop must be
   * evenly divisible by 4.
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
   * 4 complex elements will be processed in each inner loop iteration for
   * interleaved data and the 512-bit SIMD width. The number of elements to process in the inner loop
   * must be evenly divisible by 4.
   */
  protected static final int INTERLEAVED_LOOP_512 = LENGTH_512 / 2;
  /**
   * 8 complex elements will be processed in each inner loop iteration for block data
   * and a 512-bit SIMD width. The number of elements to process in the inner loop must be
   * evenly divisible by 8.
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
   * The number of FFTs to process (default = 1).
   */
  protected final int nFFTs;
  /**
   * The imaginary offset.
   */
  protected final int im;
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
   * Equal to nextInput for blocked real and imaginary arrays.
   */
  protected final int di;
  /**
   * Equal to 2 * innerLoopLimit for interleaved complex data.
   * Equal to innerLoopLimit for blocked real and imaginary arrays.
   */
  protected final int dj;
  /**
   * The twiddle factors for this pass.
   */
  protected final double[][] twiddles;
  /**
   * The real twiddle factors for this pass.
   */
  protected final double[] wr;
  /**
   * The imaginary twiddle factors for this pass.
   */
  protected final double[] wi;
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

  /**
   * Constructor for the mixed radix factor.
   *
   * @param passConstants the pass constants.
   */
  public MixedRadixFactor(PassConstants passConstants) {
    n = passConstants.n();
    nFFTs = passConstants.nFFTs();
    im = passConstants.im();
    factor = passConstants.factor();
    product = passConstants.product();
    twiddles = passConstants.twiddles();
    outerLoopLimit = n / product;
    innerLoopLimit = (product / factor) * nFFTs;
    nextInput = (n / factor) * nFFTs;
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

    int f1 = factor - 1;
    wr = new double[outerLoopLimit * f1];
    wi = new double[outerLoopLimit * f1];
    for (int k = 0; k < outerLoopLimit; k++) {
      final double[] twids = twiddles[k];
      final int index = k * f1;
      for (int j = 0; j < f1; j++) {
        wr[index + j] = twids[2 * j];
        wi[index + j] = twids[2 * j + 1];
      }
    }
  }

  /**
   * Return a string representation of the mixed radix factor.
   *
   * @return a string representation of the mixed radix factor.
   */
  public String toString() {
    int length;
    if (im == 1) {
      length = INTERLEAVED_LOOP;
    } else {
      length = BLOCK_LOOP;
    }
    return "MixedRadixFactor {" +
        "\n n =" + n +
        "\n nFFTs =" + nFFTs +
        "\n im =" + im +
        "\n factor =" + factor +
        "\n product =" + product +
        "\n outerLoopLimit =" + outerLoopLimit +
        "\n innerLoopLimit =" + innerLoopLimit +
        "\n simd length =" + length +
        "\n nextInput =" + nextInput +
        "\n di =" + di +
        "\n dj =" + dj +
        "\n ii =" + ii +
        "\n jstep =" + jstep +
        "}\n";
  }

  /**
   * Apply the mixed radix factor using scalar operations.
   *
   * @param passData the pass data.
   */
  protected abstract void passScalar(PassData passData);

  /**
   * Apply the mixed radix factor using SIMD operations.
   *
   * @param passData the pass data.
   */
  protected abstract void passSIMD(PassData passData);

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
}
