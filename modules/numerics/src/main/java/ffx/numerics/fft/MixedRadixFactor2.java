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

import static jdk.incubator.vector.DoubleVector.SPECIES_128;
import static jdk.incubator.vector.DoubleVector.SPECIES_256;
import static jdk.incubator.vector.DoubleVector.SPECIES_512;
import static jdk.incubator.vector.DoubleVector.broadcast;
import static jdk.incubator.vector.DoubleVector.fromArray;

/**
 * The MixedRadixFactor2 class handles factors of 2 in the FFT.
 */
public class MixedRadixFactor2 extends MixedRadixFactor {

  /**
   * Create a new MixedRadixFactor2 instance.
   *
   * @param passConstants The pass constants.
   */
  public MixedRadixFactor2(PassConstants passConstants) {
    super(passConstants);
  }

  /**
   * Handle factors of 2.
   *
   * @param passData The pass data.
   */
  @Override
  protected void passScalar(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0_r = data[i];
      final double z0_i = data[i + im];
      final int idi = i + di;
      final double z1_r = data[idi];
      final double z1_i = data[idi + im];
      ret[j] = z0_r + z1_r;
      ret[j + im] = z0_i + z1_i;
      final double x_r = z0_r - z1_r;
      final double x_i = z0_i - z1_i;
      final int jdj = j + dj;
      ret[jdj] = x_r;
      ret[jdj + im] = x_i;
    }
    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      final double w_r = twids[0];
      final double w_i = -sign * twids[1];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0_r = data[i];
        final double z0_i = data[i + im];
        final int idi = i + di;
        final double z1_r = data[idi];
        final double z1_i = data[idi + im];
        ret[j] = z0_r + z1_r;
        ret[j + im] = z0_i + z1_i;
        final int jdj = j + dj;
        multiplyAndStore(z0_r - z1_r, z0_i - z1_i, w_r, w_i, ret, jdj, jdj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using SIMD vectors.
   *
   * @param passData The pass data.
   */
  @Override
  protected void passSIMD(PassData passData) {
    if (im == 1) {
      interleaved(passData);
    } else {
      blocked(passData);
    }
  }

  /**
   * Available SIMD sizes for Pass 2.
   */
  private static int[] simdSizes = {8, 4, 2};

  /**
   * Handle factors of 2 using the chosen SIMD vector.
   *
   * @param passData   The pass data.
   * @param simdLength The SIMD vector length.
   */
  private void interleaved(PassData passData, int simdLength) {
    // Use the preferred SIMD vector.
    switch (simdLength) {
      case 2:
        // 1 complex number per loop iteration.
        interleaved128(passData);
        break;
      case 4:
        // 2 complex numbers per loop iteration.
        interleaved256(passData);
        break;
      case 8:
        // 4 complex numbers per loop iteration.
        interleaved512(passData);
        break;
      default:
        passScalar(passData);
    }
  }

  /**
   * Handle factors of 2 using the chosen SIMD vector.
   *
   * @param passData   The pass data.
   * @param simdLength The SIMD vector length.
   */
  private void blocked(PassData passData, int simdLength) {
    // Use the preferred SIMD vector.
    switch (simdLength) {
      case 2:
        // 2 complex numbers per loop iteration.
        blocked128(passData);
        break;
      case 4:
        // 4 complex numbers per loop iteration.
        blocked256(passData);
        break;
      case 8:
        // 8 complex numbers per loop iteration.
        blocked512(passData);
        break;
      default:
        passScalar(passData);
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   *
   * @param passData The interleaved pass data.
   */
  private void interleaved(PassData passData) {
    if (innerLoopLimit % LOOP == 0) {
      // Use the preferred SIMD vector.
      interleaved(passData, LENGTH);
    } else {
      // If the inner loop limit is odd, use the scalar method unless the inner loop limit is 1.
      if (innerLoopLimit % 2 != 0 && innerLoopLimit != 1) {
        passScalar(passData);
        return;
      }
      // Fall back to a smaller SIMD vector that fits the inner loop limit.
      for (int size : simdSizes) {
        if (size >= LENGTH) {
          // Skip anything greater than or equal to the preferred SIMD vector size (which was too big).
          continue;
        }
        // Divide the SIMD size by two because for interleaved a single SIMD vectors stores both real and imaginary parts.
        if (innerLoopLimit % (size / 2) == 0) {
          interleaved(passData, size);
        }
      }
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   *
   * @param passData The pass blocked data.
   */
  private void blocked(PassData passData) {
    if (innerLoopLimit % BLOCK_LOOP == 0) {
      // Use the preferred SIMD vector.
      blocked(passData, LENGTH);
    } else {
      // If the inner loop limit is odd, use the scalar method unless the inner loop limit is 1.
      if (innerLoopLimit % 2 != 0) {
        passScalar(passData);
        return;
      }
      // Fall back to a smaller SIMD vector that fits the inner loop limit.
      for (int size : simdSizes) {
        if (size >= LENGTH) {
          // Skip anything greater than or equal to the preferred SIMD vector size (which was too big).
          continue;
        }
        // Divide the SIMD size by two because for interleaved a single SIMD vectors stores both real and imaginary parts.
        if (innerLoopLimit % size == 0) {
          blocked(passData, size);
        }
      }
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   */
  private void blocked128(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      DoubleVector
          z0_r = fromArray(SPECIES_128, data, i),
          z1_r = fromArray(SPECIES_128, data, i + di),
          z0_i = fromArray(SPECIES_128, data, i + im),
          z1_i = fromArray(SPECIES_128, data, i + di + im);
      z0_r.add(z1_r).intoArray(ret, j);
      z0_i.add(z1_i).intoArray(ret, j + im);
      z0_r.sub(z1_r).intoArray(ret, j + dj);
      z0_i.sub(z1_i).intoArray(ret, j + dj + im);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          w_r = broadcast(SPECIES_128, twids[0]),
          w_i = broadcast(SPECIES_128, -sign * twids[1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        DoubleVector
            z0_r = fromArray(SPECIES_128, data, i),
            z1_r = fromArray(SPECIES_128, data, i + di),
            z0_i = fromArray(SPECIES_128, data, i + im),
            z1_i = fromArray(SPECIES_128, data, i + di + im);
        z0_r.add(z1_r).intoArray(ret, j);
        z0_i.add(z1_i).intoArray(ret, j + im);
        DoubleVector x_r = z0_r.sub(z1_r), x_i = z0_i.sub(z1_i);
        w_r.mul(x_r).add(w_i.neg().mul(x_i)).intoArray(ret, j + dj);
        w_r.mul(x_i).add(w_i.mul(x_r)).intoArray(ret, j + dj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using the 256-bit SIMD vectors.
   */
  private void blocked256(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      DoubleVector
          z0_r = fromArray(SPECIES_256, data, i),
          z1_r = fromArray(SPECIES_256, data, i + di),
          z0_i = fromArray(SPECIES_256, data, i + im),
          z1_i = fromArray(SPECIES_256, data, i + di + im);
      z0_r.add(z1_r).intoArray(ret, j);
      z0_i.add(z1_i).intoArray(ret, j + im);
      z0_r.sub(z1_r).intoArray(ret, j + dj);
      z0_i.sub(z1_i).intoArray(ret, j + dj + im);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          w_r = broadcast(SPECIES_256, twids[0]),
          w_i = broadcast(SPECIES_256, -sign * twids[1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        DoubleVector
            z0_r = fromArray(SPECIES_256, data, i),
            z1_r = fromArray(SPECIES_256, data, i + di),
            z0_i = fromArray(SPECIES_256, data, i + im),
            z1_i = fromArray(SPECIES_256, data, i + di + im);
        z0_r.add(z1_r).intoArray(ret, j);
        z0_i.add(z1_i).intoArray(ret, j + im);
        DoubleVector x_r = z0_r.sub(z1_r), x_i = z0_i.sub(z1_i);
        w_r.mul(x_r).add(w_i.neg().mul(x_i)).intoArray(ret, j + dj);
        w_r.mul(x_i).add(w_i.mul(x_r)).intoArray(ret, j + dj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using the 512-bit SIMD vectors.
   */
  private void blocked512(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      DoubleVector
          z0_r = fromArray(SPECIES_512, data, i),
          z1_r = fromArray(SPECIES_512, data, i + di),
          z0_i = fromArray(SPECIES_512, data, i + im),
          z1_i = fromArray(SPECIES_512, data, i + di + im);
      z0_r.add(z1_r).intoArray(ret, j);
      z0_i.add(z1_i).intoArray(ret, j + im);
      z0_r.sub(z1_r).intoArray(ret, j + dj);
      z0_i.sub(z1_i).intoArray(ret, j + dj + im);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          w_r = broadcast(SPECIES_512, twids[0]),
          w_i = broadcast(SPECIES_512, -sign * twids[1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        DoubleVector
            z0_r = fromArray(SPECIES_512, data, i),
            z1_r = fromArray(SPECIES_512, data, i + di),
            z0_i = fromArray(SPECIES_512, data, i + im),
            z1_i = fromArray(SPECIES_512, data, i + di + im);
        z0_r.add(z1_r).intoArray(ret, j);
        z0_i.add(z1_i).intoArray(ret, j + im);
        DoubleVector x_r = z0_r.sub(z1_r), x_i = z0_i.sub(z1_i);
        w_r.mul(x_r).add(w_i.neg().mul(x_i)).intoArray(ret, j + dj);
        w_r.mul(x_i).add(w_i.mul(x_r)).intoArray(ret, j + dj + im);
      }
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   */
  private void interleaved128(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      DoubleVector
          z0 = fromArray(SPECIES_128, data, i),
          z1 = fromArray(SPECIES_128, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          wr = broadcast(SPECIES_128, twids[0]),
          wi = broadcast(SPECIES_128, -sign * twids[1]).mul(NEGATE_IM_128);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        DoubleVector
            z0 = fromArray(SPECIES_128, data, i),
            z1 = fromArray(SPECIES_128, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 2 using the 256-bit SIMD vectors.
   */
  private void interleaved256(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      DoubleVector
          z0 = fromArray(SPECIES_256, data, i),
          z1 = fromArray(SPECIES_256, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          wr = broadcast(SPECIES_256, twids[0]),
          wi = broadcast(SPECIES_256, -sign * twids[1]).mul(NEGATE_IM_256);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        DoubleVector
            z0 = fromArray(SPECIES_256, data, i),
            z1 = fromArray(SPECIES_256, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 2 using the 512-bit SIMD vectors.
   */
  private void interleaved512(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      DoubleVector
          z0 = fromArray(SPECIES_512, data, i),
          z1 = fromArray(SPECIES_512, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector
          wr = broadcast(SPECIES_512, twids[0]),
          wi = broadcast(SPECIES_512, -sign * twids[1]).mul(NEGATE_IM_512);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        DoubleVector
            z0 = fromArray(SPECIES_512, data, i),
            z1 = fromArray(SPECIES_512, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj);
      }
    }
  }

}
