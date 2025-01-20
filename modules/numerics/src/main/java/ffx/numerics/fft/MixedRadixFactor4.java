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

import static jdk.incubator.vector.DoubleVector.SPECIES_128;
import static jdk.incubator.vector.DoubleVector.SPECIES_256;
import static jdk.incubator.vector.DoubleVector.SPECIES_512;
import static jdk.incubator.vector.DoubleVector.broadcast;
import static jdk.incubator.vector.DoubleVector.fromArray;

/**
 * The MixedRadixFactor4 class handles factors of 4 in the FFT.
 */
public class MixedRadixFactor4 extends MixedRadixFactor {

  private final int di2;
  private final int di3;
  private final int dj2;
  private final int dj3;

  /**
   * Construct a MixedRadixFactor4.
   * @param passConstants PassConstants.
   */
  public MixedRadixFactor4(PassConstants passConstants) {
    super(passConstants);
    di2 = 2 * di;
    di3 = 3 * di;
    dj2 = 2 * dj;
    dj3 = 3 * dj;
  }

  /**
   * Handle factors of 4.
   */
  @Override
  protected void passScalar(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0_r = data[i];
      final double z1_r = data[i + di];
      final double z2_r = data[i + di2];
      final double z3_r = data[i + di3];
      final double z0_i = data[i + im];
      final double z1_i = data[i + di + im];
      final double z2_i = data[i + di2 + im];
      final double z3_i = data[i + di3 + im];
      final double t1_r = z0_r + z2_r;
      final double t1_i = z0_i + z2_i;
      final double t2_r = z1_r + z3_r;
      final double t2_i = z1_i + z3_i;
      final double t3_r = z0_r - z2_r;
      final double t3_i = z0_i - z2_i;
      final double t4_r = sign * (z1_r - z3_r);
      final double t4_i = sign * (z1_i - z3_i);
      ret[j] = t1_r + t2_r;
      ret[j + im] = t1_i + t2_i;
      ret[j + dj] = t3_r - t4_i;
      ret[j + dj + im] = t3_i + t4_r;
      ret[j + dj2] = t1_r - t2_r;
      ret[j + dj2 + im] = t1_i - t2_i;
      ret[j + dj3] = t3_r + t4_i;
      ret[j + dj3 + im] = t3_i - t4_r;
    }
    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      final double w1_r = twids[0];
      final double w1_i = -sign * twids[1];
      final double w2_r = twids[2];
      final double w2_i = -sign * twids[3];
      final double w3_r = twids[4];
      final double w3_i = -sign * twids[5];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0_r = data[i];
        final double z1_r = data[i + di];
        final double z2_r = data[i + di2];
        final double z3_r = data[i + di3];
        final double z0_i = data[i + im];
        final double z1_i = data[i + di + im];
        final double z2_i = data[i + di2 + im];
        final double z3_i = data[i + di3 + im];
        final double t1_r = z0_r + z2_r;
        final double t1_i = z0_i + z2_i;
        final double t2_r = z1_r + z3_r;
        final double t2_i = z1_i + z3_i;
        final double t3_r = z0_r - z2_r;
        final double t3_i = z0_i - z2_i;
        final double t4_r = sign * (z1_r - z3_r);
        final double t4_i = sign * (z1_i - z3_i);
        ret[j] = t1_r + t2_r;
        ret[j + im] = t1_i + t2_i;
        multiplyAndStore(t3_r - t4_i, t3_i + t4_r, w1_r, w1_i, ret, j + dj, j + dj + im);
        multiplyAndStore(t1_r - t2_r, t1_i - t2_i, w2_r, w2_i, ret, j + dj2, j + dj2 + im);
        multiplyAndStore(t3_r + t4_i, t3_i - t4_r, w3_r, w3_i, ret, j + dj3, j + dj3 + im);
      }
    }
  }

  /**
   * Handle factors of 4 using SIMD vectors.
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
   * Available SIMD sizes for Pass 4.
   */
  private static int[] simdSizes = {8, 4, 2};

  /**
   * Handle factors of 4 using the chosen SIMD vector.
   * @param passData The pass data.
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
   * Handle factors of 4 using the chosen SIMD vector.
   * @param passData The pass data.
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
   * Handle factors of 4 using the 128-bit SIMD vectors.
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
   * Handle factors of 4 using the 128-bit SIMD vectors.
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
   * Handle factors of 4 using the 128-bit SIMD vectors.
   */
  private void blocked128(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      final DoubleVector
          z0_r = fromArray(SPECIES_128, data, i),
          z0_i = fromArray(SPECIES_128, data, i + im),
          z1_r = fromArray(SPECIES_128, data, i + di),
          z1_i = fromArray(SPECIES_128, data, i + di + im),
          z2_r = fromArray(SPECIES_128, data, i + di2),
          z2_i = fromArray(SPECIES_128, data, i + di2 + im),
          z3_r = fromArray(SPECIES_128, data, i + di3),
          z3_i = fromArray(SPECIES_128, data, i + di3 + im);
      final DoubleVector
          t1_r = z0_r.add(z2_r),
          t1_i = z0_i.add(z2_i),
          t2_r = z1_r.add(z3_r),
          t2_i = z1_i.add(z3_i),
          t3_r = z0_r.sub(z2_r),
          t3_i = z0_i.sub(z2_i),
          t4_r = z1_r.sub(z3_r).mul(sign),
          t4_i = z1_i.sub(z3_i).mul(sign);
      t1_r.add(t2_r).intoArray(ret, j);
      t1_i.add(t2_i).intoArray(ret, j + im);
      t3_r.sub(t4_i).intoArray(ret, j + dj);
      t3_i.add(t4_r).intoArray(ret, j + dj + im);
      t1_r.sub(t2_r).intoArray(ret, j + dj2);
      t1_i.sub(t2_i).intoArray(ret, j + dj2 + im);
      t3_r.add(t4_i).intoArray(ret, j + dj3);
      t3_i.sub(t4_r).intoArray(ret, j + dj3 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(SPECIES_128, twids[0]),
          w1i = broadcast(SPECIES_128, -sign * twids[1]),
          w2r = broadcast(SPECIES_128, twids[2]),
          w2i = broadcast(SPECIES_128, -sign * twids[3]),
          w3r = broadcast(SPECIES_128, twids[4]),
          w3i = broadcast(SPECIES_128, -sign * twids[5]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        final DoubleVector
            z0_r = fromArray(SPECIES_128, data, i),
            z0_i = fromArray(SPECIES_128, data, i + im),
            z1_r = fromArray(SPECIES_128, data, i + di),
            z1_i = fromArray(SPECIES_128, data, i + di + im),
            z2_r = fromArray(SPECIES_128, data, i + di2),
            z2_i = fromArray(SPECIES_128, data, i + di2 + im),
            z3_r = fromArray(SPECIES_128, data, i + di3),
            z3_i = fromArray(SPECIES_128, data, i + di3 + im);
        final DoubleVector
            t1_r = z0_r.add(z2_r),
            t1_i = z0_i.add(z2_i),
            t2_r = z1_r.add(z3_r),
            t2_i = z1_i.add(z3_i),
            t3_r = z0_r.sub(z2_r),
            t3_i = z0_i.sub(z2_i),
            t4_r = z1_r.sub(z3_r).mul(sign),
            t4_i = z1_i.sub(z3_i).mul(sign);
        t1_r.add(t2_r).intoArray(ret, j);
        t1_i.add(t2_i).intoArray(ret, j + im);
        DoubleVector
            x1_r = t3_r.sub(t4_i), x1_i = t3_i.add(t4_r),
            x2_r = t1_r.sub(t2_r), x2_i = t1_i.sub(t2_i),
            x3_r = t3_r.add(t4_i), x3_i = t3_i.sub(t4_r);
        w1r.mul(x1_r).add(w1i.neg().mul(x1_i)).intoArray(ret, j + dj);
        w2r.mul(x2_r).add(w2i.neg().mul(x2_i)).intoArray(ret, j + dj2);
        w3r.mul(x3_r).add(w3i.neg().mul(x3_i)).intoArray(ret, j + dj3);
        w1r.mul(x1_i).add(w1i.mul(x1_r)).intoArray(ret, j + dj + im);
        w2r.mul(x2_i).add(w2i.mul(x2_r)).intoArray(ret, j + dj2 + im);
        w3r.mul(x3_i).add(w3i.mul(x3_r)).intoArray(ret, j + dj3 + im);
      }
    }
  }

  /**
   * Handle factors of 4 using the 256-bit SIMD vectors.
   */
  private void blocked256(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      final DoubleVector
          z0_r = fromArray(SPECIES_256, data, i),
          z1_r = fromArray(SPECIES_256, data, i + di),
          z2_r = fromArray(SPECIES_256, data, i + di2),
          z3_r = fromArray(SPECIES_256, data, i + di3),
          z0_i = fromArray(SPECIES_256, data, i + im),
          z1_i = fromArray(SPECIES_256, data, i + di + im),
          z2_i = fromArray(SPECIES_256, data, i + di2 + im),
          z3_i = fromArray(SPECIES_256, data, i + di3 + im);
      final DoubleVector
          t1_r = z0_r.add(z2_r),
          t1_i = z0_i.add(z2_i),
          t2_r = z1_r.add(z3_r),
          t2_i = z1_i.add(z3_i),
          t3_r = z0_r.sub(z2_r),
          t3_i = z0_i.sub(z2_i),
          t4_r = z1_r.sub(z3_r).mul(sign),
          t4_i = z1_i.sub(z3_i).mul(sign);
      t1_r.add(t2_r).intoArray(ret, j);
      t3_r.sub(t4_i).intoArray(ret, j + dj);
      t1_r.sub(t2_r).intoArray(ret, j + dj2);
      t3_r.add(t4_i).intoArray(ret, j + dj3);
      t1_i.add(t2_i).intoArray(ret, j + im);
      t3_i.add(t4_r).intoArray(ret, j + dj + im);
      t1_i.sub(t2_i).intoArray(ret, j + dj2 + im);
      t3_i.sub(t4_r).intoArray(ret, j + dj3 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(SPECIES_256, twids[0]),
          w1i = broadcast(SPECIES_256, -sign * twids[1]),
          w2r = broadcast(SPECIES_256, twids[2]),
          w2i = broadcast(SPECIES_256, -sign * twids[3]),
          w3r = broadcast(SPECIES_256, twids[4]),
          w3i = broadcast(SPECIES_256, -sign * twids[5]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        final DoubleVector
            z0_r = fromArray(SPECIES_256, data, i),
            z1_r = fromArray(SPECIES_256, data, i + di),
            z2_r = fromArray(SPECIES_256, data, i + di2),
            z3_r = fromArray(SPECIES_256, data, i + di3),
            z0_i = fromArray(SPECIES_256, data, i + im),
            z1_i = fromArray(SPECIES_256, data, i + di + im),
            z2_i = fromArray(SPECIES_256, data, i + di2 + im),
            z3_i = fromArray(SPECIES_256, data, i + di3 + im);
        final DoubleVector
            t1_r = z0_r.add(z2_r),
            t1_i = z0_i.add(z2_i),
            t2_r = z1_r.add(z3_r),
            t2_i = z1_i.add(z3_i),
            t3_r = z0_r.sub(z2_r),
            t3_i = z0_i.sub(z2_i),
            t4_r = z1_r.sub(z3_r).mul(sign),
            t4_i = z1_i.sub(z3_i).mul(sign);
        t1_r.add(t2_r).intoArray(ret, j);
        t1_i.add(t2_i).intoArray(ret, j + im);
        DoubleVector
            x1_r = t3_r.sub(t4_i), x1_i = t3_i.add(t4_r),
            x2_r = t1_r.sub(t2_r), x2_i = t1_i.sub(t2_i),
            x3_r = t3_r.add(t4_i), x3_i = t3_i.sub(t4_r);
        w1r.mul(x1_r).add(w1i.neg().mul(x1_i)).intoArray(ret, j + dj);
        w2r.mul(x2_r).add(w2i.neg().mul(x2_i)).intoArray(ret, j + dj2);
        w3r.mul(x3_r).add(w3i.neg().mul(x3_i)).intoArray(ret, j + dj3);
        w1r.mul(x1_i).add(w1i.mul(x1_r)).intoArray(ret, j + dj + im);
        w2r.mul(x2_i).add(w2i.mul(x2_r)).intoArray(ret, j + dj2 + im);
        w3r.mul(x3_i).add(w3i.mul(x3_r)).intoArray(ret, j + dj3 + im);
      }
    }
  }

  /**
   * Handle factors of 4 using the 512-bit SIMD vectors.
   */
  private void blocked512(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      final DoubleVector
          z0_r = fromArray(SPECIES_512, data, i),
          z1_r = fromArray(SPECIES_512, data, i + di),
          z2_r = fromArray(SPECIES_512, data, i + di2),
          z3_r = fromArray(SPECIES_512, data, i + di3),
          z0_i = fromArray(SPECIES_512, data, i + im),
          z1_i = fromArray(SPECIES_512, data, i + di + im),
          z2_i = fromArray(SPECIES_512, data, i + di2 + im),
          z3_i = fromArray(SPECIES_512, data, i + di3 + im);
      final DoubleVector
          t1_r = z0_r.add(z2_r),
          t1_i = z0_i.add(z2_i),
          t2_r = z1_r.add(z3_r),
          t2_i = z1_i.add(z3_i),
          t3_r = z0_r.sub(z2_r),
          t3_i = z0_i.sub(z2_i),
          t4_r = z1_r.sub(z3_r).mul(sign),
          t4_i = z1_i.sub(z3_i).mul(sign);
      t1_r.add(t2_r).intoArray(ret, j);
      t3_r.sub(t4_i).intoArray(ret, j + dj);
      t1_r.sub(t2_r).intoArray(ret, j + dj2);
      t3_r.add(t4_i).intoArray(ret, j + dj3);
      t1_i.add(t2_i).intoArray(ret, j + im);
      t3_i.add(t4_r).intoArray(ret, j + dj + im);
      t1_i.sub(t2_i).intoArray(ret, j + dj2 + im);
      t3_i.sub(t4_r).intoArray(ret, j + dj3 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(SPECIES_512, twids[0]),
          w1i = broadcast(SPECIES_512, -sign * twids[1]),
          w2r = broadcast(SPECIES_512, twids[2]),
          w2i = broadcast(SPECIES_512, -sign * twids[3]),
          w3r = broadcast(SPECIES_512, twids[4]),
          w3i = broadcast(SPECIES_512, -sign * twids[5]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        final DoubleVector
            z0_r = fromArray(SPECIES_512, data, i),
            z1_r = fromArray(SPECIES_512, data, i + di),
            z2_r = fromArray(SPECIES_512, data, i + di2),
            z3_r = fromArray(SPECIES_512, data, i + di3),
            z0_i = fromArray(SPECIES_512, data, i + im),
            z1_i = fromArray(SPECIES_512, data, i + di + im),
            z2_i = fromArray(SPECIES_512, data, i + di2 + im),
            z3_i = fromArray(SPECIES_512, data, i + di3 + im);
        final DoubleVector
            t1_r = z0_r.add(z2_r),
            t1_i = z0_i.add(z2_i),
            t2_r = z1_r.add(z3_r),
            t2_i = z1_i.add(z3_i),
            t3_r = z0_r.sub(z2_r),
            t3_i = z0_i.sub(z2_i),
            t4_r = z1_r.sub(z3_r).mul(sign),
            t4_i = z1_i.sub(z3_i).mul(sign);
        t1_r.add(t2_r).intoArray(ret, j);
        t1_i.add(t2_i).intoArray(ret, j + im);
        DoubleVector
            x1_r = t3_r.sub(t4_i), x1_i = t3_i.add(t4_r),
            x2_r = t1_r.sub(t2_r), x2_i = t1_i.sub(t2_i),
            x3_r = t3_r.add(t4_i), x3_i = t3_i.sub(t4_r);
        w1r.mul(x1_r).add(w1i.neg().mul(x1_i)).intoArray(ret, j + dj);
        w2r.mul(x2_r).add(w2i.neg().mul(x2_i)).intoArray(ret, j + dj2);
        w3r.mul(x3_r).add(w3i.neg().mul(x3_i)).intoArray(ret, j + dj3);
        w1r.mul(x1_i).add(w1i.mul(x1_r)).intoArray(ret, j + dj + im);
        w2r.mul(x2_i).add(w2i.mul(x2_r)).intoArray(ret, j + dj2 + im);
        w3r.mul(x3_i).add(w3i.mul(x3_r)).intoArray(ret, j + dj3 + im);
      }
    }
  }

  /**
   * Handle factors of 4 using the 128-bit SIMD vectors.
   */
  private void interleaved128(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      DoubleVector
          z0 = fromArray(SPECIES_128, data, i),
          z1 = fromArray(SPECIES_128, data, i + di),
          z2 = fromArray(SPECIES_128, data, i + di2),
          z3 = fromArray(SPECIES_128, data, i + di3);
      DoubleVector
          t1 = z0.add(z2),
          t2 = z1.add(z3),
          t3 = z0.sub(z2),
          t4 = z1.sub(z3).mul(sign).rearrange(SHUFFLE_RE_IM_128);
      t1.add(t2).intoArray(ret, j);
      t3.add(t4.mul(NEGATE_RE_128)).intoArray(ret, j + dj);
      t1.sub(t2).intoArray(ret, j + dj2);
      t3.add(t4.mul(NEGATE_IM_128)).intoArray(ret, j + dj3);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(SPECIES_128, twids[0]),
          w1i = broadcast(SPECIES_128, -sign * twids[1]).mul(NEGATE_IM_128),
          w2r = broadcast(SPECIES_128, twids[2]),
          w2i = broadcast(SPECIES_128, -sign * twids[3]).mul(NEGATE_IM_128),
          w3r = broadcast(SPECIES_128, twids[4]),
          w3i = broadcast(SPECIES_128, -sign * twids[5]).mul(NEGATE_IM_128);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        DoubleVector
            z0 = fromArray(SPECIES_128, data, i),
            z1 = fromArray(SPECIES_128, data, i + di),
            z2 = fromArray(SPECIES_128, data, i + di2),
            z3 = fromArray(SPECIES_128, data, i + di3);
        DoubleVector
            t1 = z0.add(z2),
            t2 = z1.add(z3),
            t3 = z0.sub(z2),
            t4 = z1.sub(z3).mul(sign).rearrange(SHUFFLE_RE_IM_128);
        t1.add(t2).intoArray(ret, j);
        DoubleVector
            x1 = t3.add(t4.mul(NEGATE_RE_128)),
            x2 = t1.sub(t2),
            x3 = t3.add(t4.mul(NEGATE_IM_128));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj2);
        w3r.mul(x3).add(w3i.mul(x3).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj3);
      }
    }
  }

  /**
   * Handle factors of 4 using the 256-bit SIMD vectors.
   */
  private void interleaved256(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      DoubleVector
          z0 = fromArray(SPECIES_256, data, i),
          z1 = fromArray(SPECIES_256, data, i + di),
          z2 = fromArray(SPECIES_256, data, i + di2),
          z3 = fromArray(SPECIES_256, data, i + di3);
      DoubleVector
          t1 = z0.add(z2),
          t2 = z1.add(z3),
          t3 = z0.sub(z2),
          t4 = z1.sub(z3).mul(sign).rearrange(SHUFFLE_RE_IM_256);
      t1.add(t2).intoArray(ret, j);
      t3.add(t4.mul(NEGATE_RE_256)).intoArray(ret, j + dj);
      t1.sub(t2).intoArray(ret, j + dj2);
      t3.add(t4.mul(NEGATE_IM_256)).intoArray(ret, j + dj3);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(SPECIES_256, twids[0]),
          w1i = broadcast(SPECIES_256, -sign * twids[1]).mul(NEGATE_IM_256),
          w2r = broadcast(SPECIES_256, twids[2]),
          w2i = broadcast(SPECIES_256, -sign * twids[3]).mul(NEGATE_IM_256),
          w3r = broadcast(SPECIES_256, twids[4]),
          w3i = broadcast(SPECIES_256, -sign * twids[5]).mul(NEGATE_IM_256);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        DoubleVector
            z0 = fromArray(SPECIES_256, data, i),
            z1 = fromArray(SPECIES_256, data, i + di),
            z2 = fromArray(SPECIES_256, data, i + di2),
            z3 = fromArray(SPECIES_256, data, i + di3);
        DoubleVector
            t1 = z0.add(z2),
            t2 = z1.add(z3),
            t3 = z0.sub(z2),
            t4 = z1.sub(z3).mul(sign).rearrange(SHUFFLE_RE_IM_256);
        t1.add(t2).intoArray(ret, j);
        DoubleVector
            x1 = t3.add(t4.mul(NEGATE_RE_256)),
            x2 = t1.sub(t2),
            x3 = t3.add(t4.mul(NEGATE_IM_256));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj2);
        w3r.mul(x3).add(w3i.mul(x3).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj3);
      }
    }
  }

  /**
   * Handle factors of 4 using the 512-bit SIMD vectors.
   */
  private void interleaved512(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    // First pass of the 4-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      DoubleVector
          z0 = fromArray(SPECIES_512, data, i),
          z1 = fromArray(SPECIES_512, data, i + di),
          z2 = fromArray(SPECIES_512, data, i + di2),
          z3 = fromArray(SPECIES_512, data, i + di3);
      DoubleVector
          t1 = z0.add(z2),
          t2 = z1.add(z3),
          t3 = z0.sub(z2),
          t4 = z1.sub(z3).mul(sign).rearrange(SHUFFLE_RE_IM_512);
      t1.add(t2).intoArray(ret, j);
      t3.add(t4.mul(NEGATE_RE_512)).intoArray(ret, j + dj);
      t1.sub(t2).intoArray(ret, j + dj2);
      t3.add(t4.mul(NEGATE_IM_512)).intoArray(ret, j + dj3);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(SPECIES_512, twids[0]),
          w1i = broadcast(SPECIES_512, -sign * twids[1]).mul(NEGATE_IM_512),
          w2r = broadcast(SPECIES_512, twids[2]),
          w2i = broadcast(SPECIES_512, -sign * twids[3]).mul(NEGATE_IM_512),
          w3r = broadcast(SPECIES_512, twids[4]),
          w3i = broadcast(SPECIES_512, -sign * twids[5]).mul(NEGATE_IM_512);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        DoubleVector
            z0 = fromArray(SPECIES_512, data, i),
            z1 = fromArray(SPECIES_512, data, i + di),
            z2 = fromArray(SPECIES_512, data, i + di2),
            z3 = fromArray(SPECIES_512, data, i + di3);
        DoubleVector
            t1 = z0.add(z2),
            t2 = z1.add(z3),
            t3 = z0.sub(z2),
            t4 = z1.sub(z3).mul(sign).rearrange(SHUFFLE_RE_IM_512);
        t1.add(t2).intoArray(ret, j);
        DoubleVector
            x1 = t3.add(t4.mul(NEGATE_RE_512)),
            x2 = t1.sub(t2),
            x3 = t3.add(t4.mul(NEGATE_IM_512));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj2);
        w3r.mul(x3).add(w3i.mul(x3).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj3);
      }
    }
  }
}
