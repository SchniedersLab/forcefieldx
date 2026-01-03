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

import jdk.incubator.vector.DoubleVector;

import static java.lang.Math.fma;
import static jdk.incubator.vector.DoubleVector.SPECIES_128;
import static jdk.incubator.vector.DoubleVector.SPECIES_256;
import static jdk.incubator.vector.DoubleVector.SPECIES_512;
import static jdk.incubator.vector.DoubleVector.broadcast;
import static jdk.incubator.vector.DoubleVector.fromArray;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The MixedRadixFactor3 class handles factors of 3 in the FFT.
 */
public class MixedRadixFactor3 extends MixedRadixFactor {

  private static final double sqrt3_2 = sqrt(3.0) / 2.0;
  private final int di2;
  private final int dj2;

  /**
   * Available SIMD sizes for Pass 3.
   */
  private static int[] simdSizes = {8, 4, 2};

  /**
   * Create a new MixedRadixFactor3 instance.
   *
   * @param passConstants The pass constants.
   */
  public MixedRadixFactor3(PassConstants passConstants) {
    super(passConstants);
    di2 = 2 * di;
    dj2 = 2 * dj;
  }

  /**
   * Check if the requested SIMD length is valid.
   * @param width Requested SIMD species width.
   * @return True if this width is supported.
   */
  public boolean isValidSIMDWidth(int width) {
    if (width <= 0 || width > LENGTH) {
      return false;
    }
    // Must be a supported width.
    if (width != 2 && width != 4 && width != 8) {
      return false;
    }
    if (im == 1) {
      // Interleaved
      return innerLoopLimit % (width / 2) == 0;
    } else {
      // Blocked
      return innerLoopLimit % width == 0;
    }
  }

  /**
   * Determine the optimal SIMD width. Currently supported widths are 2, 4 and 8.
   * If no SIMD width is valid, return 0 to indicate use of the scalar path.
   * @return The optimal SIMD width.
   */
  public int getOptimalSIMDWidth() {
    // Check the platform specific preferred width.
    if (isValidSIMDWidth(LENGTH)) {
      return LENGTH;
    }
    // Fall back to a smaller SIMD vector that fits the inner loop limit.
    for (int size : simdSizes) {
      if (size >= LENGTH) {
        // Skip anything greater than or equal to the preferred SIMD vector size (which was too big).
        continue;
      }
      if (isValidSIMDWidth(size)) {
        return size;
      }
    }
    // No valid SIMD width is found.
    return 0;
  }

  /**
   * Handle factors of 3.
   */
  @Override
  protected void passScalar(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;

    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0_r = data[i];
      final double z1_r = data[i + di];
      final double z2_r = data[i + di2];
      final double z0_i = data[i + im];
      final double z1_i = data[i + di + im];
      final double z2_i = data[i + di2 + im];
      final double t1_r = z1_r + z2_r;
      final double t1_i = z1_i + z2_i;
      final double t2_r = fma(-0.5, t1_r, z0_r);
      final double t2_i = fma(-0.5, t1_i, z0_i);
      final double t3_r = tau * (z1_r - z2_r);
      final double t3_i = tau * (z1_i - z2_i);
      ret[j] = z0_r + t1_r;
      ret[j + im] = z0_i + t1_i;
      ret[j + dj] = t2_r - t3_i;
      ret[j + dj + im] = t2_i + t3_r;
      ret[j + dj2] = t2_r + t3_i;
      ret[j + dj2 + im] = t2_i - t3_r;
    }
    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final double w1_r = wr[index];
      final double w2_r = wr[index + 1];
      final double w1_i = -sign * wi[index];
      final double w2_i = -sign * wi[index + 1];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0_r = data[i];
        final double z1_r = data[i + di];
        final double z2_r = data[i + di2];
        final double z0_i = data[i + im];
        final double z1_i = data[i + di + im];
        final double z2_i = data[i + di2 + im];
        final double t1_r = z1_r + z2_r;
        final double t1_i = z1_i + z2_i;
        final double t2_r = fma(-0.5, t1_r, z0_r);
        final double t2_i = fma(-0.5, t1_i, z0_i);
        final double t3_r = tau * (z1_r - z2_r);
        final double t3_i = tau * (z1_i - z2_i);
        ret[j] = z0_r + t1_r;
        ret[j + im] = z0_i + t1_i;
        multiplyAndStore(t2_r - t3_i, t2_i + t3_r, w1_r, w1_i, ret, j + dj, j + dj + im);
        multiplyAndStore(t2_r + t3_i, t2_i - t3_r, w2_r, w2_i, ret, j + dj2, j + dj2 + im);
      }
    }
  }

  /**
   * Handle factors of 3 using SIMD vectors.
   *
   * @param passData The pass data.
   */
  @Override
  protected void passSIMD(PassData passData) {
    if (!isValidSIMDWidth(simdWidth)) {
      passScalar(passData);
    } else {
      if (im == 1) {
        interleaved(passData, simdWidth);
      } else {
        blocked(passData, simdWidth);
      }
    }
  }

  /**
   * Handle factors of 3 using the chosen SIMD vector.
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
   * Handle factors of 3 using the chosen SIMD vector.
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
   * Handle factors of 3 using the 128-bit SIMD vectors.
   */
  private void blocked128(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;

    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      final DoubleVector
          z0_r = fromArray(SPECIES_128, data, i),
          z1_r = fromArray(SPECIES_128, data, i + di),
          z2_r = fromArray(SPECIES_128, data, i + di2),
          z0_i = fromArray(SPECIES_128, data, i + im),
          z1_i = fromArray(SPECIES_128, data, i + di + im),
          z2_i = fromArray(SPECIES_128, data, i + di2 + im);
      final DoubleVector
          t1_r = z1_r.add(z2_r),
          t1_i = z1_i.add(z2_i),
          t2_r = t1_r.mul(-0.5).add(z0_r),
          t2_i = t1_i.mul(-0.5).add(z0_i),
          t3_r = z1_r.sub(z2_r).mul(tau),
          t3_i = z1_i.sub(z2_i).mul(tau);
      z0_r.add(t1_r).intoArray(ret, j);
      z0_i.add(t1_i).intoArray(ret, j + im);
      t2_r.sub(t3_i).intoArray(ret, j + dj);
      t2_i.add(t3_r).intoArray(ret, j + dj + im);
      t2_r.add(t3_i).intoArray(ret, j + dj2);
      t2_i.sub(t3_r).intoArray(ret, j + dj2 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final DoubleVector
          w1_r = broadcast(SPECIES_128, wr[index]),
          w2_r = broadcast(SPECIES_128, wr[index + 1]),
          w1_i = broadcast(SPECIES_128, -sign * wi[index]),
          w2_i = broadcast(SPECIES_128, -sign * wi[index + 1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        final DoubleVector
            z0_r = fromArray(SPECIES_128, data, i),
            z1_r = fromArray(SPECIES_128, data, i + di),
            z2_r = fromArray(SPECIES_128, data, i + di2),
            z0_i = fromArray(SPECIES_128, data, i + im),
            z1_i = fromArray(SPECIES_128, data, i + di + im),
            z2_i = fromArray(SPECIES_128, data, i + di2 + im);
        final DoubleVector
            t1_r = z1_r.add(z2_r),
            t1_i = z1_i.add(z2_i),
            t2_r = t1_r.mul(-0.5).add(z0_r),
            t2_i = t1_i.mul(-0.5).add(z0_i),
            t3_r = z1_r.sub(z2_r).mul(tau),
            t3_i = z1_i.sub(z2_i).mul(tau);
        z0_r.add(t1_r).intoArray(ret, j);
        z0_i.add(t1_i).intoArray(ret, j + im);
        final DoubleVector
            x1_r = t2_r.sub(t3_i), x1_i = t2_i.add(t3_r),
            x2_r = t2_r.add(t3_i), x2_i = t2_i.sub(t3_r);
        w1_r.mul(x1_r).add(w1_i.neg().mul(x1_i)).intoArray(ret, j + dj);
        w2_r.mul(x2_r).add(w2_i.neg().mul(x2_i)).intoArray(ret, j + dj2);
        w1_r.mul(x1_i).add(w1_i.mul(x1_r)).intoArray(ret, j + dj + im);
        w2_r.mul(x2_i).add(w2_i.mul(x2_r)).intoArray(ret, j + dj2 + im);
      }
    }
  }

  /**
   * Handle factors of 3 using the 256-bit SIMD vectors.
   */
  private void blocked256(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;

    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      final DoubleVector
          z0_r = fromArray(SPECIES_256, data, i),
          z1_r = fromArray(SPECIES_256, data, i + di),
          z2_r = fromArray(SPECIES_256, data, i + di2),
          z0_i = fromArray(SPECIES_256, data, i + im),
          z1_i = fromArray(SPECIES_256, data, i + di + im),
          z2_i = fromArray(SPECIES_256, data, i + di2 + im);
      final DoubleVector
          t1_r = z1_r.add(z2_r),
          t1_i = z1_i.add(z2_i),
          t2_r = t1_r.mul(-0.5).add(z0_r),
          t2_i = t1_i.mul(-0.5).add(z0_i),
          t3_r = z1_r.sub(z2_r).mul(tau),
          t3_i = z1_i.sub(z2_i).mul(tau);
      z0_r.add(t1_r).intoArray(ret, j);
      z0_i.add(t1_i).intoArray(ret, j + im);
      t2_r.sub(t3_i).intoArray(ret, j + dj);
      t2_i.add(t3_r).intoArray(ret, j + dj + im);
      t2_r.add(t3_i).intoArray(ret, j + dj2);
      t2_i.sub(t3_r).intoArray(ret, j + dj2 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final DoubleVector
          w1_r = broadcast(SPECIES_256, wr[index]),
          w2_r = broadcast(SPECIES_256, wr[index + 1]),
          w1_i = broadcast(SPECIES_256, -sign * wi[index]),
          w2_i = broadcast(SPECIES_256, -sign * wi[index + 1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        final DoubleVector
            z0_r = fromArray(SPECIES_256, data, i),
            z1_r = fromArray(SPECIES_256, data, i + di),
            z2_r = fromArray(SPECIES_256, data, i + di2),
            z0_i = fromArray(SPECIES_256, data, i + im),
            z1_i = fromArray(SPECIES_256, data, i + di + im),
            z2_i = fromArray(SPECIES_256, data, i + di2 + im);
        final DoubleVector
            t1_r = z1_r.add(z2_r),
            t1_i = z1_i.add(z2_i),
            t2_r = t1_r.mul(-0.5).add(z0_r),
            t2_i = t1_i.mul(-0.5).add(z0_i),
            t3_r = z1_r.sub(z2_r).mul(tau),
            t3_i = z1_i.sub(z2_i).mul(tau);
        z0_r.add(t1_r).intoArray(ret, j);
        z0_i.add(t1_i).intoArray(ret, j + im);
        final DoubleVector
            x1_r = t2_r.sub(t3_i), x1_i = t2_i.add(t3_r),
            x2_r = t2_r.add(t3_i), x2_i = t2_i.sub(t3_r);
        w1_r.mul(x1_r).add(w1_i.neg().mul(x1_i)).intoArray(ret, j + dj);
        w2_r.mul(x2_r).add(w2_i.neg().mul(x2_i)).intoArray(ret, j + dj2);
        w1_r.mul(x1_i).add(w1_i.mul(x1_r)).intoArray(ret, j + dj + im);
        w2_r.mul(x2_i).add(w2_i.mul(x2_r)).intoArray(ret, j + dj2 + im);
      }
    }
  }

  /**
   * Handle factors of 3 using the 512-bit SIMD vectors.
   */
  private void blocked512(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;

    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      final DoubleVector
          z0_r = fromArray(SPECIES_512, data, i),
          z1_r = fromArray(SPECIES_512, data, i + di),
          z2_r = fromArray(SPECIES_512, data, i + di2),
          z0_i = fromArray(SPECIES_512, data, i + im),
          z1_i = fromArray(SPECIES_512, data, i + di + im),
          z2_i = fromArray(SPECIES_512, data, i + di2 + im);
      final DoubleVector
          t1_r = z1_r.add(z2_r),
          t1_i = z1_i.add(z2_i),
          t2_r = t1_r.mul(-0.5).add(z0_r),
          t2_i = t1_i.mul(-0.5).add(z0_i),
          t3_r = z1_r.sub(z2_r).mul(tau),
          t3_i = z1_i.sub(z2_i).mul(tau);
      z0_r.add(t1_r).intoArray(ret, j);
      z0_i.add(t1_i).intoArray(ret, j + im);
      t2_r.sub(t3_i).intoArray(ret, j + dj);
      t2_i.add(t3_r).intoArray(ret, j + dj + im);
      t2_r.add(t3_i).intoArray(ret, j + dj2);
      t2_i.sub(t3_r).intoArray(ret, j + dj2 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final DoubleVector
          w1_r = broadcast(SPECIES_512, wr[index]),
          w2_r = broadcast(SPECIES_512, wr[index + 1]),
          w1_i = broadcast(SPECIES_512, -sign * wi[index]),
          w2_i = broadcast(SPECIES_512, -sign * wi[index + 1]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        final DoubleVector
            z0_r = fromArray(SPECIES_512, data, i),
            z1_r = fromArray(SPECIES_512, data, i + di),
            z2_r = fromArray(SPECIES_512, data, i + di2),
            z0_i = fromArray(SPECIES_512, data, i + im),
            z1_i = fromArray(SPECIES_512, data, i + di + im),
            z2_i = fromArray(SPECIES_512, data, i + di2 + im);
        final DoubleVector
            t1_r = z1_r.add(z2_r),
            t1_i = z1_i.add(z2_i),
            t2_r = t1_r.mul(-0.5).add(z0_r),
            t2_i = t1_i.mul(-0.5).add(z0_i),
            t3_r = z1_r.sub(z2_r).mul(tau),
            t3_i = z1_i.sub(z2_i).mul(tau);
        z0_r.add(t1_r).intoArray(ret, j);
        z0_i.add(t1_i).intoArray(ret, j + im);
        final DoubleVector
            x1_r = t2_r.sub(t3_i), x1_i = t2_i.add(t3_r),
            x2_r = t2_r.add(t3_i), x2_i = t2_i.sub(t3_r);
        w1_r.mul(x1_r).add(w1_i.neg().mul(x1_i)).intoArray(ret, j + dj);
        w2_r.mul(x2_r).add(w2_i.neg().mul(x2_i)).intoArray(ret, j + dj2);
        w1_r.mul(x1_i).add(w1_i.mul(x1_r)).intoArray(ret, j + dj + im);
        w2_r.mul(x2_i).add(w2_i.mul(x2_r)).intoArray(ret, j + dj2 + im);
      }
    }
  }

  /**
   * Handle factors of 3 using the 128-bit SIMD vectors.
   */
  private void interleaved128(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;
    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
      final DoubleVector
          z0 = fromArray(SPECIES_128, data, i),
          z1 = fromArray(SPECIES_128, data, i + di),
          z2 = fromArray(SPECIES_128, data, i + di2);
      final DoubleVector
          t1 = z1.add(z2),
          t2 = t1.mul(-0.5).add(z0),
          t3 = z1.sub(z2).mul(tau).rearrange(SHUFFLE_RE_IM_128);
      z0.add(t1).intoArray(ret, j);
      t2.add(t3.mul(NEGATE_RE_128)).intoArray(ret, j + dj);
      t2.add(t3.mul(NEGATE_IM_128)).intoArray(ret, j + dj2);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final DoubleVector
          w1r = broadcast(SPECIES_128, wr[index]),
          w2r = broadcast(SPECIES_128, wr[index + 1]),
          w1i = broadcast(SPECIES_128, -sign * wi[index]).mul(NEGATE_IM_128),
          w2i = broadcast(SPECIES_128, -sign * wi[index + 1]).mul(NEGATE_IM_128);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP_128, i += LENGTH_128, j += LENGTH_128) {
        final DoubleVector
            z0 = fromArray(SPECIES_128, data, i),
            z1 = fromArray(SPECIES_128, data, i + di),
            z2 = fromArray(SPECIES_128, data, i + di2);
        final DoubleVector
            t1 = z1.add(z2),
            t2 = t1.mul(-0.5).add(z0),
            t3 = z1.sub(z2).mul(tau).rearrange(SHUFFLE_RE_IM_128);
        z0.add(t1).intoArray(ret, j);
        z0.add(t1).intoArray(ret, j);
        final DoubleVector
            x1 = t2.add(t3.mul(NEGATE_RE_128)),
            x2 = t2.add(t3.mul(NEGATE_IM_128));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj2);
      }
    }
  }

  /**
   * Handle factors of 3 using the 256-bit SIMD vectors.
   */
  private void interleaved256(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;
    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
      final DoubleVector
          z0 = fromArray(SPECIES_256, data, i),
          z1 = fromArray(SPECIES_256, data, i + di),
          z2 = fromArray(SPECIES_256, data, i + di2);
      final DoubleVector
          t1 = z1.add(z2),
          t2 = t1.mul(-0.5).add(z0),
          t3 = z1.sub(z2).mul(tau).rearrange(SHUFFLE_RE_IM_256);
      z0.add(t1).intoArray(ret, j);
      t2.add(t3.mul(NEGATE_RE_256)).intoArray(ret, j + dj);
      t2.add(t3.mul(NEGATE_IM_256)).intoArray(ret, j + dj2);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final DoubleVector
          w1r = broadcast(SPECIES_256, wr[index]),
          w2r = broadcast(SPECIES_256, wr[index + 1]),
          w1i = broadcast(SPECIES_256, -sign * wi[index]).mul(NEGATE_IM_256),
          w2i = broadcast(SPECIES_256, -sign * wi[index + 1]).mul(NEGATE_IM_256);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP_256, i += LENGTH_256, j += LENGTH_256) {
        final DoubleVector
            z0 = fromArray(SPECIES_256, data, i),
            z1 = fromArray(SPECIES_256, data, i + di),
            z2 = fromArray(SPECIES_256, data, i + di2);
        final DoubleVector
            t1 = z1.add(z2),
            t2 = t1.mul(-0.5).add(z0),
            t3 = z1.sub(z2).mul(tau).rearrange(SHUFFLE_RE_IM_256);
        z0.add(t1).intoArray(ret, j);
        final DoubleVector
            x1 = t2.add(t3.mul(NEGATE_RE_256)),
            x2 = t2.add(t3.mul(NEGATE_IM_256));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj2);
      }
    }
  }

  /**
   * Handle factors of 3 using the 512-bit SIMD vectors.
   */
  private void interleaved512(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;
    // First pass of the 3-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
      final DoubleVector
          z0 = fromArray(SPECIES_512, data, i),
          z1 = fromArray(SPECIES_512, data, i + di),
          z2 = fromArray(SPECIES_512, data, i + di2);
      final DoubleVector
          t1 = z1.add(z2),
          t2 = t1.mul(-0.5).add(z0),
          t3 = z1.sub(z2).mul(tau).rearrange(SHUFFLE_RE_IM_512);
      z0.add(t1).intoArray(ret, j);
      t2.add(t3.mul(NEGATE_RE_512)).intoArray(ret, j + dj);
      t2.add(t3.mul(NEGATE_IM_512)).intoArray(ret, j + dj2);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 2;
      final DoubleVector
          w1r = broadcast(SPECIES_512, wr[index]),
          w2r = broadcast(SPECIES_512, wr[index + 1]),
          w1i = broadcast(SPECIES_512, -sign * wi[index]).mul(NEGATE_IM_512),
          w2i = broadcast(SPECIES_512, -sign * wi[index + 1]).mul(NEGATE_IM_512);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP_512, i += LENGTH_512, j += LENGTH_512) {
        final DoubleVector
            z0 = fromArray(SPECIES_512, data, i),
            z1 = fromArray(SPECIES_512, data, i + di),
            z2 = fromArray(SPECIES_512, data, i + di2);
        final DoubleVector
            t1 = z1.add(z2),
            t2 = t1.mul(-0.5).add(z0),
            t3 = z1.sub(z2).mul(tau).rearrange(SHUFFLE_RE_IM_512);
        z0.add(t1).intoArray(ret, j);
        final DoubleVector
            x1 = t2.add(t3.mul(NEGATE_RE_512)),
            x2 = t2.add(t3.mul(NEGATE_IM_512));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj2);
      }
    }
  }
}
