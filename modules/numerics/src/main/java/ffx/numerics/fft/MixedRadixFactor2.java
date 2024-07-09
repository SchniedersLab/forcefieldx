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

import static java.lang.Math.fma;

public class MixedRadixFactor2 extends MixedRadixFactor {

  public MixedRadixFactor2(PassConstants passConstants) {
    super(passConstants);
  }

  /**
   * Handle factors of 2.
   *
   * @param passData the data.
   */
  @Override
  protected void pass(PassData passData) {
    final int sign = passData.sign();
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();
    // First pass of the 2-point FFT has no twiddle factors.
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
      ret[jdj] = x_r;
      ret[jdj + 1] = x_i;
    }
    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
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
   * Handle factors of 2 using SIMD vectors.
   *
   * @param passData the data.
   */
  @Override
  protected void passSIMD(PassData passData) {
    if (innerLoopLimit % LOOP_INCREMENT == 0) {
      switch (SPECIES_LENGTH) {
        case 2:
          passSIMD_128(passData);
          break;
        case 4:
          passSIMD_256(passData);
          break;
        case 8:
          passSIMD_512(passData);
          break;
      }
    } else {
      // If the inner loop limit is not divisible by the loop increment, use largest SIMD vector that fits.
      switch (innerLoopLimit) {
        case 1:
          // Use the scalar method.
          pass(passData);
          break;
        case 2:
          passSIMD_128(passData);
          break;
        case 4:
          passSIMD_256(passData);
          break;
        case 8:
          passSIMD_512(passData);
          break;
        default:
          throw new IllegalStateException(" Unsupported inner loop limit: " + innerLoopLimit);
      }
    }
  }

  /**
   * Handle factors of 2 using the 128-bit SIMD vectors.
   *
   * @param passData the data.
   */
  private void passSIMD_128(PassData passData) {
    final int sign = passData.sign();
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT_128, i += SPECIES_LENGTH_128, j += SPECIES_LENGTH_128) {
      DoubleVector
          z0 = DoubleVector.fromArray(DOUBLE_SPECIES_128, data, i),
          z1 = DoubleVector.fromArray(DOUBLE_SPECIES_128, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(DOUBLE_SPECIES_128, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(DOUBLE_SPECIES_128, -sign * twids[1]).mul(NEGATE_IM_128);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT_128, i += SPECIES_LENGTH_128, j += SPECIES_LENGTH_128) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES_128, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES_128, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(SHUFFLE_RE_IM_128)).intoArray(ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 2 using the 256-bit SIMD vectors.
   *
   * @param passData the data.
   */
  private void passSIMD_256(PassData passData) {
    final int sign = passData.sign();
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT_256, i += SPECIES_LENGTH_256, j += SPECIES_LENGTH_256) {
      DoubleVector
          z0 = DoubleVector.fromArray(DOUBLE_SPECIES_256, data, i),
          z1 = DoubleVector.fromArray(DOUBLE_SPECIES_256, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(DOUBLE_SPECIES_256, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(DOUBLE_SPECIES_256, -sign * twids[1]).mul(NEGATE_IM_256);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT_256, i += SPECIES_LENGTH_256, j += SPECIES_LENGTH_256) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES_256, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES_256, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(SHUFFLE_RE_IM_256)).intoArray(ret, j + dj);
      }
    }
  }

  /**
   * Handle factors of 2 using the 512-bit SIMD vectors.
   *
   * @param passData the data.
   */
  private void passSIMD_512(PassData passData) {
    final int sign = passData.sign();
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();
    // First pass of the 2-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT_512, i += SPECIES_LENGTH_512, j += SPECIES_LENGTH_512) {
      DoubleVector
          z0 = DoubleVector.fromArray(DOUBLE_SPECIES_512, data, i),
          z1 = DoubleVector.fromArray(DOUBLE_SPECIES_512, data, i + di);
      z0.add(z1).intoArray(ret, j);
      z0.sub(z1).intoArray(ret, j + dj);
    }

    j += dj;
    for (int k = 1; k < outerLoopLimit; k++, j += dj) {
      final double[] twids = twiddles[k];
      DoubleVector wr = DoubleVector.broadcast(DOUBLE_SPECIES_512, twids[0]);
      DoubleVector wi = DoubleVector.broadcast(DOUBLE_SPECIES_512, -sign * twids[1]).mul(NEGATE_IM_512);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT_512, i += SPECIES_LENGTH_512, j += SPECIES_LENGTH_512) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES_512, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES_512, data, i + di);
        z0.add(z1).intoArray(ret, j);
        DoubleVector x = z0.sub(z1);
        x.mul(wr).add(x.mul(wi).rearrange(SHUFFLE_RE_IM_512)).intoArray(ret, j + dj);
      }
    }
  }

}
