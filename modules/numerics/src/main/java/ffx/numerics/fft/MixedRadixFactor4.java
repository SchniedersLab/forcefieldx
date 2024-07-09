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

public class MixedRadixFactor4 extends MixedRadixFactor {

  private final int di2;
  private final int di3;
  private final int dj2;
  private final int dj3;

  public MixedRadixFactor4(PassConstants passConstants) {
    super(passConstants);
    di2 = 2 * di;
    di3 = 3 * di;
    dj2 = 2 * dj;
    dj3 = 3 * dj;
  }

  /**
   * Handle factors of 4.
   *
   * @param passData the data.
   */
  protected void pass(PassData passData) {
    final int sign = passData.sign();
    int i = passData.inOffset();
    int j = passData.outOffset();
    final double[] data = passData.in();
    final double[] ret = passData.out();

    // First pass of the 4-point FFT has no twiddle factors.
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
      int jdj = j + dj;
      ret[jdj] = t3_r - t4_i;
      ret[jdj + 1] = t3_i + t4_r;
      jdj += dj;
      ret[jdj] = t1_r - t2_r;
      ret[jdj + 1] = t1_i - t2_i;
      jdj += dj;
      ret[jdj] = t3_r + t4_i;
      ret[jdj + 1] = t3_i - t4_r;
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
   * @param passData the data.
   */
  protected void passSIMD(PassData passData) {
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      // System.out.printf("Scalar %d product=%d innerLoopLimit=%d increment=%d%n",
      // factor, product, innerLoopLimit, LOOP_INCREMENT);
      pass(passData);
      return;
    }
    final int sign = passData.sign();
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();

    // First pass of the 4-point FFT has no twiddle factors.
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
      t3.add(t4.mul(negateRe)).intoArray(ret, j + dj);
      t1.sub(t2).intoArray(ret, j + dj2);
      t3.add(t4.mul(negateIm)).intoArray(ret, j + dj3);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
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

}
