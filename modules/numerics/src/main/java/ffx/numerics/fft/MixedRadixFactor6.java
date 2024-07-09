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
import static org.apache.commons.math3.util.FastMath.sqrt;

public class MixedRadixFactor6 extends MixedRadixFactor {

  private static final double sqrt3_2 = sqrt(3.0) / 2.0;

  private final int di2;
  private final int di3;
  private final int di4;
  private final int di5;
  private final int dj2;
  private final int dj3;
  private final int dj4;
  private final int dj5;

  public MixedRadixFactor6(PassConstants passConstants) {
    super(passConstants);
    di2 = 2 * di;
    di3 = 3 * di;
    di4 = 4 * di;
    di5 = 5 * di;
    dj2 = 2 * dj;
    dj3 = 3 * dj;
    dj4 = 4 * dj;
    dj5 = 5 * dj;
  }

  /**
   * Handle factors of 6.
   *
   * @param passData the data.
   */
  protected void pass(PassData passData) {
    final int sign = passData.sign();
    final double tau = sign * sqrt3_2;
    int i = passData.inOffset();
    int j = passData.outOffset();
    final double[] data = passData.in();
    final double[] ret = passData.out();

    // First pass of the 6-point FFT has no twiddle factors.
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
      int jdj = j + dj;
      ret[jdj] = a1r - b1r;
      ret[jdj + 1] = a1i - b1i;
      jdj += dj;
      ret[jdj] = a2r + b2r;
      ret[jdj + 1] = a2i + b2i;
      jdj += dj;
      ret[jdj] = a0r - b0r;
      ret[jdj + 1] = a0i - b0i;
      jdj += dj;
      ret[jdj] = a1r + b1r;
      ret[jdj + 1] = a1i + b1i;
      jdj += dj;
      ret[jdj] = a2r - b2r;
      ret[jdj + 1] = a2i - b2i;
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
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
   * @param passData the data.
   */
  protected void passSIMD(PassData passData) {
    // If the inner loop limit is not divisible by the loop increment, use the scalar method.
    if (innerLoopLimit % LOOP_INCREMENT != 0) {
      pass(passData);
      return;
    }
    final int sign = passData.sign();
    final double tau = sign * sqrt3_2;
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();

    // First pass of the 6-point FFT has no twiddle factors.
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
      a1.sub(b1).intoArray(ret, j + dj);
      a2.add(b2).intoArray(ret, j + dj2);
      a0.sub(b0).intoArray(ret, j + dj3);
      a1.add(b1).intoArray(ret, j + dj4);
      a2.sub(b2).intoArray(ret, j + dj5);
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

}
