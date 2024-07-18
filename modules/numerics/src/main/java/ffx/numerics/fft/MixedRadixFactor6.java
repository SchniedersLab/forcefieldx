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
import static jdk.incubator.vector.DoubleVector.broadcast;
import static jdk.incubator.vector.DoubleVector.fromArray;
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
   */
  @Override
  protected void passScalar(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;
    // First pass of the 6-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0r = data[i];
      final double z1r = data[i + di];
      final double z2r = data[i + di2];
      final double z3r = data[i + di3];
      final double z4r = data[i + di4];
      final double z5r = data[i + di5];
      final double z0i = data[i + im];
      final double z1i = data[i + di + im];
      final double z2i = data[i + di2 + im];
      final double z3i = data[i + di3 + im];
      final double z4i = data[i + di4 + im];
      final double z5i = data[i + di5 + im];
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
      ret[j + im] = a0i + b0i;
      ret[j + dj] = a1r - b1r;
      ret[j + dj + im] = a1i - b1i;
      ret[j + dj2] = a2r + b2r;
      ret[j + dj2 + im] = a2i + b2i;
      ret[j + dj3] = a0r - b0r;
      ret[j + dj3 + im] = a0i - b0i;
      ret[j + dj4] = a1r + b1r;
      ret[j + dj4 + im] = a1i + b1i;
      ret[j + dj5] = a2r - b2r;
      ret[j + dj5 + im] = a2i - b2i;
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
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0r = data[i];
        final double z1r = data[i + di];
        final double z2r = data[i + di2];
        final double z3r = data[i + di3];
        final double z4r = data[i + di4];
        final double z5r = data[i + di5];
        final double z0i = data[i + im];
        final double z1i = data[i + di + im];
        final double z2i = data[i + di2 + im];
        final double z3i = data[i + di3 + im];
        final double z4i = data[i + di4 + im];
        final double z5i = data[i + di5 + im];
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
        ret[j + im] = a0i + b0i;
        multiplyAndStore(a1r - b1r, a1i - b1i, w1r, w1i, ret, j + dj, j + dj + im);
        multiplyAndStore(a2r + b2r, a2i + b2i, w2r, w2i, ret, j + dj2, j + dj2 + im);
        multiplyAndStore(a0r - b0r, a0i - b0i, w3r, w3i, ret, j + dj3, j + dj3 + im);
        multiplyAndStore(a1r + b1r, a1i + b1i, w4r, w4i, ret, j + dj4, j + dj4 + im);
        multiplyAndStore(a2r - b2r, a2i - b2i, w5r, w5i, ret, j + dj5, j + dj5 + im);
      }
    }
  }

  /**
   * Handle factors of 6 using SIMD vectors.
   */
  @Override
  protected void passSIMD(PassData passData) {
    // Interleaved.
    if (im == 1) {
      // If the inner loop limit is not divisible by the loop increment, use the scalar method.
      if (innerLoopLimit % LOOP != 0) {
        passScalar(passData);
      } else {
        interleaved(passData);
      }
      // Blocked.
    } else {
      // If the inner loop limit is not divisible by the loop increment, use the scalar method.
      if (innerLoopLimit % BLOCK_LOOP != 0) {
        passScalar(passData);
      } else {
        blocked(passData);
      }
    }
  }

  /**
   * Handle factors of 6.
   */
  private void interleaved(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;
    // First pass of the 6-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP, i += LENGTH, j += LENGTH) {
      DoubleVector
          z0 = fromArray(DOUBLE_SPECIES, data, i),
          z1 = fromArray(DOUBLE_SPECIES, data, i + di),
          z2 = fromArray(DOUBLE_SPECIES, data, i + di2),
          z3 = fromArray(DOUBLE_SPECIES, data, i + di3),
          z4 = fromArray(DOUBLE_SPECIES, data, i + di4),
          z5 = fromArray(DOUBLE_SPECIES, data, i + di5);
      DoubleVector
          ta1 = z2.add(z4),
          ta2 = ta1.mul(-0.5).add(z0),
          ta3 = z2.sub(z4).mul(tau).rearrange(SHUFFLE_RE_IM),
          a0 = z0.add(ta1),
          a1 = ta2.add(ta3.mul(NEGATE_RE)),
          a2 = ta2.add(ta3.mul(NEGATE_IM)),
          tb1 = z5.add(z1),
          tb2 = tb1.mul(-0.5).add(z3),
          tb3 = z5.sub(z1).mul(tau).rearrange(SHUFFLE_RE_IM),
          b0 = z3.add(tb1),
          b1 = tb2.add(tb3.mul(NEGATE_RE)),
          b2 = tb2.add(tb3.mul(NEGATE_IM));
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
          w1r = broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(NEGATE_IM),
          w2r = broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(NEGATE_IM),
          w3r = broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(NEGATE_IM),
          w4r = broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = broadcast(DOUBLE_SPECIES, -sign * twids[7]).mul(NEGATE_IM),
          w5r = broadcast(DOUBLE_SPECIES, twids[8]),
          w5i = broadcast(DOUBLE_SPECIES, -sign * twids[9]).mul(NEGATE_IM);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP, i += LENGTH, j += LENGTH) {
        DoubleVector
            z0 = fromArray(DOUBLE_SPECIES, data, i),
            z1 = fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = fromArray(DOUBLE_SPECIES, data, i + di4),
            z5 = fromArray(DOUBLE_SPECIES, data, i + di5);
        DoubleVector
            ta1 = z2.add(z4),
            ta2 = ta1.mul(-0.5).add(z0),
            ta3 = z2.sub(z4).mul(tau).rearrange(SHUFFLE_RE_IM),
            a0 = z0.add(ta1),
            a1 = ta2.add(ta3.mul(NEGATE_RE)),
            a2 = ta2.add(ta3.mul(NEGATE_IM)),
            tb1 = z5.add(z1),
            tb2 = tb1.mul(-0.5).add(z3),
            tb3 = z5.sub(z1).mul(tau).rearrange(SHUFFLE_RE_IM),
            b0 = z3.add(tb1),
            b1 = tb2.add(tb3.mul(NEGATE_RE)),
            b2 = tb2.add(tb3.mul(NEGATE_IM));
        a0.add(b0).intoArray(ret, j);
        DoubleVector
            x1 = a1.sub(b1),
            x2 = a2.add(b2),
            x3 = a0.sub(b0),
            x4 = a1.add(b1),
            x5 = a2.sub(b2);
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj2);
        w3r.mul(x3).add(w3i.mul(x3).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj3);
        w4r.mul(x4).add(w4i.mul(x4).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj4);
        w5r.mul(x5).add(w5i.mul(x5).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj5);
      }
    }
  }

  /**
   * Handle factors of 6.
   */
  private void blocked(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double tau = sign * sqrt3_2;
    // First pass of the 6-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP, i += LENGTH, j += LENGTH) {
      final DoubleVector
          z0r = fromArray(DOUBLE_SPECIES, data, i),
          z1r = fromArray(DOUBLE_SPECIES, data, i + di),
          z2r = fromArray(DOUBLE_SPECIES, data, i + di2),
          z3r = fromArray(DOUBLE_SPECIES, data, i + di3),
          z4r = fromArray(DOUBLE_SPECIES, data, i + di4),
          z5r = fromArray(DOUBLE_SPECIES, data, i + di5),
          z0i = fromArray(DOUBLE_SPECIES, data, i + im),
          z1i = fromArray(DOUBLE_SPECIES, data, i + di + im),
          z2i = fromArray(DOUBLE_SPECIES, data, i + di2 + im),
          z3i = fromArray(DOUBLE_SPECIES, data, i + di3 + im),
          z4i = fromArray(DOUBLE_SPECIES, data, i + di4 + im),
          z5i = fromArray(DOUBLE_SPECIES, data, i + di5 + im);
      final DoubleVector
          ta1r = z2r.add(z4r),
          ta1i = z2i.add(z4i),
          ta2r = ta1r.mul(-0.5).add(z0r),
          ta2i = ta1i.mul(-0.5).add(z0i),
          ta3r = z2r.sub(z4r).mul(tau),
          ta3i = z2i.sub(z4i).mul(tau),
          a0r = z0r.add(ta1r),
          a0i = z0i.add(ta1i),
          a1r = ta2r.sub(ta3i),
          a1i = ta2i.add(ta3r),
          a2r = ta2r.add(ta3i),
          a2i = ta2i.sub(ta3r),
          tb1r = z5r.add(z1r),
          tb1i = z5i.add(z1i),
          tb2r = tb1r.mul(-0.5).add(z3r),
          tb2i = tb1i.mul(-0.5).add(z3i),
          tb3r = z5r.sub(z1r).mul(tau),
          tb3i = z5i.sub(z1i).mul(tau),
          b0r = z3r.add(tb1r),
          b0i = z3i.add(tb1i),
          b1r = tb2r.sub(tb3i),
          b1i = tb2i.add(tb3r),
          b2r = tb2r.add(tb3i),
          b2i = tb2i.sub(tb3r);
      a0r.add(b0r).intoArray(ret, j);
      a0i.add(b0i).intoArray(ret, j + im);
      a1r.sub(b1r).intoArray(ret, j + dj);
      a1i.sub(b1i).intoArray(ret, j + dj + im);
      a2r.add(b2r).intoArray(ret, j + dj2);
      a2i.add(b2i).intoArray(ret, j + dj2 + im);
      a0r.sub(b0r).intoArray(ret, j + dj3);
      a0i.sub(b0i).intoArray(ret, j + dj3 + im);
      a1r.add(b1r).intoArray(ret, j + dj4);
      a1i.add(b1i).intoArray(ret, j + dj4 + im);
      a2r.sub(b2r).intoArray(ret, j + dj5);
      a2i.sub(b2i).intoArray(ret, j + dj5 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = broadcast(DOUBLE_SPECIES, -sign * twids[1]),
          w2r = broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = broadcast(DOUBLE_SPECIES, -sign * twids[3]),
          w3r = broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = broadcast(DOUBLE_SPECIES, -sign * twids[5]),
          w4r = broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = broadcast(DOUBLE_SPECIES, -sign * twids[7]),
          w5r = broadcast(DOUBLE_SPECIES, twids[8]),
          w5i = broadcast(DOUBLE_SPECIES, -sign * twids[9]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP, i += LENGTH, j += LENGTH) {
        final DoubleVector
            z0r = fromArray(DOUBLE_SPECIES, data, i),
            z1r = fromArray(DOUBLE_SPECIES, data, i + di),
            z2r = fromArray(DOUBLE_SPECIES, data, i + di2),
            z3r = fromArray(DOUBLE_SPECIES, data, i + di3),
            z4r = fromArray(DOUBLE_SPECIES, data, i + di4),
            z5r = fromArray(DOUBLE_SPECIES, data, i + di5),
            z0i = fromArray(DOUBLE_SPECIES, data, i + im),
            z1i = fromArray(DOUBLE_SPECIES, data, i + di + im),
            z2i = fromArray(DOUBLE_SPECIES, data, i + di2 + im),
            z3i = fromArray(DOUBLE_SPECIES, data, i + di3 + im),
            z4i = fromArray(DOUBLE_SPECIES, data, i + di4 + im),
            z5i = fromArray(DOUBLE_SPECIES, data, i + di5 + im);
        final DoubleVector
            ta1r = z2r.add(z4r),
            ta1i = z2i.add(z4i),
            ta2r = ta1r.mul(-0.5).add(z0r),
            ta2i = ta1i.mul(-0.5).add(z0i),
            ta3r = z2r.sub(z4r).mul(tau),
            ta3i = z2i.sub(z4i).mul(tau),
            a0r = z0r.add(ta1r),
            a0i = z0i.add(ta1i),
            a1r = ta2r.sub(ta3i),
            a1i = ta2i.add(ta3r),
            a2r = ta2r.add(ta3i),
            a2i = ta2i.sub(ta3r),
            tb1r = z5r.add(z1r),
            tb1i = z5i.add(z1i),
            tb2r = tb1r.mul(-0.5).add(z3r),
            tb2i = tb1i.mul(-0.5).add(z3i),
            tb3r = z5r.sub(z1r).mul(tau),
            tb3i = z5i.sub(z1i).mul(tau),
            b0r = z3r.add(tb1r),
            b0i = z3i.add(tb1i),
            b1r = tb2r.sub(tb3i),
            b1i = tb2i.add(tb3r),
            b2r = tb2r.add(tb3i),
            b2i = tb2i.sub(tb3r);
        a0r.add(b0r).intoArray(ret, j);
        a0i.add(b0i).intoArray(ret, j + im);
        DoubleVector
            x1r = a1r.sub(b1r), x1i = a1i.sub(b1i),
            x2r = a2r.add(b2r), x2i = a2i.add(b2i),
            x3r = a0r.sub(b0r), x3i = a0i.sub(b0i),
            x4r = a1r.add(b1r), x4i = a1i.add(b1i),
            x5r = a2r.sub(b2r), x5i = a2i.sub(b2i);
        w1r.mul(x1r).add(w1i.neg().mul(x1i)).intoArray(ret, j + dj);
        w2r.mul(x2r).add(w2i.neg().mul(x2i)).intoArray(ret, j + dj2);
        w3r.mul(x3r).add(w3i.neg().mul(x3i)).intoArray(ret, j + dj3);
        w4r.mul(x4r).add(w4i.neg().mul(x4i)).intoArray(ret, j + dj4);
        w5r.mul(x5r).add(w5i.neg().mul(x5i)).intoArray(ret, j + dj5);
        w1r.mul(x1i).add(w1i.mul(x1r)).intoArray(ret, j + dj + im);
        w2r.mul(x2i).add(w2i.mul(x2r)).intoArray(ret, j + dj2 + im);
        w3r.mul(x3i).add(w3i.mul(x3r)).intoArray(ret, j + dj3 + im);
        w4r.mul(x4i).add(w4i.mul(x4r)).intoArray(ret, j + dj4 + im);
        w5r.mul(x5i).add(w5i.mul(x5r)).intoArray(ret, j + dj5 + im);
      }
    }
  }
}
