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
import static jdk.incubator.vector.DoubleVector.broadcast;
import static jdk.incubator.vector.DoubleVector.fromArray;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The MixedRadixFactor5 class handles factors of 5 in the FFT.
 */
public class MixedRadixFactor5 extends MixedRadixFactor {

  private static final double sqrt5_4 = sqrt(5.0) / 4.0;
  private static final double sinPI_5 = sin(PI / 5.0);
  private static final double sin2PI_5 = sin(2.0 * PI / 5.0);

  private final int di2;
  private final int di3;
  private final int di4;
  private final int dj2;
  private final int dj3;
  private final int dj4;
  private final double tau = sqrt5_4;

  /**
   * Construct a MixedRadixFactor5.
   *
   * @param passConstants PassConstants.
   */
  public MixedRadixFactor5(PassConstants passConstants) {
    super(passConstants);
    di2 = 2 * di;
    di3 = 3 * di;
    di4 = 4 * di;
    dj2 = 2 * dj;
    dj3 = 3 * dj;
    dj4 = 4 * dj;
  }

  /**
   * Handle factors of 5.
   *
   * @param passData PassData.
   */
  @Override
  protected void passScalar(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    // First pass of the 5-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0r = data[i];
      final double z1r = data[i + di];
      final double z2r = data[i + di2];
      final double z3r = data[i + di3];
      final double z4r = data[i + di4];
      final double z0i = data[i + im];
      final double z1i = data[i + di + im];
      final double z2i = data[i + di2 + im];
      final double z3i = data[i + di3 + im];
      final double z4i = data[i + di4 + im];
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
      ret[j + im] = z0i + t5i;
      ret[j + dj] = t8r - t10i;
      ret[j + dj + im] = t8i + t10r;
      ret[j + dj2] = t9r - t11i;
      ret[j + dj2 + im] = t9i + t11r;
      ret[j + dj3] = t9r + t11i;
      ret[j + dj3 + im] = t9i - t11r;
      ret[j + dj4] = t8r + t10i;
      ret[j + dj4 + im] = t8i - t10r;
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
//      final double[] twids = twiddles[k];
//      final double w1r = twids[0];
//      final double w1i = -sign * twids[1];
//      final double w2r = twids[2];
//      final double w2i = -sign * twids[3];
//      final double w3r = twids[4];
//      final double w3i = -sign * twids[5];
//      final double w4r = twids[6];
//      final double w4i = -sign * twids[7];
      final int index = k * 4;
      final double w1r = wr[index];
      final double w2r = wr[index + 1];
      final double w3r = wr[index + 2];
      final double w4r = wr[index + 3];
      final double w1i = -sign * wi[index];
      final double w2i = -sign * wi[index + 1];
      final double w3i = -sign * wi[index + 2];
      final double w4i = -sign * wi[index + 3];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0r = data[i];
        final double z1r = data[i + di];
        final double z2r = data[i + di2];
        final double z3r = data[i + di3];
        final double z4r = data[i + di4];
        final double z0i = data[i + im];
        final double z1i = data[i + di + im];
        final double z2i = data[i + di2 + im];
        final double z3i = data[i + di3 + im];
        final double z4i = data[i + di4 + im];
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
        ret[j + im] = z0i + t5i;

        multiplyAndStore(t8r - t10i, t8i + t10r, w1r, w1i, ret, j + dj, j + dj + im);
        multiplyAndStore(t9r - t11i, t9i + t11r, w2r, w2i, ret, j + dj2, j + dj2 + im);
        multiplyAndStore(t9r + t11i, t9i - t11r, w3r, w3i, ret, j + dj3, j + dj3 + im);
        multiplyAndStore(t8r + t10i, t8i - t10r, w4r, w4i, ret, j + dj4, j + dj4 + im);
      }
    }
  }

  /**
   * Handle factors of 5 using SIMD vectors.
   *
   * @param passData PassData.
   */
  @Override
  protected void passSIMD(PassData passData) {
    if (!isValidSIMDWidth(simdWidth)) {
      passScalar(passData);
    } else {
      if (im == 1) {
        interleaved(passData);
      } else {
        blocked(passData);
      }
    }
  }

  /**
   * Handle factors of 5.
   *
   * @param passData PassData.
   */
  protected void interleaved(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    // First pass of the 5-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP, i += LENGTH, j += LENGTH) {
      final DoubleVector
          z0 = fromArray(DOUBLE_SPECIES, data, i),
          z1 = fromArray(DOUBLE_SPECIES, data, i + di),
          z2 = fromArray(DOUBLE_SPECIES, data, i + di2),
          z3 = fromArray(DOUBLE_SPECIES, data, i + di3),
          z4 = fromArray(DOUBLE_SPECIES, data, i + di4);
      final DoubleVector
          t1 = z1.add(z4),
          t2 = z2.add(z3),
          t3 = z1.sub(z4),
          t4 = z2.sub(z3),
          t5 = t1.add(t2),
          t6 = t1.sub(t2).mul(tau),
          t7 = t5.mul(-0.25).add(z0),
          t8 = t7.add(t6),
          t9 = t7.sub(t6),
          t10 = t3.mul(sin2PI_5s).add(t4.mul(sinPI_5s)).rearrange(SHUFFLE_RE_IM),
          t11 = t4.mul(-sin2PI_5s).add(t3.mul(sinPI_5s)).rearrange(SHUFFLE_RE_IM);
      z0.add(t5).intoArray(ret, j);
      t8.add(t10.mul(NEGATE_RE)).intoArray(ret, j + dj);
      t9.add(t11.mul(NEGATE_RE)).intoArray(ret, j + dj2);
      t9.add(t11.mul(NEGATE_IM)).intoArray(ret, j + dj3);
      t8.add(t10.mul(NEGATE_IM)).intoArray(ret, j + dj4);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 4;
      final DoubleVector
          w1r = broadcast(DOUBLE_SPECIES, wr[index]),
          w2r = broadcast(DOUBLE_SPECIES, wr[index + 1]),
          w3r = broadcast(DOUBLE_SPECIES, wr[index + 2]),
          w4r = broadcast(DOUBLE_SPECIES, wr[index + 3]),
          w1i = broadcast(DOUBLE_SPECIES, -sign * wi[index]).mul(NEGATE_IM),
          w2i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 1]).mul(NEGATE_IM),
          w3i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 2]).mul(NEGATE_IM),
          w4i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 3]).mul(NEGATE_IM);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP, i += LENGTH, j += LENGTH) {
        final DoubleVector
            z0 = fromArray(DOUBLE_SPECIES, data, i),
            z1 = fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = fromArray(DOUBLE_SPECIES, data, i + di4);
        final DoubleVector
            t1 = z1.add(z4),
            t2 = z2.add(z3),
            t3 = z1.sub(z4),
            t4 = z2.sub(z3),
            t5 = t1.add(t2),
            t6 = t1.sub(t2).mul(tau),
            t7 = t5.mul(-0.25).add(z0),
            t8 = t7.add(t6),
            t9 = t7.sub(t6),
            t10 = t3.mul(sin2PI_5s).add(t4.mul(sinPI_5s)).rearrange(SHUFFLE_RE_IM),
            t11 = t4.mul(-sin2PI_5s).add(t3.mul(sinPI_5s)).rearrange(SHUFFLE_RE_IM);
        z0.add(t5).intoArray(ret, j);
        final DoubleVector
            x1 = t8.add(t10.mul(NEGATE_RE)),
            x2 = t9.add(t11.mul(NEGATE_RE)),
            x3 = t9.add(t11.mul(NEGATE_IM)),
            x4 = t8.add(t10.mul(NEGATE_IM));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj2);
        w3r.mul(x3).add(w3i.mul(x3).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj3);
        w4r.mul(x4).add(w4i.mul(x4).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj4);
      }
    }
  }

  /**
   * Handle factors of 5.
   *
   * @param passData PassData.
   */
  protected void blocked(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    // First pass of the 5-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP, i += LENGTH, j += LENGTH) {
      final DoubleVector
          z0r = fromArray(DOUBLE_SPECIES, data, i),
          z1r = fromArray(DOUBLE_SPECIES, data, i + di),
          z2r = fromArray(DOUBLE_SPECIES, data, i + di2),
          z3r = fromArray(DOUBLE_SPECIES, data, i + di3),
          z4r = fromArray(DOUBLE_SPECIES, data, i + di4),
          z0i = fromArray(DOUBLE_SPECIES, data, i + im),
          z1i = fromArray(DOUBLE_SPECIES, data, i + di + im),
          z2i = fromArray(DOUBLE_SPECIES, data, i + di2 + im),
          z3i = fromArray(DOUBLE_SPECIES, data, i + di3 + im),
          z4i = fromArray(DOUBLE_SPECIES, data, i + di4 + im);
      final DoubleVector
          t1r = z1r.add(z4r),
          t1i = z1i.add(z4i),
          t2r = z2r.add(z3r),
          t2i = z2i.add(z3i),
          t3r = z1r.sub(z4r),
          t3i = z1i.sub(z4i),
          t4r = z2r.sub(z3r),
          t4i = z2i.sub(z3i),
          t5r = t1r.add(t2r),
          t5i = t1i.add(t2i),
          t6r = t1r.sub(t2r).mul(tau),
          t6i = t1i.sub(t2i).mul(tau),
          t7r = t5r.mul(-0.25).add(z0r),
          t7i = t5i.mul(-0.25).add(z0i),
          t8r = t7r.add(t6r),
          t8i = t7i.add(t6i),
          t9r = t7r.sub(t6r),
          t9i = t7i.sub(t6i),
          t10r = t3r.mul(sin2PI_5s).add(t4r.mul(sinPI_5s)),
          t10i = t3i.mul(sin2PI_5s).add(t4i.mul(sinPI_5s)),
          t11r = t4r.mul(-sin2PI_5s).add(t3r.mul(sinPI_5s)),
          t11i = t4i.mul(-sin2PI_5s).add(t3i.mul(sinPI_5s));
      z0r.add(t5r).intoArray(ret, j);
      z0i.add(t5i).intoArray(ret, j + im);
      t8r.sub(t10i).intoArray(ret, j + dj);
      t8i.add(t10r).intoArray(ret, j + dj + im);
      t9r.sub(t11i).intoArray(ret, j + dj2);
      t9i.add(t11r).intoArray(ret, j + dj2 + im);
      t9r.add(t11i).intoArray(ret, j + dj3);
      t9i.sub(t11r).intoArray(ret, j + dj3 + im);
      t8r.add(t10i).intoArray(ret, j + dj4);
      t8i.sub(t10r).intoArray(ret, j + dj4 + im);
    }
    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 4;
      final DoubleVector
          w1r = broadcast(DOUBLE_SPECIES, wr[index]),
          w2r = broadcast(DOUBLE_SPECIES, wr[index + 1]),
          w3r = broadcast(DOUBLE_SPECIES, wr[index + 2]),
          w4r = broadcast(DOUBLE_SPECIES, wr[index + 3]),
          w1i = broadcast(DOUBLE_SPECIES, -sign * wi[index]),
          w2i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 1]),
          w3i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 2]),
          w4i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 3]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP, i += LENGTH, j += LENGTH) {
        final DoubleVector
            z0r = fromArray(DOUBLE_SPECIES, data, i),
            z1r = fromArray(DOUBLE_SPECIES, data, i + di),
            z2r = fromArray(DOUBLE_SPECIES, data, i + di2),
            z3r = fromArray(DOUBLE_SPECIES, data, i + di3),
            z4r = fromArray(DOUBLE_SPECIES, data, i + di4),
            z0i = fromArray(DOUBLE_SPECIES, data, i + im),
            z1i = fromArray(DOUBLE_SPECIES, data, i + di + im),
            z2i = fromArray(DOUBLE_SPECIES, data, i + di2 + im),
            z3i = fromArray(DOUBLE_SPECIES, data, i + di3 + im),
            z4i = fromArray(DOUBLE_SPECIES, data, i + di4 + im);
        final DoubleVector
            t1r = z1r.add(z4r),
            t1i = z1i.add(z4i),
            t2r = z2r.add(z3r),
            t2i = z2i.add(z3i),
            t3r = z1r.sub(z4r),
            t3i = z1i.sub(z4i),
            t4r = z2r.sub(z3r),
            t4i = z2i.sub(z3i),
            t5r = t1r.add(t2r),
            t5i = t1i.add(t2i),
            t6r = t1r.sub(t2r).mul(tau),
            t6i = t1i.sub(t2i).mul(tau),
            t7r = t5r.mul(-0.25).add(z0r),
            t7i = t5i.mul(-0.25).add(z0i),
            t8r = t7r.add(t6r),
            t8i = t7i.add(t6i),
            t9r = t7r.sub(t6r),
            t9i = t7i.sub(t6i),
            t10r = t3r.mul(sin2PI_5s).add(t4r.mul(sinPI_5s)),
            t10i = t3i.mul(sin2PI_5s).add(t4i.mul(sinPI_5s)),
            t11r = t4r.mul(-sin2PI_5s).add(t3r.mul(sinPI_5s)),
            t11i = t4i.mul(-sin2PI_5s).add(t3i.mul(sinPI_5s));
        z0r.add(t5r).intoArray(ret, j);
        z0i.add(t5i).intoArray(ret, j + im);
        final DoubleVector
            x1r = t8r.sub(t10i), x1i = t8i.add(t10r),
            x2r = t9r.sub(t11i), x2i = t9i.add(t11r),
            x3r = t9r.add(t11i), x3i = t9i.sub(t11r),
            x4r = t8r.add(t10i), x4i = t8i.sub(t10r);
        w1r.mul(x1r).add(w1i.neg().mul(x1i)).intoArray(ret, j + dj);
        w2r.mul(x2r).add(w2i.neg().mul(x2i)).intoArray(ret, j + dj2);
        w3r.mul(x3r).add(w3i.neg().mul(x3i)).intoArray(ret, j + dj3);
        w4r.mul(x4r).add(w4i.neg().mul(x4i)).intoArray(ret, j + dj4);
        w1r.mul(x1i).add(w1i.mul(x1r)).intoArray(ret, j + dj + im);
        w2r.mul(x2i).add(w2i.mul(x2r)).intoArray(ret, j + dj2 + im);
        w3r.mul(x3i).add(w3i.mul(x3r)).intoArray(ret, j + dj3 + im);
        w4r.mul(x4i).add(w4i.mul(x4r)).intoArray(ret, j + dj4 + im);
      }
    }
  }
}
