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

import static jdk.incubator.vector.DoubleVector.broadcast;
import static jdk.incubator.vector.DoubleVector.fromArray;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;

/**
 * The MixedRadixFactor7 class handles factors of 7 in the FFT.
 */
public class MixedRadixFactor7 extends MixedRadixFactor {

  private static final double oneThird = 1.0 / 3.0;
  private static final double sin2PI_7 = sin(2.0 * PI / 7.0);
  private static final double sin4PI_7 = sin(4.0 * PI / 7.0);
  private static final double sin6PI_7 = sin(6.0 * PI / 7.0);
  private static final double cos2PI_7 = cos(2.0 * PI / 7.0);
  private static final double cos4PI_7 = cos(4.0 * PI / 7.0);
  private static final double cos6PI_7 = cos(6.0 * PI / 7.0);
  private static final double c1 = cos2PI_7;
  private static final double c2 = cos4PI_7;
  private static final double c3 = cos6PI_7;
  private static final double v1 = (c1 + c2 + c3) * oneThird - 1.0;
  private static final double v2 = (2.0 * c1 - c2 - c3) * oneThird;
  private static final double v3 = (c1 - 2.0 * c2 + c3) * oneThird;
  private static final double v4 = (c1 + c2 - 2.0 * c3) * oneThird;

  private final int di2;
  private final int di3;
  private final int di4;
  private final int di5;
  private final int di6;
  private final int dj2;
  private final int dj3;
  private final int dj4;
  private final int dj5;
  private final int dj6;

  /**
   * Construct a MixedRadixFactor7.
   *
   * @param passConstants PassConstants.
   */
  public MixedRadixFactor7(PassConstants passConstants) {
    super(passConstants);
    di2 = 2 * di;
    di3 = 3 * di;
    di4 = 4 * di;
    di5 = 5 * di;
    di6 = 6 * di;
    dj2 = 2 * dj;
    dj3 = 3 * dj;
    dj4 = 4 * dj;
    dj5 = 5 * dj;
    dj6 = 6 * dj;
  }

  /**
   * Handle factors of 7.
   *
   * @param passData PassData.
   */
  protected void passScalar(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    // First pass of the 7-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
      final double z0r = data[i];
      final double z1r = data[i + di];
      final double z2r = data[i + di2];
      final double z3r = data[i + di3];
      final double z4r = data[i + di4];
      final double z5r = data[i + di5];
      final double z6r = data[i + di6];
      final double z0i = data[i + im];
      final double z1i = data[i + di + im];
      final double z2i = data[i + di2 + im];
      final double z3i = data[i + di3 + im];
      final double z4i = data[i + di4 + im];
      final double z5i = data[i + di5 + im];
      final double z6i = data[i + di6 + im];
      final double t0r = z1r + z6r;
      final double t0i = z1i + z6i;
      final double t1r = z1r - z6r;
      final double t1i = z1i - z6i;
      final double t2r = z2r + z5r;
      final double t2i = z2i + z5i;
      final double t3r = z2r - z5r;
      final double t3i = z2i - z5i;
      final double t4r = z4r + z3r;
      final double t4i = z4i + z3i;
      final double t5r = z4r - z3r;
      final double t5i = z4i - z3i;
      final double t6r = t2r + t0r;
      final double t6i = t2i + t0i;
      final double t7r = t5r + t3r;
      final double t7i = t5i + t3i;
      final double b0r = z0r + t6r + t4r;
      final double b0i = z0i + t6i + t4i;
      final double b1r = v1 * (t6r + t4r);
      final double b1i = v1 * (t6i + t4i);
      final double b2r = v2 * (t0r - t4r);
      final double b2i = v2 * (t0i - t4i);
      final double b3r = v3 * (t4r - t2r);
      final double b3i = v3 * (t4i - t2i);
      final double b4r = v4 * (t2r - t0r);
      final double b4i = v4 * (t2i - t0i);
      final double b5r = v5 * (t7r + t1r);
      final double b5i = v5 * (t7i + t1i);
      final double b6r = v6 * (t1r - t5r);
      final double b6i = v6 * (t1i - t5i);
      final double b7r = v7 * (t5r - t3r);
      final double b7i = v7 * (t5i - t3i);
      final double b8r = v8 * (t3r - t1r);
      final double b8i = v8 * (t3i - t1i);
      final double u0r = b0r + b1r;
      final double u0i = b0i + b1i;
      final double u1r = b2r + b3r;
      final double u1i = b2i + b3i;
      final double u2r = b4r - b3r;
      final double u2i = b4i - b3i;
      final double u3r = -b2r - b4r;
      final double u3i = -b2i - b4i;
      final double u4r = b6r + b7r;
      final double u4i = b6i + b7i;
      final double u5r = b8r - b7r;
      final double u5i = b8i - b7i;
      final double u6r = -b8r - b6r;
      final double u6i = -b8i - b6i;
      final double u7r = u0r + u1r;
      final double u7i = u0i + u1i;
      final double u8r = u0r + u2r;
      final double u8i = u0i + u2i;
      final double u9r = u0r + u3r;
      final double u9i = u0i + u3i;
      final double u10r = u4r + b5r;
      final double u10i = u4i + b5i;
      final double u11r = u5r + b5r;
      final double u11i = u5i + b5i;
      final double u12r = u6r + b5r;
      final double u12i = u6i + b5i;
      ret[j] = b0r;
      ret[j + im] = b0i;
      ret[j + dj] = u7r + u10i;
      ret[j + dj + im] = u7i - u10r;
      ret[j + dj2] = u9r + u12i;
      ret[j + dj2 + im] = u9i - u12r;
      ret[j + dj3] = u8r - u11i;
      ret[j + dj3 + im] = u8i + u11r;
      ret[j + dj4] = u8r + u11i;
      ret[j + dj4 + im] = u8i - u11r;
      ret[j + dj5] = u9r - u12i;
      ret[j + dj5 + im] = u9i + u12r;
      ret[j + dj6] = u7r - u10i;
      ret[j + dj6 + im] = u7i + u10r;
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 6;
      final double w1r = wr[index];
      final double w2r = wr[index + 1];
      final double w3r = wr[index + 2];
      final double w4r = wr[index + 3];
      final double w5r = wr[index + 4];
      final double w6r = wr[index + 5];
      final double w1i = -sign * wi[index];
      final double w2i = -sign * wi[index + 1];
      final double w3i = -sign * wi[index + 2];
      final double w4i = -sign * wi[index + 3];
      final double w5i = -sign * wi[index + 4];
      final double w6i = -sign * wi[index + 5];
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += ii, j += ii) {
        final double z0r = data[i];
        final double z1r = data[i + di];
        final double z2r = data[i + di2];
        final double z3r = data[i + di3];
        final double z4r = data[i + di4];
        final double z5r = data[i + di5];
        final double z6r = data[i + di6];
        final double z0i = data[i + im];
        final double z1i = data[i + di + im];
        final double z2i = data[i + di2 + im];
        final double z3i = data[i + di3 + im];
        final double z4i = data[i + di4 + im];
        final double z5i = data[i + di5 + im];
        final double z6i = data[i + di6 + im];
        final double t0r = z1r + z6r;
        final double t0i = z1i + z6i;
        final double t1r = z1r - z6r;
        final double t1i = z1i - z6i;
        final double t2r = z2r + z5r;
        final double t2i = z2i + z5i;
        final double t3r = z2r - z5r;
        final double t3i = z2i - z5i;
        final double t4r = z4r + z3r;
        final double t4i = z4i + z3i;
        final double t5r = z4r - z3r;
        final double t5i = z4i - z3i;
        final double t6r = t2r + t0r;
        final double t6i = t2i + t0i;
        final double t7r = t5r + t3r;
        final double t7i = t5i + t3i;
        final double b0r = z0r + t6r + t4r;
        final double b0i = z0i + t6i + t4i;
        final double b1r = v1 * (t6r + t4r);
        final double b1i = v1 * (t6i + t4i);
        final double b2r = v2 * (t0r - t4r);
        final double b2i = v2 * (t0i - t4i);
        final double b3r = v3 * (t4r - t2r);
        final double b3i = v3 * (t4i - t2i);
        final double b4r = v4 * (t2r - t0r);
        final double b4i = v4 * (t2i - t0i);
        final double b5r = v5 * (t7r + t1r);
        final double b5i = v5 * (t7i + t1i);
        final double b6r = v6 * (t1r - t5r);
        final double b6i = v6 * (t1i - t5i);
        final double b7r = v7 * (t5r - t3r);
        final double b7i = v7 * (t5i - t3i);
        final double b8r = v8 * (t3r - t1r);
        final double b8i = v8 * (t3i - t1i);
        final double u0r = b0r + b1r;
        final double u0i = b0i + b1i;
        final double u1r = b2r + b3r;
        final double u1i = b2i + b3i;
        final double u2r = b4r - b3r;
        final double u2i = b4i - b3i;
        final double u3r = -b2r - b4r;
        final double u3i = -b2i - b4i;
        final double u4r = b6r + b7r;
        final double u4i = b6i + b7i;
        final double u5r = b8r - b7r;
        final double u5i = b8i - b7i;
        final double u6r = -b8r - b6r;
        final double u6i = -b8i - b6i;
        final double u7r = u0r + u1r;
        final double u7i = u0i + u1i;
        final double u8r = u0r + u2r;
        final double u8i = u0i + u2i;
        final double u9r = u0r + u3r;
        final double u9i = u0i + u3i;
        final double u10r = u4r + b5r;
        final double u10i = u4i + b5i;
        final double u11r = u5r + b5r;
        final double u11i = u5i + b5i;
        final double u12r = u6r + b5r;
        final double u12i = u6i + b5i;
        ret[j] = b0r;
        ret[j + im] = b0i;
        multiplyAndStore(u7r + u10i, u7i - u10r, w1r, w1i, ret, j + dj, j + dj + im);
        multiplyAndStore(u9r + u12i, u9i - u12r, w2r, w2i, ret, j + dj2, j + dj2 + im);
        multiplyAndStore(u8r - u11i, u8i + u11r, w3r, w3i, ret, j + dj3, j + dj3 + im);
        multiplyAndStore(u8r + u11i, u8i - u11r, w4r, w4i, ret, j + dj4, j + dj4 + im);
        multiplyAndStore(u9r - u12i, u9i + u12r, w5r, w5i, ret, j + dj5, j + dj5 + im);
        multiplyAndStore(u7r - u10i, u7i + u10r, w6r, w6i, ret, j + dj6, j + dj6 + im);
      }
    }
  }

  /**
   * Handle factors of 6 using SIMD vectors.
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
   * Handle factors of 7.
   *
   * @param passData PassData.
   */
  private void interleaved(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    // First pass of the 7-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP, i += LENGTH, j += LENGTH) {
      final DoubleVector
          z0 = fromArray(DOUBLE_SPECIES, data, i),
          z1 = fromArray(DOUBLE_SPECIES, data, i + di),
          z2 = fromArray(DOUBLE_SPECIES, data, i + di2),
          z3 = fromArray(DOUBLE_SPECIES, data, i + di3),
          z4 = fromArray(DOUBLE_SPECIES, data, i + di4),
          z5 = fromArray(DOUBLE_SPECIES, data, i + di5),
          z6 = fromArray(DOUBLE_SPECIES, data, i + di6);
      final DoubleVector
          t0 = z1.add(z6),
          t1 = z1.sub(z6),
          t2 = z2.add(z5),
          t3 = z2.sub(z5),
          t4 = z4.add(z3),
          t5 = z4.sub(z3),
          t6 = t2.add(t0),
          t7 = t5.add(t3);
      final DoubleVector
          b0 = z0.add(t6).add(t4),
          b1 = t6.add(t4).mul(v1),
          b2 = t0.sub(t4).mul(v2),
          b3 = t4.sub(t2).mul(v3),
          b4 = t2.sub(t0).mul(v4),
          b5 = t7.add(t1).mul(v5),
          b6 = t1.sub(t5).mul(v6),
          b7 = t5.sub(t3).mul(v7),
          b8 = t3.sub(t1).mul(v8);
      final DoubleVector
          u0 = b0.add(b1),
          u1 = b2.add(b3),
          u2 = b4.sub(b3),
          u3 = b2.add(b4).neg(),
          u4 = b6.add(b7),
          u5 = b8.sub(b7),
          u6 = b8.add(b6).neg(),
          u7 = u0.add(u1),
          u8 = u0.add(u2),
          u9 = u0.add(u3),
          u10 = u4.add(b5).rearrange(SHUFFLE_RE_IM),
          u11 = u5.add(b5).rearrange(SHUFFLE_RE_IM),
          u12 = u6.add(b5).rearrange(SHUFFLE_RE_IM);
      b0.intoArray(ret, j);
      u7.add(u10.mul(NEGATE_IM)).intoArray(ret, j + dj);
      u9.add(u12.mul(NEGATE_IM)).intoArray(ret, j + dj2);
      u8.add(u11.mul(NEGATE_RE)).intoArray(ret, j + dj3);
      u8.add(u11.mul(NEGATE_IM)).intoArray(ret, j + dj4);
      u9.add(u12.mul(NEGATE_RE)).intoArray(ret, j + dj5);
      u7.add(u10.mul(NEGATE_RE)).intoArray(ret, j + dj6);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 6;
      final DoubleVector
          w1r = broadcast(DOUBLE_SPECIES, wr[index]),
          w2r = broadcast(DOUBLE_SPECIES, wr[index + 1]),
          w3r = broadcast(DOUBLE_SPECIES, wr[index + 2]),
          w4r = broadcast(DOUBLE_SPECIES, wr[index + 3]),
          w5r = broadcast(DOUBLE_SPECIES, wr[index + 4]),
          w6r = broadcast(DOUBLE_SPECIES, wr[index + 5]),
          w1i = broadcast(DOUBLE_SPECIES, -sign * wi[index]).mul(NEGATE_IM),
          w2i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 1]).mul(NEGATE_IM),
          w3i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 2]).mul(NEGATE_IM),
          w4i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 3]).mul(NEGATE_IM),
          w5i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 4]).mul(NEGATE_IM),
          w6i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 5]).mul(NEGATE_IM);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += INTERLEAVED_LOOP, i += LENGTH, j += LENGTH) {
        final DoubleVector
            z0 = fromArray(DOUBLE_SPECIES, data, i),
            z1 = fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = fromArray(DOUBLE_SPECIES, data, i + di4),
            z5 = fromArray(DOUBLE_SPECIES, data, i + di5),
            z6 = fromArray(DOUBLE_SPECIES, data, i + di6);
        final DoubleVector
            t0 = z1.add(z6),
            t1 = z1.sub(z6),
            t2 = z2.add(z5),
            t3 = z2.sub(z5),
            t4 = z4.add(z3),
            t5 = z4.sub(z3),
            t6 = t2.add(t0),
            t7 = t5.add(t3);
        final DoubleVector
            b0 = z0.add(t6).add(t4),
            b1 = t6.add(t4).mul(v1),
            b2 = t0.sub(t4).mul(v2),
            b3 = t4.sub(t2).mul(v3),
            b4 = t2.sub(t0).mul(v4),
            b5 = t7.add(t1).mul(v5),
            b6 = t1.sub(t5).mul(v6),
            b7 = t5.sub(t3).mul(v7),
            b8 = t3.sub(t1).mul(v8);
        final DoubleVector
            u0 = b0.add(b1),
            u1 = b2.add(b3),
            u2 = b4.sub(b3),
            u3 = b2.add(b4).neg(),
            u4 = b6.add(b7),
            u5 = b8.sub(b7),
            u6 = b8.add(b6).neg(),
            u7 = u0.add(u1),
            u8 = u0.add(u2),
            u9 = u0.add(u3),
            u10 = u4.add(b5).rearrange(SHUFFLE_RE_IM),
            u11 = u5.add(b5).rearrange(SHUFFLE_RE_IM),
            u12 = u6.add(b5).rearrange(SHUFFLE_RE_IM);
        b0.intoArray(ret, j);
        final DoubleVector
            x1 = u7.add(u10.mul(NEGATE_IM)),
            x2 = u9.add(u12.mul(NEGATE_IM)),
            x3 = u8.add(u11.mul(NEGATE_RE)),
            x4 = u8.add(u11.mul(NEGATE_IM)),
            x5 = u9.add(u12.mul(NEGATE_RE)),
            x6 = u7.add(u10.mul(NEGATE_RE));
        w1r.mul(x1).add(w1i.mul(x1).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj);
        w2r.mul(x2).add(w2i.mul(x2).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj2);
        w3r.mul(x3).add(w3i.mul(x3).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj3);
        w4r.mul(x4).add(w4i.mul(x4).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj4);
        w5r.mul(x5).add(w5i.mul(x5).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj5);
        w6r.mul(x6).add(w6i.mul(x6).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj6);
      }
    }
  }

  /**
   * Handle factors of 7.
   *
   * @param passData PassData.
   */
  private void blocked(PassData passData) {
    final double[] data = passData.in;
    final double[] ret = passData.out;
    int sign = passData.sign;
    int i = passData.inOffset;
    int j = passData.outOffset;
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    // First pass of the 7-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP, i += LENGTH, j += LENGTH) {
      final DoubleVector
          z0r = fromArray(DOUBLE_SPECIES, data, i),
          z1r = fromArray(DOUBLE_SPECIES, data, i + di),
          z2r = fromArray(DOUBLE_SPECIES, data, i + di2),
          z3r = fromArray(DOUBLE_SPECIES, data, i + di3),
          z4r = fromArray(DOUBLE_SPECIES, data, i + di4),
          z5r = fromArray(DOUBLE_SPECIES, data, i + di5),
          z6r = fromArray(DOUBLE_SPECIES, data, i + di6),
          z0i = fromArray(DOUBLE_SPECIES, data, i + im),
          z1i = fromArray(DOUBLE_SPECIES, data, i + di + im),
          z2i = fromArray(DOUBLE_SPECIES, data, i + di2 + im),
          z3i = fromArray(DOUBLE_SPECIES, data, i + di3 + im),
          z4i = fromArray(DOUBLE_SPECIES, data, i + di4 + im),
          z5i = fromArray(DOUBLE_SPECIES, data, i + di5 + im),
          z6i = fromArray(DOUBLE_SPECIES, data, i + di6 + im);
      final DoubleVector
          t0r = z1r.add(z6r), t0i = z1i.add(z6i),
          t1r = z1r.sub(z6r), t1i = z1i.sub(z6i),
          t2r = z2r.add(z5r), t2i = z2i.add(z5i),
          t3r = z2r.sub(z5r), t3i = z2i.sub(z5i),
          t4r = z4r.add(z3r), t4i = z4i.add(z3i),
          t5r = z4r.sub(z3r), t5i = z4i.sub(z3i),
          t6r = t2r.add(t0r), t6i = t2i.add(t0i),
          t7r = t5r.add(t3r), t7i = t5i.add(t3i);
      final DoubleVector
          b0r = z0r.add(t6r).add(t4r), b0i = z0i.add(t6i).add(t4i),
          b1r = t6r.add(t4r).mul(v1), b1i = t6i.add(t4i).mul(v1),
          b2r = t0r.sub(t4r).mul(v2), b2i = t0i.sub(t4i).mul(v2),
          b3r = t4r.sub(t2r).mul(v3), b3i = t4i.sub(t2i).mul(v3),
          b4r = t2r.sub(t0r).mul(v4), b4i = t2i.sub(t0i).mul(v4),
          b5r = t7r.add(t1r).mul(v5), b5i = t7i.add(t1i).mul(v5),
          b6r = t1r.sub(t5r).mul(v6), b6i = t1i.sub(t5i).mul(v6),
          b7r = t5r.sub(t3r).mul(v7), b7i = t5i.sub(t3i).mul(v7),
          b8r = t3r.sub(t1r).mul(v8), b8i = t3i.sub(t1i).mul(v8);
      final DoubleVector
          u0r = b0r.add(b1r), u0i = b0i.add(b1i),
          u1r = b2r.add(b3r), u1i = b2i.add(b3i),
          u2r = b4r.sub(b3r), u2i = b4i.sub(b3i),
          u3r = b2r.add(b4r).neg(), u3i = b2i.add(b4i).neg(),
          u4r = b6r.add(b7r), u4i = b6i.add(b7i),
          u5r = b8r.sub(b7r), u5i = b8i.sub(b7i),
          u6r = b8r.add(b6r).neg(), u6i = b8i.add(b6i).neg(),
          u7r = u0r.add(u1r), u7i = u0i.add(u1i),
          u8r = u0r.add(u2r), u8i = u0i.add(u2i),
          u9r = u0r.add(u3r), u9i = u0i.add(u3i),
          u10r = u4r.add(b5r), u10i = u4i.add(b5i),
          u11r = u5r.add(b5r), u11i = u5i.add(b5i),
          u12r = u6r.add(b5r), u12i = u6i.add(b5i);
      b0r.intoArray(ret, j);
      b0i.intoArray(ret, j + im);
      u7r.add(u10i).intoArray(ret, j + dj);
      u7i.sub(u10r).intoArray(ret, j + dj + im);
      u9r.add(u12i).intoArray(ret, j + dj2);
      u9i.sub(u12r).intoArray(ret, j + dj2 + im);
      u8r.sub(u11i).intoArray(ret, j + dj3);
      u8i.add(u11r).intoArray(ret, j + dj3 + im);
      u8r.add(u11i).intoArray(ret, j + dj4);
      u8i.sub(u11r).intoArray(ret, j + dj4 + im);
      u9r.sub(u12i).intoArray(ret, j + dj5);
      u9i.add(u12r).intoArray(ret, j + dj5 + im);
      u7r.sub(u10i).intoArray(ret, j + dj6);
      u7i.add(u10r).intoArray(ret, j + dj6 + im);
    }

    j += jstep;
    for (int k = 1; k < outerLoopLimit; k++, j += jstep) {
      final int index = k * 6;
      final DoubleVector
          w1r = broadcast(DOUBLE_SPECIES, wr[index]),
          w2r = broadcast(DOUBLE_SPECIES, wr[index + 1]),
          w3r = broadcast(DOUBLE_SPECIES, wr[index + 2]),
          w4r = broadcast(DOUBLE_SPECIES, wr[index + 3]),
          w5r = broadcast(DOUBLE_SPECIES, wr[index + 4]),
          w6r = broadcast(DOUBLE_SPECIES, wr[index + 5]),
          w1i = broadcast(DOUBLE_SPECIES, -sign * wi[index]),
          w2i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 1]),
          w3i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 2]),
          w4i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 3]),
          w5i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 4]),
          w6i = broadcast(DOUBLE_SPECIES, -sign * wi[index + 5]);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += BLOCK_LOOP, i += LENGTH, j += LENGTH) {
        final DoubleVector
            z0r = fromArray(DOUBLE_SPECIES, data, i),
            z1r = fromArray(DOUBLE_SPECIES, data, i + di),
            z2r = fromArray(DOUBLE_SPECIES, data, i + di2),
            z3r = fromArray(DOUBLE_SPECIES, data, i + di3),
            z4r = fromArray(DOUBLE_SPECIES, data, i + di4),
            z5r = fromArray(DOUBLE_SPECIES, data, i + di5),
            z6r = fromArray(DOUBLE_SPECIES, data, i + di6),
            z0i = fromArray(DOUBLE_SPECIES, data, i + im),
            z1i = fromArray(DOUBLE_SPECIES, data, i + di + im),
            z2i = fromArray(DOUBLE_SPECIES, data, i + di2 + im),
            z3i = fromArray(DOUBLE_SPECIES, data, i + di3 + im),
            z4i = fromArray(DOUBLE_SPECIES, data, i + di4 + im),
            z5i = fromArray(DOUBLE_SPECIES, data, i + di5 + im),
            z6i = fromArray(DOUBLE_SPECIES, data, i + di6 + im);
        final DoubleVector
            t0r = z1r.add(z6r), t0i = z1i.add(z6i),
            t1r = z1r.sub(z6r), t1i = z1i.sub(z6i),
            t2r = z2r.add(z5r), t2i = z2i.add(z5i),
            t3r = z2r.sub(z5r), t3i = z2i.sub(z5i),
            t4r = z4r.add(z3r), t4i = z4i.add(z3i),
            t5r = z4r.sub(z3r), t5i = z4i.sub(z3i),
            t6r = t2r.add(t0r), t6i = t2i.add(t0i),
            t7r = t5r.add(t3r), t7i = t5i.add(t3i);
        final DoubleVector
            b0r = z0r.add(t6r).add(t4r), b0i = z0i.add(t6i).add(t4i),
            b1r = t6r.add(t4r).mul(v1), b1i = t6i.add(t4i).mul(v1),
            b2r = t0r.sub(t4r).mul(v2), b2i = t0i.sub(t4i).mul(v2),
            b3r = t4r.sub(t2r).mul(v3), b3i = t4i.sub(t2i).mul(v3),
            b4r = t2r.sub(t0r).mul(v4), b4i = t2i.sub(t0i).mul(v4),
            b5r = t7r.add(t1r).mul(v5), b5i = t7i.add(t1i).mul(v5),
            b6r = t1r.sub(t5r).mul(v6), b6i = t1i.sub(t5i).mul(v6),
            b7r = t5r.sub(t3r).mul(v7), b7i = t5i.sub(t3i).mul(v7),
            b8r = t3r.sub(t1r).mul(v8), b8i = t3i.sub(t1i).mul(v8);
        final DoubleVector
            u0r = b0r.add(b1r), u0i = b0i.add(b1i),
            u1r = b2r.add(b3r), u1i = b2i.add(b3i),
            u2r = b4r.sub(b3r), u2i = b4i.sub(b3i),
            u3r = b2r.add(b4r).neg(), u3i = b2i.add(b4i).neg(),
            u4r = b6r.add(b7r), u4i = b6i.add(b7i),
            u5r = b8r.sub(b7r), u5i = b8i.sub(b7i),
            u6r = b8r.add(b6r).neg(), u6i = b8i.add(b6i).neg(),
            u7r = u0r.add(u1r), u7i = u0i.add(u1i),
            u8r = u0r.add(u2r), u8i = u0i.add(u2i),
            u9r = u0r.add(u3r), u9i = u0i.add(u3i),
            u10r = u4r.add(b5r), u10i = u4i.add(b5i),
            u11r = u5r.add(b5r), u11i = u5i.add(b5i),
            u12r = u6r.add(b5r), u12i = u6i.add(b5i);
        b0r.intoArray(ret, j);
        b0i.intoArray(ret, j + im);
        final DoubleVector
            x1r = u7r.add(u10i), x1i = u7i.sub(u10r),
            x2r = u9r.add(u12i), x2i = u9i.sub(u12r),
            x3r = u8r.sub(u11i), x3i = u8i.add(u11r),
            x4r = u8r.add(u11i), x4i = u8i.sub(u11r),
            x5r = u9r.sub(u12i), x5i = u9i.add(u12r),
            x6r = u7r.sub(u10i), x6i = u7i.add(u10r);
        w1r.mul(x1r).add(w1i.neg().mul(x1i)).intoArray(ret, j + dj);
        w2r.mul(x2r).add(w2i.neg().mul(x2i)).intoArray(ret, j + dj2);
        w3r.mul(x3r).add(w3i.neg().mul(x3i)).intoArray(ret, j + dj3);
        w4r.mul(x4r).add(w4i.neg().mul(x4i)).intoArray(ret, j + dj4);
        w5r.mul(x5r).add(w5i.neg().mul(x5i)).intoArray(ret, j + dj5);
        w6r.mul(x6r).add(w6i.neg().mul(x6i)).intoArray(ret, j + dj6);
        w1r.mul(x1i).add(w1i.mul(x1r)).intoArray(ret, j + dj + im);
        w2r.mul(x2i).add(w2i.mul(x2r)).intoArray(ret, j + dj2 + im);
        w3r.mul(x3i).add(w3i.mul(x3r)).intoArray(ret, j + dj3 + im);
        w4r.mul(x4i).add(w4i.mul(x4r)).intoArray(ret, j + dj4 + im);
        w5r.mul(x5i).add(w5i.mul(x5r)).intoArray(ret, j + dj5 + im);
        w6r.mul(x6i).add(w6i.mul(x6r)).intoArray(ret, j + dj6 + im);
      }
    }
  }
}
