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
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;

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
   * @param passData the data.
   */
  protected void pass(PassData passData) {
    final int sign = passData.sign();
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    int i = passData.inOffset();
    int j = passData.outOffset();
    final double[] data = passData.in();
    final double[] ret = passData.out();

    // First pass of the 7-point FFT has no twiddle factors.
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
      idi += di;
      final double z6r = data[idi];
      final double z6i = data[idi + 1];
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
      ret[j + 1] = b0i;
      int jdj = j + dj;
      ret[jdj] = u7r + u10i;
      ret[jdj + 1] = u7i - u10r;
      jdj += dj;
      ret[jdj] = u9r + u12i;
      ret[jdj + 1] = u9i - u12r;
      jdj += dj;
      ret[jdj] = u8r - u11i;
      ret[jdj + 1] = u8i + u11r;
      jdj += dj;
      ret[jdj] = u8r + u11i;
      ret[jdj + 1] = u8i - u11r;
      jdj += dj;
      ret[jdj] = u9r - u12i;
      ret[jdj + 1] = u9i + u12r;
      jdj += dj;
      ret[jdj] = u7r - u10i;
      ret[jdj + 1] = u7i + u10r;
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
      final double w6r = twids[10];
      final double w6i = -sign * twids[11];
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
        idi += di;
        final double z6r = data[idi];
        final double z6i = data[idi + 1];
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
        ret[j + 1] = b0i;
        double xr = u7r + u10i;
        double xi = u7i - u10r;
        int jdj = j + dj;
        ret[jdj] = fma(w1r, xr, -w1i * xi);
        ret[jdj + 1] = fma(w1r, xi, w1i * xr);
        xr = u9r + u12i;
        xi = u9i - u12r;
        jdj += dj;
        ret[jdj] = fma(w2r, xr, -w2i * xi);
        ret[jdj + 1] = fma(w2r, xi, w2i * xr);
        xr = u8r - u11i;
        xi = u8i + u11r;
        jdj += dj;
        ret[jdj] = fma(w3r, xr, -w3i * xi);
        ret[jdj + 1] = fma(w3r, xi, w3i * xr);
        xr = u8r + u11i;
        xi = u8i - u11r;
        jdj += dj;
        ret[jdj] = fma(w4r, xr, -w4i * xi);
        ret[jdj + 1] = fma(w4r, xi, w4i * xr);
        xr = u9r - u12i;
        xi = u9i + u12r;
        jdj += dj;
        ret[jdj] = fma(w5r, xr, -w5i * xi);
        ret[jdj + 1] = fma(w5r, xi, w5i * xr);
        xr = u7r - u10i;
        xi = u7i + u10r;
        jdj += dj;
        ret[jdj] = fma(w6r, xr, -w6i * xi);
        ret[jdj + 1] = fma(w6r, xi, w6i * xr);
      }
    }
  }

  /**
   * Handle factors of 7.
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
    final double s1 = (-sign) * sin2PI_7;
    final double s2 = (-sign) * sin4PI_7;
    final double s3 = (-sign) * sin6PI_7;
    final double v5 = (s1 + s2 - s3) * oneThird;
    final double v6 = (2.0 * s1 - s2 + s3) * oneThird;
    final double v7 = (s1 - 2.0 * s2 - s3) * oneThird;
    final double v8 = (s1 + s2 + 2.0 * s3) * oneThird;
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();

    // First pass of the 7-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
      DoubleVector
          z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
          z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
          z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
          z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
          z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4),
          z5 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di5),
          z6 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di6);
      DoubleVector
          t0 = z1.add(z6),
          t1 = z1.sub(z6),
          t2 = z2.add(z5),
          t3 = z2.sub(z5),
          t4 = z4.add(z3),
          t5 = z4.sub(z3),
          t6 = t2.add(t0),
          t7 = t5.add(t3);
      DoubleVector
          b0 = z0.add(t6).add(t4),
          b1 = t6.add(t4).mul(v1),
          b2 = t0.sub(t4).mul(v2),
          b3 = t4.sub(t2).mul(v3),
          b4 = t2.sub(t0).mul(v4),
          b5 = t7.add(t1).mul(v5),
          b6 = t1.sub(t5).mul(v6),
          b7 = t5.sub(t3).mul(v7),
          b8 = t3.sub(t1).mul(v8);
      DoubleVector
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
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(NEGATE_IM),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(NEGATE_IM),
          w3r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(NEGATE_IM),
          w4r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[7]).mul(NEGATE_IM),
          w5r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[8]),
          w5i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[9]).mul(NEGATE_IM),
          w6r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[10]),
          w6i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[11]).mul(NEGATE_IM);
      for (int k1 = 0; k1 < innerLoopLimit; k1++, i += 2, j += 2) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4),
            z5 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di5),
            z6 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di6);
        DoubleVector
            t0 = z1.add(z6),
            t1 = z1.sub(z6),
            t2 = z2.add(z5),
            t3 = z2.sub(z5),
            t4 = z4.add(z3),
            t5 = z4.sub(z3),
            t6 = t2.add(t0),
            t7 = t5.add(t3);
        DoubleVector
            b0 = z0.add(t6).add(t4),
            b1 = t6.add(t4).mul(v1),
            b2 = t0.sub(t4).mul(v2),
            b3 = t4.sub(t2).mul(v3),
            b4 = t2.sub(t0).mul(v4),
            b5 = t7.add(t1).mul(v5),
            b6 = t1.sub(t5).mul(v6),
            b7 = t5.sub(t3).mul(v7),
            b8 = t3.sub(t1).mul(v8);
        DoubleVector
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
        DoubleVector x = u7.add(u10.mul(NEGATE_IM));
        w1r.mul(x).add(w1i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj);
        x = u9.add(u12.mul(NEGATE_IM));
        w2r.mul(x).add(w2i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj2);
        x = u8.add(u11.mul(NEGATE_RE));
        w3r.mul(x).add(w3i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj3);
        x = u8.add(u11.mul(NEGATE_IM));
        w4r.mul(x).add(w4i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj4);
        x = u9.add(u12.mul(NEGATE_RE));
        w5r.mul(x).add(w5i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj5);
        x = u7.add(u10.mul(NEGATE_RE));
        w6r.mul(x).add(w6i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj6);
      }
    }
  }

}
