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
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

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
   * @param passData the data.
   */
  protected void pass(PassData passData) {
    final int sign = passData.sign();
    final double tau = sqrt5_4;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    int i = passData.inOffset();
    int j = passData.outOffset();
    final double[] data = passData.in();
    final double[] ret = passData.out();

    // First pass of the 5-point FFT has no twiddle factors.
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
      ret[j + 1] = z0i + t5i;
      int jdj = j + dj;
      ret[jdj] = t8r - t10i;
      ret[jdj + 1] = t8i + t10r;
      jdj += dj;
      ret[jdj] = t9r - t11i;
      ret[jdj + 1] = t9i + t11r;
      jdj += dj;
      ret[jdj] = t9r + t11i;
      ret[jdj + 1] = t9i - t11r;
      jdj += dj;
      ret[jdj] = t8r + t10i;
      ret[jdj + 1] = t8i - t10r;
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
        ret[j + 1] = z0i + t5i;
        double xr = t8r - t10i;
        double xi = t8i + t10r;
        int jdj = j + dj;
        ret[jdj] = fma(w1r, xr, -w1i * xi);
        ret[jdj + 1] = fma(w1r, xi, w1i * xr);
        xr = t9r - t11i;
        xi = t9i + t11r;
        jdj += dj;
        ret[jdj] = fma(w2r, xr, -w2i * xi);
        ret[jdj + 1] = fma(w2r, xi, w2i * xr);
        xr = t9r + t11i;
        xi = t9i - t11r;
        jdj += dj;
        ret[jdj] = fma(w3r, xr, -w3i * xi);
        ret[jdj + 1] = fma(w3r, xi, w3i * xr);
        xr = t8r + t10i;
        xi = t8i - t10r;
        jdj += dj;
        ret[jdj] = fma(w4r, xr, -w4i * xi);
        ret[jdj + 1] = fma(w4r, xi, w4i * xr);
      }
    }
  }

  /**
   * Handle factors of 5.
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
    final double tau = sqrt5_4;
    final double sin2PI_5s = sign * sin2PI_5;
    final double sinPI_5s = sign * sinPI_5;
    final double[] data = passData.in();
    final double[] ret = passData.out();
    int i = passData.inOffset();
    int j = passData.outOffset();

    // First pass of the 5-point FFT has no twiddle factors.
    for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT, i += SPECIES_LENGTH, j += SPECIES_LENGTH) {
      DoubleVector
          z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
          z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
          z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
          z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
          z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4);
      DoubleVector
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
      final double[] twids = twiddles[k];
      DoubleVector
          w1r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[0]),
          w1i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[1]).mul(NEGATE_IM),
          w2r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[2]),
          w2i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[3]).mul(NEGATE_IM),
          w3r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[4]),
          w3i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[5]).mul(NEGATE_IM),
          w4r = DoubleVector.broadcast(DOUBLE_SPECIES, twids[6]),
          w4i = DoubleVector.broadcast(DOUBLE_SPECIES, -sign * twids[7]).mul(NEGATE_IM);
      for (int k1 = 0; k1 < innerLoopLimit; k1 += LOOP_INCREMENT, i += SPECIES_LENGTH, j += SPECIES_LENGTH) {
        DoubleVector
            z0 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i),
            z1 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di),
            z2 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di2),
            z3 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di3),
            z4 = DoubleVector.fromArray(DOUBLE_SPECIES, data, i + di4);
        DoubleVector
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
        DoubleVector x = t8.add(t10.mul(NEGATE_RE));
        w1r.mul(x).add(w1i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj);
        x = t9.add(t11.mul(NEGATE_RE));
        w2r.mul(x).add(w2i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj2);
        x = t9.add(t11.mul(NEGATE_IM));
        w3r.mul(x).add(w3i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj3);
        x = t8.add(t10.mul(NEGATE_IM));
        w4r.mul(x).add(w4i.mul(x).rearrange(SHUFFLE_RE_IM)).intoArray(ret, j + dj4);
      }
    }
  }

}
