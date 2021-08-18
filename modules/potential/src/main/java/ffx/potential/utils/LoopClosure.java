// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.potential.utils;

import static ffx.numerics.math.DoubleMath.X;
import static ffx.numerics.math.DoubleMath.dot;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.atan2;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.tan;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.potential.bonded.SturmMethod;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * LoopClosure class.
 *
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class LoopClosure {

  private static final Logger logger = Logger.getLogger(LoopClosure.class.getName());

  private final SturmMethod sturmMethod;

  private int max_soln = 16;
  private int[] deg_pol = new int[1];
  private double[] len0 = new double[6];
  private double[] b_ang0 = new double[7];
  private double[] t_ang0 = new double[2];
  private double aa13_min_sqr;
  private double aa13_max_sqr;
  private double[][] delta = new double[4][1];
  private double[][] xi = new double[3][1];
  private double[][] eta = new double[3][1];
  private double[] alpha = new double[3];
  private double[] theta = new double[3];
  private double[] cos_alpha = new double[3];
  private double[] sin_alpha = new double[3];
  private double[] cos_theta = new double[3];
  private double[] cos_delta = new double[4];
  private double[] sin_delta = new double[4];
  private double[] cos_xi = new double[3];
  private double[] cos_eta = new double[3];
  private double[] sin_xi = new double[3];
  private double[] sin_eta = new double[3];
  private double[] r_a1a3 = new double[3];
  private double[] r_a1n1 = new double[3];
  private double[] r_a3c3 = new double[3];
  private double[] b_a1a3 = new double[3];
  private double[] b_a1n1 = new double[3];
  private double[] b_a3c3 = new double[3];
  private double[] len_na = new double[3];
  private double[] len_ac = new double[3];
  private double[] len_aa = new double[3];
  private double[][] C0 = new double[3][3];
  private double[][] C1 = new double[3][3];
  private double[][] C2 = new double[3][3];
  private double[][] Q = new double[5][17];
  private double[][] R = new double[3][17];

  /** LoopClosure Constructor. */
  public LoopClosure() {
    deg_pol[0] = 16;
    sturmMethod = new SturmMethod();
    initializeLoopClosure();
  }

  /**
   * Calculate bond angle.
   *
   * @param r1 an array of {@link double} objects.
   * @param r2 an array of {@link double} objects.
   * @param angle an array of {@link double} objects.
   */
  public void calcBondAngle(double[] r1, double[] r2, double[] angle) {
    double arg = dot(r1, r2);
    arg = sign(min(abs(arg), 1.0), arg);
    angle[0] = acos(arg);
  }

  /**
   * Calculate dihedral angle.
   *
   * @param r1 an array of {@link double} objects.
   * @param r2 an array of {@link double} objects.
   * @param r3 an array of {@link double} objects.
   * @param angle an array of {@link double} objects.
   */
  public void calcDihedralAngle(double[] r1, double[] r2, double[] r3, double[] angle) {
    double[] p = new double[3];
    double[] q = new double[3];
    double[] s = new double[3];

    X(r1, r2, p);
    X(r2, r3, q);
    X(r3, r1, s);
    double arg = dot(p, q) / sqrt(dot(p, p) * dot(q, q));
    arg = sign(min(abs(arg), 1.0), arg);
    angle[0] = sign(acos(arg), dot(s, r2));
  }

  /**
   * Calculate T1.
   *
   * @param t0 a double.
   * @param t2 a double.
   * @return return a double T1.
   */
  public double calcT1(double t0, double t2) {
    double t0_2 = t0 * t0;
    double t2_2 = t2 * t2;
    double U11 = C0[0][0] + C0[0][1] * t0 + C0[0][2] * t0_2;
    double U12 = C1[0][0] + C1[0][1] * t0 + C1[0][2] * t0_2;
    double U13 = C2[0][0] + C2[0][1] * t0 + C2[0][2] * t0_2;
    double U31 = C0[1][0] + C0[1][1] * t2 + C0[1][2] * t2_2;
    double U32 = C1[1][0] + C1[1][1] * t2 + C1[1][2] * t2_2;
    double U33 = C2[1][0] + C2[1][1] * t2 + C2[1][2] * t2_2;
    double tmp_value = (U31 * U13 - U11 * U33) / (U12 * U33 - U13 * U32);

    return tmp_value;
  }

  /**
   * Calculate T2.
   *
   * @param t0 a double.
   * @return a double T2.
   */
  public double calcT2(double t0) {
    double t0_2 = t0 * t0;
    double t0_3 = t0_2 * t0;
    double t0_4 = t0_3 * t0;

    double A0 = Q[0][0] + Q[0][1] * t0 + Q[0][2] * t0_2 + Q[0][3] * t0_3 + Q[0][4] * t0_4;
    double A1 = Q[1][0] + Q[1][1] * t0 + Q[1][2] * t0_2 + Q[1][3] * t0_3 + Q[1][4] * t0_4;
    double A2 = Q[2][0] + Q[2][1] * t0 + Q[2][2] * t0_2 + Q[2][3] * t0_3 + Q[2][4] * t0_4;
    double A3 = Q[3][0] + Q[3][1] * t0 + Q[3][2] * t0_2 + Q[3][3] * t0_3 + Q[3][4] * t0_4;
    double A4 = Q[4][0] + Q[4][1] * t0 + Q[4][2] * t0_2 + Q[4][3] * t0_3 + Q[4][4] * t0_4;

    double B0 = R[0][0] + R[0][1] * t0 + R[0][2] * t0_2;
    double B1 = R[1][0] + R[1][1] * t0 + R[1][2] * t0_2;
    double B2 = R[2][0] + R[2][1] * t0 + R[2][2] * t0_2;

    double B2_2 = B2 * B2;
    double B2_3 = B2_2 * B2;

    double K0 = A2 * B2 - A4 * B0;
    double K1 = A3 * B2 - A4 * B1;
    double K2 = A1 * B2_2 - K1 * B0;
    double K3 = K0 * B2 - K1 * B1;
    double tmp_value = (K3 * B0 - A0 * B2_3) / (K2 * B2 - K3 * B1);

    return tmp_value;
  }

  /**
   * Get coordinates from polynomial roots.
   *
   * @param n_soln an array of {@link int} objects.
   * @param roots an array of {@link double} objects.
   * @param r_n1 an array of {@link double} objects.
   * @param r_a1 an array of {@link double} objects.
   * @param r_a3 an array of {@link double} objects.
   * @param r_c3 an array of {@link double} objects.
   * @param r_soln_n an array of {@link double} objects.
   * @param r_soln_a an array of {@link double} objects.
   * @param r_soln_c an array of {@link double} objects.
   */
  public void getCoordsFromPolyRoots(
      int[] n_soln,
      double[] roots,
      double[] r_n1,
      double[] r_a1,
      double[] r_a3,
      double[] r_c3,
      double[][][] r_soln_n,
      double[][][] r_soln_a,
      double[][][] r_soln_c) {
    double[] ex = new double[3];
    double[] ey = new double[3];
    double[] ez = new double[3];
    double[] b_a1a2 = new double[3];
    double[] b_a3a2 = new double[3];
    double[] r_tmp = new double[3];
    double[][] p_s = new double[3][3];
    double[][] s1 = new double[3][3];
    double[][] s2 = new double[3][3];
    double[][] p_t = new double[3][3];
    double[][] t1 = new double[3][3];
    double[][] t2 = new double[3][3];
    double[][] p_s_c = new double[3][3];
    double[][] s1_s = new double[3][3];
    double[][] s2_s = new double[3][3];
    double[][] p_t_c = new double[3][3];
    double[][] t1_s = new double[3][3];
    double[][] t2_s = new double[3][3];
    double[] angle = new double[1];
    double[] half_tan = new double[3];
    double[] cos_tau = new double[4];
    double[] sin_tau = new double[4];
    double[] cos_sig = new double[3];
    double[] sin_sig = new double[3];
    double[] r_s = new double[3];
    double[] r_t = new double[3];
    double[] r0 = new double[3];
    double[][] r_n = new double[3][3];
    double[][] r_a = new double[3][3];
    double[][] r_c = new double[3][3];
    double[] p = new double[4];
    double[][] Us = new double[3][3];
    double[] rr_a1c1 = new double[3];
    double[] rr_c1n2 = new double[3];
    double[] rr_n2a2 = new double[3];
    double[] rr_a2c2 = new double[3];
    double[] rr_c2n3 = new double[3];
    double[] rr_n3a3 = new double[3];
    double[] rr_a1a2 = new double[3];
    double[] rr_a2a3 = new double[3];
    double[] a3a1a2 = new double[1];
    double[] a2a3a1 = new double[1];
    double[] n1a1c1 = new double[1];
    double[] n2a2c2 = new double[1];
    double[] n3a3c3 = new double[1];
    double[] a1c1n2a2 = new double[1];
    double[] a2c2n3a3 = new double[1];
    double[] ex_tmp = new double[3];
    double[] tmp_array = new double[3];
    double[] tmp_array1 = new double[3];
    double[] tmp_array2 = new double[3];
    double[] tmp_array3 = new double[3];
    double[] mat1 = new double[3];
    double[] mat2 = new double[3];
    double[] mat3 = new double[3];
    double[] mat4 = new double[3];
    double[] mat5 = new double[3];
    double[] mat11 = new double[3];
    double[] mat22 = new double[3];
    double[] mat33 = new double[3];
    double[] mat44 = new double[3];
    double[] mat55 = new double[3];

    if (n_soln[0] == 0) {
      return;
    }

    for (int i = 0; i < 3; i++) {
      ex[i] = b_a1a3[i];
    }
    X(r_a1n1, ex, ez);
    double tmp_value = sqrt(dot(ez, ez));
    for (int i = 0; i < 3; i++) {
      ez[i] = ez[i] / tmp_value;
    }
    X(ez, ex, ey);

    for (int i = 0; i < 3; i++) {
      b_a1a2[i] = -cos_alpha[0] * ex[i] + sin_alpha[0] * ey[i];
      b_a3a2[i] = cos_alpha[2] * ex[i] + sin_alpha[2] * ey[i];
    }

    for (int i = 0; i < 3; i++) {
      p_s[0][i] = -ex[i];
      s1[0][i] = ez[i];
      s2[0][i] = ey[i];
      p_t[0][i] = b_a1a2[i];
      t1[0][i] = ez[i];
      t2[0][i] = sin_alpha[0] * ex[i] + cos_alpha[0] * ey[i];
    }

    for (int i = 0; i < 3; i++) {
      p_s[1][i] = -b_a1a2[i];
      s1[1][i] = -ez[i];
      s2[1][i] = t2[0][i];
      p_t[1][i] = -b_a3a2[i];
      t1[1][i] = -ez[i];
      t2[1][i] = sin_alpha[2] * ex[i] - cos_alpha[2] * ey[i];
    }

    for (int i = 0; i < 3; i++) {
      p_s[2][i] = b_a3a2[i];
      s2[2][i] = t2[1][i];
      s1[2][i] = ez[i];
      p_t[2][i] = ex[i];
      t1[2][i] = ez[i];
      t2[2][i] = -ey[i];
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        p_s_c[i][j] = p_s[i][j] * cos_xi[i];
        s1_s[i][j] = s1[i][j] * sin_xi[i];
        s2_s[i][j] = s2[i][j] * sin_xi[i];
        p_t_c[i][j] = p_t[i][j] * cos_eta[i];
        t1_s[i][j] = t1[i][j] * sin_eta[i];
        t2_s[i][j] = t2[i][j] * sin_eta[i];
      }
    }

    for (int i = 0; i < 3; i++) {
      r_tmp[i] = (r_a1n1[i] / len_na[0] - p_s_c[0][i]) / sin_xi[0];
    }
    calcBondAngle(s1[0], r_tmp, angle);
    double sig1Init = sign(angle[0], dot(r_tmp, s2[0]));

    for (int i = 0; i < 3; i++) {
      r_a[0][i] = r_a1[i];
      r_a[1][i] = r_a1[i] + len_aa[1] * b_a1a2[i];
      r_a[2][i] = r_a3[i];
      r0[i] = r_a1[i];
    }

    for (int i_soln = 0; i_soln < n_soln[0]; i_soln++) {
      half_tan[2] = roots[i_soln];
      half_tan[1] = calcT2(half_tan[2]);
      half_tan[0] = calcT1(half_tan[2], half_tan[1]);

      for (int i = 1; i <= 3; i++) {
        double ht = half_tan[i - 1];
        double tmp = 1.0 + ht * ht;
        cos_tau[i] = (1.0 - ht * ht) / tmp;
        sin_tau[i] = 2.0 * ht / tmp;
      }

      cos_tau[0] = cos_tau[3];
      sin_tau[0] = sin_tau[3];

      for (int i = 0; i < 3; i++) {
        cos_sig[i] = cos_delta[i] * cos_tau[i] + sin_delta[i] * sin_tau[i];
        sin_sig[i] = sin_delta[i] * cos_tau[i] - cos_delta[i] * sin_tau[i];
      }

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          r_s[j] = p_s_c[i][j] + cos_sig[i] * s1_s[i][j] + sin_sig[i] * s2_s[i][j];
          r_t[j] = p_t_c[i][j] + cos_tau[i + 1] * t1_s[i][j] + sin_tau[i + 1] * t2_s[i][j];
          r_n[i][j] = r_s[j] * len_na[i] + r_a[i][j];
          r_c[i][j] = r_t[j] * len_ac[i] + r_a[i][j];
        }
      }

      double sig1 = atan2(sin_sig[0], cos_sig[0]);
      ex_tmp[0] = -ex[0];
      ex_tmp[1] = -ex[1];
      ex_tmp[2] = -ex[2];
      tmp_value = -(sig1 - sig1Init) * 0.25;
      quaternion(ex_tmp, tmp_value, p);

      rotationMatrix(p, Us);

      for (int i = 0; i < 3; i++) {
        mat11[i] = r_c[0][i] - r0[i];
        mat22[i] = r_n[1][i] - r0[i];
        mat33[i] = r_a[1][i] - r0[i];
        mat44[i] = r_c[1][i] - r0[i];
        mat55[i] = r_n[2][i] - r0[i];
      }

      matMul(Us, mat11, mat1);
      matMul(Us, mat22, mat2);
      matMul(Us, mat33, mat3);
      matMul(Us, mat44, mat4);
      matMul(Us, mat55, mat5);

      for (int i = 0; i < 3; i++) {
        r_soln_n[i_soln][0][i] = r_n1[i];
        r_soln_a[i_soln][0][i] = r_a1[i];
        r_soln_c[i_soln][0][i] = mat1[i] + r0[i];
        r_soln_n[i_soln][1][i] = mat2[i] + r0[i];
        r_soln_a[i_soln][1][i] = mat3[i] + r0[i];
        r_soln_c[i_soln][1][i] = mat4[i] + r0[i];
        r_soln_n[i_soln][2][i] = mat5[i] + r0[i];
        r_soln_a[i_soln][2][i] = r_a3[i];
        r_soln_c[i_soln][2][i] = r_c3[i];
      }

      if (logger.isLoggable(Level.FINE)) {
        StringBuilder string1 = new StringBuilder();
        string1.append(String.format("roots: t0\t\t t2\t\t t1\t\t %d\n", i_soln));
        string1.append(
            String.format("%15.6f %15.6f %15.6f\n", half_tan[2], half_tan[1], half_tan[0]));
        logger.fine(string1.toString());
      }

      for (int i = 0; i < 3; i++) {
        rr_a1c1[i] = r_soln_c[i_soln][0][i] - r_soln_a[i_soln][0][i];
        rr_c1n2[i] = r_soln_n[i_soln][1][i] - r_soln_c[i_soln][0][i];
        rr_n2a2[i] = r_soln_a[i_soln][1][i] - r_soln_n[i_soln][1][i];
        rr_a2c2[i] = r_soln_c[i_soln][1][i] - r_soln_a[i_soln][1][i];
        rr_c2n3[i] = r_soln_n[i_soln][2][i] - r_soln_c[i_soln][1][i];
        rr_n3a3[i] = r_soln_a[i_soln][2][i] - r_soln_n[i_soln][2][i];
        rr_a1a2[i] = r_soln_a[i_soln][1][i] - r_soln_a[i_soln][0][i];
        rr_a2a3[i] = r_soln_a[i_soln][2][i] - r_soln_a[i_soln][1][i];
      }

      double a1c1 = sqrt(dot(rr_a1c1, rr_a1c1));
      double c1n2 = sqrt(dot(rr_c1n2, rr_c1n2));
      double n2a2 = sqrt(dot(rr_n2a2, rr_n2a2));
      double a2c2 = sqrt(dot(rr_a2c2, rr_a2c2));
      double c2n3 = sqrt(dot(rr_c2n3, rr_c2n3));
      double n3a3 = sqrt(dot(rr_n3a3, rr_n3a3));
      double a1a2 = sqrt(dot(rr_a1a2, rr_a1a2));
      double a2a3 = sqrt(dot(rr_a2a3, rr_a2a3));

      if (logger.isLoggable(Level.FINE)) {
        StringBuilder string2 = new StringBuilder();
        string2.append(
            String.format("na: n2a2, n3a3 = %9.3f%9.3f%9.3f%9.3f\n", len0[2], n2a2, len0[5], n3a3));
        string2.append(
            String.format("ac: a1c1, a2c2 = %9.3f%9.3f%9.3f%9.3f\n", len0[0], a1c1, len0[3], a2c2));
        string2.append(
            String.format("cn: c1n2, c2n3 = %9.3f%9.3f%9.3f%9.3f\n", len0[1], c1n2, len0[4], c2n3));
        string2.append(
            String.format(
                "aa: a1a2, a2a3 = %9.3f%9.3f%9.3f%9.3f\n", len_aa[1], a1a2, len_aa[2], a2a3));
        logger.info(string2.toString());
      }

      for (int i = 0; i < 3; i++) {
        tmp_array[i] = rr_a1a2[i] / a1a2;
      }
      calcBondAngle(b_a1a3, tmp_array, a3a1a2);

      for (int i = 0; i < 3; i++) {
        tmp_array[i] = rr_a2a3[i] / a2a3;
      }
      calcBondAngle(tmp_array, b_a1a3, a2a3a1);

      for (int i = 0; i < 3; i++) {
        tmp_array[i] = rr_a1c1[i] / a1c1;
      }
      calcBondAngle(b_a1n1, tmp_array, n1a1c1);

      for (int i = 0; i < 3; i++) {
        tmp_array[i] = -rr_n3a3[i] / n3a3;
      }
      calcBondAngle(b_a3c3, tmp_array, n3a3c3);

      for (int i = 0; i < 3; i++) {
        tmp_array1[i] = rr_a2c2[i] / a2c2;
        tmp_array2[i] = -rr_n2a2[i] / n2a2;
      }
      calcBondAngle(tmp_array1, tmp_array2, n2a2c2);

      if (logger.isLoggable(Level.FINE)) {
        StringBuilder string3 = new StringBuilder();
        string3.append(
            String.format("ang_nac = %9.3f%9.3f\n", toDegrees(b_ang0[0]), toDegrees(n1a1c1[0])));
        string3.append(
            String.format("ang_nac = %9.3f%9.3f\n", toDegrees(b_ang0[3]), toDegrees(n2a2c2[0])));
        string3.append(
            String.format("ang_nac = %9.3f%9.3f\n", toDegrees(b_ang0[6]), toDegrees(n3a3c3[0])));
        logger.fine(string3.toString());
      }

      for (int i = 0; i < 3; i++) {
        tmp_array1[i] = rr_a1c1[i] / a1c1;
        tmp_array2[i] = rr_c1n2[i] / c1n2;
        tmp_array3[i] = rr_n2a2[i] / n2a2;
      }
      calcDihedralAngle(tmp_array1, tmp_array2, tmp_array3, a1c1n2a2);

      for (int i = 0; i < 3; i++) {
        tmp_array1[i] = rr_a2c2[i] / a2c2;
        tmp_array2[i] = rr_c2n3[i] / c2n3;
        tmp_array3[i] = rr_n3a3[i] / n3a3;
      }
      calcDihedralAngle(tmp_array1, tmp_array2, tmp_array3, a2c2n3a3);

      if (logger.isLoggable(Level.FINE)) {
        StringBuilder string4 = new StringBuilder();
        string4.append(
            String.format("t_ang1 = %9.3f%9.3f\n", toDegrees(t_ang0[0]), toDegrees(a1c1n2a2[0])));
        string4.append(
            String.format("t_ang2 = %9.3f%9.3f\n", toDegrees(t_ang0[1]), toDegrees(a2c2n3a3[0])));
        logger.fine(string4.toString());
      }
    }
  }

  /**
   * Get the input angles.
   *
   * @param n_soln an array of {@link int} objects.
   * @param r_n1 an array of {@link double} objects.
   * @param r_a1 an array of {@link double} objects.
   * @param r_a3 an array of {@link double} objects.
   * @param r_c3 an array of {@link double} objects.
   */
  public void getInputAngles(
      int[] n_soln, double[] r_n1, double[] r_a1, double[] r_a3, double[] r_c3) {
    double[] tmp_val = new double[3];
    char[] cone_type = new char[2];

    n_soln[0] = max_soln;

    for (int i = 0; i < 3; i++) {
      r_a1a3[i] = r_a3[i] - r_a1[i];
    }

    double dr_sqr = dot(r_a1a3, r_a1a3);
    len_aa[0] = sqrt(dr_sqr);

    if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr)) {
      n_soln[0] = 0;
      return;
    }

    for (int i = 0; i < 3; i++) {
      r_a1n1[i] = r_n1[i] - r_a1[i];
    }

    len_na[0] = sqrt(dot(r_a1n1, r_a1n1));
    len_na[1] = len0[2];
    len_na[2] = len0[5];

    for (int i = 0; i < 3; i++) {
      r_a3c3[i] = r_c3[i] - r_a3[i];
    }

    len_ac[0] = len0[0];
    len_ac[1] = len0[3];
    len_ac[2] = sqrt(dot(r_a3c3, r_a3c3));

    for (int i = 0; i < 3; i++) {
      b_a1n1[i] = r_a1n1[i] / len_na[0];
      b_a3c3[i] = r_a3c3[i] / len_ac[2];
      b_a1a3[i] = r_a1a3[i] / len_aa[0];
    }

    for (int i = 0; i < 3; i++) {
      tmp_val[i] = -b_a1n1[i];
    }

    calcDihedralAngle(tmp_val, b_a1a3, b_a3c3, delta[3]);

    delta[0][0] = delta[3][0];

    for (int i = 0; i < 3; i++) {
      tmp_val[i] = -b_a1a3[i];
    }

    calcBondAngle(tmp_val, b_a1n1, xi[0]);

    calcBondAngle(b_a1a3, b_a3c3, eta[2]);

    for (int i = 0; i < 3; i++) {
      cos_delta[i + 1] = cos(delta[i + 1][0]);
      sin_delta[i + 1] = sin(delta[i + 1][0]);
      cos_xi[i] = cos(xi[i][0]);
      sin_xi[i] = sin(xi[i][0]);
      sin_xi[i] = sin(xi[i][0]);
      cos_eta[i] = cos(eta[i][0]);
      cos_eta[i] = cos(eta[i][0]);
      sin_eta[i] = sin(eta[i][0]);
      sin_eta[i] = sin(eta[i][0]);
    }

    cos_delta[0] = cos_delta[3];
    sin_delta[0] = sin_delta[3];
    theta[0] = b_ang0[0];
    theta[1] = b_ang0[3];
    theta[2] = b_ang0[6];

    for (int i = 0; i < 3; i++) {
      cos_theta[i] = cos(theta[i]);
    }

    cos_alpha[0] =
        -((len_aa[0] * len_aa[0]) + (len_aa[1] * len_aa[1]) - (len_aa[2] * len_aa[2]))
            / (2.0 * len_aa[0] * len_aa[1]);
    alpha[0] = acos(cos_alpha[0]);
    sin_alpha[0] = sin(alpha[0]);
    cos_alpha[1] =
        ((len_aa[1] * len_aa[1]) + (len_aa[2] * len_aa[2]) - (len_aa[0] * len_aa[0]))
            / (2.0 * len_aa[1] * len_aa[2]);
    alpha[1] = acos(cos_alpha[1]);
    sin_alpha[1] = sin(alpha[1]);
    alpha[2] = PI - alpha[0] + alpha[1];
    cos_alpha[2] = cos(alpha[2]);
    sin_alpha[2] = sin(alpha[2]);

    if (logger.isLoggable(Level.FINE)) {
      StringBuilder sb = new StringBuilder();
      sb.append(
          String.format(
              " xi = %9.4f%9.4f%9.4f\n",
              toDegrees(xi[0][0]), toDegrees(xi[1][0]), toDegrees(xi[2][0])));
      sb.append(
          String.format(
              " eta = %9.4f%9.4f%9.4f\n",
              toDegrees(eta[0][0]), toDegrees(eta[1][0]), toDegrees(eta[2][0])));
      sb.append(
          String.format(
              " delta = %9.4f%9.4f%9.4f\n",
              toDegrees(delta[1][0]), toDegrees(delta[2][0]), toDegrees(delta[3][0])));
      sb.append(
          String.format(
              " theta = %9.4f%9.4f%9.4f\n",
              toDegrees(theta[0]), toDegrees(theta[1]), toDegrees(theta[2])));
      sb.append(
          String.format(
              " alpha = %9.4f%9.4f%9.4f\n",
              toDegrees(alpha[0]), toDegrees(alpha[1]), toDegrees(alpha[2])));
      logger.fine(sb.toString());
    }

    for (int i = 0; i < 3; i++) {
      testTwoConeExistenceSoln(theta[i], xi[i][0], eta[i][0], alpha[i], n_soln, cone_type);
      if (n_soln[0] == 0) {
        return;
      }
    }
  }

  /**
   * Get Polynomial Coefficient.
   *
   * @param polyCoeff an array of {@link double} objects.
   */
  public void getPolyCoeff(double[] polyCoeff) {
    double[] B0 = new double[3];
    double[] B1 = new double[3];
    double[] B2 = new double[3];
    double[] B3 = new double[3];
    double[] B4 = new double[3];
    double[] B5 = new double[3];
    double[] B6 = new double[3];
    double[] B7 = new double[3];
    double[] B8 = new double[3];
    double[][] u11 = new double[5][5];
    double[][] u12 = new double[5][5];
    double[][] u13 = new double[5][5];
    double[][] u31 = new double[5][5];
    double[][] u32 = new double[5][5];
    double[][] u33 = new double[5][5];
    double[][] um1 = new double[5][5];
    double[][] um2 = new double[5][5];
    double[][] um3 = new double[5][5];
    double[][] um4 = new double[5][5];
    double[][] um5 = new double[5][5];
    double[][] um6 = new double[5][5];
    double[][] q_tmp = new double[5][5];
    int[] p1 = new int[2];
    int[] p3 = new int[2];
    int[] p_um1 = new int[2];
    int[] p_um2 = new int[2];
    int[] p_um3 = new int[2];
    int[] p_um4 = new int[2];
    int[] p_um5 = new int[2];
    int[] p_um6 = new int[2];
    int[] p_Q = new int[2];
    int[] p_f1 = new int[1];
    int[] p_f2 = new int[1];
    int[] p_f3 = new int[1];
    int[] p_f4 = new int[1];
    int[] p_f5 = new int[1];
    int[] p_f6 = new int[1];
    int[] p_f7 = new int[1];
    int[] p_f8 = new int[1];
    int[] p_f9 = new int[1];
    int[] p_f10 = new int[1];
    int[] p_f11 = new int[1];
    int[] p_f12 = new int[1];
    int[] p_f13 = new int[1];
    int[] p_f14 = new int[1];
    int[] p_f15 = new int[1];
    int[] p_f16 = new int[1];
    int[] p_f17 = new int[1];
    int[] p_f18 = new int[1];
    int[] p_f19 = new int[1];
    int[] p_f20 = new int[1];
    int[] p_f21 = new int[1];
    int[] p_f22 = new int[1];
    int[] p_f23 = new int[1];
    int[] p_f24 = new int[1];
    int[] p_f25 = new int[1];
    int[] p_f26 = new int[1];
    int[] p_final = new int[1];
    double[] f1 = new double[17];
    double[] f2 = new double[17];
    double[] f3 = new double[17];
    double[] f4 = new double[17];
    double[] f5 = new double[17];
    double[] f6 = new double[17];
    double[] f7 = new double[17];
    double[] f8 = new double[17];
    double[] f9 = new double[17];
    double[] f10 = new double[17];
    double[] f11 = new double[17];
    double[] f12 = new double[17];
    double[] f13 = new double[17];
    double[] f14 = new double[17];
    double[] f15 = new double[17];
    double[] f16 = new double[17];
    double[] f17 = new double[17];
    double[] f18 = new double[17];
    double[] f19 = new double[17];
    double[] f20 = new double[17];
    double[] f21 = new double[17];
    double[] f22 = new double[17];
    double[] f23 = new double[17];
    double[] f24 = new double[17];
    double[] f25 = new double[17];
    double[] f26 = new double[17];

    for (int i = 0; i < 3; i++) {
      double A0 = cos_alpha[i] * cos_xi[i] * cos_eta[i] - cos_theta[i];
      double A1 = -sin_alpha[i] * cos_xi[i] * sin_eta[i];
      double A2 = sin_alpha[i] * sin_xi[i] * cos_eta[i];
      double A3 = sin_xi[i] * sin_eta[i];
      double A4 = A3 * cos_alpha[i];
      int j = i;
      double A21 = A2 * cos_delta[j];
      double A22 = A2 * sin_delta[j];
      double A31 = A3 * cos_delta[j];
      double A32 = A3 * sin_delta[j];
      double A41 = A4 * cos_delta[j];
      double A42 = A4 * sin_delta[j];
      B0[i] = A0 + A22 + A31;
      B1[i] = 2.0 * (A1 + A42);
      B2[i] = 2.0 * (A32 - A21);
      B3[i] = -4.0 * A41;
      B4[i] = A0 + A22 - A31;
      B5[i] = A0 - A22 - A31;
      B6[i] = -2.0 * (A21 + A32);
      B7[i] = 2.0 * (A1 - A42);
      B8[i] = A0 - A22 + A31;
    }

    C0[0][0] = B0[0];
    C0[0][1] = B2[0];
    C0[0][2] = B5[0];
    C1[0][0] = B1[0];
    C1[0][1] = B3[0];
    C1[0][2] = B7[0];
    C2[0][0] = B4[0];
    C2[0][1] = B6[0];
    C2[0][2] = B8[0];

    for (int i = 1; i < 3; i++) {
      C0[i][0] = B0[i];
      C0[i][1] = B1[i];
      C0[i][2] = B4[i];
      C1[i][0] = B2[i];
      C1[i][1] = B3[i];
      C1[i][2] = B6[i];
      C2[i][0] = B5[i];
      C2[i][1] = B7[i];
      C2[i][2] = B8[i];
    }

    for (int i = 0; i < 3; i++) {
      u11[0][i] = C0[0][i];
      u12[0][i] = C1[0][i];
      u13[0][i] = C2[0][i];
      u31[i][0] = C0[1][i];
      u32[i][0] = C1[1][i];
      u33[i][0] = C2[1][i];
    }

    p1[0] = 2;
    p1[1] = 0;
    p3[0] = 0;
    p3[1] = 2;

    polyMulSub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
    polyMulSub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
    polyMulSub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
    polyMulSub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
    polyMulSub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
    polyMulSub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
    polyMulSub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);

    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        Q[i][j] = q_tmp[i][j];
      }
    }

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 17; j++) {
        R[i][j] = 0.0;
      }
    }

    for (int i = 0; i < 3; i++) {
      R[0][i] = C0[2][i];
      R[1][i] = C1[2][i];
      R[2][i] = C2[2][i];
    }

    int p2 = 2;
    int p4 = 4;

    polyMulSub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, p_f1);
    polyMul1(R[1], R[2], p2, p2, f2, p_f2);
    polyMulSub1(R[1], f1, R[0], f2, p2, p_f1[0], p2, p_f2[0], f3, p_f3);
    polyMul1(R[2], f1, p2, p_f1[0], f4, p_f4);
    polyMulSub1(R[1], f3, R[0], f4, p2, p_f3[0], p2, p_f4[0], f5, p_f5);
    polyMulSub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, p_f6);
    polyMulSub1(Q[2], f1, R[2], f6, p4, p_f1[0], p2, p_f6[0], f7, p_f7);
    polyMulSub1(Q[3], f3, R[2], f7, p4, p_f3[0], p2, p_f7[0], f8, p_f8);
    polyMulSub1(Q[4], f5, R[2], f8, p4, p_f5[0], p2, p_f8[0], f9, p_f9);
    polyMulSub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, p_f10);
    polyMulSub1(Q[2], f1, R[0], f10, p4, p_f1[0], p2, p_f10[0], f11, p_f11);
    polyMulSub1(Q[1], f3, R[0], f11, p4, p_f3[0], p2, p_f11[0], f12, p_f12);
    polyMulSub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, p_f13);
    polyMulSub1(Q[3], f1, R[2], f13, p4, p_f1[0], p2, p_f13[0], f14, p_f14);
    polyMulSub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, p_f15);
    polyMulSub1(Q[4], f1, R[2], f15, p4, p_f1[0], p2, p_f15[0], f16, p_f16);
    polyMulSub1(Q[1], f14, Q[0], f16, p4, p_f14[0], p4, p_f16[0], f17, p_f17);
    polyMulSub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, p_f18);
    polyMulSub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, p_f19);
    polyMulSub1(Q[3], f19, Q[2], f18, p4, p_f19[0], p4, p_f18[0], f20, p_f20);
    polyMulSub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, p_f21);
    polyMul1(Q[4], f21, p4, p_f21[0], f22, p_f22);
    polySub1(f20, f22, p_f20, p_f22, f23, p_f23);
    polyMul1(R[0], f23, p2, p_f23[0], f24, p_f24);
    polySub1(f17, f24, p_f17, p_f24, f25, p_f25);
    polyMulSub1(Q[4], f12, R[2], f25, p4, p_f12[0], p2, p_f25[0], f26, p_f26);
    polyMulSub1(Q[0], f9, R[0], f26, p4, p_f9[0], p2, p_f26[0], polyCoeff, p_final);

    if (p_final[0] != 16) {
      logger.info(String.format("Error. Degree of polynomial is not 16!"));
      return;
    }

    if (polyCoeff[16] < 0.0e0) {
      for (int i = 0; i < 17; i++) {
        polyCoeff[i] *= -1.0;
      }
    }

    if (logger.isLoggable(Level.FINE)) {
      StringBuilder string = new StringBuilder();
      string.append(String.format(" Polynomial Coefficients\n"));
      for (int i = 0; i < 17; i++) {
        string.append(String.format(" %5d %15.6f\n", i, polyCoeff[i]));
      }
      string.append(String.format("\n"));
      logger.fine(string.toString());
    }
  }

  /**
   * A Basic 2D matrix multiplication. This function multiplies matrices ma and mb, storing the
   * answer in mc.
   *
   * @param ma an array of {@link double} objects.
   * @param mb an array of {@link double} objects.
   * @param mc an array of {@link double} objects.
   */
  public void matMul(double[][] ma, double[] mb, double[] mc) {
    for (int i = 0; i < 3; i++) {
      mc[i] = 0.0;
      for (int j = 0; j < 3; j++) {
        mc[i] += ma[i][j] * mb[j];
      }
    }
  }

  /**
   * Polynomial multiply 1.
   *
   * @param u1 an array of {@link double} objects.
   * @param u2 an array of {@link double} objects.
   * @param p1 a int.
   * @param p2 a int.
   * @param u3 an array of {@link double} objects.
   * @param p3 an array of {@link int} objects.
   */
  public void polyMul1(double[] u1, double[] u2, int p1, int p2, double[] u3, int[] p3) {
    p3[0] = p1 + p2;
    for (int i = 0; i < 17; i++) {
      u3[i] = 0.0;
    }
    for (int i1 = 0; i1 <= p1; i1++) {
      double u1i = u1[i1];
      for (int i2 = 0; i2 <= p2; i2++) {
        int i3 = i1 + i2;
        u3[i3] = u3[i3] + u1i * u2[i2];
      }
    }
  }

  /**
   * Polynomial multiplication 2.
   *
   * @param u1 an array of {@link double} objects.
   * @param u2 an array of {@link double} objects.
   * @param p1 an array of {@link int} objects.
   * @param p2 an array of {@link int} objects.
   * @param u3 an array of {@link double} objects.
   * @param p3 an array of {@link int} objects.
   */
  public void polyMul2(double[][] u1, double[][] u2, int[] p1, int[] p2, double[][] u3, int[] p3) {
    for (int i = 0; i < 2; i++) {
      p3[i] = p1[i] + p2[i];
    }
    for (int i = 0; i < 5; i++) {
      for (int j = 0; j < 5; j++) {
        u3[i][j] = 0.0;
      }
    }
    int p11 = p1[0];
    int p12 = p1[1];
    int p21 = p2[0];
    int p22 = p2[1];
    for (int i1 = 0; i1 <= p12; i1++) {
      for (int j1 = 0; j1 <= p11; j1++) {
        double u1ij = u1[i1][j1];
        for (int i2 = 0; i2 <= p22; i2++) {
          int i3 = i1 + i2;
          for (int j2 = 0; j2 <= p21; j2++) {
            int j3 = j1 + j2;
            u3[i3][j3] = u3[i3][j3] + u1ij * u2[i2][j2];
          }
        }
      }
    }
  }

  /**
   * Polynomial multiply subtraction 1.
   *
   * @param u1 an array of {@link double} objects.
   * @param u2 an array of {@link double} objects.
   * @param u3 an array of {@link double} objects.
   * @param u4 an array of {@link double} objects.
   * @param p1 a int.
   * @param p2 a int.
   * @param p3 a int.
   * @param p4 a int.
   * @param u5 an array of {@link double} objects.
   * @param p5 an array of {@link int} objects.
   */
  public void polyMulSub1(
      double[] u1,
      double[] u2,
      double[] u3,
      double[] u4,
      int p1,
      int p2,
      int p3,
      int p4,
      double[] u5,
      int[] p5) {

    double[] d1 = new double[17];
    double[] d2 = new double[17];
    int[] pd1 = new int[1];
    int[] pd2 = new int[1];

    polyMul1(u1, u2, p1, p2, d1, pd1);
    polyMul1(u3, u4, p3, p4, d2, pd2);
    polySub1(d1, d2, pd1, pd2, u5, p5);
  }

  /**
   * Polynomial multiplication subtraction 2.
   *
   * @param u1 an array of {@link double} objects.
   * @param u2 an array of {@link double} objects.
   * @param u3 an array of {@link double} objects.
   * @param u4 an array of {@link double} objects.
   * @param p1 an array of {@link int} objects.
   * @param p2 an array of {@link int} objects.
   * @param p3 an array of {@link int} objects.
   * @param p4 an array of {@link int} objects.
   * @param u5 an array of {@link double} objects.
   * @param p5 an array of {@link int} objects.
   */
  public void polyMulSub2(
      double[][] u1,
      double[][] u2,
      double[][] u3,
      double[][] u4,
      int[] p1,
      int[] p2,
      int[] p3,
      int[] p4,
      double[][] u5,
      int[] p5) {
    double[][] d1 = new double[5][5];
    double[][] d2 = new double[5][5];
    int[] pd1 = new int[2];
    int[] pd2 = new int[2];

    polyMul2(u1, u2, p1, p2, d1, pd1);
    polyMul2(u3, u4, p3, p4, d2, pd2);
    polySub2(d1, d2, pd1, pd2, u5, p5);
  }

  /**
   * Polynomial subtraction 1.
   *
   * @param u1 an array of {@link double} objects.
   * @param u2 an array of {@link double} objects.
   * @param p1 an array of {@link int} objects.
   * @param p2 an array of {@link int} objects.
   * @param u3 an array of {@link double} objects.
   * @param p3 an array of {@link int} objects.
   */
  public void polySub1(double[] u1, double[] u2, int[] p1, int[] p2, double[] u3, int[] p3) {
    p3[0] = max(p1[0], p2[0]);
    for (int i = 0; i <= p3[0]; i++) {
      if (i > p2[0]) {
        u3[i] = u1[i];
      } else if (i > p1[0]) {
        u3[i] = -u2[i];
      } else {
        u3[i] = u1[i] - u2[i];
      }
    }
  }

  /**
   * Polynomial subtraction 2.
   *
   * @param u1 an array of {@link double} objects.
   * @param u2 an array of {@link double} objects.
   * @param p1 an array of {@link int} objects.
   * @param p2 an array of {@link int} objects.
   * @param u3 an array of {@link double} objects.
   * @param p3 an array of {@link int} objects.
   */
  public void polySub2(double[][] u1, double[][] u2, int[] p1, int[] p2, double[][] u3, int[] p3) {
    int p11 = p1[0];
    int p12 = p1[1];
    int p21 = p2[0];
    int p22 = p2[1];
    p3[0] = max(p11, p21);
    p3[1] = max(p12, p22);
    for (int i = 0; i <= p3[1]; i++) {
      boolean i1_ok = (i > p12);
      boolean i2_ok = (i > p22);
      for (int j = 0; j <= p3[0]; j++) {
        if (i2_ok || (j > p21)) {
          u3[i][j] = u1[i][j];
        } else if (i1_ok || (j > p11)) {
          u3[i][j] = -u2[i][j];
        } else {
          u3[i][j] = u1[i][j] - u2[i][j];
        }
      }
    }
  }

  /**
   * Quotient of two vectors in three dimensional space.
   *
   * @param axis an array of {@link double} objects.
   * @param quarter_ang a double.
   * @param p an array of {@link double} objects.
   */
  public void quaternion(double[] axis, double quarter_ang, double[] p) {
    double tan_w = tan(quarter_ang);
    double tan_sqr = tan_w * tan_w;
    double tan1 = 1.0 + tan_sqr;
    double cosine = (1.0 - tan_sqr) / tan1;
    double sine = 2.0 * tan_w / tan1;
    p[0] = cosine;
    p[1] = axis[0] * sine;
    p[2] = axis[1] * sine;
    p[3] = axis[2] * sine;
  }

  /**
   * Means of representing a rotation of an axis about an origin with no angular specification.
   *
   * @param q input Axis angle.
   * @param U output rotation matrix.
   */
  public void rotationMatrix(double[] q, double[][] U) {

    double q0 = q[0];
    double q1 = q[1];
    double q2 = q[2];
    double q3 = q[3];
    double b0 = 2.0 * q0;
    double b1 = 2.0 * q1;
    double q00 = b0 * q0 - 1.0;
    double q02 = b0 * q2;
    double q03 = b0 * q3;
    double q11 = b1 * q1;
    double q12 = b1 * q2;
    double q13 = b1 * q3;
    double b2 = 2.0 * q2;
    double b3 = 2.0 * q3;
    double q01 = b0 * q1;
    double q22 = b2 * q2;
    double q23 = b2 * q3;
    double q33 = b3 * q3;

    U[0][0] = q00 + q11;
    U[0][1] = q12 - q03;
    U[0][2] = q13 + q02;
    U[1][0] = q12 + q03;
    U[1][1] = q00 + q22;
    U[1][2] = q23 - q01;
    U[2][0] = q13 - q02;
    U[2][1] = q23 + q01;
    U[2][2] = q00 + q33;
  }

  /**
   * Returns the sign (positive or negative) of a variable.
   *
   * @param a a double.
   * @param b a double.
   * @return If b is positive or zero, return abs(a), else return -abs(a).
   */
  public double sign(double a, double b) {
    if (b >= 0.0) {
      return FastMath.abs(a);
    } else {
      return -FastMath.abs(a);
    }
  }

  /**
   * solve3PepPoly.
   *
   * @param r_n1 an array of {@link double} objects.
   * @param r_a1 an array of {@link double} objects.
   * @param r_a3 an array of {@link double} objects.
   * @param r_c3 an array of {@link double} objects.
   * @param r_soln_n an array of {@link double} objects.
   * @param r_soln_a an array of {@link double} objects.
   * @param r_soln_c an array of {@link double} objects.
   * @param n_soln an array of {@link int} objects.
   */
  public void solve3PepPoly(
      double[] r_n1,
      double[] r_a1,
      double[] r_a3,
      double[] r_c3,
      double[][][] r_soln_n,
      double[][][] r_soln_a,
      double[][][] r_soln_c,
      int[] n_soln) {
    double[] polyCoeff = new double[deg_pol[0] + 1];
    double[] roots = new double[max_soln];

    getInputAngles(n_soln, r_n1, r_a1, r_a3, r_c3);

    if (n_soln[0] == 0) {
      return;
    }

    getPolyCoeff(polyCoeff);
    sturmMethod.solveSturm(deg_pol, n_soln, polyCoeff, roots);
    if (n_soln[0] == 0) {
      logger.info("Could not find alternative loop solutions using KIC.");
    }

    getCoordsFromPolyRoots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);
  }

  /**
   * testTwoConeExistenceSoln.
   *
   * @param tt a double.
   * @param kx a double.
   * @param et a double.
   * @param ap a double.
   * @param n_soln an array of {@link int} objects.
   * @param cone_type an array of {@link char} objects.
   */
  public void testTwoConeExistenceSoln(
      double tt, double kx, double et, double ap, int[] n_soln, char[] cone_type) {

    n_soln[0] = max_soln;
    double ap1 = ap;
    double kx1 = kx;
    double et1 = et;
    double at = ap1 - tt;
    double ex = kx1 + et1;
    double abs_at = abs(at);
    if (abs_at > ex) {
      n_soln[0] = 0;
      return;
    }

    boolean complicated = false;
    if (complicated) {
      double cos_tx1 = cos(tt + kx1);
      double cos_tx2 = cos(tt - kx1);
      double cos_te1 = cos(tt + et1);
      double cos_te2 = cos(tt - et1);
      double cos_ea1 = cos(et1 + ap1);
      double cos_ea2 = cos(et1 - ap1);
      double cos_xa1 = cos(kx1 + ap1);
      double cos_xa2 = cos(kx1 - ap1);

      double s1 = 0;
      double s2 = 0;
      double t1 = 0;
      double t2 = 0;

      if ((cos_te1 - cos_xa2) * (cos_te1 - cos_xa1) <= 0.0e0) {
        s1 = 0;
      }

      if ((cos_te2 - cos_xa2) * (cos_te2 - cos_xa1) <= 0.0e0) {
        s2 = 0;
      }

      if ((cos_tx1 - cos_ea2) * (cos_tx1 - cos_ea1) <= 0.0e0) {
        t1 = 0;
      }

      if ((cos_tx2 - cos_ea2) * (cos_tx2 - cos_ea1) <= 0.0e0) {
        t2 = 0;
      }
    }
  }

  /** Initialize Loop Closure. */
  private void initializeLoopClosure() {
    double[] axis = new double[3];
    double[] rr_a1 = new double[3];
    double[] rr_c1 = new double[3];
    double[] rr_n2 = new double[3];
    double[] rr_a2 = new double[3];
    double[] rr_n2a2_ref = new double[3];
    double[] rr_c1a1 = new double[3];
    double[] rr_a1a2 = new double[3];
    double[] dr = new double[3];
    double[] bb_c1a1 = new double[3];
    double[] bb_a1a2 = new double[3];
    double[] bb_a2n2 = new double[3];
    double[] p = new double[4];
    double[][] Us = new double[3][3];
    double[] mulpro = new double[3];
    double[] tmp_val = new double[3];

    // initializing initial length, angle, and torsion arrays
    len0[0] = 1.52;
    len0[1] = 1.33;
    len0[2] = 1.45;
    len0[3] = 1.52;
    len0[4] = 1.33;
    len0[5] = 1.45;

    b_ang0[0] = Math.toRadians(111.6);
    b_ang0[1] = Math.toRadians(117.5);
    b_ang0[2] = Math.toRadians(120.0);
    b_ang0[3] = Math.toRadians(111.6);
    b_ang0[4] = Math.toRadians(117.5);
    b_ang0[5] = Math.toRadians(120.0);
    b_ang0[6] = Math.toRadians(111.6);

    t_ang0[0] = Math.PI;
    t_ang0[1] = Math.PI;

    for (int i = 0; i < 3; i++) {
      rr_c1[i] = 0.0;
    }
    // initializing axis
    axis[0] = 1.0;
    axis[1] = 0.0;
    axis[2] = 0.0;
    for (int i = 0; i < 2; i++) {
      // iniitalizing rr_a1 array
      rr_a1[0] = cos(b_ang0[3 * i + 1]) * len0[3 * i];
      rr_a1[1] = sin(b_ang0[3 * i + 1]) * len0[3 * i];
      rr_a1[2] = 0.0e0;
      // initializing rr_n2 array
      rr_n2[0] = len0[3 * i + 1];
      rr_n2[1] = 0.0e0;
      rr_n2[2] = 0.0e0;
      // initializing rr_c1a1 array
      for (int j = 0; j < 3; j++) {
        rr_c1a1[j] = rr_a1[j] - rr_c1[j];
      }
      // initializing rr_n2a2_ref array
      rr_n2a2_ref[0] = -cos(b_ang0[3 * i + 2]) * len0[3 * i + 2];
      rr_n2a2_ref[1] = sin(b_ang0[3 * i + 2]) * len0[3 * i + 2];
      rr_n2a2_ref[2] = 0.0e0;
      // quaternion is the quotient of two vectors in 3D space
      quaternion(axis, t_ang0[i] * 0.25e0, p);
      // means of representing a rotation of an axis about an origin
      //      with no angular specification
      rotationMatrix(p, Us);
      // basic matrix multiplication
      matMul(Us, rr_n2a2_ref, mulpro);
      for (int j = 0; j < 3; j++) {
        rr_a2[j] = mulpro[j] + rr_n2[j];
        rr_a1a2[j] = rr_a2[j] - rr_a1[j];
        dr[j] = rr_a1a2[j];
      }
      double len2 = dot(dr, dr);
      double len1 = sqrt(len2);
      len_aa[i + 1] = len1;
      for (int j = 0; j < 3; j++) {
        bb_c1a1[j] = rr_c1a1[j] / len0[3 * i];
        bb_a1a2[j] = rr_a1a2[j] / len1;
        bb_a2n2[j] = (rr_n2[j] - rr_a2[j]) / len0[3 * i + 2];
      }
      for (int j = 0; j < 3; j++) {
        tmp_val[j] = -bb_a1a2[j];
      }
      calcBondAngle(tmp_val, bb_a2n2, xi[i + 1]);
      for (int j = 0; j < 3; j++) {
        tmp_val[j] = -bb_c1a1[j];
      }
      calcBondAngle(bb_a1a2, tmp_val, eta[i]);
      calcDihedralAngle(bb_c1a1, bb_a1a2, bb_a2n2, delta[i + 1]);
      delta[i + 1][0] = PI - delta[i + 1][0];
    }

    double a_min = b_ang0[3] - (xi[1][0] + eta[1][0]);
    double a_max = min((b_ang0[3] + (xi[1][0] + eta[1][0])), PI);
    aa13_min_sqr = pow(len_aa[1], 2) + pow(len_aa[2], 2) - 2.0 * len_aa[1] * len_aa[2] * cos(a_min);
    aa13_max_sqr = pow(len_aa[1], 2) + pow(len_aa[2], 2) - 2.0 * len_aa[1] * len_aa[2] * cos(a_max);
  }
}
