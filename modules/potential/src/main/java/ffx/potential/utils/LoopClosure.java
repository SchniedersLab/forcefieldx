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
package ffx.potential.utils;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.SturmMethod;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.DoubleMath.X;
import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static ffx.numerics.math.DoubleMath.dot;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
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

/**
 * LoopClosure class. The program can be used to find loop structures involving six backbone torsion angles given the
 * position of the two atoms before the loop and two atoms after the loop. Possible structures for a three residue gap
 * in a protein can be found given the coordinates of the N and CA atoms of the first residue and the CA ans C atoms of
 * the third residue. Multiple conformations are suggested by this algorithm when multiple solutions exist.
 *
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class LoopClosure {

    private static final Logger logger = Logger.getLogger(LoopClosure.class.getName());

    private static final int MAX_SOLUTIONS = 16;

    private static final double[] BOND_LENS;
    private static final double[] BOND_ANGLES;
    private static final double AA_LEN;
    private static final double AA13_MIN_SQR;
    private static final double AA13_MAX_SQR;
    private static final double[] XI = new double[3];
    private static final double[] DELTA = new double[4];
    private static final double[] ETA = new double[3];

    private final double[] alpha = new double[3];
    private final double[] theta = new double[3];
    private final double[] cos_alpha = new double[3];
    private final double[] sin_alpha = new double[3];
    private final double[] cos_theta = new double[3];
    private final double[] cos_delta = new double[4];
    private final double[] sin_delta = new double[4];
    private final double[] cos_xi = new double[3];
    private final double[] cos_eta = new double[3];
    private final double[] sin_xi = new double[3];
    private final double[] sin_eta = new double[3];
    private final double[] r_a1a3 = new double[3];
    private final double[] r_a1n1 = new double[3];
    private final double[] r_a3c3 = new double[3];
    private final double[] b_a1a3 = new double[3];
    private final double[] b_a1n1 = new double[3];
    private final double[] b_a3c3 = new double[3];
    private final double[] len_na = new double[3];
    private final double[] len_ac = new double[3];
    private final double[][] C0 = new double[3][3];
    private final double[][] C1 = new double[3][3];
    private final double[][] C2 = new double[3][3];
    private final double[][] Q = new double[5][17];
    private final double[][] R = new double[3][17];

    static {
        // Initialize bond length and angle arrays.
        BOND_LENS = new double[3];
        BOND_LENS[0] = 1.52;  // Ca - C
        BOND_LENS[1] = 1.33;  // C - N
        BOND_LENS[2] = 1.45;  // N - Ca
        BOND_ANGLES = new double[3];
        BOND_ANGLES[0] = Math.toRadians(111.6);
        BOND_ANGLES[1] = Math.toRadians(117.5);
        BOND_ANGLES[2] = Math.toRadians(120.0);

        // Calculate the pi/4 rotation matrix about the x axis.
        var rotMatrix = rotationMatrix(new double[]{1.0, 0.0, 0.0}, Math.PI * 0.25e0);

        // Relative cartesian coordinates for the carboxyl carbon, nitrogen, and other alpha carbon. The base alpha
        // carbon is assumed to be at the origin.
        var relCarboxyl = new double[3];
        relCarboxyl[0] = cos(BOND_ANGLES[1]) * BOND_LENS[0];
        relCarboxyl[1] = sin(BOND_ANGLES[1]) * BOND_LENS[0];
        relCarboxyl[2] = 0.0e0;
        var relNitrogen = new double[3];
        relNitrogen[0] = BOND_LENS[1];
        relNitrogen[1] = 0.0e0;
        relNitrogen[2] = 0.0e0;
        var relAlpha2 = new double[3];
        relAlpha2[0] = -cos(BOND_ANGLES[2]) * BOND_LENS[2];
        relAlpha2[1] = sin(BOND_ANGLES[2]) * BOND_LENS[2];
        relAlpha2[2] = 0.0e0;

        // Calculate the relative position of the second alpha carbon.
        relAlpha2 = matrixMultiplication(rotMatrix, relAlpha2);
        var alpha1Alpha2 = new double[3];
        var alpha2N = new double[3];
        for (int i = 0; i < 3; i++) {
            alpha1Alpha2[i] = relAlpha2[i] + relNitrogen[i] - relCarboxyl[i];
            alpha2N[i] = -relAlpha2[i];
        }

        // Calculate the space between the two alpha carbons.
        AA_LEN = sqrt(dot(alpha1Alpha2, alpha1Alpha2));
        double[] tmp_val = new double[3];
        for (int i = 0; i < 2; i++) {
            // Temporary variables to calculate relevant angles.
            for (int j = 0; j < 3; j++) {
                tmp_val[j] = -alpha1Alpha2[j];
            }

            XI[i + 1] = DoubleMath.angle(tmp_val, alpha2N);
            for (int j = 0; j < 3; j++) {
                tmp_val[j] = -relCarboxyl[j];
            }

            ETA[i] = DoubleMath.angle(alpha1Alpha2, tmp_val);
            DELTA[i + 1] = PI - dihedralAngle(relCarboxyl, alpha1Alpha2, alpha2N);
        }

        double a_min = BOND_ANGLES[0] - (XI[1] + ETA[1]);
        double a_max = min((BOND_ANGLES[0] + (XI[1] + ETA[1])), PI);
        AA13_MIN_SQR = pow(AA_LEN, 2) + pow(AA_LEN, 2) - 2.0 * AA_LEN * AA_LEN * cos(a_min);
        AA13_MAX_SQR = pow(AA_LEN, 2) + pow(AA_LEN, 2) - 2.0 * AA_LEN * AA_LEN * cos(a_max);
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

        return (U31 * U13 - U11 * U33) / (U12 * U33 - U13 * U32);
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

        return (K3 * B0 - A0 * B2_3) / (K2 * B2 - K3 * B1);
    }

    /**
     * Get coordinates from polynomial roots.
     *
     * @param numSolutions the number of potential solutions to be tested.
     * @param roots    a double array.
     * @param r_n1     a double array.
     * @param r_a1     a double array.
     * @param r_a3     a double array.
     * @param r_c3     a double array.
     * @param r_soln_n a double array.
     * @param r_soln_a a double array.
     * @param r_soln_c a double array.
     */
    private void getCoordsFromPolyRoots(int numSolutions, double[] roots, double[] r_n1, double[] r_a1,
                                        double[] r_a3, double[] r_c3, double[][][] r_soln_n, double[][][] r_soln_a,
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
        double[] rr_a1c1 = new double[3];
        double[] rr_c1n2 = new double[3];
        double[] rr_n2a2 = new double[3];
        double[] rr_a2c2 = new double[3];
        double[] rr_c2n3 = new double[3];
        double[] rr_n3a3 = new double[3];
        double[] rr_a1a2 = new double[3];
        double[] rr_a2a3 = new double[3];
        double[] ex_tmp = new double[3];
        double[] tmp_array = new double[3];
        double[] tmp_array1 = new double[3];
        double[] tmp_array2 = new double[3];
        double[] tmp_array3 = new double[3];
        double[] mat11 = new double[3];
        double[] mat22 = new double[3];
        double[] mat33 = new double[3];
        double[] mat44 = new double[3];
        double[] mat55 = new double[3];

        arraycopy(b_a1a3, 0, ex, 0, 3);
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

        // Determine sign of the angle based on the dot product.
        double sig1Init = DoubleMath.angle(s1[0], r_tmp);
        if (dot(r_tmp, s2[0]) >= 0.0) {
            sig1Init = FastMath.abs(sig1Init);
        } else {
            sig1Init = -FastMath.abs(sig1Init);
        }

        for (int i = 0; i < 3; i++) {
            r_a[0][i] = r_a1[i];
            r_a[1][i] = r_a1[i] + AA_LEN * b_a1a2[i];
            r_a[2][i] = r_a3[i];
            r0[i] = r_a1[i];
        }

        for (int i_soln = 0; i_soln < numSolutions; i_soln++) {
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
            var rotMatrix = rotationMatrix(ex_tmp, tmp_value);

            for (int i = 0; i < 3; i++) {
                mat11[i] = r_c[0][i] - r0[i];
                mat22[i] = r_n[1][i] - r0[i];
                mat33[i] = r_a[1][i] - r0[i];
                mat44[i] = r_c[1][i] - r0[i];
                mat55[i] = r_n[2][i] - r0[i];
            }

            var mat1 = matrixMultiplication(rotMatrix, mat11);
            var mat2 = matrixMultiplication(rotMatrix, mat22);
            var mat3 = matrixMultiplication(rotMatrix, mat33);
            var mat4 = matrixMultiplication(rotMatrix, mat44);
            var mat5 = matrixMultiplication(rotMatrix, mat55);

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
                logger.fine(format("roots: t0\t\t t2\t\t t1\t\t %d\n%15.6f %15.6f %15.6f\n",
                        i_soln, half_tan[2], half_tan[1], half_tan[0]));
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

            if (logger.isLoggable(Level.FINE)) {
                double a1c1 = sqrt(dot(rr_a1c1, rr_a1c1));
                double c1n2 = sqrt(dot(rr_c1n2, rr_c1n2));
                double n2a2 = sqrt(dot(rr_n2a2, rr_n2a2));
                double a2c2 = sqrt(dot(rr_a2c2, rr_a2c2));
                double c2n3 = sqrt(dot(rr_c2n3, rr_c2n3));
                double n3a3 = sqrt(dot(rr_n3a3, rr_n3a3));
                double a1a2 = sqrt(dot(rr_a1a2, rr_a1a2));
                double a2a3 = sqrt(dot(rr_a2a3, rr_a2a3));

                StringBuilder sb = new StringBuilder();
                sb.append(format("na: n2a2, n3a3 = %9.3f%9.3f%9.3f\n", BOND_LENS[2], n2a2, n3a3));
                sb.append(format("ac: a1c1, a2c2 = %9.3f%9.3f%9.3f\n", BOND_LENS[0], a1c1, a2c2));
                sb.append(format("cn: c1n2, c2n3 = %9.3f%9.3f%9.3f\n", BOND_LENS[1], c1n2, c2n3));
                sb.append(
                        format("aa: a1a2, a2a3 = %9.3f%9.3f%9.3f%9.3f\n", AA_LEN, a1a2, AA_LEN, a2a3));
                logger.fine(sb.toString());
                sb.setLength(0);

                for (int i = 0; i < 3; i++) {
                    tmp_array[i] = rr_a1a2[i] / a1a2;
                }

                DoubleMath.angle(b_a1a3, tmp_array);

                for (int i = 0; i < 3; i++) {
                    tmp_array[i] = rr_a2a3[i] / a2a3;
                }

                DoubleMath.angle(tmp_array, b_a1a3);

                for (int i = 0; i < 3; i++) {
                    tmp_array[i] = rr_a1c1[i] / a1c1;
                }

                DoubleMath.angle(b_a1n1, tmp_array);

                for (int i = 0; i < 3; i++) {
                    tmp_array[i] = -rr_n3a3[i] / n3a3;
                }

                DoubleMath.angle(b_a3c3, tmp_array);

                for (int i = 0; i < 3; i++) {
                    tmp_array1[i] = rr_a2c2[i] / a2c2;
                    tmp_array2[i] = -rr_n2a2[i] / n2a2;
                }

                DoubleMath.angle(tmp_array1, tmp_array2);

                for (int i = 0; i < 3; i++) {
                    tmp_array1[i] = rr_a1c1[i] / a1c1;
                    tmp_array2[i] = rr_c1n2[i] / c1n2;
                    tmp_array3[i] = rr_n2a2[i] / n2a2;
                }

                var a1c1n2a2 = dihedralAngle(tmp_array1, tmp_array2, tmp_array3);

                for (int i = 0; i < 3; i++) {
                    tmp_array1[i] = rr_a2c2[i] / a2c2;
                    tmp_array2[i] = rr_c2n3[i] / c2n3;
                    tmp_array3[i] = rr_n3a3[i] / n3a3;
                }

                var a2c2n3a3 = dihedralAngle(tmp_array1, tmp_array2, tmp_array3);
                sb.append(format("t_ang1 = %9.3f\n", toDegrees(a1c1n2a2)));
                sb.append(format("t_ang2 = %9.3f\n", toDegrees(a2c2n3a3)));
                logger.fine(sb.toString());
            }
        }
    }

    /**
     * Get the input angles.
     *
     * @param r_n1   a double array.
     * @param r_a1   a double array.
     * @param r_a3   a double array.
     * @param r_c3   a double array.
     */
    public boolean getInputAngles(double[] r_n1, double[] r_a1, double[] r_a3, double[] r_c3) {
        double[] tmp_val = new double[3];
        for (int i = 0; i < 3; i++) {
            r_a1a3[i] = r_a3[i] - r_a1[i];
        }

        var dr_sqr = dot(r_a1a3, r_a1a3);
        if ((dr_sqr < AA13_MIN_SQR) || (dr_sqr > AA13_MAX_SQR)) {
            return false;
        }

        var a1A3Len = sqrt(dr_sqr);
        for (int i = 0; i < 3; i++) {
            r_a1n1[i] = r_n1[i] - r_a1[i];
        }

        len_na[0] = sqrt(dot(r_a1n1, r_a1n1));
        len_na[1] = BOND_LENS[2];
        len_na[2] = BOND_LENS[2];

        for (int i = 0; i < 3; i++) {
            r_a3c3[i] = r_c3[i] - r_a3[i];
        }

        len_ac[0] = BOND_LENS[0];
        len_ac[1] = BOND_LENS[0];
        len_ac[2] = sqrt(dot(r_a3c3, r_a3c3));

        for (int i = 0; i < 3; i++) {
            b_a1n1[i] = r_a1n1[i] / len_na[0];
            b_a3c3[i] = r_a3c3[i] / len_ac[2];
            b_a1a3[i] = r_a1a3[i] / a1A3Len;
        }

        for (int i = 0; i < 3; i++) {
            tmp_val[i] = -b_a1n1[i];
        }

        DELTA[3] = dihedralAngle(r_n1, r_a1, r_a3, r_c3);
        DELTA[0] = DELTA[3];

        for (int i = 0; i < 3; i++) {
            tmp_val[i] = -b_a1a3[i];
        }

        XI[0] = DoubleMath.angle(tmp_val, b_a1n1);
        ETA[2] = DoubleMath.angle(b_a1a3, b_a3c3);

        for (int i = 0; i < 3; i++) {
            cos_delta[i + 1] = cos(DELTA[i + 1]);
            sin_delta[i + 1] = sin(DELTA[i + 1]);
            cos_xi[i] = cos(XI[i]);
            sin_xi[i] = sin(XI[i]);
            sin_xi[i] = sin(XI[i]);
            cos_eta[i] = cos(ETA[i]);
            sin_eta[i] = sin(ETA[i]);
        }

        cos_delta[0] = cos_delta[3];
        sin_delta[0] = sin_delta[3];
        theta[0] = BOND_ANGLES[0];
        theta[1] = BOND_ANGLES[0];
        theta[2] = BOND_ANGLES[0];

        for (int i = 0; i < 3; i++) {
            cos_theta[i] = cos(theta[i]);
        }

        cos_alpha[0] =
                -((a1A3Len * a1A3Len) + (AA_LEN * AA_LEN) - (AA_LEN * AA_LEN)) / (2.0
                        * a1A3Len * AA_LEN);
        alpha[0] = acos(cos_alpha[0]);
        sin_alpha[0] = sin(alpha[0]);
        cos_alpha[1] =
                ((AA_LEN * AA_LEN) + (AA_LEN * AA_LEN) - (a1A3Len * a1A3Len)) / (2.0
                        * AA_LEN * AA_LEN);
        alpha[1] = acos(cos_alpha[1]);
        sin_alpha[1] = sin(alpha[1]);
        alpha[2] = PI - alpha[0] + alpha[1];
        cos_alpha[2] = cos(alpha[2]);
        sin_alpha[2] = sin(alpha[2]);

        if (logger.isLoggable(Level.FINE)) {
            StringBuilder sb = new StringBuilder();
            sb.append(format(" xi = %9.4f%9.4f%9.4f\n", toDegrees(XI[0]), toDegrees(XI[1]),
                    toDegrees(XI[2])));
            sb.append(format(" eta = %9.4f%9.4f%9.4f\n", toDegrees(ETA[0]), toDegrees(ETA[1]),
                    toDegrees(ETA[2])));
            sb.append(format(" delta = %9.4f%9.4f%9.4f\n", toDegrees(DELTA[1]), toDegrees(DELTA[2]),
                    toDegrees(DELTA[3])));
            sb.append(format(" theta = %9.4f%9.4f%9.4f\n", toDegrees(theta[0]), toDegrees(theta[1]),
                    toDegrees(theta[2])));
            sb.append(format(" alpha = %9.4f%9.4f%9.4f\n", toDegrees(alpha[0]), toDegrees(alpha[1]),
                    toDegrees(alpha[2])));
            logger.fine(sb.toString());
        }

        for (int i = 0; i < 3; i++) {
            if (!testTwoConeExistenceSoln(theta[i], XI[i], ETA[i], alpha[i])) {
                return false;
            }
        }

        // No errors, return true;
        return true;
    }

    /**
     * Get Polynomial Coefficient.
     *
     * @param polyCoeff a double array.
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
            double A21 = A2 * cos_delta[i];
            double A22 = A2 * sin_delta[i];
            double A31 = A3 * cos_delta[i];
            double A32 = A3 * sin_delta[i];
            double A41 = A4 * cos_delta[i];
            double A42 = A4 * sin_delta[i];
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
            arraycopy(q_tmp[i], 0, Q[i], 0, 5);
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
            logger.info("Error. Degree of polynomial is not 16!");
            return;
        }

        if (polyCoeff[16] < 0.0e0) {
            for (int i = 0; i < 17; i++) {
                polyCoeff[i] *= -1.0;
            }
        }

        if (logger.isLoggable(Level.FINE)) {
            StringBuilder string = new StringBuilder();
            string.append(" Polynomial Coefficients\n");
            for (int i = 0; i < 17; i++) {
                string.append(format(" %5d %15.6f\n", i, polyCoeff[i]));
            }
            string.append("\n");
            logger.fine(string.toString());
        }
    }

    /**
     * Multiplication of a matrix and a vector using the RealMatrix and RealVector libraries.
     *
     * @param ma Matrix A, two dimensional matrix.
     * @param mb Matrix B, a one dimensional vector.
     * @return Returns an array containing the product of ma and mb.
     */
    private static double[] matrixMultiplication(double[][] ma, double[] mb) {
        RealMatrix matrix = new Array2DRowRealMatrix(ma);
        matrix.multiply(matrix);
        RealVector vector = new ArrayRealVector(mb);
        RealVector result = matrix.operate(vector);
        return result.toArray();
    }

    /**
     * Polynomial multiply 1.
     *
     * @param u1 a double array.
     * @param u2 a double array.
     * @param p1 an int.
     * @param p2 an int.
     * @param u3 a double array.
     * @param p3 an int array.
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
     * @param u1 a double array.
     * @param u2 a double array.
     * @param p1 an int array.
     * @param p2 an int array.
     * @param u3 a double array.
     * @param p3 an int array.
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
     * @param u1 a double array.
     * @param u2 a double array.
     * @param u3 a double array.
     * @param u4 a double array.
     * @param p1 an int.
     * @param p2 an int.
     * @param p3 an int.
     * @param p4 an int.
     * @param u5 a double array.
     * @param p5 an int array.
     */
    public void polyMulSub1(double[] u1, double[] u2, double[] u3, double[] u4, int p1, int p2, int p3,
                            int p4, double[] u5, int[] p5) {

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
     * @param u1 a double array.
     * @param u2 a double array.
     * @param u3 a double array.
     * @param u4 a double array.
     * @param p1 an int array.
     * @param p2 an int array.
     * @param p3 an int array.
     * @param p4 an int array.
     * @param u5 a double array.
     * @param p5 an int array.
     */
    public void polyMulSub2(double[][] u1, double[][] u2, double[][] u3, double[][] u4, int[] p1,
                            int[] p2, int[] p3, int[] p4, double[][] u5, int[] p5) {
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
     * @param u1 a double array.
     * @param u2 a double array.
     * @param p1 an int array.
     * @param p2 an int array.
     * @param u3 a double array.
     * @param p3 an int array.
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
     * @param u1 a double array.
     * @param u2 a double array.
     * @param p1 an int array.
     * @param p2 an int array.
     * @param u3 a double array.
     * @param p3 an int array.
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
     * Calculates a rotation matrix about an origin with no angular specification using the quaternion.
     *
     * @param axis Quaternion input position vector.
     * @param angle Angle to calculate the quaternion of given axis.
     * @return Returns an output rotation matrix.
     */
    private static double[][] rotationMatrix(double[] axis, double angle) {
        // Calculate the quaternion.
        var tan_w = tan(angle);
        var tan_sqr = tan_w * tan_w;
        var tan1 = 1.0 + tan_sqr;
        var cosine = (1.0 - tan_sqr) / tan1;
        var sine = 2.0 * tan_w / tan1;
        var quaternion = new double[]{cosine, axis[0] * sine, axis[1] * sine, axis[2] * sine};

        // Calculate intermediate values.
        var b0 = 2.0 * quaternion[0];
        var b1 = 2.0 * quaternion[1];
        var q00 = b0 * quaternion[0] - 1.0;
        var q02 = b0 * quaternion[2];
        var q03 = b0 * quaternion[3];
        var q11 = b1 * quaternion[1];
        var q12 = b1 * quaternion[2];
        var q13 = b1 * quaternion[3];
        var b2 = 2.0 * quaternion[2];
        var b3 = 2.0 * quaternion[3];
        var q01 = b0 * quaternion[1];
        var q22 = b2 * quaternion[2];
        var q23 = b2 * quaternion[3];
        var q33 = b3 * quaternion[3];

        // Finish calculating the result.
        var result = new double[3][3];
        result[0][0] = q00 + q11;
        result[0][1] = q12 - q03;
        result[0][2] = q13 + q02;
        result[1][0] = q12 + q03;
        result[1][1] = q00 + q22;
        result[1][2] = q23 - q01;
        result[2][0] = q13 - q02;
        result[2][1] = q23 + q01;
        result[2][2] = q00 + q33;

        return result;
    }

    /**
     * Close a 3-residue loop by filling in the backbone atom coordinates and return the possible solution
     * set. This can include up to {@value MAX_SOLUTIONS} solutions.
     *
     * @param r_n1     a double array.
     * @param r_a1     a double array.
     * @param r_a3     a double array.
     * @param r_c3     a double array.
     * @param r_soln_n a double array.
     * @param r_soln_a a double array.
     * @param r_soln_c a double array.
     */
    public int solve3PepPoly(double[] r_n1, double[] r_a1, double[] r_a3, double[] r_c3,
                             double[][][] r_soln_n, double[][][] r_soln_a, double[][][] r_soln_c) {
        var polyCoeff = new double[MAX_SOLUTIONS + 1];
        var roots = new double[MAX_SOLUTIONS];

        if (!getInputAngles(r_n1, r_a1, r_a3, r_c3)) {
            // No solutions are available in this case.
            return 0;
        }

        getPolyCoeff(polyCoeff);

        var sturmMethod = new SturmMethod();
        var solutions = sturmMethod.solveSturm(MAX_SOLUTIONS, polyCoeff, roots);
        if (solutions > 0) {
            getCoordsFromPolyRoots(solutions, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c);
        } else {
            logger.info("Could not find alternative loop solutions using KIC.");
        }

        return solutions;
    }

    /**
     * testTwoConeExistenceSoln.
     *
     * @param tt     a double.
     * @param kx     a double.
     * @param et     a double.
     * @param ap     a double.
     */
    private boolean testTwoConeExistenceSoln(double tt, double kx, double et, double ap) {
        double at = ap - tt;
        double ex = kx + et;
        double abs_at = abs(at);
        return (abs_at <= ex);
    }
}
