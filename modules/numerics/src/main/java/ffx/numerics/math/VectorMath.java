/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.numerics.math;

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The VectorMath class is a simple math library that operates on 3-coordinate
 * double and float arrays. The design objectives are speed and no memory
 * consumption.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class VectorMath {

    private static final Logger logger = Logger.getLogger(VectorMath.class.getName());
    /**
     * Internal float work arrays.
     */
    private static final float fa[] = new float[3];
    private static final float fb[] = new float[3];
    private static final float fc[] = new float[3];
    private static final float fd[] = new float[3];
    private static final float fba[] = new float[3];
    private static final float fcb[] = new float[3];
    private static final float fdc[] = new float[3];
    private static final float ft[] = new float[3];
    private static final float fu[] = new float[3];
    private static final float ftu[] = new float[3];
    /**
     * Internal double work arrays.
     */
    private static final double da[] = new double[3];
    private static final double db[] = new double[3];
    private static final double dc[] = new double[3];
    private static final double dd[] = new double[3];
    private static final double dba[] = new double[3];
    private static final double dcb[] = new double[3];
    private static final double ddc[] = new double[3];
    private static final double dt[] = new double[3];
    private static final double du[] = new double[3];
    private static final double dtu[] = new double[3];
    private static final double eightpi2 = 8.0 * Math.PI * Math.PI;

    /**
     * <p>log.</p>
     *
     * @param v     an array of {@link double} objects.
     * @param label a {@link java.lang.String} object.
     */
    public static void log(double[] v, String label) {
        if (v == null) {
            return;
        }
        StringBuilder sb = null;
        if (label != null) {
            sb = new StringBuilder(String.format(" %16s = [", label));
        } else {
            sb = new StringBuilder(String.format(" %16s = [", "v"));
        }

        for (int i = 0; i < v.length; i++) {
            sb.append(String.format(" %16.8f", v[i]));
        }
        sb.append(" ]");
        logger.info(sb.toString());
    }

    /**
     * <p>
     * angle</p>
     *
     * @param i an array of double.
     * @param j an array of double.
     * @return a double.
     */
    public static double angle(double[] i, double[] j) {
        double x;
        norm(i, da);
        norm(j, db);
        x = dot(da, db);
        if (abs(x) > 1) {
            if (x > 0) {
                x = 1;
            } else {
                x = -1;
            }
        }
        return acos(x);
    }

    /**
     * <p>
     * angle</p>
     *
     * @param i an array of float.
     * @param j an array of float.
     * @return a float.
     */
    public static float angle(float[] i, float[] j) {
        float x;
        norm(i, fa);
        norm(j, fb);
        x = dot(fa, fb);
        if (abs(x) > 1) {
            logger.warning(String.format(" Angle: abs(dot) > 1 %10.6f", x));
            if (x > 0) {
                x = 1;
            } else {
                x = -1;
            }
        }
        return (float) acos(x);
    }

    /**
     * <p>
     * binomial</p>
     *
     * @param n a long.
     * @param k a long.
     * @return a long.
     */
    public static long binomial(long n, long k) {
        return factorial(n) / (factorial(n - k) * factorial(k));
    }

    /**
     * Finds the angle formed by three atoms
     *
     * @param i Atom position vector
     * @param j Atom position vector (central atom)
     * @param k Atom position vector
     * @return The angle in the range [ -pi, pi ]
     */
    public static double bondAngle(double[] i, double[] j, double[] k) {
        double x;
        diff(i, j, da);
        diff(k, j, db);
        norm(da, dc);
        norm(db, dd);
        x = dot(dc, dd);
        if (abs(x) > 1) {
            if (x > 0) {
                x = 1;
            } else {
                x = -1;
            }
        }
        return acos(x);
    }

    /**
     * Finds the angle formed by three atoms
     *
     * @param i Atom position vector
     * @param j Atom position vector (central atom)
     * @param k Atom position vector
     * @return The angle in the range [ -pi, pi ]
     */
    public static float bondAngle(float[] i, float[] j, float[] k) {
        float x;
        diff(i, j, fa);
        diff(k, j, fb);
        norm(fa, fc);
        norm(fb, fd);
        x = dot(fc, fd);
        if (abs(x) > 1) {
            logger.warning(String.format(" Bond Angle: abs(dot) > 1 %10.6f", x));
            if (x > 0) {
                x = 1;
            } else {
                x = -1;
            }
        }
        return (float) acos((double) x);
    }

    /**
     * Finds the cross-product between two vectors
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret The cross-product a x b
     */
    public static void cross(double[] a, double[] b, double ret[]) {
        ret[0] = a[1] * b[2] - a[2] * b[1];
        ret[1] = a[2] * b[0] - a[0] * b[2];
        ret[2] = a[0] * b[1] - a[1] * b[0];
    }

    /**
     * Finds the cross-product between two vectors
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret The cross-product a x b
     */
    public static void cross(float[] a, float[] b, float ret[]) {
        ret[0] = a[1] * b[2] - a[2] * b[1];
        ret[1] = a[2] * b[0] - a[0] * b[2];
        ret[2] = a[0] * b[1] - a[1] * b[0];
    }

    /**
     * returns the determinant for a 3x3 matrix
     *
     * @param m input matrix
     * @return determinant
     */
    public static double determinant3(double m[][]) {
        return (m[0][0] * m[1][1] * m[2][2]
                - m[0][0] * m[1][2] * m[2][1]
                + m[0][1] * m[1][2] * m[2][0]
                - m[0][1] * m[1][0] * m[2][2]
                + m[0][2] * m[1][0] * m[2][1]
                - m[0][2] * m[1][1] * m[2][0]);
    }

    /**
     * <p>
     * determinant3</p>
     *
     * @param m an array of double.
     * @return a double.
     */
    public static double determinant3(double m[]) {
        return (m[0] * m[1] * m[2]
                - m[0] * m[5] * m[5]
                + m[3] * m[5] * m[4]
                - m[3] * m[3] * m[2]
                + m[4] * m[3] * m[5]
                - m[4] * m[1] * m[4]);
    }

    /**
     * Finds the difference between two vectors
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret Return Values
     */
    public static void diff(double[] a, double[] b, double[] ret) {
        ret[0] = a[0] - b[0];
        ret[1] = a[1] - b[1];
        ret[2] = a[2] - b[2];
    }

    /**
     * Finds the difference between two vectors
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret Return Values
     */
    public static void diff(float[] a, float[] b, float[] ret) {
        ret[0] = a[0] - b[0];
        ret[1] = a[1] - b[1];
        ret[2] = a[2] - b[2];
    }

    /**
     * Finds the dihedral angle formed between 4 atoms
     *
     * @param a Atom position vector
     * @param b Atom position vector
     * @param c Atom position vector
     * @param d Atom position vector
     * @return The dihedral angle in the range [ -pi, pi ]
     */
    public static double dihedralAngle(double[] a, double[] b, double[] c,
                                       double[] d) {
        diff(b, a, dba);
        diff(c, b, dcb);
        diff(d, c, ddc);
        cross(dba, dcb, dt);
        cross(dcb, ddc, du);
        cross(dt, du, dtu);
        double rt = dot(dt, dt);
        double ru = dot(du, du);
        double rtu = sqrt(rt * ru);
        if (rtu != 0.0) {
            double rcb = r(dcb);
            double cosine = dot(dt, du) / rtu;
            double sine = dot(dcb, dtu) / (rcb * rtu);
            cosine = min(1.0, max(-1.0, cosine));
            double angle = acos(cosine);
            if (sine < 0.0) {
                angle = -angle;
            }
            return angle;
        }
        return 0;
    }

    /**
     * Finds the dihedral angle formed between 4 atoms
     *
     * @param a Atom position vector
     * @param b Atom position vector
     * @param c Atom position vector
     * @param d Atom position vector
     * @return The dihedral angle in the range [ -pi, pi ]
     */
    public static float dihedralAngle(float[] a, float[] b, float[] c, float[] d) {
        diff(b, a, fba);
        diff(c, b, fcb);
        diff(d, c, fdc);
        cross(fba, fcb, ft);
        cross(fcb, fdc, fu);
        cross(ft, fu, ftu);
        float rt = dot(ft, ft);
        float ru = dot(fu, fu);
        float rtu = (float) sqrt(rt * ru);
        if (rtu != 0.0) {
            float rcb = r(fcb);
            float cosine = dot(ft, fu) / rtu;
            float sine = dot(fcb, ftu) / (rcb * rtu);
            cosine = min(1.0f, max(-1.0f, cosine));
            float angle = (float) acos((double) cosine);
            if (sine < 0.0) {
                angle = -angle;
            }
            return angle;
        }
        return 0;
    }

    /**
     * Finds the squared distance between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return The squared distance between vectors a and b.
     */
    public static double dist2(double[] a, double[] b) {
        double dx = a[0] - b[0];
        double dy = a[1] - b[1];
        double dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * Finds the distance between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return The distance between vectors a and b
     */
    public static double dist(double[] a, double[] b) {
        return sqrt(dist2(a, b));
    }

    /**
     * Finds the squared distance between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return The squared distance between vectors a and b
     */
    public static float dist2(float[] a, float[] b) {
        float dx = a[0] - b[0];
        float dy = a[1] - b[1];
        float dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * Finds the distance between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return The distance between vectors a and b
     */
    public static float dist(float[] a, float[] b) {
        return (float) sqrt((double) dist2(a, b));
    }

    /**
     * Finds the dot product between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return The dot product of a and b
     */
    public static double dot(double[] a, double[] b) {
        return ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]));
    }

    /**
     * Finds the dot product between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return The dot product of a and b
     */
    public static float dot(float[] a, float[] b) {
        return ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]));
    }

    /**
     * Returns n!! Precondition: n .GE. -1 Returning 1 for -1 input is analogous
     * to Maple behavior.
     *
     * @param n long
     * @return long
     */
    public static long doubleFactorial(long n) {
        if (n < -1) {
            throw new RuntimeException("Underflow error in doublefactorial");
        } else if (n == 0 || n == 1 || n == -1) {
            return 1;
        } else {
            return n * doubleFactorial(n - 2);
        }
    }

    /**
     * Returns n!
     * <br>
     * Precondition: n .GE. 0 and n .LE. 20
     * <br>
     * Max long = 9223372036854775807
     * <br>
     * 20! = 2432902008176640000 is ok.
     * <br>
     * 21! returns an overflow: -4249290049419214848
     *
     * @param n long
     * @return long
     */
    public static long factorial(long n) {
        if (n < 0) {
            throw new RuntimeException("Underflow error in factorial");
        } else if (n > 20) {
            throw new RuntimeException("Overflow error in factorial");
        } else if (n == 0) {
            return 1;
        } else {
            return n * factorial(n - 1);
        }
    }

    /**
     * inverse of a 3x3 matrix
     *
     * @param m input matrix
     * @return matrix inverse
     */
    public static double[][] mat3Inverse(double m[][]) {
        double res[][] = new double[3][3];
        mat3Inverse(m, res);
        return res;
    }

    /**
     * <p>
     * mat3inverse</p>
     *
     * @param m   an array of double.
     * @param res an array of double.
     */
    public static void mat3Inverse(double m[][], double res[][]) {
        double det = determinant3(m);
        res[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
        res[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
        res[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
        res[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
        res[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
        res[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
        res[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
        res[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
        res[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
    }

    /**
     * vector times a matrix
     *
     * @param v input vector
     * @param m input matrix
     * @return vector product
     */
    public static double[] vec3Mat3(double v[], double m[][]) {
        double res[] = new double[3];
        vec3Mat3(v, m, res);
        return res;
    }

    /**
     * <p>
     * vec3mat3</p>
     *
     * @param v   an array of double.
     * @param m   an array of double.
     * @param res an array of double.
     */
    public static void vec3Mat3(double v[], double m[][], double res[]) {
        res[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
        res[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
        res[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];
    }

    /**
     * matrix times a vector
     *
     * @param m input matrix
     * @param v input vector
     * @return vector product
     */
    public static double[] mat3Vec3(double v[], double m[][]) {
        double res[] = new double[3];
        mat3Vec3(v, m, res);
        return res;
    }

    /**
     * <p>
     * mat3vec3</p>
     *
     * @param v   an array of double.
     * @param m   an array of double.
     * @param res an array of double.
     */
    public static void mat3Vec3(double v[], double m[][], double res[]) {
        res[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
        res[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
        res[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
    }

    /**
     * vector representation of a symmetric 3x3 matrix times a matrix
     *
     * @param v input vector of the form 11, 22, 33, 12, 13, 23
     * @param m input matrix
     * @return matrix product
     */
    public static double[][] symVec6Mat3(double v[], double m[][]) {
        double res[][] = new double[3][3];
        symVec6Mat3(v, m, res);
        return res;
    }

    /**
     * <p>
     * symvec6mat3</p>
     *
     * @param v   an array of double.
     * @param m   an array of double.
     * @param res an array of double.
     */
    public static void symVec6Mat3(double v[], double m[][], double res[][]) {
        res[0][0] = v[0] * m[0][0] + v[3] * m[1][0] + v[4] * m[2][0];
        res[0][1] = v[0] * m[0][1] + v[3] * m[1][1] + v[4] * m[2][1];
        res[0][2] = v[0] * m[0][2] + v[3] * m[1][2] + v[4] * m[2][2];
        res[1][0] = v[3] * m[0][0] + v[1] * m[1][0] + v[5] * m[2][0];
        res[1][1] = v[3] * m[0][1] + v[1] * m[1][1] + v[5] * m[2][1];
        res[1][2] = v[3] * m[0][2] + v[1] * m[1][2] + v[5] * m[2][2];
        res[2][0] = v[4] * m[0][0] + v[5] * m[1][0] + v[2] * m[2][0];
        res[2][1] = v[4] * m[0][1] + v[5] * m[1][1] + v[2] * m[2][1];
        res[2][2] = v[4] * m[0][2] + v[5] * m[1][2] + v[2] * m[2][2];
    }

    /**
     * matrix times a vector representation of a symmetric 3x3 matrix
     *
     * @param m input matrix
     * @param v input vector of the form 11, 22, 33, 12, 13, 23
     * @return matrix product
     */
    public static double[][] mat3SymVec6(double m[][], double v[]) {
        double res[][] = new double[3][3];
        VectorMath.mat3SymVec6(m, v, res);
        return res;
    }

    /**
     * <p>
     * mat3symvec6</p>
     *
     * @param m   an array of double.
     * @param v   an array of double.
     * @param res an array of double.
     */
    public static void mat3SymVec6(double m[][], double v[], double res[][]) {
        res[0][0] = m[0][0] * v[0] + m[0][1] * v[3] + m[0][2] * v[4];
        res[0][1] = m[0][0] * v[3] + m[0][1] * v[1] + m[0][2] * v[5];
        res[0][2] = m[0][0] * v[4] + m[0][1] * v[5] + m[0][2] * v[2];
        res[1][0] = m[1][0] * v[0] + m[1][1] * v[3] + m[1][2] * v[4];
        res[1][1] = m[1][0] * v[3] + m[1][1] * v[1] + m[1][2] * v[5];
        res[1][2] = m[1][0] * v[4] + m[1][1] * v[5] + m[1][2] * v[2];
        res[2][0] = m[2][0] * v[0] + m[2][1] * v[3] + m[2][2] * v[4];
        res[2][1] = m[2][0] * v[3] + m[2][1] * v[1] + m[2][2] * v[5];
        res[2][2] = m[2][0] * v[4] + m[2][1] * v[5] + m[2][2] * v[2];
    }

    /**
     * matrix times a matrix
     *
     * @param m1 first input matrix
     * @param m2 second input matrix
     * @return matrix product
     */
    public static double[][] mat3Mat3(double m1[][], double m2[][]) {
        double res[][] = new double[3][3];
        mat3Mat3(m1, m2, res);
        return res;
    }

    /**
     * <p>
     * mat3mat3</p>
     *
     * @param m1  an array of double.
     * @param m2  an array of double.
     * @param res an array of double.
     */
    public static void mat3Mat3(double m1[][], double m2[][], double res[][]) {
        res[0][0] = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0];
        res[0][1] = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1];
        res[0][2] = m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2];
        res[1][0] = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0];
        res[1][1] = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1];
        res[1][2] = m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2];
        res[2][0] = m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0];
        res[2][1] = m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1];
        res[2][2] = m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2];
    }

    /**
     * scalar times a matrix times a matrix
     *
     * @param scalar input scalar
     * @param m1     first input matrix
     * @param m2     second input matrix
     * @return matrix product
     */
    public static double[][] scalarMat3Mat3(double scalar, double m1[][], double m2[][]) {
        double res[][] = new double[3][3];
        scalarMat3Mat3(scalar, m1, m2, res);
        return res;
    }

    /**
     * <p>
     * scalarmat3mat3</p>
     *
     * @param scalar a double.
     * @param m1     an array of double.
     * @param m2     an array of double.
     * @param res    an array of double.
     */
    public static void scalarMat3Mat3(double scalar, double m1[][], double m2[][], double res[][]) {
        res[0][0] = (scalar * m1[0][0]) * m2[0][0] + (scalar * m1[0][1]) * m2[1][0] + (scalar * m1[0][2]) * m2[2][0];
        res[0][1] = (scalar * m1[0][0]) * m2[0][1] + (scalar * m1[0][1]) * m2[1][1] + (scalar * m1[0][2]) * m2[2][1];
        res[0][2] = (scalar * m1[0][0]) * m2[0][2] + (scalar * m1[0][1]) * m2[1][2] + (scalar * m1[0][2]) * m2[2][2];
        res[1][0] = (scalar * m1[1][0]) * m2[0][0] + (scalar * m1[1][1]) * m2[1][0] + (scalar * m1[1][2]) * m2[2][0];
        res[1][1] = (scalar * m1[1][0]) * m2[0][1] + (scalar * m1[1][1]) * m2[1][1] + (scalar * m1[1][2]) * m2[2][1];
        res[1][2] = (scalar * m1[1][0]) * m2[0][2] + (scalar * m1[1][1]) * m2[1][2] + (scalar * m1[1][2]) * m2[2][2];
        res[2][0] = (scalar * m1[2][0]) * m2[0][0] + (scalar * m1[2][1]) * m2[1][0] + (scalar * m1[2][2]) * m2[2][0];
        res[2][1] = (scalar * m1[2][0]) * m2[0][1] + (scalar * m1[2][1]) * m2[1][1] + (scalar * m1[2][2]) * m2[2][1];
        res[2][2] = (scalar * m1[2][0]) * m2[0][2] + (scalar * m1[2][1]) * m2[1][2] + (scalar * m1[2][2]) * m2[2][2];
    }

    /**
     * Normalizes a vector
     *
     * @param n   A vector to be normalized.
     * @param ret The normalized vector.
     */
    public static void norm(double[] n, double[] ret) {
        double length;
        length = r(n);
        ret[0] = n[0] / length;
        ret[1] = n[1] / length;
        ret[2] = n[2] / length;
    }

    /**
     * Normalizes a vector
     *
     * @param n   A vector to be normalized.
     * @param ret The normalized vector.
     */
    public static void norm(float[] n, float[] ret) {
        float length;
        length = r(n);
        ret[0] = n[0] / length;
        ret[1] = n[1] / length;
        ret[2] = n[2] / length;
    }

    /**
     * <p>
     * printVector</p>
     *
     * @param v an array of double.
     */
    public static void printVector(double v[]) {
        StringBuilder sb = new StringBuilder("Vector ( ");
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%g ", v[i]));
        }
        sb.append(")");
        logger.info(sb.toString());
    }

    /**
     * Finds the length of a vector
     *
     * @param d A vector to find the length of.
     * @return Length of vector d.
     */
    public static double r(double[] d) {
        return sqrt((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
    }

    /**
     * Finds the length of a vector
     *
     * @param d A vector to find the length of.
     * @return Length of vector d.
     */
    public static float r(float[] d) {
        return (float) sqrt((double) (d[0] * d[0] + d[1] * d[1] + d[2]
                * d[2]));
    }

    /**
     * Finds the length^2 of a vector
     *
     * @param d A vector to find the length of.
     * @return Length^2 of vector d.
     */
    public static double rsq(double[] d) {
        return (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    }

    /**
     * Scales a vector
     *
     * @param n   A vector to be scaled
     * @param a   A scaler value
     * @param ret The scaled vector
     */
    public static void scalar(double[] n, double a, double[] ret) {
        ret[0] = n[0] * a;
        ret[1] = n[1] * a;
        ret[2] = n[2] * a;
    }

    /**
     * Scales a vector
     *
     * @param n   A vector to be scaled
     * @param a   A scaler value
     * @param ret The scaled Vector
     */
    public static void scalar(float[] n, float a, float[] ret) {
        ret[0] = n[0] * a;
        ret[1] = n[1] * a;
        ret[2] = n[2] * a;
    }

    /**
     * <p>
     * sum</p>
     *
     * @param a   an array of double.
     * @param b   an array of double.
     * @param ret an array of double.
     */
    public static void sum(double[] a, double[] b, double[] ret) {
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
    }

    /**
     * <p>
     * sum</p>
     *
     * @param a   an array of float.
     * @param b   an array of float.
     * @param ret an array of float.
     */
    public static void sum(float[] a, float[] b, float[] ret) {
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
    }

    /**
     * <p>
     * transpose3</p>
     *
     * @param m an array of double.
     * @return an array of double.
     */
    public static double[][] transpose3(double m[][]) {
        double t[][] = new double[3][3];
        transpose3(m, t);
        return t;
    }

    /**
     * <p>
     * transpose3</p>
     *
     * @param m an array of double.
     * @param t an array of double.
     */
    public static void transpose3(double m[][], double t[][]) {
        t[0][0] = m[0][0];
        t[0][1] = m[1][0];
        t[0][2] = m[2][0];
        t[1][0] = m[0][1];
        t[1][1] = m[1][1];
        t[1][2] = m[2][1];
        t[2][0] = m[0][2];
        t[2][1] = m[1][2];
        t[2][2] = m[2][2];
    }

    /**
     * <p>
     * b2u</p>
     *
     * @param b a double.
     * @return a double.
     */
    public static double b2u(double b) {
        return b / eightpi2;
    }

    /**
     * <p>
     * u2b</p>
     *
     * @param u a double.
     * @return a double.
     */
    public static double u2b(double u) {
        return u * eightpi2;
    }
}
