/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.numerics;

import static java.lang.Math.*;

import java.util.logging.Logger;

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
    private static final float fa[] = new float[3];
    private static final float fb[] = new float[3];
    private static final float fc[] = new float[3];
    private static final float fd[] = new float[3];
    private static final float[] fba = new float[3];
    private static final float[] fcb = new float[3];
    private static final float[] fdc = new float[3];
    private static final float[] ft = new float[3];
    private static final float[] fu = new float[3];
    private static final float[] ftu = new float[3];
    private static final double da[] = new double[3];
    private static final double db[] = new double[3];
    private static final double dc[] = new double[3];
    private static final double dd[] = new double[3];
    private static final double[] dba = new double[3];
    private static final double[] dcb = new double[3];
    private static final double[] ddc = new double[3];
    private static final double[] dt = new double[3];
    private static final double[] du = new double[3];
    private static final double[] dtu = new double[3];

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

    public static float angle(float[] i, float[] j) {
        float x;
        norm(i, fa);
        norm(j, fb);
        x = dot(fa, fb);
        if (abs(x) > 1) {
            logger.warning("angle: abs(dot) > 1 " + x);
            if (x > 0) {
                x = 1;
            } else {
                x = -1;
            }
        }
        return (float) acos(x);
    }

    public static long binomial(long n, long k) {
        return factorial(n) / (factorial(n - k) * factorial(k));
    }

    /**
     * Finds the angle formed by three atoms
     *
     * @param i
     *            Atom position vector
     * @param j
     *            Atom position vector (central atom)
     * @param k
     *            Atom position vector
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
     * @param i
     *            Atom position vector
     * @param j
     *            Atom position vector (central atom)
     * @param k
     *            Atom position vector
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
            logger.warning("bondAngle: abs(dot) > 1 " + x);
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
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @param ret
     *            The cross-product a x b
     */
    public static void cross(double[] a, double[] b, double ret[]) {
        ret[0] = a[1] * b[2] - a[2] * b[1];
        ret[1] = a[2] * b[0] - a[0] * b[2];
        ret[2] = a[0] * b[1] - a[1] * b[0];
    }

    /**
     * Finds the cross-product between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @param ret
     *            The cross-product a x b
     */
    public static void cross(float[] a, float[] b, float ret[]) {
        ret[0] = a[1] * b[2] - a[2] * b[1];
        ret[1] = a[2] * b[0] - a[0] * b[2];
        ret[2] = a[0] * b[1] - a[1] * b[0];
    }

    /**
     * returns the determinant for a 3x3 matrix
     *
     * @param m  input matrix
     *
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
     * Finds the difference between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @param ret
     *            Return Values
     */
    public static void diff(double[] a, double[] b, double[] ret) {
        ret[0] = a[0] - b[0];
        ret[1] = a[1] - b[1];
        ret[2] = a[2] - b[2];
    }

    /**
     * Finds the difference between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @param ret
     *            Return Values
     */
    public static void diff(float[] a, float[] b, float[] ret) {
        ret[0] = a[0] - b[0];
        ret[1] = a[1] - b[1];
        ret[2] = a[2] - b[2];
    }

    /**
     * Finds the dihedral angle formed between 4 atoms
     *
     * @param a
     *            Atom position vector
     * @param b
     *            Atom position vector
     * @param c
     *            Atom position vector
     * @param d
     *            Atom position vector
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
     * @param a
     *            Atom position vector
     * @param b
     *            Atom position vector
     * @param c
     *            Atom position vector
     * @param d
     *            Atom position vector
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
     * Finds the distance between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @return The distance between vectors a and b
     */
    public static double dist(double[] a, double[] b) {
        double x;
        x = (a[0] - b[0]) * (a[0] - b[0]);
        x += (a[1] - b[1]) * (a[1] - b[1]);
        x += (a[2] - b[2]) * (a[2] - b[2]);
        return sqrt(x);
    }

    /**
     * Finds the distance between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @return The distance between vectors a and b
     */
    public static float dist(float[] a, float[] b) {
        float x;
        x = (a[0] - b[0]) * (a[0] - b[0]);
        x += (a[1] - b[1]) * (a[1] - b[1]);
        x += (a[2] - b[2]) * (a[2] - b[2]);
        return (float) sqrt((double) x);
    }

    /**
     * Finds the dot product between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @return The dot product of a and b
     */
    public static double dot(double[] a, double[] b) {
        return ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]));
    }

    /**
     * Finds the dot product between two vectors
     *
     * @param a
     *            First vector
     * @param b
     *            Second vector
     * @return The dot product of a and b
     */
    public static float dot(float[] a, float[] b) {
        return ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]));
    }

    /**
     * Returns n!! Precondition: n >= -1 Returning 1 for -1 input is analogous
     * to Maple behavior.
     *
     * @param n
     *            long
     * @return long
     */
    public static long doublefactorial(long n) {
        if (n < -1) {
            throw new RuntimeException("Underflow error in doublefactorial");
        } else if (n == 0 || n == 1 || n == -1) {
            return 1;
        } else {
            return n * doublefactorial(n - 2);
        }
    }

    /**
     * Returns n! Precondition: n >= 0 and n <= 20 Max long =
     * 9223372036854775807 20! = 2432902008176640000 is ok. 21! returns an
     * overflow: -4249290049419214848
     *
     * @param n
     *            long
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
     * vector times a matrix
     *
     * @param v input vector
     *
     * @param m input matrix
     *
     * @return vector product
     */
    public static double[] vec3mat3(double v[], double m[][]) {
        double res[] = new double[3];
        res[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
        res[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
        res[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];

        return res;
    }

    /**
     * matrix times a vector
     *
     * @param m input matrix
     *
     * @param v input vector
     *
     * @return vector product
     */
    public static double[] mat3vec3(double v[], double m[][]) {
        double res[] = new double[3];
        res[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
        res[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
        res[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];

        return res;
    }

    /**
     * matrix times a matrix
     *
     * @param m1 first input matrix
     *
     * @param m2 second input matrix
     *
     * @return matrix product
     */
    public static double[][] mat3mat3(double m1[][], double m2[][]) {
        double res[][] = new double[3][3];
        res[0][0] = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0];
        res[0][1] = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1];
        res[0][2] = m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2];
        res[1][0] = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0];
        res[1][1] = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1];
        res[1][2] = m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2];
        res[2][0] = m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0];
        res[2][1] = m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1];
        res[2][2] = m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2];

        return res;
    }

    /**
     * scalar times a matrix times a matrix
     *
     * @param scalar input scalar
     * @param m1 first input matrix
     *
     * @param m2 second input matrix
     *
     * @return matrix product
     */
    public static double[][] scalarmat3mat3(double scalar, double m1[][], double m2[][]) {
        double res[][] = new double[3][3];
        res[0][0] = scalar * (m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0]);
        res[0][1] = scalar * (m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1]);
        res[0][2] = scalar * (m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2]);
        res[1][0] = scalar * (m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0]);
        res[1][1] = scalar * (m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1]);
        res[1][2] = scalar * (m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2]);
        res[2][0] = scalar * (m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0]);
        res[2][1] = scalar * (m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1]);
        res[2][2] = scalar * (m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2]);

        return res;
    }

    /**
     * Normalizes a vector
     *
     * @param n
     *            A vector to be normalized.
     * @param ret
     *            The normalized vector.
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
     * @param n
     *            A vector to be normalized.
     * @param ret
     *            The normalized vector.
     */
    public static void norm(float[] n, float[] ret) {
        float length;
        length = r(n);
        ret[0] = n[0] / length;
        ret[1] = n[1] / length;
        ret[2] = n[2] / length;
    }

    public static void printVector(double v[]) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < v.length; i++) {
            sb.append(String.format("%g ", v[i]));
        }
        logger.info("Vector ( " + sb.toString() + ")");
    }

    /**
     * Finds the length of a vector
     *
     * @param d
     *            A vector to find the length of.
     * @return Length of vector d.
     */
    public static double r(double[] d) {
        return sqrt((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
    }

    /**
     * Finds the length of a vector
     *
     * @param d
     *            A vector to find the length of.
     * @return Length of vector d.
     */
    public static float r(float[] d) {
        return (float) sqrt((double) (d[0] * d[0] + d[1] * d[1] + d[2]
                * d[2]));
    }

    /**
     * Finds the length^2 of a vector
     *
     * @param d
     *            A vector to find the length of.
     * @return Length^2 of vector d.
     */
    public static double rsq(double[] d) {
        return (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    }

    /**
     * Scales a vector
     *
     * @param n
     *            A vector to be scaled
     * @param a
     *            A scaler value
     * @param ret
     *            The scaled vector
     */
    public static void scalar(double[] n, double a, double[] ret) {
        ret[0] = n[0] * a;
        ret[1] = n[1] * a;
        ret[2] = n[2] * a;
    }

    /**
     * Scales a vector
     *
     * @param n
     *            A vector to be scaled
     * @param a
     *            A scaler value
     * @param ret
     *            The scaled Vector
     */
    public static void scalar(float[] n, float a, float[] ret) {
        ret[0] = n[0] * a;
        ret[1] = n[1] * a;
        ret[2] = n[2] * a;
    }

    public static void sum(double[] a, double[] b, double[] ret) {
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
    }

    public static void sum(float[] a, float[] b, float[] ret) {
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
    }
}
