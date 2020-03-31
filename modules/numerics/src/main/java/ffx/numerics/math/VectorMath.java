//******************************************************************************
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
//******************************************************************************
package ffx.numerics.math;

import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The VectorMath class is a simple math library that operates mainly on 3-coordinate
 * double and float arrays.
 * <p>
 * All methods are thread-safe.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class VectorMath {

    private static final Logger logger = Logger.getLogger(VectorMath.class.getName());
    private static final double eightPi2 = 8.0 * PI * PI;

    /**
     * <p>
     * angle</p>
     *
     * @param i an array of double.
     * @param j an array of double.
     * @return Returns the angle.
     */
    public static double angle(double[] i, double[] j) {
        var x = dot(norm(i), norm(j));
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
     * @return Returns the angle.
     */
    public static float angle(float[] i, float[] j) {
        var x = dot(norm(i), norm(j));
        if (abs(x) > 1) {
            logger.warning(format(" Angle: abs(dot) > 1 %10.6f", x));
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
     * b2u</p>
     *
     * @param b a double.
     * @return a double.
     */
    public static double b2u(double b) {
        return b / eightPi2;
    }

    /**
     * <p>
     * binomial</p>
     *
     * @param n a long.
     * @param k a long.
     * @return Returns the binomial of (n,k).
     */
    public static long binomial(long n, long k) {
        return factorial(n) / (factorial(n - k) * factorial(k));
    }

    /**
     * Finds the angle formed by three atoms.
     *
     * @param i Atom position vector.
     * @param j Atom position vector (central atom).
     * @param k Atom position vector.
     * @return Return the angle in the range [ -pi, pi ].
     */
    public static double bondAngle(double[] i, double[] j, double[] k) {
        return angle(diff(i, j), diff(k, j));
    }

    /**
     * Finds the angle formed by three atoms
     *
     * @param i Atom position vector.
     * @param j Atom position vector (central atom).
     * @param k Atom position vector.
     * @return Returns the angle in the range [ -pi, pi ].
     */
    public static float bondAngle(float[] i, float[] j, float[] k) {
        return angle(diff(i, j), diff(k, j));
    }

    /**
     * Finds the cross-product between two vectors
     *
     * @param a First vector
     * @param b Second vector
     * @return Returns the cross-product.
     */
    public static double[] cross(double[] a, double[] b) {
        return cross(a, b, new double[3]);
    }

    /**
     * Finds the cross-product between two vectors
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret The cross-product a x b.
     * @return Returns the cross-product ret.
     */
    public static double[] cross(double[] a, double[] b, double[] ret) {
        ret[0] = a[1] * b[2] - a[2] * b[1];
        ret[1] = a[2] * b[0] - a[0] * b[2];
        ret[2] = a[0] * b[1] - a[1] * b[0];
        return ret;
    }

    /**
     * Finds the cross-product between two vectors.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the cross-product.
     */
    public static float[] cross(float[] a, float[] b) {
        return cross(a, b, new float[3]);
    }

    /**
     * Finds the cross-product between two vectors.
     *
     * @param a   First vector.
     * @param b   Second vector.
     * @param ret The cross-product a x b.
     * @return Returns the cross-product ret.
     */
    public static float[] cross(float[] a, float[] b, float[] ret) {
        ret[0] = a[1] * b[2] - a[2] * b[1];
        ret[1] = a[2] * b[0] - a[0] * b[2];
        ret[2] = a[0] * b[1] - a[1] * b[0];
        return ret;
    }

    /**
     * Returns the determinant for a 3x3 matrix.
     *
     * @param m input matrix.
     * @return Returns the determinant.
     */
    public static double determinant3(double[][] m) {
        return m[0][0] * m[1][1] * m[2][2]
                - m[0][0] * m[1][2] * m[2][1]
                + m[0][1] * m[1][2] * m[2][0]
                - m[0][1] * m[1][0] * m[2][2]
                + m[0][2] * m[1][0] * m[2][1]
                - m[0][2] * m[1][1] * m[2][0];
    }

    /**
     * <p>
     * determinant3</p>
     *
     * @param m an array of double.
     * @return Returns the determinant.
     */
    public static double determinant3(double[] m) {
        return m[0] * m[1] * m[2]
                - m[0] * m[5] * m[5]
                + m[3] * m[5] * m[4]
                - m[3] * m[3] * m[2]
                + m[4] * m[3] * m[5]
                - m[4] * m[1] * m[4];
    }

    /**
     * Finds the difference between two vectors.
     *
     * @param a First vector
     * @param b Second vector
     * @return Returns the difference.
     */
    public static double[] diff(double[] a, double[] b) {
        return diff(a, b, new double[3]);
    }

    /**
     * Finds the difference between two vectors.
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret Return Values
     * @return Returns the difference ret.
     */
    public static double[] diff(double[] a, double[] b, double[] ret) {
        ret[0] = a[0] - b[0];
        ret[1] = a[1] - b[1];
        ret[2] = a[2] - b[2];
        return ret;
    }

    /**
     * Finds the difference between two vectors.
     *
     * @param a First vector
     * @param b Second vector
     * @return Returns the difference ret.
     */
    public static float[] diff(float[] a, float[] b) {
        return diff(a, b, new float[3]);
    }

    /**
     * Finds the difference between two vectors
     *
     * @param a   First vector
     * @param b   Second vector
     * @param ret Return Values
     * @return Returns the difference ret.
     */
    public static float[] diff(float[] a, float[] b, float[] ret) {
        ret[0] = a[0] - b[0];
        ret[1] = a[1] - b[1];
        ret[2] = a[2] - b[2];
        return ret;
    }

    /**
     * Finds the dihedral angle formed between 4 atoms.
     *
     * @param a Atom position vector.
     * @param b Atom position vector.
     * @param c Atom position vector.
     * @param d Atom position vector.
     * @return The dihedral angle in the range [ -pi, pi ].
     */
    public static double dihedralAngle(double[] a, double[] b, double[] c, double[] d) {
        var ba = diff(b, a);
        var cb = diff(c, b);
        var dc = diff(d, c);
        var t = cross(ba, cb);
        var u = cross(cb, dc);
        var rt = dot(t, t);
        var ru = dot(u, u);
        var rtu = sqrt(rt * ru);
        if (rtu != 0.0) {
            var rcb = r(cb);
            var cosine = dot(t, u) / rtu;
            var tu = cross(t, u);
            var sine = dot(cb, tu) / (rcb * rtu);
            cosine = min(1.0, max(-1.0, cosine));
            var angle = acos(cosine);
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
     * @param a Atom position vector.
     * @param b Atom position vector.
     * @param c Atom position vector.
     * @param d Atom position vector.
     * @return The dihedral angle in the range [ -pi, pi ].
     */
    public static float dihedralAngle(float[] a, float[] b, float[] c, float[] d) {
        var ba = diff(b, a);
        var cb = diff(c, b);
        var dc = diff(d, c);
        var t = cross(ba, cb);
        var u = cross(cb, dc);
        var rt = dot(t, t);
        var ru = dot(u, u);
        var rtu = (float) sqrt(rt * ru);
        if (rtu != 0.0) {
            var rcb = r(cb);
            var cosine = dot(t, u) / rtu;
            var tu = cross(t, u);
            var sine = dot(cb, tu) / (rcb * rtu);
            cosine = min(1.0f, max(-1.0f, cosine));
            var angle = (float) acos(cosine);
            if (sine < 0.0) {
                angle = -angle;
            }
            return angle;
        }
        return 0;
    }

    /**
     * Finds the distance between two vectors.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the distance between vectors a and b.
     */
    public static double dist(double[] a, double[] b) {
        return sqrt(dist2(a, b));
    }

    /**
     * Finds the distance between two vectors.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the distance between vectors a and b.
     */
    public static float dist(float[] a, float[] b) {
        return (float) sqrt(dist2(a, b));
    }

    /**
     * Finds the squared distance between two vectors
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the squared distance between vectors a and b.
     */
    public static double dist2(double[] a, double[] b) {
        var dx = a[0] - b[0];
        var dy = a[1] - b[1];
        var dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * Finds the squared distance between two vectors.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the squared distance between vectors a and b.
     */
    public static float dist2(float[] a, float[] b) {
        var dx = a[0] - b[0];
        var dy = a[1] - b[1];
        var dz = a[2] - b[2];
        return dx * dx + dy * dy + dz * dz;
    }

    /**
     * Finds the dot product between two vectors.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the dot product of a and b.
     */
    public static double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    /**
     * Finds the dot product between two vectors.
     *
     * @param a First vector.
     * @param b Second vector.
     * @return Returns the dot product of a and b.
     */
    public static float dot(float[] a, float[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    /**
     * Returns n!! Precondition: n .GE. -1 Returning 1 for -1 input is analogous
     * to Maple behavior.
     *
     * @param n long.
     * @return Returns the n!!.
     */
    public static long doubleFactorial(long n) {
        if (n < -1) {
            throw new RuntimeException("Underflow error in doubleFactorial");
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
     * @param n long.
     * @return Returns n!.
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
     * <p>
     * logVector</p>
     *
     * @param v an array of double.
     */
    public static void logVector(double[] v) {
        logger.info(vectorToString(v));
    }

    /**
     * <p>logVector.</p>
     *
     * @param v     an array of {@link double} objects.
     * @param label a {@link java.lang.String} object.
     */
    public static void logVector(double[] v, String label) {
        logger.info(vectorToString(v, label));
    }

    /**
     * inverse of a 3x3 matrix.
     *
     * @param m input matrix.
     * @return Returns the matrix inverse.
     */
    public static double[][] mat3Inverse(double[][] m) {
        return mat3Inverse(m, new double[3][3]);
    }

    /**
     * <p>
     * mat3inverse</p>
     *
     * @param m   an array of double.
     * @param res an array of double.
     * @return Returns the matrix res.
     */
    public static double[][] mat3Inverse(double[][] m, double[][] res) {
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
        return res;
    }

    /**
     * Matrix times a matrix.
     *
     * @param m1 First input matrix.
     * @param m2 Second input matrix.
     * @return Returns the matrix product.
     */
    public static double[][] mat3Mat3(double[][] m1, double[][] m2) {
        return mat3Mat3(m1, m2, new double[3][3]);
    }

    /**
     * <p>
     * mat3mat3</p>
     *
     * @param m1  an array of double.
     * @param m2  an array of double.
     * @param res an array of double.
     * @return Returns the matrix res.
     */
    public static double[][] mat3Mat3(double[][] m1, double[][] m2, double[][] res) {
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
     * matrix times a vector representation of a symmetric 3x3 matrix
     *
     * @param m input matrix
     * @param v input vector of the form 11, 22, 33, 12, 13, 23
     * @return matrix product
     */
    public static double[][] mat3SymVec6(double[][] m, double[] v) {
        return mat3SymVec6(m, v, new double[3][3]);
    }

    /**
     * <p>
     * mat3SymVec6</p>
     *
     * @param m   an array of double.
     * @param v   an array of double.
     * @param res an array of double.
     * @return Returns the matrix res.
     */
    public static double[][] mat3SymVec6(double[][] m, double[] v, double[][] res) {
        res[0][0] = m[0][0] * v[0] + m[0][1] * v[3] + m[0][2] * v[4];
        res[0][1] = m[0][0] * v[3] + m[0][1] * v[1] + m[0][2] * v[5];
        res[0][2] = m[0][0] * v[4] + m[0][1] * v[5] + m[0][2] * v[2];
        res[1][0] = m[1][0] * v[0] + m[1][1] * v[3] + m[1][2] * v[4];
        res[1][1] = m[1][0] * v[3] + m[1][1] * v[1] + m[1][2] * v[5];
        res[1][2] = m[1][0] * v[4] + m[1][1] * v[5] + m[1][2] * v[2];
        res[2][0] = m[2][0] * v[0] + m[2][1] * v[3] + m[2][2] * v[4];
        res[2][1] = m[2][0] * v[3] + m[2][1] * v[1] + m[2][2] * v[5];
        res[2][2] = m[2][0] * v[4] + m[2][1] * v[5] + m[2][2] * v[2];
        return res;
    }

    /**
     * matrix times a vector
     *
     * @param m input matrix.
     * @param v input vector.
     * @return Returns the vector product.
     */
    public static double[] mat3Vec3(double[] v, double[][] m) {
        return mat3Vec3(v, m, new double[3]);
    }

    /**
     * <p>
     * mat3vec3</p>
     *
     * @param v   an array of double.
     * @param m   an array of double.
     * @param res an array of double.
     * @return Returns the matrix res.
     */
    public static double[] mat3Vec3(double[] v, double[][] m, double[] res) {
        res[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
        res[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
        res[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
        return res;
    }

    /**
     * This is an atypical mod function used by crystallography methods.
     *
     * <p>
     * mod</p>
     *
     * @param a Value to mod.
     * @param b Value to mod by.
     * @return Positive a % b.
     */
    public static double mod(double a, double b) {
        var res = a % b;
        if (res < 0.0) {
            res += b;
        }
        return res;
    }

    /**
     * Atypical mod function used to move a value into the range lb &lt;= value &lt; ub,
     * assuming the domain is periodic with a period of (ub - lb).
     *
     * @param value Value to move between bounds.
     * @param lb    Lower bound.
     * @param ub    Upper bound.
     * @return Returns periodic copy of value, in the range lb &lt;= value &lt; ub.
     */
    public static double modToRange(double value, double lb, double ub) {
        value -= lb;
        var range = ub - lb;
        value = mod(value, range);
        value += lb;
        return value;
    }

    /**
     * Normalizes a vector.
     *
     * @param n A vector to be normalized.
     * @return Returns the normalized vector.
     */
    public static double[] norm(double[] n) {
        return scalar(n, 1.0 / r(n), new double[3]);
    }

    /**
     * Normalizes a vector.
     *
     * @param n   A vector to be normalized.
     * @param ret The normalized vector.
     * @return Returns the normalized vector.
     */
    public static double[] norm(double[] n, double[] ret) {
        return scalar(n, 1.0 / r(n), ret);
    }

    /**
     * Normalizes a vector.
     *
     * @param n A vector to be normalized.
     * @return Returns the normalized vector.
     */
    public static float[] norm(float[] n) {
        return scalar(n, (float) 1.0 / r(n), new float[3]);
    }

    /**
     * Normalizes a vector.
     *
     * @param n   A vector to be normalized.
     * @param ret The normalized vector.
     * @return Returns the normalized vector.
     */
    public static float[] norm(float[] n, float[] ret) {
        return scalar(n, (float) 1.0 / r(n), ret);
    }

    /**
     * Finds the length of a vector.
     *
     * @param d A vector to find the length of.
     * @return Length of vector d.
     */
    public static double r(double[] d) {
        return sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    }

    /**
     * Finds the length of a vector.
     *
     * @param d A vector to find the length of.
     * @return Returns the length of vector d.
     */
    public static float r(float[] d) {
        return (float) sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    }

    /**
     * Finds the length^2 of a vector.
     *
     * @param d A vector to find the length of.
     * @return Returns the length^2 of vector d.
     */
    public static double rsq(double[] d) {
        return d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    }

    /**
     * Scales a vector.
     *
     * @param n A vector to be scaled.
     * @param a A scalar value.
     * @return Returns the scaled vector.
     */
    public static double[] scalar(double[] n, double a) {
        return scalar(n, a, new double[3]);
    }

    /**
     * Scales a vector.
     *
     * @param n   A vector to be scaled.
     * @param a   A scalar value.
     * @param ret The scaled vector.
     * @return Returns the array ret.
     */
    public static double[] scalar(double[] n, double a, double[] ret) {
        ret[0] = n[0] * a;
        ret[1] = n[1] * a;
        ret[2] = n[2] * a;
        return ret;
    }

    /**
     * Scales a vector.
     *
     * @param n A vector to be scaled.
     * @param a A scalar value.
     * @return Returns the scaled vector.
     */
    public static float[] scalar(float[] n, float a) {
        return scalar(n, a, new float[3]);
    }

    /**
     * Scales a vector.
     *
     * @param n   A vector to be scaled.
     * @param a   A scalar value.
     * @param ret The scaled Vector.
     * @return Returns the array ret.
     */
    public static float[] scalar(float[] n, float a, float[] ret) {
        ret[0] = n[0] * a;
        ret[1] = n[1] * a;
        ret[2] = n[2] * a;
        return ret;
    }

    /**
     * scalar times a matrix times a matrix.
     *
     * @param scalar Input scalar.
     * @param m1     First input matrix.
     * @param m2     Second input matrix.
     * @return Returns the matrix product.
     */
    public static double[][] scalarMat3Mat3(double scalar, double[][] m1, double[][] m2) {
        return scalarMat3Mat3(scalar, m1, m2, new double[3][3]);
    }

    /**
     * <p>
     * scalarMat3mat3</p>
     *
     * @param scalar a double.
     * @param m1     an array of double.
     * @param m2     an array of double.
     * @param res    an array of double.
     * @return Returns the matrix res.
     */
    public static double[][] scalarMat3Mat3(double scalar, double[][] m1, double[][] m2, double[][] res) {
        res[0][0] = (scalar * m1[0][0]) * m2[0][0] + (scalar * m1[0][1]) * m2[1][0] + (scalar * m1[0][2]) * m2[2][0];
        res[0][1] = (scalar * m1[0][0]) * m2[0][1] + (scalar * m1[0][1]) * m2[1][1] + (scalar * m1[0][2]) * m2[2][1];
        res[0][2] = (scalar * m1[0][0]) * m2[0][2] + (scalar * m1[0][1]) * m2[1][2] + (scalar * m1[0][2]) * m2[2][2];
        res[1][0] = (scalar * m1[1][0]) * m2[0][0] + (scalar * m1[1][1]) * m2[1][0] + (scalar * m1[1][2]) * m2[2][0];
        res[1][1] = (scalar * m1[1][0]) * m2[0][1] + (scalar * m1[1][1]) * m2[1][1] + (scalar * m1[1][2]) * m2[2][1];
        res[1][2] = (scalar * m1[1][0]) * m2[0][2] + (scalar * m1[1][1]) * m2[1][2] + (scalar * m1[1][2]) * m2[2][2];
        res[2][0] = (scalar * m1[2][0]) * m2[0][0] + (scalar * m1[2][1]) * m2[1][0] + (scalar * m1[2][2]) * m2[2][0];
        res[2][1] = (scalar * m1[2][0]) * m2[0][1] + (scalar * m1[2][1]) * m2[1][1] + (scalar * m1[2][2]) * m2[2][1];
        res[2][2] = (scalar * m1[2][0]) * m2[0][2] + (scalar * m1[2][1]) * m2[1][2] + (scalar * m1[2][2]) * m2[2][2];
        return res;
    }

    /**
     * <p>
     * sum</p>
     *
     * @param a an array of double.
     * @param b an array of double.
     * @return Returns the sum array.
     */
    public static double[] sum(double[] a, double[] b) {
        return sum(a, b, new double[3]);
    }

    /**
     * <p>
     * sum</p>
     *
     * @param a   an array of double.
     * @param b   an array of double.
     * @param ret an array of double.
     * @return Returns the array ret.
     */
    public static double[] sum(double[] a, double[] b, double[] ret) {
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
        return ret;
    }

    /**
     * <p>
     * sum</p>
     *
     * @param a an array of float.
     * @param b an array of float.
     * @return Returns the array ret.
     */
    public static float[] sum(float[] a, float[] b) {
        return sum(a, b, new float[3]);
    }

    /**
     * <p>
     * sum</p>
     *
     * @param a   an array of float.
     * @param b   an array of float.
     * @param ret an array of float.
     * @return Returns the array ret.
     */
    public static float[] sum(float[] a, float[] b, float[] ret) {
        ret[0] = a[0] + b[0];
        ret[1] = a[1] + b[1];
        ret[2] = a[2] + b[2];
        return ret;
    }

    /**
     * vector representation of a symmetric 3x3 matrix times a matrix
     *
     * @param v input vector of the form 11, 22, 33, 12, 13, 23.
     * @param m input matrix.
     * @return Returns matrix product.
     */
    public static double[][] symVec6Mat3(double[] v, double[][] m) {
        double[][] res = new double[3][3];
        symVec6Mat3(v, m, res);
        return res;
    }

    /**
     * <p>
     * symVec6mat3</p>
     *
     * @param v   an array of double.
     * @param m   an array of double.
     * @param res an array of double.
     * @return Returns the matrix res.
     */
    public static double[][] symVec6Mat3(double[] v, double[][] m, double[][] res) {
        res[0][0] = v[0] * m[0][0] + v[3] * m[1][0] + v[4] * m[2][0];
        res[0][1] = v[0] * m[0][1] + v[3] * m[1][1] + v[4] * m[2][1];
        res[0][2] = v[0] * m[0][2] + v[3] * m[1][2] + v[4] * m[2][2];
        res[1][0] = v[3] * m[0][0] + v[1] * m[1][0] + v[5] * m[2][0];
        res[1][1] = v[3] * m[0][1] + v[1] * m[1][1] + v[5] * m[2][1];
        res[1][2] = v[3] * m[0][2] + v[1] * m[1][2] + v[5] * m[2][2];
        res[2][0] = v[4] * m[0][0] + v[5] * m[1][0] + v[2] * m[2][0];
        res[2][1] = v[4] * m[0][1] + v[5] * m[1][1] + v[2] * m[2][1];
        res[2][2] = v[4] * m[0][2] + v[5] * m[1][2] + v[2] * m[2][2];
        return res;
    }

    /**
     * <p>
     * transpose3</p>
     *
     * @param m Matrix m.
     * @return Returns the transposed matrix in a new array.
     */
    public static double[][] transpose3(double[][] m) {
        return transpose3(m, new double[3][3]);
    }

    /**
     * <p>
     * transpose3</p>
     *
     * @param m Matrix m.
     * @param t Matrix t.
     * @return Returns the matrix t.
     */
    public static double[][] transpose3(double[][] m, double[][] t) {
        t[0][0] = m[0][0];
        t[0][1] = m[1][0];
        t[0][2] = m[2][0];
        t[1][0] = m[0][1];
        t[1][1] = m[1][1];
        t[1][2] = m[2][1];
        t[2][0] = m[0][2];
        t[2][1] = m[1][2];
        t[2][2] = m[2][2];
        return t;
    }

    /**
     * <p>
     * u2b</p>
     *
     * @param u a double.
     * @return a double.
     */
    public static double u2b(double u) {
        return u * eightPi2;
    }

    /**
     * Vector times a matrix.
     *
     * @param v input vector.
     * @param m input matrix.
     * @return Returns the vector product.
     */
    public static double[] vec3Mat3(double[] v, double[][] m) {
        return vec3Mat3(v, m, new double[3]);
    }

    /**
     * <p>
     * vec3mat3</p>
     *
     * @param v   an array of double.
     * @param m   an array of double.
     * @param res an array of double.
     * @return Returns the array res.
     */
    public static double[] vec3Mat3(double[] v, double[][] m, double[] res) {
        res[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
        res[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
        res[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];
        return res;
    }

    /**
     * <p>
     * logVector.</p>
     *
     * @param v an array of double.
     * @return Returns the a String description of the vector.
     */
    public static String vectorToString(double[] v) {
        StringBuilder sb = new StringBuilder("Vector ( ");
        for (double d : v) {
            sb.append(format("%g ", d));
        }
        sb.append(")");
        return sb.toString();
    }

    /**
     * <p>vectorToString.</p>
     *
     * @param v     an array of {@link double} objects.
     * @param label a {@link java.lang.String} object.
     * @return Returns the a String description of the vector.
     */
    public static String vectorToString(double[] v, String label) {
        if (v == null) {
            return null;
        }
        StringBuilder sb;
        if (label != null) {
            sb = new StringBuilder(format(" %16s = [", label));
        } else {
            sb = new StringBuilder(format(" %16s = [", "v"));
        }

        for (double value : v) {
            sb.append(format(" %16.8f", value));
        }
        sb.append(" ]");
        return sb.toString();
    }

}
