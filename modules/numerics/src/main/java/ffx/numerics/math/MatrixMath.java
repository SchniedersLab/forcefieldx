// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.numerics.math;

/**
 * The MatrixMath class is a simple matrix math library used mainly by the X-ray package.
 * <p>
 * All methods are thread-safe and static.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class MatrixMath {

    private MatrixMath() {
        // Prevent instantiation.
    }

    /**
     * Calculated the determinant for a 3x3 matrix.
     *
     * @param m input matrix.
     * @return The determinant.
     */
    public static double mat3Determinant(double[][] m) {
        return m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1]
                + m[0][1] * m[1][2] * m[2][0] - m[0][1] * m[1][0] * m[2][2]
                + m[0][2] * m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0];
    }

    /**
     * Returns the determinant of a symmetric 3x3 matrix stored
     * in vector form as 11, 22, 33, 12, 13, 23.
     *
     * @param m symmetric input matrix in vector format.
     * @return The determinant.
     */
    public static double mat3Determinant(double[] m) {
        return m[0] * m[1] * m[2] - m[0] * m[5] * m[5]
                + m[3] * m[5] * m[4] - m[3] * m[3] * m[2]
                + m[4] * m[3] * m[5] - m[4] * m[1] * m[4];
    }

    /**
     * Compute the inverse of 3x3 matrix. The result is returned in a newly
     * allocated matrix.
     *
     * @param m The input 3x3 matrix.
     * @return Returns the inverse of the input matrix m.
     */
    public static double[][] mat3Inverse(double[][] m) {
        return mat3Inverse(m, new double[3][3]);
    }

    /**
     * Compute the inverse of 3x3 matrix. If the input matrix m can be specified
     * for the output matrix to store the result in place.
     *
     * @param m      The input 3x3 matrix.
     * @param output The output 3x3 matrix.
     * @return Returns the inverse of the input matrix m.
     */
    public static double[][] mat3Inverse(double[][] m, double[][] output) {
        double inverseDet = 1.0 / mat3Determinant(m);
        double r00 = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inverseDet;
        double r01 = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inverseDet;
        double r02 = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inverseDet;
        double r10 = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inverseDet;
        double r11 = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inverseDet;
        double r12 = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inverseDet;
        double r20 = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inverseDet;
        double r21 = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inverseDet;
        double r22 = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inverseDet;
        output[0][0] = r00;
        output[0][1] = r01;
        output[0][2] = r02;
        output[1][0] = r10;
        output[1][1] = r11;
        output[1][2] = r12;
        output[2][0] = r20;
        output[2][1] = r21;
        output[2][2] = r22;
        return output;
    }

    /**
     * Multiple a 3x3 matrix m and a 3x3 matrix n. The output is returned in a newly allocated
     * 3x3 matrix.
     *
     * @param m an input 3x3 matrix.
     * @param n an input 3x3 matrix
     * @return Returns the 3x3 matrix result.
     */
    public static double[][] mat3Mat3Multiply(double[][] m, double[][] n) {
        return mat3Mat3Multiply(m, n, new double[3][3]);
    }

    /**
     * Multiple a 3x3 matrix m and a 3x3 matrix n. The output is stored in result. The matrix m or n
     * can be used for the result matrix.
     *
     * @param m      an input 3x3 matrix.
     * @param n      an input 3x3 matrix
     * @param result a 3x3 matrix to store the result.
     * @return Returns the matrix result.
     */
    public static double[][] mat3Mat3Multiply(double[][] m, double[][] n, double[][] result) {
        double r00 = m[0][0] * n[0][0] + m[0][1] * n[1][0] + m[0][2] * n[2][0];
        double r01 = m[0][0] * n[0][1] + m[0][1] * n[1][1] + m[0][2] * n[2][1];
        double r02 = m[0][0] * n[0][2] + m[0][1] * n[1][2] + m[0][2] * n[2][2];
        double r10 = m[1][0] * n[0][0] + m[1][1] * n[1][0] + m[1][2] * n[2][0];
        double r11 = m[1][0] * n[0][1] + m[1][1] * n[1][1] + m[1][2] * n[2][1];
        double r12 = m[1][0] * n[0][2] + m[1][1] * n[1][2] + m[1][2] * n[2][2];
        double r20 = m[2][0] * n[0][0] + m[2][1] * n[1][0] + m[2][2] * n[2][0];
        double r21 = m[2][0] * n[0][1] + m[2][1] * n[1][1] + m[2][2] * n[2][1];
        double r22 = m[2][0] * n[0][2] + m[2][1] * n[1][2] + m[2][2] * n[2][2];
        result[0][0] = r00;
        result[0][1] = r01;
        result[0][2] = r02;
        result[1][0] = r10;
        result[1][1] = r11;
        result[1][2] = r12;
        result[2][0] = r20;
        result[2][1] = r21;
        result[2][2] = r22;
        return result;
    }

    /**
     * Matrix times a matrix (both 4x4).
     *
     * @param m1 First input matrix.
     * @param m2 Second input matrix.
     * @return Returns the matrix product.
     */
    public static double[][] mat4Mat4(double[][] m1, double[][] m2) {
        return mat4Mat4(m1, m2, new double[4][4]);
    }

    /**
     * Multiply two 4x4 matrices.
     *
     * @param m1  an array of double (first matrix).
     * @param m2  an array of double (second matrix).
     * @param res Resultant matrix.
     * @return Returns the matrix res.
     */
    public static double[][] mat4Mat4(double[][] m1, double[][] m2, double[][] res) {
        double r00 = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0] + m1[0][2] * m2[2][0] + m1[0][3] * m2[3][0];
        double r01 = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1] + m1[0][2] * m2[2][1] + m1[0][3] * m2[3][1];
        double r02 = m1[0][0] * m2[0][2] + m1[0][1] * m2[1][2] + m1[0][2] * m2[2][2] + m1[0][3] * m2[3][2];
        double r03 = m1[0][0] * m2[0][3] + m1[0][1] * m2[1][3] + m1[0][2] * m2[2][3] + m1[0][3] * m2[3][3];
        double r10 = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0] + m1[1][2] * m2[2][0] + m1[1][3] * m2[3][0];
        double r11 = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1] + m1[1][2] * m2[2][1] + m1[1][3] * m2[3][1];
        double r12 = m1[1][0] * m2[0][2] + m1[1][1] * m2[1][2] + m1[1][2] * m2[2][2] + m1[1][3] * m2[3][2];
        double r13 = m1[1][0] * m2[0][3] + m1[1][1] * m2[1][3] + m1[1][2] * m2[2][3] + m1[1][3] * m2[3][3];
        double r20 = m1[2][0] * m2[0][0] + m1[2][1] * m2[1][0] + m1[2][2] * m2[2][0] + m1[2][3] * m2[3][0];
        double r21 = m1[2][0] * m2[0][1] + m1[2][1] * m2[1][1] + m1[2][2] * m2[2][1] + m1[2][3] * m2[3][1];
        double r22 = m1[2][0] * m2[0][2] + m1[2][1] * m2[1][2] + m1[2][2] * m2[2][2] + m1[2][3] * m2[3][2];
        double r23 = m1[2][0] * m2[0][3] + m1[2][1] * m2[1][3] + m1[2][2] * m2[2][3] + m1[2][3] * m2[3][3];
        double r30 = m1[3][0] * m2[0][0] + m1[3][1] * m2[1][0] + m1[3][2] * m2[2][0] + m1[3][3] * m2[3][0];
        double r31 = m1[3][0] * m2[0][1] + m1[3][1] * m2[1][1] + m1[3][2] * m2[2][1] + m1[3][3] * m2[3][1];
        double r32 = m1[3][0] * m2[0][2] + m1[3][1] * m2[1][2] + m1[3][2] * m2[2][2] + m1[3][3] * m2[3][2];
        double r33 = m1[3][0] * m2[0][3] + m1[3][1] * m2[1][3] + m1[3][2] * m2[2][3] + m1[3][3] * m2[3][3];
        res[0][0] = r00;
        res[0][1] = r01;
        res[0][2] = r02;
        res[0][3] = r03;
        res[1][0] = r10;
        res[1][1] = r11;
        res[1][2] = r12;
        res[1][3] = r13;
        res[2][0] = r20;
        res[2][1] = r21;
        res[2][2] = r22;
        res[2][3] = r23;
        res[3][0] = r30;
        res[3][1] = r31;
        res[3][2] = r32;
        res[3][3] = r33;
        return res;
    }

    /**
     * Matrix times a vector representation of a symmetric 3x3 matrix
     *
     * @param m input 3x3 matrix.
     * @param v input vector of the form 11, 22, 33, 12, 13, 23
     * @return matrix product in newly allocated space.
     */
    public static double[][] mat3SymVec6(double[][] m, double[] v) {
        return mat3SymVec6(m, v, new double[3][3]);
    }

    /**
     * Matrix times a vector representation of a symmetric 3x3 matrix. The input matrix m
     * can also be used for the output, which overwrites m with the output.
     *
     * @param m      input 3x3 matrix.
     * @param v      input vector of the form 11, 22, 33, 12, 13, 23
     * @param output output 3x3 matrix.
     * @return matrix product in newly allocated space.
     */
    public static double[][] mat3SymVec6(double[][] m, double[] v, double[][] output) {
        double r00 = m[0][0] * v[0] + m[0][1] * v[3] + m[0][2] * v[4];
        double r01 = m[0][0] * v[3] + m[0][1] * v[1] + m[0][2] * v[5];
        double r02 = m[0][0] * v[4] + m[0][1] * v[5] + m[0][2] * v[2];
        double r10 = m[1][0] * v[0] + m[1][1] * v[3] + m[1][2] * v[4];
        double r11 = m[1][0] * v[3] + m[1][1] * v[1] + m[1][2] * v[5];
        double r12 = m[1][0] * v[4] + m[1][1] * v[5] + m[1][2] * v[2];
        double r20 = m[2][0] * v[0] + m[2][1] * v[3] + m[2][2] * v[4];
        double r21 = m[2][0] * v[3] + m[2][1] * v[1] + m[2][2] * v[5];
        double r22 = m[2][0] * v[4] + m[2][1] * v[5] + m[2][2] * v[2];
        output[0][0] = r00;
        output[0][1] = r01;
        output[0][2] = r02;
        output[1][0] = r10;
        output[1][1] = r11;
        output[1][2] = r12;
        output[2][0] = r20;
        output[2][1] = r21;
        output[2][2] = r22;
        return output;
    }

    /**
     * Returns the transpose of a 3x3 Matrix m in newly allocated memory.
     *
     * @param m The input matrix.
     * @return An allocated a 3x3 matrix with transpose of m.
     */
    public static double[][] mat3Transpose(double[][] m) {
        return mat3Transpose(m, new double[3][3]);
    }

    /**
     * Transpose a 3x3 Matrix m and store the result in output. The input matrix m and
     * output matrix output can be the same variable to store the transpose in-place.
     *
     * @param m      The input matrix.
     * @param output The output transposed matrix.
     * @return Returns the transposed matrix output.
     */
    public static double[][] mat3Transpose(double[][] m, double[][] output) {
        double m00 = m[0][0];
        double m01 = m[1][0];
        double m02 = m[2][0];
        double m10 = m[0][1];
        double m11 = m[1][1];
        double m12 = m[2][1];
        double m20 = m[0][2];
        double m21 = m[1][2];
        double m22 = m[2][2];
        output[0][0] = m00;
        output[0][1] = m10;
        output[0][2] = m20;
        output[1][0] = m01;
        output[1][1] = m11;
        output[1][2] = m21;
        output[2][0] = m02;
        output[2][1] = m12;
        output[2][2] = m22;
        return output;
    }

    /**
     * Multiply a 1x3 vector and 3x3 matrix. The output is returned in a newly allocated 1x3 vector.
     *
     * @param v input 1x3 vector.
     * @param m input 3x3 matrix.
     * @return Returns the output vector.
     */
    public static double[] vec3Mat3(double[] v, double[][] m) {
        return vec3Mat3(v, m, new double[3]);
    }

    /**
     * Multiply a 1x3 vector and 3x3 matrix. If the vector v is also used for the output, the
     * output result will overwrite the input values of v.
     *
     * @param v      input 1x3 vector.
     * @param m      input 3x3 matrix.
     * @param output output 1x3 vector.
     * @return Returns the output vector.
     */
    public static double[] vec3Mat3(double[] v, double[][] m, double[] output) {
        double r0 = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
        double r1 = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
        double r2 = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];
        output[0] = r0;
        output[1] = r1;
        output[2] = r2;
        return output;
    }

}
