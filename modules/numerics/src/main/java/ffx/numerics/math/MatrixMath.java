// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
 *
 * <p>All methods are thread-safe and static.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class MatrixMath {

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
   * determinant3
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
   * inverse of a 3x3 matrix.
   *
   * @param m input matrix.
   * @return Returns the matrix inverse.
   */
  public static double[][] mat3Inverse(double[][] m) {
    return mat3Inverse(m, new double[3][3]);
  }

  /**
   * mat3inverse
   *
   * @param m an array of double.
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
   * mat3mat3
   *
   * @param m1 an array of double.
   * @param m2 an array of double.
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
   * mat3SymVec6
   *
   * @param m an array of double.
   * @param v an array of double.
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
   * mat3vec3
   *
   * @param v an array of double.
   * @param m an array of double.
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
   * scalar times a matrix times a matrix.
   *
   * @param scalar Input scalar.
   * @param m1 First input matrix.
   * @param m2 Second input matrix.
   * @return Returns the matrix product.
   */
  public static double[][] scalarMat3Mat3(double scalar, double[][] m1, double[][] m2) {
    return scalarMat3Mat3(scalar, m1, m2, new double[3][3]);
  }

  /**
   * scalarMat3mat3
   *
   * @param scalar a double.
   * @param m1 an array of double.
   * @param m2 an array of double.
   * @param res an array of double.
   * @return Returns the matrix res.
   */
  public static double[][] scalarMat3Mat3(
      double scalar, double[][] m1, double[][] m2, double[][] res) {
    res[0][0] =
        (scalar * m1[0][0]) * m2[0][0]
            + (scalar * m1[0][1]) * m2[1][0]
            + (scalar * m1[0][2]) * m2[2][0];
    res[0][1] =
        (scalar * m1[0][0]) * m2[0][1]
            + (scalar * m1[0][1]) * m2[1][1]
            + (scalar * m1[0][2]) * m2[2][1];
    res[0][2] =
        (scalar * m1[0][0]) * m2[0][2]
            + (scalar * m1[0][1]) * m2[1][2]
            + (scalar * m1[0][2]) * m2[2][2];
    res[1][0] =
        (scalar * m1[1][0]) * m2[0][0]
            + (scalar * m1[1][1]) * m2[1][0]
            + (scalar * m1[1][2]) * m2[2][0];
    res[1][1] =
        (scalar * m1[1][0]) * m2[0][1]
            + (scalar * m1[1][1]) * m2[1][1]
            + (scalar * m1[1][2]) * m2[2][1];
    res[1][2] =
        (scalar * m1[1][0]) * m2[0][2]
            + (scalar * m1[1][1]) * m2[1][2]
            + (scalar * m1[1][2]) * m2[2][2];
    res[2][0] =
        (scalar * m1[2][0]) * m2[0][0]
            + (scalar * m1[2][1]) * m2[1][0]
            + (scalar * m1[2][2]) * m2[2][0];
    res[2][1] =
        (scalar * m1[2][0]) * m2[0][1]
            + (scalar * m1[2][1]) * m2[1][1]
            + (scalar * m1[2][2]) * m2[2][1];
    res[2][2] =
        (scalar * m1[2][0]) * m2[0][2]
            + (scalar * m1[2][1]) * m2[1][2]
            + (scalar * m1[2][2]) * m2[2][2];
    return res;
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
   * symVec6mat3
   *
   * @param v an array of double.
   * @param m an array of double.
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
   * transpose3
   *
   * @param m Matrix m.
   * @return Returns the transposed matrix in a new array.
   */
  public static double[][] transpose3(double[][] m) {
    return transpose3(m, new double[3][3]);
  }

  /**
   * transpose3
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
   * vec3mat3
   *
   * @param v an array of double.
   * @param m an array of double.
   * @param res an array of double.
   * @return Returns the array res.
   */
  public static double[] vec3Mat3(double[] v, double[][] m, double[] res) {
    res[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
    res[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
    res[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];
    return res;
  }
}
