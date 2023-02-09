// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import java.util.logging.Logger;

/**
 * The FloatMath class is a simple math library that operates on 3-coordinate float arrays.
 *
 * <p>All methods are static and thread-safe.
 *
 * <p>Use instances of Float3 for convenience.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class FloatMath {

  private static final Logger logger = Logger.getLogger(FloatMath.class.getName());

  /**
   * Compute a * b + c and return the result in a new array.
   *
   * @param a First vector.
   * @param b Scalar.
   * @param c Second vector.
   * @return Returns a * b + c.
   */
  public static float[] fma(float[] a, float b, float[] c) {
    return fma(a, b, c, new float[3]);
  }

  /**
   * Compute a * b + c.
   *
   * @param a First vector.
   * @param b Scalar.
   * @param c Second vector.
   * @param ret Result vector.
   * @return Returns a * b + c in the vector ret.
   */
  public static float[] fma(float[] a, float b, float[] c, float[] ret) {
    ret[0] = Math.fma(a[0], b, c[0]);
    ret[1] = Math.fma(a[1], b, c[1]);
    ret[2] = Math.fma(a[2], b, c[2]);
    return ret;
  }

  /**
   * Finds the cross-product between two vectors.
   *
   * @param a First vector.
   * @param b Second vector.
   * @return Returns the cross-product.
   */
  public static float[] X(float[] a, float[] b) {
    return X(a, b, new float[3]);
  }

  /**
   * Finds the cross-product between two vectors.
   *
   * @param a First vector.
   * @param b Second vector.
   * @param ret The cross-product of a x b.
   * @return Returns the cross-product ret.
   */
  public static float[] X(float[] a, float[] b, float[] ret) {
    ret[0] = a[1] * b[2] - a[2] * b[1];
    ret[1] = a[2] * b[0] - a[0] * b[2];
    ret[2] = a[0] * b[1] - a[1] * b[0];
    return ret;
  }

  /**
   * sum
   *
   * @param a an array of float.
   * @param b an array of float.
   * @return Returns the array ret.
   */
  public static float[] add(float[] a, float[] b) {
    return add(a, b, new float[3]);
  }

  /**
   * sum
   *
   * @param a an array of float.
   * @param b an array of float.
   * @param ret an array of float.
   * @return Returns the array ret.
   */
  public static float[] add(float[] a, float[] b, float[] ret) {
    ret[0] = a[0] + b[0];
    ret[1] = a[1] + b[1];
    ret[2] = a[2] + b[2];
    return ret;
  }

  /**
   * angle
   *
   * @param i an array of float.
   * @param j an array of float.
   * @return Returns the angle.
   */
  public static float angle(float[] i, float[] j) {
    var x = dot(normalize(i), normalize(j));
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
   * Finds the angle formed by three atoms
   *
   * @param i Atom position vector.
   * @param j Atom position vector (central atom).
   * @param k Atom position vector.
   * @return Returns the angle in the range [ -pi, pi ].
   */
  public static float bondAngle(float[] i, float[] j, float[] k) {
    return angle(sub(i, j), sub(k, j));
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
    var ba = sub(b, a);
    var cb = sub(c, b);
    var dc = sub(d, c);
    var t = X(ba, cb);
    var u = X(cb, dc);
    var rt = dot(t, t);
    var ru = dot(u, u);
    var rtu = (float) sqrt(rt * ru);
    if (rtu != 0.0) {
      var rcb = length(cb);
      var cosine = dot(t, u) / rtu;
      var tu = X(t, u);
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
  public static float dist(float[] a, float[] b) {
    return (float) sqrt(dist2(a, b));
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
    return Math.fma(dx, dx, Math.fma(dy, dy, dz * dz));
  }

  /**
   * Finds the dot product between two vectors.
   *
   * @param a First vector.
   * @param b Second vector.
   * @return Returns the dot product of a and b.
   */
  public static float dot(float[] a, float[] b) {
    return Math.fma(a[0], b[0], Math.fma(a[1], b[1], a[2] * b[2]));
  }

  /**
   * Finds the length of a vector.
   *
   * @param d A vector to find the length of.
   * @return Returns the length of vector d.
   */
  public static float length(float[] d) {
    return (float) sqrt(length2(d));
  }

  /**
   * Finds the length of a vector squared.
   *
   * @param d A vector to find the length of squared.
   * @return Returns the length of vector d squared.
   */
  public static float length2(float[] d) {
    return Math.fma(d[0], d[0], Math.fma(d[1], d[1], d[2] * d[2]));
  }

  /**
   * logVector
   *
   * @param v an array of float.
   */
  public static void log(float[] v) {
    logger.info(toString(v));
  }

  /**
   * logVector.
   *
   * @param v an array of {@link float} objects.
   * @param label a {@link String} object.
   */
  public static void log(float[] v, String label) {
    logger.info(toString(v, label));
  }

  /**
   * Normalizes a vector.
   *
   * @param n A vector to be normalized.
   * @return Returns the normalized vector.
   */
  public static float[] normalize(float[] n) {
    return scale(n, (float) 1.0 / length(n), new float[3]);
  }

  /**
   * Normalizes a vector.
   *
   * @param n A vector to be normalized.
   * @param ret The normalized vector.
   * @return Returns the normalized vector.
   */
  public static float[] normalize(float[] n, float[] ret) {
    return scale(n, (float) 1.0 / length(n), ret);
  }

  /**
   * Scales a vector.
   *
   * @param n A vector to be scaled.
   * @param a A scalar value.
   * @return Returns the scaled vector.
   */
  public static float[] scale(float[] n, float a) {
    return scale(n, a, new float[3]);
  }

  /**
   * Scales a vector.
   *
   * @param n A vector to be scaled.
   * @param a A scalar value.
   * @param ret The scaled Vector.
   * @return Returns the array ret.
   */
  public static float[] scale(float[] n, float a, float[] ret) {
    ret[0] = n[0] * a;
    ret[1] = n[1] * a;
    ret[2] = n[2] * a;
    return ret;
  }


  /**
   * Squares values of a vector.
   *
   * @param n A vector to be squared.
   * @return Returns the squared vector.
   */
  public static float[] square(float[] n) {return square(n, new float[3]);}

  /**
   * Squares values of a vector.
   *
   * @param n A vector to be squared.
   * @param ret The squared vector.
   * @return Returns the array ret.
   */
  public static float[] square(float[] n, float[] ret) {
    ret[0] = n[0] * n[0];
    ret[1] = n[1] * n[1];
    ret[2] = n[2] * n[2];
    return ret;
  }

  /**
   * Square root values of a vector.
   *
   * @param n A vector to determine square root.
   * @return Returns the rooted vector.
   */
  public static float[] squareRoot(float[] n) {return squareRoot(n, new float[3]);}

  /**
   * Square root values of a vector.
   *
   * @param n A vector to determine square root.
   * @param ret The rooted vector.
   * @return Returns the array ret.
   */
  public static float[] squareRoot(float[] n, float[] ret) {
    ret[0] = (float) sqrt(n[0]);
    ret[1] = (float) sqrt(n[1]);
    ret[2] = (float) sqrt(n[2]);
    return ret;
  }

  /**
   * Finds the difference between two vectors.
   *
   * @param a First vector
   * @param b Second vector
   * @return Returns the difference ret.
   */
  public static float[] sub(float[] a, float[] b) {
    return sub(a, b, new float[3]);
  }

  /**
   * Finds the difference between two vectors
   *
   * @param a First vector
   * @param b Second vector
   * @param ret Return Values
   * @return Returns the difference ret.
   */
  public static float[] sub(float[] a, float[] b, float[] ret) {
    ret[0] = a[0] - b[0];
    ret[1] = a[1] - b[1];
    ret[2] = a[2] - b[2];
    return ret;
  }

  /**
   * logVector.
   *
   * @param v an array of float.
   * @return Returns a String description of the vector.
   */
  public static String toString(float[] v) {
    StringBuilder sb = new StringBuilder("Vector ( ");
    for (float d : v) {
      sb.append(format("%g ", d));
    }
    sb.append(")");
    return sb.toString();
  }

  /**
   * vectorToString.
   *
   * @param v an array of {@link float} objects.
   * @param label a {@link String} object.
   * @return Returns a String description of the vector.
   */
  public static String toString(float[] v, String label) {
    if (v == null) {
      return null;
    }
    if (label == null) {
      label = "v";
    }
    StringBuilder sb = new StringBuilder(format(" %16s = [", label));
    for (float value : v) {
      sb.append(format(" %16.8f", value));
    }
    sb.append(" ]");
    return sb.toString();
  }
}
