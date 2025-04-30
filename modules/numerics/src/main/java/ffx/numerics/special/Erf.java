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
package ffx.numerics.special;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Static methods to evaluate erf(x) and erfc(x) for a real argument x. Rational functions are used
 * that approximate erf(x) and erfc(x) to machine precision (approximately 15 decimal digits).
 *
 * <p>The error function erf(x) is defined as:
 * erf(x) = (2/sqrt(pi)) * integral from 0 to x of exp(-t^2) dt
 *
 * <p>The complementary error function erfc(x) is defined as:
 * erfc(x) = 1 - erf(x) = (2/sqrt(pi)) * integral from x to infinity of exp(-t^2) dt
 *
 * <p>This implementation uses different approximation formulas for different ranges of the input:
 * <ul>
 *   <li>|x| &lt;= 0.46875: Uses a rational function approximation.</li>
 *   <li>0.46875 &lt; |x| &lt;= 4.0: Uses a rational function approximation for erfc.</li>
 *   <li>|x| &gt; 4.0: Uses a rational function approximation for erfc with additional scaling.</li>
 * </ul>
 *
 * <p>Adapted from an original program written by W. J. Cody, Mathematics and Computer Science
 * Division, Argonne National Laboratory, Argonne, IL 60439
 *
 * <p>Reference: W. J. Cody, "Rational Chebyshev Approximations for the Error Function,"
 * Mathematics of Computation, Vol. 23, No. 107 (Jul., 1969), pp. 631-637.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Erf {

  // Mathematical and machine-dependent constants.

  /**
   * The reciprocal of the square root of PI, used in the error function calculations.
   */
  private static final double sqrtPI = 1.0 / sqrt(PI);

  /**
   * One-sixteenth (1/16), used in the computation of exp(-y*y) for large y.
   */
  private static final double oneSixteenth = 1.0 / 16.0;

  /**
   * The threshold value that determines which approximation formula to use.
   * For |x| <= thresh, a direct approximation of erf(x) is used.
   * For |x| > thresh, the complementary error function erfc(x) is computed first.
   */
  private static final double thresh = 0.46875;

  /**
   * The xSmall argument below which erf(x) may be represented by 2*x/sqrt(pi) and above which x*x won't underflow.
   * <p>
   * A conservative value is the largest machine number X such that 1.0 + X = 1.0 to machine precision.
   * This is approximately the square root of the double precision machine epsilon.
   */
  private static final double xSmall = 1.11e-16;

  /**
   * xBig is the largest argument acceptable for erfc.
   * <p>
   * Solution to the equation:
   * <p>
   * W(x) * (1-0.5/x**2) = xMin, where
   * <p>
   * W(x) = exp(-x*x)/[x*sqrt(pi)]
   * <p>
   * For x > xBig, erfc(x) is effectively 0 in double precision.
   */
  private static final double xBig = 26.543;

  private Erf() {
  }

  /**
   * Evaluates erf(x) for a real argument x.
   * <p>
   * The error function erf(x) is defined as:
   * erf(x) = (2/sqrt(pi)) * integral from 0 to x of exp(-t^2) dt
   * <p>
   * Special cases:
   * <ul>
   *   <li>If arg is NaN, then the result is NaN.</li>
   *   <li>If arg is +infinity, then the result is 1.0.</li>
   *   <li>If arg is -infinity, then the result is -1.0.</li>
   *   <li>If arg is 0, then the result is 0.</li>
   * </ul>
   *
   * @param arg the value to evaluate erf at.
   * @return erf of the argument.
   * @since 1.0
   */
  public static double erf(final double arg) {
    // Handle special cases
    if (Double.isNaN(arg)) {
      return Double.NaN;
    }
    if (Double.isInfinite(arg)) {
      return arg > 0 ? 1.0 : -1.0;
    }

    return erfCore(arg, false);
  }

  /**
   * Evaluate erfc(x) for a real argument x.
   * <p>
   * The complementary error function erfc(x) is defined as:
   * erfc(x) = 1 - erf(x) = (2/sqrt(pi)) * integral from x to infinity of exp(-t^2) dt
   * <p>
   * Special cases:
   * <ul>
   *   <li>If arg is NaN, then the result is NaN.</li>
   *   <li>If arg is +infinity, then the result is 0.0.</li>
   *   <li>If arg is -infinity, then the result is 2.0.</li>
   *   <li>If arg is 0, then the result is 1.0.</li>
   * </ul>
   *
   * @param arg the value to evaluate erfc at.
   * @return erfc of the argument.
   * @since 1.0
   */
  public static double erfc(final double arg) {
    // Handle special cases
    if (Double.isNaN(arg)) {
      return Double.NaN;
    }
    if (Double.isInfinite(arg)) {
      return arg > 0 ? 0.0 : 2.0;
    }

    return erfCore(arg, true);
  }

  /**
   * Evaluates erf(x) or erfc(x) for a real argument x. When called with mode = false, erf is
   * returned, while with mode = true, erfc is returned.
   * <p>
   * This method implements the core algorithm for both erf and erfc functions using rational
   * function approximations for different ranges of the input value. The implementation is based
   * on the algorithm developed by W. J. Cody.
   * <p>
   * The algorithm uses three different approximation formulas:
   * <ul>
   *   <li>For |x| <= 0.46875: Direct rational approximation of erf(x)</li>
   *   <li>For 0.46875 < |x| <= 4.0: Rational approximation of erfc(x)</li>
   *   <li>For |x| > 4.0: Continued fraction approximation of erfc(x)</li>
   * </ul>
   *
   * @param x    the value to evaluate erf or erfc at.
   * @param mode if mode is true, evaluate erfc, otherwise evaluate erf.
   * @return if (!mode) erf(arg), else erfc(arg)
   * @since 1.0
   */
  private static double erfCore(final double x, final boolean mode) {

    // Store the argument and its absolute value.
    final double y = abs(x);
    double result = 0.0;

    // Evaluate error function for |x| less than 0.46875.
    if (y <= thresh) {
      // For very small values, avoid underflow in y^2 calculation
      double ysq = 0.0;
      if (y > xSmall) {
        ysq = y * y;
      }

      // Calculate rational approximation for erf(x)
      // P(y^2)/Q(y^2) where P and Q are polynomials
      // These coefficients are from Cody's Chebyshev approximation
      double xNum = 1.85777706184603153e-1 * ysq;
      double xDen = ysq;

      // Build up the polynomials term by term
      xNum = (xNum + 3.16112374387056560e0) * ysq;
      xDen = (xDen + 2.36012909523441209e1) * ysq;

      xNum = (xNum + 1.13864154151050156e2) * ysq;
      xDen = (xDen + 2.44024637934444173e2) * ysq;

      xNum = (xNum + 3.77485237685302021e2) * ysq;
      xDen = (xDen + 1.28261652607737228e3) * ysq;

      // Final term and calculation
      // Multiply by x to get erf(x)
      result = x * (xNum + 3.20937758913846947e3) / (xDen + 2.84423683343917062e3);

      // If calculating erfc, return 1 - erf(x)
      if (mode) {
        result = 1.0 - result;
      }
    } else if (y <= 4.0) {
      // Get complementary error function for 0.46875 <= |x| <= 4.0.
      // For this range, we use a rational approximation for erfc(y)
      // These coefficients are from Cody's Chebyshev approximation

      // Calculate rational approximation for erfc(y)
      // R(y) = P(y)/Q(y) where P and Q are polynomials in y
      double xNum = 2.15311535474403846e-8 * y;
      double xDen = y;

      // Build up the polynomials term by term
      xNum = (xNum + 5.64188496988670089e-1) * y;
      xDen = (xDen + 1.57449261107098347e1) * y;

      xNum = (xNum + 8.88314979438837594e0) * y;
      xDen = (xDen + 1.17693950891312499e2) * y;

      xNum = (xNum + 6.61191906371416295e1) * y;
      xDen = (xDen + 5.37181101862009858e2) * y;

      xNum = (xNum + 2.98635138197400131e2) * y;
      xDen = (xDen + 1.62138957456669019e3) * y;

      xNum = (xNum + 8.81952221241769090e2) * y;
      xDen = (xDen + 3.29079923573345963e3) * y;

      xNum = (xNum + 1.71204761263407058e3) * y;
      xDen = (xDen + 4.36261909014324716e3) * y;

      xNum = (xNum + 2.05107837782607147e3) * y;
      xDen = (xDen + 3.43936767414372164e3) * y;

      // Final term and calculation
      result = (xNum + 1.23033935479799725e3) / (xDen + 1.23033935480374942e3);

      // Compute exp(-y*y) efficiently for large y
      // Break y into integer and fractional parts to avoid overflow
      double ysq = floor(16.0 * y) * oneSixteenth;  // Integer part divided by 16
      double del = (y - ysq) * (y + ysq);              // Compute remainder using difference of squares

      // Multiply by exp(-y*y) to get erfc(y)
      result = exp(-ysq * ysq - del) * result;
      if (!mode) {
        result = 1.0 - result;
        if (x < 0.0) {
          result = -result;
        }
      } else if (x < 0.0) {
        result = 2.0 - result;
      }
    } else {
      // Get complementary error function for |x| greater than 4.0.
      // For large |x|, we use a continued fraction approximation for erfc
      if (y < xBig) {
        // For very large y, erfc(y) approaches 0, but we need to compute it accurately
        // Use a rational approximation in 1/y and 1/y^2

        double iy = 1.0 / y;  // Inverse of y
        double ysq = iy * iy; // 1/y^2

        // Calculate rational approximation for erfc(y) * exp(y^2) * y * sqrt(pi)
        // These coefficients are from Cody's Chebyshev approximation
        double xNum = 1.63153871373020978e-2 * ysq;
        double xDen = ysq;

        // Build up the polynomials term by term
        xNum = (xNum + 3.05326634961232344e-1) * ysq;
        xDen = (xDen + 2.56852019228982242e0) * ysq;

        xNum = (xNum + 3.60344899949804439e-1) * ysq;
        xDen = (xDen + 1.87295284992346047e0) * ysq;

        xNum = (xNum + 1.25781726111229246e-1) * ysq;
        xDen = (xDen + 5.27905102951428412e-1) * ysq;

        xNum = (xNum + 1.60837851487422766e-2) * ysq;
        xDen = (xDen + 6.05183413124413191e-2) * ysq;

        // Final term and calculation
        result = ysq * (xNum + 6.58749161529837803e-4) / (xDen + 2.33520497626869185e-3);

        // Complete the calculation of erfc(y)
        result = (sqrtPI - result) * iy;

        // Compute exp(-y*y) efficiently for large y
        // Break y into integer and fractional parts to avoid overflow
        ysq = floor(16.0 * y) * oneSixteenth;  // Integer part divided by 16
        double del = (y - ysq) * (y + ysq);       // Compute remainder using difference of squares

        // Multiply by exp(-y*y) to get erfc(y)
        result = exp(-ysq * ysq - del) * result;
      }
      // For y >= xBig, erfc(y) is effectively 0 in double precision
      // result is already initialized to 0.0
      // Convert erfc result to erf result if needed
      if (!mode) {
        // For erf, use the relationship erf(x) = 1 - erfc(x)
        result = 1.0 - result;

        // erf is an odd function: erf(-x) = -erf(x)
        if (x < 0.0) {
          result = -result;
        }
      } else {
        // erfc is neither odd nor even: erfc(-x) = 2 - erfc(x)
        if (x < 0.0) {
          result = 2.0 - result;
        }
      }
    }
    return result;
  }
}
