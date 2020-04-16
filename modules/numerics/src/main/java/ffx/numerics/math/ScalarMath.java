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
package ffx.numerics.math;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * The ScalarMath class is a simple math library that operates on single variables
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ScalarMath {

  private static final double eightPi2 = 8.0 * PI * PI;

  /**
   * b2u
   *
   * @param b a double.
   * @return a double.
   */
  public static double b2u(double b) {
    return b / eightPi2;
  }

  /**
   * binomial
   *
   * @param n a long.
   * @param k a long.
   * @return Returns the binomial of (n,k).
   */
  public static long binomial(long n, long k) {
    return factorial(n) / (factorial(n - k) * factorial(k));
  }

  /**
   * Returns n!! Precondition: n .GE. -1 Returning 1 for -1 input is analogous to Maple behavior.
   *
   * @param n long.
   * @return Returns the n!!.
   */
  public static long doubleFactorial(long n) {
    if (n < -1) {
      throw new RuntimeException("Underflow error in doubleFactorial");
    } else if (n == -1 || n == 0 || n == 1) {
      return 1;
    } else {
      return n * doubleFactorial(n - 2);
    }
  }

  /**
   * Returns n! <br>
   * Precondition: n .GE. 0 and n .LE. 20 <br>
   * Max long = 9223372036854775807 <br>
   * 20! = 2432902008176640000 is ok. <br>
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
   * Compute 1.0 / (1.0 + exp(x)).
   *
   * @param x Input.
   * @return Returns 1.0 / (1.0 + exp(x)).
   */
  public static double fermiFunction(double x) {
    return 1.0 / (1.0 + exp(x));
  }

  /**
   * This is an atypical mod function used by crystallography methods.
   *
   * <p>mod
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
   * Atypical mod function used to move a value into the range lb &lt;= value &lt; ub, assuming the
   * domain is periodic with a period of (ub - lb).
   *
   * @param value Value to move between bounds.
   * @param lb Lower bound.
   * @param ub Upper bound.
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
   * u2b
   *
   * @param u a double.
   * @return a double.
   */
  public static double u2b(double u) {
    return u * eightPi2;
  }
}
