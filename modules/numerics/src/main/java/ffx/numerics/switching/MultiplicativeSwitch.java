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
package ffx.numerics.switching;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The 6 coefficients of the multiplicative polynomial switch are unique given the distances "a" and
 * "b". They are found by solving a system of 6 equations, which define the boundary conditions of
 * the switch. <br>
 * f(a) = 1 <br>
 * f'(a) = f"(a) = 0 <br>
 * f(b) = f'(b) = f"(b) = 0
 *
 * @author Michael J. Schnieders
 */
public class MultiplicativeSwitch implements UnivariateSwitchingFunction {

  private final double b;
  private final double a;
  private final double c0;
  private final double c1;
  private final double c2;
  private final double c3;
  private final double c4;
  private final double c5;
  private final double twoC2;
  private final double threeC3;
  private final double fourC4;
  private final double fiveC5;

  /**
   * Constructs a MultiplicativeSwitch that starts at f(0)=1 and ends at f(1)=0. The switch smoothly
   * interpolates from 1 to 0 across that range, with zero first and second derivatives at off and
   * cut.
   */
  public MultiplicativeSwitch() {
    this(0.0, 1.0);
  }

  /**
   * Constructs a MultiplicativeSwitch that starts at f(a)=1 and ends at f(b)=0. The switch smoothly
   * interpolates from 1 to 0 across that range, with zero first and second derivatives at off and
   * cut.
   *
   * @param a f(a)=1
   * @param b f(b)=0
   */
  public MultiplicativeSwitch(double a, double b) {

    // f(a) = 1.0
    this.a = a;
    // f(b) = 0.0
    this.b = b;

    double a2 = a * a;
    double b2 = b * b;

    double denom = pow(b - a, 5.0);
    c0 = b * b2 * (b2 - 5.0 * a * b + 10.0 * a2) / denom;
    c1 = -30.0 * a2 * b2 / denom;
    c2 = 30.0 * b * a * (b + a) / denom;
    c3 = -10.0 * (a2 + 4.0 * a * b + b2) / denom;
    c4 = 15.0 * (a + b) / denom;
    c5 = -6.0 / denom;
    twoC2 = 2.0 * c2;
    threeC3 = 3.0 * c3;
    fourC4 = 4.0 * c4;
    fiveC5 = 5.0 * c5;
  }

  /** {@inheritDoc} */
  @Override
  public boolean constantOutsideBounds() {
    return false;
  }

  /**
   * First derivative of the switching function at r.
   *
   * @param r r
   * @param r2 r^2
   * @param r3 r^3
   * @param r4 r^4
   * @return First derivative of switch at r
   */
  public double dtaper(double r, double r2, double r3, double r4) {
    return fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
  }

  /**
   * First derivative of the switching function at r.
   *
   * @param r r
   * @return First derivative of switch at r
   */
  public double dtaper(double r) {
    // Minimize number of multiply operations by storing r^2.
    double r2 = r * r;
    return dtaper(r, r2, r2 * r, r2 * r2);
  }

  /** {@inheritDoc} */
  @Override
  public double firstDerivative(double x) {
    return dtaper(x);
  }

  /** {@inheritDoc} */
  @Override
  public int getHighestOrderZeroDerivative() {
    return 2;
  }

  /** {@inheritDoc} */
  @Override
  public double getOneBound() {
    return a > b ? a : b;
  }

  /**
   * Get the value where the switch starts.
   *
   * @return Switch start.
   */
  public double getSwitchEnd() {
    return b;
  }

  /**
   * Get the value where the switch starts.
   *
   * @return Switch start.
   */
  public double getSwitchStart() {
    return a;
  }

  /** {@inheritDoc} */
  @Override
  public double getZeroBound() {
    return b < a ? b : a;
  }

  /** {@inheritDoc} */
  @Override
  public double nthDerivative(double x, int order) throws IllegalArgumentException {
    if (order < 1) {
      throw new IllegalArgumentException("Order must be >= 1");
    }
    switch (order) {
      case 1:
        return dtaper(x);
      case 2:
        return secondDerivative(x);
      case 3:
        double val = 60.0 * c5 * x * x;
        val += 24.0 * c4 * x;
        val += 6.0 * c3;
        return val;
      case 4:
        val = 120.0 * c5 * x;
        val += 24.0 * c4;
        return val;
      case 5:
        return 120.0 * c5;
      default:
        return 0;
    }
  }

  /** {@inheritDoc} */
  @Override
  public double secondDerivative(double x) {
    double x2 = x * x;
    double val = 20.0 * c5 * x2 * x;
    val += 12.0 * c4 * x2;
    val += 6.0 * c3 * x;
    val += 2.0 * c2;
    return val;
  }

  /** {@inheritDoc} */
  @Override
  public boolean symmetricToUnity() {
    return true;
  }

  /**
   * Value of the switching function at r.
   *
   * @param r r
   * @param r2 r^2
   * @param r3 r^3
   * @param r4 r^4
   * @param r5 r^5
   * @return Value of switch at r
   */
  public double taper(double r, double r2, double r3, double r4, double r5) {
    return c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
  }

  /**
   * Value of the switching function at r.
   *
   * @param r r
   * @return Value of switch at r
   */
  public double taper(double r) {
    // Minimize number of multiply operations by storing r^2, r^3.
    double r2 = r * r;
    double r3 = r2 * r;
    return taper(r, r2, r3, r2 * r2, r3 * r2);
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return format(
        "Multiplicative switch of form f(x) = %8.4g*x^5 + "
            + "%8.4g*x^4 + %8.4g*x^3 + %8.4g*x^2 + %8.4g*x + %8.4g",
        c5, c4, c3, c2, c1, c0);
  }

  /** {@inheritDoc} */
  @Override
  public boolean validOutsideBounds() {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public double valueAt(double x) throws IllegalArgumentException {
    return taper(x);
  }
}
