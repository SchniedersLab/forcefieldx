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
package ffx.numerics.func1d;

import static ffx.numerics.math.ScalarMath.modToRange;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;

/**
 * A QuasiLinearThetaMap implements a map of theta[-pi, +pi] to lambda[0,1] in a mostly-linear
 * fashion (i.e. rectangular sampling of theta produces roughly rectangular sampling of lambda).
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class QuasiLinearThetaMap implements UnivariateDiffFunction {
  private final double theta0;
  private final double r;
  private final double a;
  private final double b;
  private final double c;
  private final double piMinusTheta0;
  private final double negPiPlusTheta0;
  private final double halfR;

  /** Constructs a QuasiLinearThetaMap with a theta0 of 0.1. */
  public QuasiLinearThetaMap() {
    this(0.1);
  }

  /**
   * Constructs a QuasiLinearThetaMap which is roughly V-shaped from [-pi,+pi], is periodic, and
   * uses trigonometric functions to spline between the linear ranges (theta0-pi, -theta0),
   * (+theta0, pi-theta0) and the trigonometric interpolating regions [-pi, theta0-pi],
   * [-theta0,+theta0] and [pi-theta0, pi].
   *
   * @param theta0 Defines the width of the trigonometric interpolating regions.
   */
  public QuasiLinearThetaMap(double theta0) {
    if (theta0 <= 0 || theta0 >= PI) {
      throw new IllegalArgumentException(
          String.format(
              " QuasiLinearThetaMap " + "must receive theta0 from (0 to +pi), received %11.5g",
              theta0));
    }

    double sinT = sin(theta0);
    double cosT = cos(theta0);
    r = 1.0 / (1 - cosT + (0.5 * sinT * (PI - 2 * theta0)));

    this.theta0 = theta0;
    piMinusTheta0 = PI - theta0;
    negPiPlusTheta0 = 0.1 - PI;

    b = r * 0.5 * (1 - cosT - (theta0 * sinT));
    halfR = 0.5 * r;
    a = r * sinT * 0.5;
    double temp = sin(piMinusTheta0 * 0.5);
    c = (r * 0.5 * sinT * piMinusTheta0) + b - (r * (temp * temp));
  }

  @Override
  public double firstDerivative(double x) throws IllegalArgumentException {
    return fd(modToRange(x, -PI, PI));
  }

  @Override
  public double nthDerivative(double x, int order) throws IllegalArgumentException {
    x = modToRange(x, -PI, PI);
    switch (order) {
      case 0:
        return val(x);
      case 1:
        return fd(x);
      case 2:
        return sd(x);
      default:
        return nd(x, order);
    }
  }

  @Override
  public double secondDerivative(double x) throws IllegalArgumentException {
    return sd(modToRange(x, -PI, PI));
  }

  @Override
  public double valueAt(double x) throws IllegalArgumentException {
    return val(modToRange(x, -PI, PI));
  }

  final double[] getConstants() {
    return new double[] {r, a, b, c};
  }

  Branch getBranch(double x) {
    assert x >= -1.0 * PI && x <= PI;
    if (x < negPiPlusTheta0) {
      return Branch.D;
    } else if (x < -1.0 * theta0) {
      return Branch.C;
    } else if (x <= theta0) {
      return Branch.A;
    } else if (x < piMinusTheta0) {
      return Branch.B;
    } else {
      return Branch.D;
    }
  }

  private double val(double x) throws IllegalArgumentException {
    return val(x, getBranch(x));
  }

  double val(double x, Branch branch) throws IllegalArgumentException {
    switch (branch) {
      case A:
        {
          double sinT = sin(x * 0.5);
          return r * sinT * sinT;
        }
      case B:
        return b + a * x;
      case C:
        return b - a * x;
      case D:
        {
          double sinT = sin(x * 0.5);
          return (r * sinT * sinT) + c;
        }
      default:
        throw new IllegalArgumentException("Could not pick a branch! Should be impossible!");
    }
  }

  private double fd(double x) throws IllegalArgumentException {
    return fd(x, getBranch(x));
  }

  double fd(double x, Branch branch) throws IllegalArgumentException {
    switch (branch) {
      case A:
      case D:
        {
          return halfR * sin(x);
        }
      case B:
        {
          return a;
        }
      case C:
        {
          return -1.0 * a;
        }
      default:
        throw new IllegalArgumentException("Could not pick a branch! Should be impossible!");
    }
  }

  private double sd(double x) throws IllegalArgumentException {
    return sd(x, getBranch(x));
  }

  double sd(double x, Branch branch) throws IllegalArgumentException {
    switch (branch) {
      case A:
      case D:
        {
          return halfR * cos(x);
        }
      case B:
      case C:
        {
          return 0;
        }
      default:
        throw new IllegalArgumentException("Could not pick a branch! Should be impossible!");
    }
  }

  private double nd(double x, int order) throws IllegalArgumentException {
    Branch br = getBranch(x);
    if (order < 3) {
      throw new IllegalArgumentException(" Order was " + order + ", must be > 2!");
    }
    switch (br) {
      case B:
      case C:
        return 0;
    }

    switch (order % 4) {
      case 0:
        {
          return halfR * sin(x);
        }
      case 1:
        {
          return halfR * cos(x);
        }
      case 2:
        {
          return -1.0 * halfR * sin(x);
        }
      case 3:
        {
          return -1.0 * halfR * cos(x);
        }
      default:
        {
          throw new ArithmeticException(
              String.format(" Value %d modulo 4 somehow not 0-3!", order));
        }
    }
  }

  enum Branch {
    A,
    B,
    C,
    D
  }
}
