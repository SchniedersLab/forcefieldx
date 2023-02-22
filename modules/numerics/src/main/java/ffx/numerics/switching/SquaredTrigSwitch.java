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
package ffx.numerics.switching;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;

import java.util.function.DoubleUnaryOperator;

/**
 * A SquaredTrigSwitch implements a 0-1 switch of form f(x) = sin^2(ax) or of form f(x) = cos^2(ax).
 * Cosine implementation is achieved by phase-shifting a sine wave.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class SquaredTrigSwitch implements UnivariateSwitchingFunction {

  private static final double PI_OVER_TWO = Math.PI * 0.5;
  private final double multiplier;
  private final double halfPeriod;
  private final DoubleUnaryOperator xTransform;
  private final boolean cosine;

  /**
   * Default constructor, creating a switch sin^2(pi*x/2), or cos^2 if flag set. The sine form will
   * switch 0-1, and the cosine form switch 1-0.
   *
   * @param cosine Use a cos^2(pi*x/2) transform instead of sin^2(pi*x/2).
   */
  public SquaredTrigSwitch(boolean cosine) {
    this(PI_OVER_TWO, cosine);
  }

  /**
   * Constructor permitting a custom frequency "a" in the form sin^2(a*x) or cos^2(a*x). The sine
   * form will switch 0-1, and the cosine form switch 1-0.
   *
   * @param coefficient Value of a
   * @param cosine Use a cos^2(ax) transform instead of sin^2(ax).
   */
  public SquaredTrigSwitch(double coefficient, boolean cosine) {
    multiplier = coefficient;
    halfPeriod = PI_OVER_TWO / multiplier;
    xTransform = cosine ? this::cosineTransform : this::sineTransform;
    this.cosine = cosine;
  }

  /** {@inheritDoc} */
  @Override
  public boolean constantOutsideBounds() {
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public double firstDerivative(double x) throws IllegalArgumentException {
    x = xTransform.applyAsDouble(x);
    return 2.0 * sin(x) * cos(x) * multiplier;
  }

  /** {@inheritDoc} */
  @Override
  public int getHighestOrderZeroDerivative() {
    return 1;
  }

  /** {@inheritDoc} */
  @Override
  public double getOneBound() {
    return cosine ? 0 : halfPeriod;
  }

  /**
   * Get the repeating period of this switch.
   *
   * @return Period
   */
  public double getPeriod() {
    return 2.0 * halfPeriod;
  }

  /** {@inheritDoc} */
  @Override
  public double getZeroBound() {
    return cosine ? halfPeriod : 0;
  }

  /**
   * Return true if a cos^2(ax) transform, false if a sin^2(ax) transform.
   *
   * @return a boolean.
   */
  public boolean isCosine() {
    return cosine;
  }

  /** {@inheritDoc} */
  @Override
  public double nthDerivative(double x, int order) throws IllegalArgumentException {
    if (order < 1) {
      throw new IllegalArgumentException("Order must be >= 1");
    }
    x = xTransform.applyAsDouble(x);
    double sinVal = sin(x);
    double cosVal = cos(x);
    double multPow = pow(multiplier, order);

    return switch (order % 4) {
      case 0 -> pow(2.0, order - 1) * multPow * ((sinVal * sinVal) - (cosVal * cosVal));
      case 1 -> pow(2.0, order) * multPow * sinVal * cosVal;
      case 2 -> pow(2.0, order - 1) * multPow * ((cosVal * cosVal) - (sinVal * sinVal));
      case 3 -> -1.0 * pow(2.0, order) * multPow * sinVal * cosVal;
      default -> throw new ArithmeticException("A positive number modulo 4 was not 0-3");
    };
  }

  /** {@inheritDoc} */
  @Override
  public double secondDerivative(double x) throws IllegalArgumentException {
    double val = 2.0 * multiplier * multiplier;
    x = xTransform.applyAsDouble(x);
    double cosTerm = cos(x);
    double sinTerm = sin(x);

    return val * ((cosTerm * cosTerm) - (sinTerm * sinTerm));
  }

  /** {@inheritDoc} */
  @Override
  public boolean symmetricToUnity() {
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    if (cosine) {
      return format(
          "Cosine-squared switching function of form f(x) = cos^2(%8.4g * x)", multiplier);
    } else {
      return format("Sine-squared switching function of form f(x) = sin^2(%8.4g * x)", multiplier);
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean validOutsideBounds() {
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public double valueAt(double x) throws IllegalArgumentException {
    x = xTransform.applyAsDouble(x);
    x = sin(x);
    x *= x;
    return x;
  }

  /**
   * Operate on a*x.
   *
   * @param x Input
   * @return Position along sine wave.
   */
  private double sineTransform(double x) {
    return x * multiplier;
  }

  /**
   * Operate on (pi/2) + a*x, period-shifting a sine wave into the desired cosine wave.
   *
   * @param x Input
   * @return Position along cosine wave.
   */
  private double cosineTransform(double x) {
    return PI_OVER_TWO + (x * multiplier);
  }
}
