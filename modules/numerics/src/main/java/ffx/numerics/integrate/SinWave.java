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
package ffx.numerics.integrate;

import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;

/**
 * A SinWave describes points along a sine wave of f(x) = a*sin(jx).
 *
 * @author Jacob M. Litman
 */
public class SinWave extends FunctionDataCurve {

  /** Magnitude. */
  private final double a;
  /** Periodicity. */
  private final double n;
  /** Inverse periodicity. */
  private final double nInverse;

  /**
   * Constructs f(x) = a*sin(nx).
   *
   * @param x an array of {@link double} objects.
   * @param a magnitude.
   * @param n periodicity.
   */
  public SinWave(double[] x, double a, double n) {
    this(x, false, a, n);
  }

  /**
   * Constructs f(x) = a*sin(nx).
   *
   * @param x an array of {@link double} objects.
   * @param halfWidthEnds Use half-width start and end bins.
   * @param a magnitude.
   * @param n periodicity.
   */
  public SinWave(double[] x, boolean halfWidthEnds, double a, double n) {
    int npoints = x.length;
    points = new double[npoints];
    this.a = a;
    this.n = n;
    nInverse = 1.0 / n;
    this.halfWidthEnd = halfWidthEnds;

    for (int i = 0; i < points.length; i++) {
      points[i] = sinAt(x[i]);
    }
    lb = x[0];
    ub = x[npoints - 1];
    assertXIntegrity(x);
    this.x = new double[x.length];
    arraycopy(x, 0, this.x, 0, x.length);
  }

  /** {@inheritDoc} */
  @Override
  public double fX(double x) {
    return sinAt(x);
  }

  /** {@inheritDoc} */
  @Override
  public double integralAt(double x) {
    return -1 * a * nInverse * cos(n * x);
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder("Sine wave f(x) = ");
    sb.append(a).append("*sin(").append(n).append("x)");
    sb.append(
        format(
            " with %d points from lower bound %9.3g and upper bound %9.3g", points.length, lb, ub));
    if (halfWidthEnd) {
      sb.append(" and half-width start/end bins");
    }

    return sb.toString();
  }

  // Private, non-overrideable method for use in the constructor.
  private double sinAt(double x) {
    return a * sin(n * x);
  }
}
