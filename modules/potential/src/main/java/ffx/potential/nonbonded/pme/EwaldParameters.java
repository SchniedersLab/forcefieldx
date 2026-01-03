// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.nonbonded.pme;

import ffx.utilities.FFXProperty;
import ffx.utilities.PropertyGroup;

import static ffx.numerics.special.Erf.erfc;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Mutable Particle Mesh Ewald constants.
 */
public class EwaldParameters {

  /**
   * Default cutoff values for PME under periodic boundary conditions.
   */
  public static final double DEFAULT_EWALD_CUTOFF = 7.0;
  /**
   * The sqrt of PI.
   */
  private static final double SQRT_PI = sqrt(Math.PI);

  /**
   * The default Ewald coefficient.
   * <br>
   * For charged systems, the converged Ewald electrostatic energy is a function of the
   * Ewald coefficient. For this reason, we've chosen to use a default value of 0.545,
   * for all real space Ewald cutoff values.
   * <br>
   * In this way, systemically more accurate values for the real space cutoff,
   * b-spline order and reciprocal space grid will converge the total electrostatic energy.
   */
  public static final double DEFAULT_EWALD_COEFFICIENT = 0.545;

  @FFXProperty(name = "ewald-alpha", propertyGroup = PropertyGroup.ParticleMeshEwald, defaultValue = "0.545",
      description = """
          Sets the value of the Ewald coefficient, which controls the width of the Gaussian screening charges during
          particle mesh Ewald summation for multipole electrostatics. In the absence of the ewald-alpha keyword,
          the default value is 0.545, which is appropriate for most applications.
          """)
  public double aewald;
  public double aewald3;
  public double an0;
  public double an1;
  public double an2;
  public double an3;
  public double an4;
  public double an5;
  public double off;
  public double off2;

  public EwaldParameters(double cutoff, double aewald) {
    setEwaldParameters(cutoff, aewald);
  }

  /**
   * Determine the real space Ewald parameters and permanent multipole self energy.
   *
   * @param off    Real space cutoff.
   * @param aewald Ewald convergence parameter (0.0 turns off reciprocal space).
   */
  public void setEwaldParameters(double off, double aewald) {
    this.off = off;
    this.aewald = aewald;
    off2 = off * off;
    double alsq2 = 2.0 * aewald * aewald;
    double piEwald = Double.POSITIVE_INFINITY;
    if (aewald > 0.0) {
      piEwald = 1.0 / (SQRT_PI * aewald);
    }
    aewald3 = 4.0 / 3.0 * pow(aewald, 3.0) / SQRT_PI;
    if (aewald > 0.0) {
      an0 = alsq2 * piEwald;
      an1 = alsq2 * an0;
      an2 = alsq2 * an1;
      an3 = alsq2 * an2;
      an4 = alsq2 * an3;
      an5 = alsq2 * an4;
    } else {
      an0 = 0.0;
      an1 = 0.0;
      an2 = 0.0;
      an3 = 0.0;
      an4 = 0.0;
      an5 = 0.0;
    }
  }

  /**
   * A precision of 1.0e-8 results in an Ewald coefficient that ensures continuity in the real space
   * gradient, but at the cost of increased amplitudes for high frequency reciprocal space structure
   * factors.
   */
  private double ewaldCoefficient(double cutoff, double precision) {

    double eps = 1.0e-8;
    if (precision < 1.0e-1) {
      eps = precision;
    }

    /*
     * Get an approximate value from cutoff and tolerance.
     */
    double ratio = eps + 1.0;
    double x = 0.5;
    int i = 0;
    // Larger values lead to a more "delta-function-like" Gaussian
    while (ratio >= eps) {
      i++;
      x *= 2.0;
      ratio = erfc(x * cutoff) / cutoff;
    }
    /*
     * Use a binary search to refine the coefficient.
     */
    int k = i + 60;
    double xlo = 0.0;
    double xhi = x;
    for (int j = 0; j < k; j++) {
      x = (xlo + xhi) / 2.0;
      ratio = erfc(x * cutoff) / cutoff;
      if (ratio >= eps) {
        xlo = x;
      } else {
        xhi = x;
      }
    }

    return x;
  }

  /**
   * Determine the Ewald real space cutoff given the Ewald coefficient and a target precision.
   *
   * @param coeff     The Ewald coefficient in use.
   * @param maxCutoff The maximum cutoff.
   * @param eps       The target precision.
   * @return The determined real space Ewald cutoff.
   */
  public static double ewaldCutoff(double coeff, double maxCutoff, double eps) {
    // Set the tolerance value; use of 1.0d-8 requires strict convergence of the real Space sum.
    double ratio = erfc(coeff * maxCutoff) / maxCutoff;

    if (ratio > eps) {
      return maxCutoff;
    }

    // Use a binary search to refine the coefficient.
    double xlo = 0.0;
    double xhi = maxCutoff;
    double cutoff = 0.0;
    for (int j = 0; j < 100; j++) {
      cutoff = (xlo + xhi) / 2.0;
      ratio = erfc(coeff * cutoff) / cutoff;
      if (ratio >= eps) {
        xlo = cutoff;
      } else {
        xhi = cutoff;
      }
    }
    return cutoff;
  }
}
