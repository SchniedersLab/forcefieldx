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
package ffx.crystal;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.random;

/**
 * Enumeration of the 7 lattice systems.
 * <p>
 * Currently, the SpaceGroup class uses the HEXAGONAL_LATTICE in all cases where its also possible to
 * use a RHOMBOHEDRAL_LATTICE.
 * <p>
 * This includes space groups 146, 148, 155, 160, 161, 166 and 167.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public enum LatticeSystem {
  TRICLINIC_LATTICE,
  MONOCLINIC_LATTICE,
  ORTHORHOMBIC_LATTICE,
  TETRAGONAL_LATTICE,
  RHOMBOHEDRAL_LATTICE,
  HEXAGONAL_LATTICE,
  CUBIC_LATTICE;

  /**
   * Tolerance for checking if the lattice system restrictions are satisfied.
   * <p>
   * Set this to 0.0 for strict checking of lattice parameters.
   * <p>
   * For an acetamide crystal minimization, 1.0e-15 was too small a tolerance for equivalent lattice
   * parameters to equate as equal.
   */
  private static final double tolerance = 1.0e-15;

  /**
   * If the two passed values are the same, within the tolerance, return true.
   *
   * @param x1 First value.
   * @param x2 Second value.
   * @return Return true if the two values are the same within specified tolerance.
   */
  public static boolean check(double x1, double x2) {
    return abs(x1 - x2) < tolerance;
  }

  /**
   * Reset lattice parameters for the given lattice systems.
   *
   * @return New unit cell parameters.
   */
  public double[] resetUnitCellParams() {
    double alpha = 60.0 + random() * 60.0;
    double beta = 60.0 + random() * 60.0;
    double gamma = 60.0 + random() * 60.0;
    double[] params = {0.25 + random(), 0.25 + random(), 0.25 + random(), alpha, beta, gamma};
    double ab = 0.5 * (params[0] + params[1]);
    double abc = (params[0] + params[1] + params[2]) / 3.0;
    switch (this) {
      case TRICLINIC_LATTICE:
        break;
      case MONOCLINIC_LATTICE:
        // alpha = gamma = 90
        params[3] = 90.0;
        params[5] = 90.0;
        break;
      case ORTHORHOMBIC_LATTICE:
        // alpha = beta = gamma = 90
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
      case TETRAGONAL_LATTICE:
        // a = b, alpha = beta = gamma = 90
        params[0] = ab;
        params[1] = ab;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
      case RHOMBOHEDRAL_LATTICE:
        // a = b = c, alpha = beta = gamma.
        double angles = (params[3] + params[4] + params[5]) / 3.0;
        params[0] = abc;
        params[1] = abc;
        params[2] = abc;
        params[3] = angles;
        params[4] = angles;
        params[5] = angles;
        break;
      case HEXAGONAL_LATTICE:
        // a = b, alpha = beta = 90, gamma = 120
        params[0] = ab;
        params[1] = ab;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 120.0;
        break;
      case CUBIC_LATTICE:
      default:
        // a = b = c, alpha = beta = gamma = 90
        params[0] = abc;
        params[1] = abc;
        params[2] = abc;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
    }
    return params;
  }

  /**
   * Check that the lattice parameters satisfy the restrictions of the lattice systems.
   *
   * @param a the a-axis length.
   * @param b the b-axis length.
   * @param c the c-axis length.
   * @param alpha the alpha angle.
   * @param beta the beta angle.
   * @param gamma the gamma angle.
   * @return True if the restrictions are satisfied, false otherwise.
   */
  public boolean validParameters(double a, double b, double c, double alpha, double beta,
      double gamma) {
    switch (this) {
      case TRICLINIC_LATTICE -> {
        // No restrictions.
        return true;
      }
      case MONOCLINIC_LATTICE -> {
        // alpha = gamma = 90
        return check(alpha, 90.0) && check(gamma, 90.0);
      }
      case ORTHORHOMBIC_LATTICE -> {
        // alpha = beta = gamma = 90
        return check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0);
      }
      case TETRAGONAL_LATTICE -> {
        // a = b, alpha = beta = gamma = 90
        return check(a, b) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0);
      }
      case RHOMBOHEDRAL_LATTICE -> {
        // a = b = c, alpha = beta = gamma.
        return check(a, b) && check(b, c) && check(alpha, beta) && check(beta, gamma);
      }
      case HEXAGONAL_LATTICE -> {
        // a = b, alpha = beta = 90, gamma = 120
        return check(a, b) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 120.0);
      }
      case CUBIC_LATTICE -> {
        // a = b = c; alpha = beta = gamma = 90
        return check(a, b) && check(b, c) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma,
            90.0);
      }
      default -> {
        assert (2 != 2);
        return false;
      }
    }
  }

  /**
   * Change the lattice parameters to satisfy the restrictions of the lattice system.
   *
   * @param a the proposed a-axis length.
   * @param b the proposed b-axis length.
   * @param c the proposed c-axis length.
   * @param alpha the proposed alpha angle.
   * @param beta the proposed beta angle.
   * @param gamma the proposed gamma angle.
   * @return Adjusted parameters if the restrictions are satisfied, original parameters otherwise.
   */
  public double[] fixParameters(double a, double b, double c, double alpha, double beta,
      double gamma) {
    double[] parameters = {a, b, c, alpha, beta, gamma};
    double ab = (parameters[0] + parameters[1]) / 2;
    double abc = (parameters[0] + parameters[1] + parameters[2]) / 3;
    switch (this) {
      case TRICLINIC_LATTICE -> {
        // No restrictions.
        return parameters;
      }
      case MONOCLINIC_LATTICE -> {
        // alpha = gamma = 90
        parameters[3] = 90.0;
        parameters[5] = 90.0;
        return parameters;
      }
      case ORTHORHOMBIC_LATTICE -> {
        // alpha = beta = gamma = 90
        parameters[3] = 90.0;
        parameters[4] = 90.0;
        parameters[5] = 90.0;
        return parameters;
      }
      case TETRAGONAL_LATTICE -> {
        // a = b, alpha = beta = gamma = 90
        parameters[0] = ab;
        parameters[1] = ab;
        parameters[3] = 90.0;
        parameters[4] = 90.0;
        parameters[5] = 90.0;
        return parameters;
      }
      case RHOMBOHEDRAL_LATTICE -> {
        // a = b = c, alpha = beta = gamma.
        double angles = (parameters[3] + parameters[4] + parameters[5]) / 3;
        parameters[0] = abc;
        parameters[1] = abc;
        parameters[2] = abc;
        parameters[3] = angles;
        parameters[4] = angles;
        parameters[5] = angles;
        return parameters;
      }
      case HEXAGONAL_LATTICE -> {
        // a = b, alpha = beta = 90, gamma = 120
        parameters[0] = ab;
        parameters[1] = ab;
        parameters[3] = 90.0;
        parameters[4] = 90.0;
        parameters[5] = 120.0;
        return parameters;
      }
      case CUBIC_LATTICE -> {
        // a = b = c; alpha = beta = gamma = 90
        parameters[0] = abc;
        parameters[1] = abc;
        parameters[2] = abc;
        parameters[3] = 90.0;
        parameters[4] = 90.0;
        parameters[5] = 90.0;
        return parameters;
      }
      default -> {
        assert (2 != 2);
        return parameters;
      }
    }
  }

  /**
   * Returns the default b-axis for the lattice system.
   *
   * @param aaxis the a-axis length is the best guess for b-axis.
   * @return default b-axis value
   */
  public double getDefaultBAxis(double aaxis) {
    return aaxis;
  }

  /**
   * Returns the default c-axis for the lattice system.
   *
   * @return default c-axis value
   */
  public double getDefaultCAxis(double aaxis, double baxis) {
    return (aaxis + baxis) / 2;
  }

  /**
   * Returns the default alpha for the lattice system.
   *
   * @return default alpha value
   */
  public double getDefaultAlpha() {
    return 90.0;
  }

  /**
   * Returns the default beta for the lattice system.
   *
   * @return default beta value
   */
  public double getDefaultBeta() {
    return 90.0;
  }

  /**
   * Returns the default gamma for the lattice system.
   *
   * @return default gamma value
   */
  public double getDefaultGamma() {
    double gamma = 90.0;
    if (this == HEXAGONAL_LATTICE) {
      gamma = 120.0;
    }
    return gamma;
  }
}
