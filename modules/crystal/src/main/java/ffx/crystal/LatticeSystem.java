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
package ffx.crystal;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.random;

/**
 * Enumeration of the 7 lattice systems.
 * <p>
 * Currently the SpaceGroup class uses the HEXAGONAL_LATTICE in all cases where its also possible to
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
    double[] params = {0.1 + random(), 0.1 + random(), 0.1 + random(), alpha, beta, gamma};
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
        double ab = 0.5 * (params[0] + params[1]);
        params[0] = ab;
        params[1] = ab;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 90.0;
        break;
      case RHOMBOHEDRAL_LATTICE:
        // a = b = c, alpha = beta = gamma.
        double abc = (params[0] + params[1] + params[2]) / 3.0;
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
        ab = 0.5 * (params[0] + params[1]);
        params[0] = ab;
        params[1] = ab;
        params[3] = 90.0;
        params[4] = 90.0;
        params[5] = 120.0;
        break;
      case CUBIC_LATTICE:
      default:
        // a = b = c, alpha = beta = gamma = 90
        abc = (params[0] + params[1] + params[2]) / 3.0;
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
      case TRICLINIC_LATTICE:
        // No restrictions.
        return true;
      case MONOCLINIC_LATTICE:
        // alpha = gamma = 90
        return check(alpha, 90.0) && check(gamma, 90.0);
      case ORTHORHOMBIC_LATTICE:
        // alpha = beta = gamma = 90
        return check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0);
      case TETRAGONAL_LATTICE:
        // a = b, alpha = beta = gamma = 90
        return check(a, b) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0);
      case RHOMBOHEDRAL_LATTICE:
        // a = b = c, alpha = beta = gamma.
        return check(a, b) && check(b, c) && check(alpha, beta) && check(beta, gamma);
      case HEXAGONAL_LATTICE:
        // a = b, alpha = beta = 90, gamma = 120
        return check(a, b) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 120.0);
      case CUBIC_LATTICE:
        // a = b = c; alpha = beta = gamma = 90
        return check(a, b) && check(b, c) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma,
            90.0);
      default:
        assert (2 != 2);
        return false;
    }
  }
}
