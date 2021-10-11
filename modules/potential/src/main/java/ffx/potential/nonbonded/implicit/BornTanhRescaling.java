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
package ffx.potential.nonbonded.implicit;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.tanh;

/**
 * Rescale the Born radius integral to account for interstitial spaces.
 * <p>
 * Ri^-1 = [rhoi^-3 - (4*pi/3(rhoi^-3 - 50^-3)) *
 *         tanh(beta0*Psi*rhoi^3 - beta1*(Psi*rhoi^3)^2 + beta2*(Psi*rhoi^3)^3)]^(1/3)
 * <p>
 * Citations:
 *   Aguilar, B.; Shadrach, R.; Onufriev, A. V. Reducing the
 *   secondary structure bias in the generalized Born model via R6 effective
 *   radii. J. Chem. Theory Comput. 2010, 6, 3613−3630.
 *
 *   Onufriev, A.; Bashford, D.; Case, D. Exploring protein native
 *   states and large-scale conformational changes with a modified
 *   generalized born model. Proteins 2004, 55, 383−394.
 *
 * @author Rae A. Corrigan
 * @since 1.0
 */
public class BornTanhRescaling {

  /**
   * Maximum Born radius.
   */
  public static final double MAX_BORN_RADIUS = 50.0;
  /**
   * 1/50^3 where 50 Angstroms is the maximum Born radius
   */
  private static final double RECIP_MAX_BORN_RADIUS3 = pow(MAX_BORN_RADIUS, -3.0);
  private static final double PI4_3 = 4.0 / 3.0 * PI;
  /**
   * Tanh coefficients from Corrigan et al.
   */
  private static double beta0 = 0.4694;
  private static double beta1 = 0.0391;
  private static double beta2 = 0.0008;

  /**
   * Rescale the Born radius integral to account for interstitial spaces.
   *
   * @param Ii The total integral of 1/r^6 over vdW spheres and pairwise neck integrals.
   * @param rhoi The base radius of for the atom being descreened.
   * @return The rescaled integral, which is greater than or equal to the input integral.
   */
  public static double tanhRescaling(double Ii, double rhoi) {
    // Set up tanh function components
    double rhoi3 = rhoi * rhoi * rhoi;
    double rhoi3Psi = rhoi3 * Ii;
    double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
    double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
    // If the output of the tanh function is 1.0, then the Born radius will be MaxBornRadius
    double tanh_constant = PI4_3 * ((1.0 / rhoi3) - RECIP_MAX_BORN_RADIUS3);
    return tanh_constant * tanh(beta0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3);
  }

  /**
   * The chain rule derivative for rescaling the Born radius integral to account for interstitial
   * spaces.
   *
   * @param Ii The total integral of 1/r^6 over vdW spheres and pairwise neck integrals.
   * @param rhoi The base radius of for the atom being descreened.
   * @return The chain rule derivative.
   */
  public static double tanhRescalingChainRule(double Ii, double rhoi) {
    double rhoi3 = rhoi * rhoi * rhoi;
    double rhoi3Psi = rhoi3 * Ii;
    double rhoi6Psi2 = rhoi3Psi * rhoi3Psi;
    double rhoi9Psi3 = rhoi6Psi2 * rhoi3Psi;
    double rhoi6Psi = rhoi3 * rhoi3 * Ii;
    double rhoi9Psi2 = rhoi6Psi2 * rhoi3;
    double tanhTerm = tanh(beta0 * rhoi3Psi - beta1 * rhoi6Psi2 + beta2 * rhoi9Psi3);
    double tanh2 = tanhTerm * tanhTerm;
    double chainRuleTerm = beta0 * rhoi3 - 2.0 * beta1 * rhoi6Psi + 3.0 * beta2 * rhoi9Psi2;
    double tanh_constant = PI4_3 * ((1.0 / rhoi3) - RECIP_MAX_BORN_RADIUS3);
    return tanh_constant * chainRuleTerm * (1.0 - tanh2);
  }

  // Getters and setters for optimizing constants
  public static double getBeta0() {
    return beta0;
  }

  public static double getBeta1() {
    return beta1;
  }

  public static double getBeta2() {
    return beta2;
  }

  public static void setBeta0(double beta0) {
    BornTanhRescaling.beta0 = beta0;
  }

  public static void setBeta1(double beta1) {
    BornTanhRescaling.beta1 = beta1;
  }

  public static void setBeta2(double beta2) {
    BornTanhRescaling.beta2 = beta2;
  }

}
