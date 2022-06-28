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
package ffx.numerics.multipole;

import static ffx.numerics.special.Erf.erfc;
import static java.lang.Math.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The EwaldTensorQI class computes derivatives of erfc(<b>r</b>)/|<b>r</b>| via recursion to
 * arbitrary order for Cartesian multipoles in a quasi-internal frame.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 *     Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 *     computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 *     Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @since 1.0
 */
public class EwaldTensorQI extends CoulombTensorQI {

  /** Constant <code>sqrtPI = sqrt(PI)</code> */
  private static final double sqrtPI = sqrt(PI);

  /**
   * These are the "source" terms for the recursion for the screened Coulomb operator erfc(R)/R.
   */
  private final double[] ewaldSource;

  /**
   * Ewald convergence parameter.
   */
  private final double beta;

  /**
   * Constructor for EwaldTensorQI.
   *
   * @param order Tensor order.
   * @param beta Ewald covergence parameter.
   */
  public EwaldTensorQI(int order, double beta) {
    super(order);
    this.beta = beta;
    operator = OPERATOR.SCREENED_COULOMB;

    // Auxillary terms for screened Coulomb (Sagui et al. Eq. 2.28)
    ewaldSource = new double[o1];
    double prefactor = 2.0 * beta / sqrtPI;
    double twoBeta2 = -2.0 * beta * beta;
    for (int n = 0; n <= order; n++) {
      ewaldSource[n] = prefactor * pow(twoBeta2, n);
    }

  }

  /**
   * Generate source terms for the Ewald Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  protected void source(double[] T000) {
    // Generate source terms for real space Ewald summation.
    if (beta > 0.0) {
      // Sagui et al. Eq. 2.22
      double betaR = beta * R;
      double betaR2 = betaR * betaR;
      double iBetaR2 = 1.0 / (2.0 * betaR2);
      double expBR2 = exp(-betaR2);
      // Fnc(x^2) = Sqrt(PI) * erfc(x) / (2*x)
      // where x = Beta*R
      double Fnc = sqrtPI * erfc(betaR) / (2.0 * betaR);
      for (int n = 0; n < o1; n++) {
        T000[n] = ewaldSource[n] * Fnc;
        // Generate F(n+1)c from Fnc (Eq. 2.24 in Sagui et al.)
        // F(n+1)c = [(2*n+1) Fnc(x) + exp(-x)] / 2x
        // where x = (Beta*R)^2
        Fnc = ((2.0 * n + 1.0) * Fnc + expBR2) * iBetaR2;
      }
    } else {
      // For beta = 0, generate tensors for the Coulomb operator.
      super.source(T000);
    }
  }

}
