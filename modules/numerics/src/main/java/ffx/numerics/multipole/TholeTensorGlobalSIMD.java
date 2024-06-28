// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;

/**
 * The TholeTensorGlobal class computes derivatives of Thole damping via recursion to order &lt;= 4 for
 * Cartesian multipoles in either a global frame.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 * Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 * computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 * Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @since 1.0
 */
public class TholeTensorGlobalSIMD extends CoulombTensorGlobalSIMD {

  /**
   * Constant <code>threeFifths=3.0 / 5.0</code>
   */
  private static final double threeFifths = 3.0 / 5.0;

  /**
   * Constant <code>oneThirtyFifth=1.0 / 35.0</code>
   */
  private static final double oneThirtyFifth = 1.0 / 35.0;

  /**
   * Thole damping parameter is set to min(pti,ptk)).
   */
  private DoubleVector thole;

  /**
   * AiAk parameter = 1/(alphaI^6*alphaK^6) where alpha is polarizability.
   */
  private DoubleVector AiAk;

  /**
   * Constructor for EwaldMultipoleTensorGlobal.
   *
   * @param order Tensor order.
   * @param thole Thole damping parameter is set to min(pti,ptk)).
   * @param AiAk  parameter = 1/(alphaI^6*alphaK^6) where alpha is polarizability.
   */
  public TholeTensorGlobalSIMD(int order, DoubleVector thole, DoubleVector AiAk) {
    super(order);
    this.thole = thole;
    this.AiAk = AiAk;
    this.operator = Operator.THOLE_FIELD;

    // Source terms are currently defined up to order 4.
    assert (order <= 4);
  }

  /**
   * Set Thole damping parameters
   *
   * @param thole a double.
   * @param AiAk  a double.
   */
  public void setThole(DoubleVector thole, DoubleVector AiAk) {
    this.thole = thole;
    this.AiAk = AiAk;
  }

  /**
   * Check if the Thole damping is exponential is greater than zero (or the interaction can be
   * neglected).
   *
   * @param r The separation distance.
   * @return True if -thole*u^3 is greater than -50.0.
   */
  public boolean checkThole(DoubleVector r) {
    return checkThole(thole, AiAk, r);
  }

  /**
   * Check if the Thole damping is exponential is greater than zero (or the interaction can be neglected).
   *
   * @param thole Thole damping parameter is set to min(pti,ptk)).
   * @param AiAk  parameter = 1/(alphaI^6*alphaK^6) where alpha is polarizability.
   * @param r     The separation distance.
   * @return True if -thole*u^3 is greater than -50.0.
   */
  protected static boolean checkThole(DoubleVector thole, DoubleVector AiAk, DoubleVector r) {
    DoubleVector rAiAk = r.mul(AiAk);
    VectorMask<Double> check = thole.mul(rAiAk).mul(rAiAk).mul(rAiAk).lt(50);
    return check.anyTrue();
  }

  /**
   * Generate source terms for the Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  @Override
  protected void source(DoubleVector[] T000) {
    // Compute the normal Coulomb auxiliary term.
    super.source(T000);

    // Add the Thole damping terms: edamp = exp(-thole*u^3).
    tholeSource(thole, AiAk, R, T000);
  }

  /**
   * Generate source terms for the Challacombe et al. recursion.
   *
   * @param thole Thole damping parameter is set to min(pti,ptk)).
   * @param AiAk  parameter = 1/(alphaI^6*alphaK^6) where alpha is polarizability.
   * @param R     The separation distance.
   * @param T000  Location to store the source terms.
   */
  protected static void tholeSource(DoubleVector thole, DoubleVector AiAk, DoubleVector R, DoubleVector[] T000) {
    // Add the Thole damping terms: edamp = exp(-thole*u^3).
    DoubleVector u = R.mul(AiAk);
    DoubleVector u3 = thole.mul(u.mul(u).mul(u));
    DoubleVector u6 = u3.mul(u3);
    DoubleVector u9 = u6.mul(u3);
    DoubleVector expU3 = u3.neg().lanewise(VectorOperators.EXP);

    // The zeroth order term is not calculated for Thole damping.
    T000[0] = DoubleVector.zero(R.species());
    T000[1] = T000[1].mul(expU3);
    T000[2] = T000[2].mul((u3.add(1.0)).mul(expU3));
    T000[3] = T000[3].mul((u3.add(1.0).add(u6.mul(threeFifths)).mul(expU3)));
    T000[4] = T000[4].mul((u3.add(1.0).add(u6.mul(18.0).add(u9.mul(9.0)).mul(oneThirtyFifth)).mul(expU3)));
  }

}
