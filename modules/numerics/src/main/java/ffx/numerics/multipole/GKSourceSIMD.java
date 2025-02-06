// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import static ffx.numerics.math.ScalarMath.doubleFactorial;
import static ffx.numerics.multipole.GKTensorMode.POTENTIAL;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static jdk.incubator.vector.DoubleVector.SPECIES_PREFERRED;
import static jdk.incubator.vector.VectorOperators.EXP;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The GKSource class generates the source terms for the Generalized Kirkwood version of the tensor recursion.
 */
public class GKSourceSIMD {

  private GKTensorMode mode = POTENTIAL;

  /**
   * Born radius of atom i multiplied by the Born radius of atom j.
   */
  private DoubleVector rb2;

  /**
   * The GK term exp(-r2 / (gc * ai * aj)).
   */
  private DoubleVector expTerm;

  /**
   * The GK effective separation distance.
   */
  private DoubleVector f;

  /**
   * f1 = 1.0 - expTerm * igc;
   */
  private DoubleVector f1;

  /**
   * f2 = 2.0 * expTerm / (gc * gcAiAj);
   */
  private DoubleVector f2;

  /**
   * Generalized Kirkwood constant.
   */
  private final double gc;

  /**
   * Inverse generalized Kirkwood constant.
   */
  private final double igc;

  /**
   * The product: gc * Ai * Aj.
   */
  private DoubleVector gcAiAj;

  /**
   * The ratio -r2 / (gc * Ai * Aj).
   */
  private DoubleVector ratio;

  /**
   * The quantity: -2.0 / (gc * Ai * Aj);
   */
  private DoubleVector fr;

  /**
   * Separation distance squared.
   */
  private DoubleVector r2;

  /**
   * Recursion order.
   */
  private final int order;

  /**
   * Compute the "source" terms for the recursion.
   */
  private final double[] kirkwoodSource;

  /**
   * Coefficients needed when taking derivatives of auxiliary functions.
   */
  private final double[][] anmc;

  /**
   * Chain rule terms from differentiating zeroth order auxiliary functions (an0) with respect to x,
   * y or z.
   */
  private final DoubleVector[] fn;

  /**
   * Chain rule terms from differentiating zeroth order born radii auxiliary functions (bn0) with respect to Ai
   * or Aj.
   */
  protected final DoubleVector[] bn;

  /**
   * Source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  private final DoubleVector[][] anm;

  /**
   * Source terms for the Kirkwood version of the Challacombe et al. recursion for Born-chain rule
   * derivatives.
   */
  private final DoubleVector[][] bnm;

  private final DoubleVector zero = DoubleVector.zero(SPECIES_PREFERRED);
  private final DoubleVector one = zero.add(1.0);

  /**
   * Construct a new GKSource object.
   *
   * @param order Recursion order.
   * @param gc    Generalized Kirkwood constant.
   */
  public GKSourceSIMD(int order, double gc) {
    this.order = order;
    this.gc = gc;
    this.igc = 1.0 / gc;

    // Auxiliary terms for Generalized Kirkwood (equivalent to Coulomb and Thole Screening).
    kirkwoodSource = new double[order + 1];
    for (short n = 0; n <= order; n++) {
      kirkwoodSource[n] = pow(-1, n) * doubleFactorial(2 * n - 1);
    }

    anmc = new double[order + 1][];
    for (int n = 0; n <= order; n++) {
      anmc[n] = GKSource.anmc(n);
    }

    fn = new DoubleVector[order + 1];
    anm = new DoubleVector[order + 1][order + 1];
    bn = new DoubleVector[order + 1];
    bnm = new DoubleVector[order + 1][order + 1];
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   *
   * @param work           The array to store the source terms.
   * @param multipoleOrder The multipole order.
   */
  protected void source(DoubleVector[] work, GKMultipoleOrder multipoleOrder) {
    int mpoleOrder = multipoleOrder.getOrder();
    fill(work, 0, mpoleOrder, zero);
    if (mode == POTENTIAL) {
      // Max derivatives.
      int derivatives = order - mpoleOrder;
      arraycopy(anm[mpoleOrder], 0, work, mpoleOrder, derivatives + 1);
    } else {
      // Max derivatives.
      int derivatives = (order - 1) - mpoleOrder;
      arraycopy(bnm[mpoleOrder], 0, work, mpoleOrder, derivatives + 1);
    }
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   *
   * @param mode      The tensor mode.
   * @param multipole The multipole order.
   * @param r2        Separation distance squared.
   * @param ai        Born radius of atom i.
   * @param aj        Born radius of atom j.
   */
  public void generateSource(GKTensorMode mode, GKMultipoleOrder multipole,
                             DoubleVector r2, DoubleVector ai, DoubleVector aj) {
    int multipoleOrder = multipole.getOrder();
    this.mode = mode;

    this.r2 = r2;
    rb2 = ai.mul(aj);
    gcAiAj = rb2.mul(gc);
    ratio = r2.neg().div(gcAiAj);
    expTerm = ratio.lanewise(EXP);
    fr = zero.sub(2.0).div(gcAiAj);
    f = r2.add(rb2.mul(expTerm)).sqrt();
    f1 = one.sub(expTerm.mul(igc));
    f2 = zero.add(2.0).mul(expTerm).div(gcAiAj.mul(gc));
    if (mode == POTENTIAL) {
      // Prepare the GK Potential tensor.
      anm(order, order - multipoleOrder);
    } else {
      // Prepare the GK Born-chain rule tensor.
      bnm(order - 1, (order - 1) - multipoleOrder);
    }
  }

  /**
   * Fill the GK auxiliary matrix.
   * <p>
   * The first row is the GK potential and derivatives for monopoles. The second row is the GK
   * potential and derivatives for dipoles. The third row is the GK potential and derivatives for
   * quadrupoles.
   *
   * @param n           Order.
   * @param derivatives Number of derivatives.
   */
  private void anm(int n, int derivatives) {
    // Fill the fn chain rule terms.
    fn(n);

    // Fill the auxiliary potential.
    an0(n);

    // Derivative loop over columns.
    for (int d = 1; d <= derivatives; d++) {
      // The coefficients are the same for each order.
      var coef = anmc[d];
      // The current derivative can only be computed for order n - d.
      int limit = n - d;
      // Order loop over rows.
      for (int order = 0; order <= limit; order++) {
        // Compute this term from previous terms for 1 higher order.
        var terms = anm[order + 1];
        var sum = zero;
        for (int i = 1; i <= d; i++) {
          sum = sum.add(fn[i].mul(terms[d - i].mul(coef[i - 1])));
        }
        anm[order][d] = sum;
      }
    }
  }

  /**
   * Fill the GK auxiliary matrix derivatives with respect to Born radii.
   * <p>
   * The first row are derivatives for the monopole potential. The second row are derivatives for the
   * dipole potential. The third row are derivatives for the quadrupole potential.
   *
   * @param n           Order.
   * @param derivatives Number of derivatives.
   */
  private void bnm(int n, int derivatives) {
    // Fill the bn chain rule terms.
    bn(n);

    // Fill the auxiliary potential derivatives.
    bn0(n);

    // Derivative loop over columns.
    for (int d = 1; d <= derivatives; d++) {
      // The coefficients are the same for each order.
      var coef = anmc[d];
      // The current derivative can only be computed for order n - d.
      int limit = n - d;
      // Order loop over rows.
      for (int order = 0; order <= limit; order++) {
        // Compute this term from previous terms for 1 higher order.
        var terma = anm[order + 1];
        var termb = bnm[order + 1];
        var sum = zero;
        for (int i = 1; i <= d; i++) {
          sum = sum.add(bn[i].mul(terma[d - i].mul(coef[i - 1])));
          sum = sum.add(fn[i].mul(termb[d - i].mul(coef[i - 1])));
        }
        bnm[order][d] = sum;
      }
    }
  }

  /**
   * Sets the function f, which are chain rule terms from differentiating zeroth order auxiliary
   * functions (an0) with respect to x, y or z.
   *
   * @param n Multipole order.
   */
  private void fn(int n) {
    switch (n) {
      case 0 -> fn[0] = f;
      case 1 -> {
        fn[0] = f;
        fn[1] = f1;
      }
      case 2 -> {
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
      }
      case 3 -> {
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr.mul(f2);
      }
      case 4 -> {
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr.mul(f2);
        fn[4] = fr.mul(fn[3]);
      }
      case 5 -> {
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr.mul(f2);
        fn[4] = fr.mul(fn[3]);
        fn[5] = fr.mul(fn[4]);
      }
      case 6 -> {
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr.mul(f2);
        fn[4] = fr.mul(fn[3]);
        fn[5] = fr.mul(fn[4]);
        fn[6] = fr.mul(fn[5]);
      }
      default -> {
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        for (int i = 3; i <= n; i++) {
          fn[i] = fr.mul(fn[i - 1]);
        }
      }
    }
  }

  /**
   * Compute the function b, which are chain rule terms from differentiating zeroth order auxiliary
   * functions (an0) with respect to Ai or Aj.
   *
   * @param n Multipole order.
   */
  protected void bn(int n) {
    var b2 = expTerm.mul(2.0).div(gcAiAj.mul(gcAiAj)).mul(ratio.neg().sub(1.0));
    switch (n) {
      case 0 -> bn[0] = expTerm.mul(0.5).mul(ratio.neg().add(1.0));
      case 1 -> {
        bn[0] = expTerm.mul(0.5).mul(ratio.neg().add(1.0));
        bn[1] = r2.mul(expTerm).div(gcAiAj.mul(gcAiAj)).neg();
      }
      case 2 -> {
        bn[0] = expTerm.mul(0.5).mul(ratio.neg().add(1.0));
        bn[1] = r2.mul(expTerm).div(gcAiAj.mul(gcAiAj)).neg();
        bn[2] = b2;
      }
      default -> {
        bn[0] = expTerm.mul(0.5).mul(ratio.neg().add(1.0));
        bn[1] = r2.mul(expTerm).div(gcAiAj.mul(gcAiAj)).neg();
        bn[2] = b2;
        var br = zero.add(2.0).div(gcAiAj.mul(rb2));
        var f2 = zero.add(2.0).div(gcAiAj.mul(gc)).mul(expTerm);
        var frA = one;
        var frB = fr;
        for (int i = 3; i < n; i++) {
          // bn[i] = (i - 2) * pow(fr, i - 3) * br * f2 + pow(fr, i - 2) * b2;
          bn[i] = frA.mul(i - 2).mul(br).mul(f2).add(frB.mul(b2));
          frA = frA.mul(fr);
          frB = frB.mul(fr);
        }
      }
    }
  }

  /**
   * Compute the potential auxiliary function up to order n.
   *
   * @param n order.
   */
  private void an0(int n) {
    DoubleVector inverseF = one.div(f);
    DoubleVector inverseF2 = inverseF.mul(inverseF);
    for (int i = 0; i <= n; i++) {
      anm[i][0] = inverseF.mul(kirkwoodSource[i]);
      inverseF = inverseF.mul(inverseF2);
    }
  }

  /**
   * Compute the Born chain-rule auxiliary function up to order n.
   *
   * @param n order.
   */
  private void bn0(int n) {
    for (int i = 0; i <= n; i++) {
      bnm[i][0] = bn[0].mul(anm[i + 1][0]);
    }
  }

}
