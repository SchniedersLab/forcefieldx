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

import static ffx.numerics.math.ScalarMath.doubleFactorial;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

public class GKSource {

  /**
   * The "mode" for the tensor (either POTENTIAL or BORN).
   */
  public enum GK_TENSOR_MODE {POTENTIAL, BORN}

  /**
   * The GK tensor can be constructed for a monopole potential (GB), a dipole potential or a
   * quadrupole potential.
   */
  public enum GK_MULTIPOLE_ORDER {
    MONOPOLE(0), DIPOLE(1), QUADRUPOLE(2);

    private final int order;

    GK_MULTIPOLE_ORDER(int order) {
      this.order = order;
    }

    public int getOrder() {
      return order;
    }
  }

  private GK_TENSOR_MODE mode = GK_TENSOR_MODE.POTENTIAL;

  /**
   * Born radius of atom i multiplied by the Born radius of atom j.
   */
  private double rb2;

  /**
   * The GK term exp(-r2 / (gc * ai * aj)).
   */
  private double expTerm;

  /**
   * The GK effective separation distance.
   */
  private double f;

  /**
   * The GK effective separation distance.
   */
  private double f1;

  /**
   * The GK effective separation distance.
   */
  private double f2;

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
  private double gcAiAj;

  /**
   * The ratio -r2 / (gc * Ai * Aj).
   */
  private double ratio;

  /**
   * The quantity: -2.0 / (gc * Ai * Aj);
   */
  private double fr;

  /**
   * Separation distance squared.
   */
  private double r2;

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
  private final double[] fn;

  /**
   * Chain rule terms from differentiating zeroth order auxiliary functions (an0) with respect to Ai
   * or Aj.
   */
  protected final double[] bn;

  /**
   * Source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  private final double[][] anm;

  /**
   * Source terms for the Kirkwood version of the Challacombe et al. recursion for Born-chain rule
   * derivatives.
   */
  private final double[][] bnm;

  public GKSource(int order, double gc) {
    this.order = order;
    this.gc = gc;
    this.igc = 1.0 / gc;

    // Auxiliary terms for Generalized Kirkwood (equivalent to Coulomb and Thole Screening).
    kirkwoodSource = new double[order + 1];
    for (int n = 0; n <= order; n++) {
      kirkwoodSource[n] = pow(-1, n) * doubleFactorial(2 * n - 1);
    }

    anmc = new double[order + 1][];
    for (int n = 0; n <= order; n++) {
      anmc[n] = anmc(n);
    }

    fn = new double[order + 1];
    anm = new double[order + 1][order + 1];
    bn = new double[order + 1];
    bnm = new double[order + 1][order + 1];
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  protected void source(double[] work, GK_MULTIPOLE_ORDER multipoleOrder) {
    int mpoleOrder = multipoleOrder.getOrder();
    fill(work, 0, mpoleOrder, 0.0);
    if (mode == GK_TENSOR_MODE.POTENTIAL) {
      // Max derivatives.
      int derivatives = order - mpoleOrder;
      if (derivatives + 1 >= 0) {
        arraycopy(anm[mpoleOrder], 0, work, mpoleOrder, derivatives + 1);
      }
    } else {
      // Max derivatives.
      int derivatives = (order - 1) - mpoleOrder;
      if (derivatives + 1 >= 0) {
        arraycopy(bnm[mpoleOrder], 0, work, mpoleOrder, derivatives + 1);
      }
    }
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  public void generateSource(GK_TENSOR_MODE mode, GK_MULTIPOLE_ORDER multipole, double r2,
      double ai, double aj) {
    int multipoleOrder = multipole.getOrder();
    this.mode = mode;

    this.r2 = r2;
    rb2 = ai * aj;
    gcAiAj = gc * rb2;
    ratio = -r2 / gcAiAj;
    expTerm = exp(ratio);
    fr = -2.0 / gcAiAj;

    f = sqrt(r2 + rb2 * expTerm);
    f1 = 1.0 - expTerm * igc;
    f2 = 2.0 * expTerm / (gc * gcAiAj);

    if (mode == GK_TENSOR_MODE.POTENTIAL) {
      // Prepare the GK Potential tensor.
      anm(order, order - multipoleOrder);
    } else {
      // Prepare the GK Born-chain rule tensor.
      bnm(order - 1, (order - 1) - multipoleOrder);
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
      case 0:
        fn[0] = f;
        break;
      case 1:
        fn[0] = f;
        fn[1] = f1;
        break;
      case 2:
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        break;
      case 3:
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr * f2;
        break;
      case 4:
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr * f2;
        fn[4] = fr * fn[3];
        break;
      case 5:
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr * f2;
        fn[4] = fr * fn[3];
        fn[5] = fr * fn[4];
        break;
      case 6:
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        fn[3] = fr * f2;
        fn[4] = fr * fn[3];
        fn[5] = fr * fn[4];
        fn[6] = fr * fn[5];
        break;
      default:
        fn[0] = f;
        fn[1] = f1;
        fn[2] = f2;
        for (int i = 3; i <= n; i++) {
          fn[i] = fr * fn[i - 1];
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
    var b2 = 2.0 * expTerm / (gcAiAj * gcAiAj) * (-ratio - 1.0);
    switch (n) {
      case 0:
        bn[0] = 0.5 * expTerm * (1.0 - ratio);
        break;
      case 1:
        bn[0] = 0.5 * expTerm * (1.0 - ratio);
        bn[1] = -r2 * expTerm / (gcAiAj * gcAiAj);
        break;
      case 2:
        bn[0] = 0.5 * expTerm * (1.0 - ratio);
        bn[1] = -r2 * expTerm / (gcAiAj * gcAiAj);
        bn[2] = b2;
        break;
      default:
        bn[0] = 0.5 * expTerm * (1.0 - ratio);
        bn[1] = -r2 * expTerm / (gcAiAj * gcAiAj);
        bn[2] = b2;
        var br = 2.0 / (gcAiAj * rb2);
        var f2 = 2.0 / (gc * gcAiAj) * expTerm;
        var frA = 1.0;
        var frB = fr;
        // (n - 2) * pow(fr, n - 3) * br * f2 + pow(fr, n - 2) * b2;
        for (int i = 3; i < n; i++) {
          bn[i] = (i - 2) * frA * br * f2 + frB * b2;
          frA *= fr;
          frB *= fr;
        }
    }
  }

  /**
   * Compute the potential auxiliary function up to order n.
   *
   * @param n order.
   */
  private void an0(int n) {
    double inverseF = 1.0 / f;
    double inverseF2 = inverseF * inverseF;
    for (int i = 0; i <= n; i++) {
      anm[i][0] = kirkwoodSource[i] * inverseF;
      inverseF *= inverseF2;
    }
  }

  /**
   * Compute the Born chain-rule auxiliary function up to order n.
   *
   * @param n order.
   */
  private void bn0(int n) {
    for (int i = 0; i <= n; i++) {
      bnm[i][0] = bn[0] * anm[i + 1][0];
    }
  }

  /**
   * Fill the GK auxiliary matrix.
   * <p>
   * The first row is the GK potential and derivatives for monpoles. The second row is the GK
   * potential and derivatives for dipoles. The third row is the GK potential and derivatives for
   * quadrupoles.
   *
   * @param n Order.
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
        var sum = 0.0;
        for (int i = 1; i <= d; i++) {
          sum += coef[i - 1] * fn[i] * terms[d - i];
        }
        anm[order][d] = sum;
      }
    }
  }

  /**
   * Fill the GK auxiliary matrix derivatives with respect to Born radii.
   * <p>
   * The first row are derivatives for the monpole potential. The second row are derivatives for the
   * dipole potential. The third row are derivatives for the quadrupole potential.
   *
   * @param n Order.
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
        var sum = 0.0;
        for (int i = 1; i <= d; i++) {
          sum += coef[i - 1] * bn[i] * terma[d - i];
          sum += coef[i - 1] * fn[i] * termb[d - i];
        }
        bnm[order][d] = sum;
      }
    }
  }

  /**
   * Return coefficients needed when taking derivatives of auxiliary functions.
   *
   * @param n Multipole order.
   * @return Returns coefficients needed when taking derivatives of auxiliary functions.
   */
  protected static double[] anmc(int n) {
    double[] ret = new double[n + 1];
    ret[0] = 1.0;
    switch (n) {
      case 0:
        return ret;
      case 1:
        ret[1] = 1.0;
        return ret;
      default:
        ret[1] = 1.0;
        double[] prev = new double[n];
        prev[0] = 1.0;
        prev[1] = 1.0;
        for (int i = 3; i <= n; i++) {
          for (int j = 2; j <= i - 1; j++) {
            ret[j - 1] = prev[j - 2] + prev[j - 1];
          }
          ret[i - 1] = 1.0;
          arraycopy(ret, 0, prev, 0, i);
        }
        return ret;
    }
  }

  public static double selfEnergy(PolarizableMultipole polarizableMultipole, double ai, double aj,
      double Eh, double Es) {
    double q2 = polarizableMultipole.q * polarizableMultipole.q;
    double dx = polarizableMultipole.dx;
    double dy = polarizableMultipole.dy;
    double dz = polarizableMultipole.dz;
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dz2 = dz * dz;
    double ux = polarizableMultipole.ux;
    double uy = polarizableMultipole.uy;
    double uz = polarizableMultipole.uz;
    double qxy2 = polarizableMultipole.qxy * polarizableMultipole.qxy;
    double qxz2 = polarizableMultipole.qxz * polarizableMultipole.qxz;
    double qyz2 = polarizableMultipole.qyz * polarizableMultipole.qyz;
    double qxx2 = polarizableMultipole.qxx * polarizableMultipole.qxx;
    double qyy2 = polarizableMultipole.qyy * polarizableMultipole.qyy;
    double qzz2 = polarizableMultipole.qzz * polarizableMultipole.qzz;

    double a = sqrt(ai * aj);
    double a2 = a * a;
    double a3 = a * a2;
    double a5 = a2 * a3;

    // Born partial charge
    double e0 = cn(0, Eh, Es) * q2 / a;
    // Permanent Dipole
    double e1 = cn(1, Eh, Es) * (dx2 + dy2 + dz2) / a3;
    // Permanent Quadrupole
    double e2 = cn(2, Eh, Es) * (3.0 * (qxy2 + qxz2 + qyz2) + 6.0 * (qxx2 + qyy2 + qzz2)) / a5;

    // Induced self-energy
    double ei = cn(1, Eh, Es) * (dx * ux + dy * uy + dz * uz) / a3;

    return 0.5 * (e0 + e1 + e2 + ei);
  }

  /**
   * Compute the Kirkwood dielectric function for a multipole of order n.
   *
   * @param n Multipole order.
   * @param Eh Homogeneous dielectric.
   * @param Es Solvent dielectric.
   * @return Returns (n+1)*(Eh-Es)/((n+1)*Es + n*Eh))
   */
  public static double cn(int n, double Eh, double Es) {
    return (n + 1) * (Eh - Es) / ((n + 1) * Es + n * Eh);
  }

}
