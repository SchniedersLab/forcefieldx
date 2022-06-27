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
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The GeneralizedKirkwoodTensor class contains utilities for generated Generalized Kirkwood
 * interaction tensors.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class GKTensorGlobal extends CoulombTensorGlobal {

  /**
   * Generalized Kirkwood constant.
   */
  private final double gc;

  /**
   * Homogeneous dielectric constant.
   */
  private final double Eh;

  /**
   * Solvent dielectric constant.
   */
  private final double Es;

  /**
   * Order of the tensor recursion (5th is needed for AMOEBA forces).
   */
  private final int order;

  /**
   * Multipole order (2nd is needed for AMOEBA forces).
   */
  protected final int multipoleOrder;

  /**
   * The Kirkwood dielectric function for the given multipole order.
   */
  private final double c;

  /**
   * Compute the "source" terms for the recursion.
   */
  protected final double[] kirkwoodSource;

  /**
   * Coefficients needed when taking derivatives of auxiliary functions.
   */
  private final double[][] anmc;

  /**
   * Source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  private final double[] an0;

  /**
   * Chain rule terms from differentiating zeroth order auxiliary functions (an0) with respect to x,
   * y or z.
   */
  private final double[] fn;

  /**
   * Chain rule terms from differentiating zeroth order auxiliary functions (an0) with respect to Ai
   * or Aj.
   */
  private final double[] bn;

  /**
   * Born radius of atom i.
   */
  private double ai;

  /**
   * Born radius of atom j.
   */
  private double aj;

  /**
   * The GK term exp(-r2 / (gc * ai * aj)).
   */
  private double expTerm;

  /**
   * The GK effective separation distance.
   */
  private double f;

  /**
   * @param multipoleOrder The multipole order.
   * @param order The number of derivatives to complete.
   * @param gc Generalized Kirkwood constant.
   * @param Eh Homogeneous dielectric constant.
   * @param Es Solvent dielectric constant.
   */
  public GKTensorGlobal(int multipoleOrder, int order, double gc, double Eh, double Es) {
    super(order);
    this.multipoleOrder = multipoleOrder;
    this.order = order;
    this.gc = gc;
    this.Eh = Eh;
    this.Es = Es;

    // Load the dielectric function
    c = cn(multipoleOrder, Eh, Es);

    // Auxiliary terms for Generalized Kirkwood (equivalent to Coulomb and Thole Screening).
    kirkwoodSource = new double[order + 1];
    for (int n = 0; n <= order; n++) {
      kirkwoodSource[n] = c * pow(-1, n) * doubleFactorial(2 * n - 1);
    }

    anmc = new double[order + 1][];
    for (int n = 0; n <= order; n++) {
      anmc[n] = anmc(n);
    }

    an0 = new double[order + 1];
    fn = new double[order + 1];
    bn = new double[order + 1];
  }

  public double selfEnergy(PolarizableMultipole polarizableMultipole) {
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
   * Set the separation vector.
   *
   * @param r Separation vector.
   * @param ai Born radius for Atom i.
   * @param aj Born radius for Atom j.
   */
  protected void setR(double[] r, double ai, double aj) {
    setR(r[0], r[1], r[2], ai, aj);
  }

  /**
   * Set the separation vector.
   *
   * @param dx Separation along the X-axis.
   * @param dy Separation along the Y-axis.
   * @param dz Separation along the Z-axis.
   * @param ai Born radius for Atom i.
   * @param aj Born radius for Atom j.
   */
  protected void setR(double dx, double dy, double dz, double ai, double aj) {
    setR(dx, dy, dz);
    this.ai = ai;
    this.aj = aj;
    expTerm = exp(-r2 / (gc * ai * aj));
    f = sqrt(r2 + ai * aj * expTerm);
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  protected void source(double[] work) {
    for (int n = 0; n <= order; n++) {
      an0[n] = an0(n);
      fn[n] = fn(n);
      bn[n] = bn(n);
    }

    for (int n = 0; n <= order; n++) {
      if (n < multipoleOrder) {
        work[n] = 0.0;
      } else {
        work[n] = anm(multipoleOrder, n - multipoleOrder);
      }
    }
  }

  /**
   * Compute the potential auxiliary function for a multipole of order n.
   *
   * @param n Multipole order.
   * @return The potential auxiliary function for a multipole of order n.
   */
  public double an0(int n) {
    return kirkwoodSource[n] / pow(f, 2 * n + 1);
  }

  /**
   * Compute the mth potential gradient auxiliary function for a multipole of order n.
   *
   * @param n Multipole order.
   * @param m Mth potential gradient auxiliary function.
   * @return Returns the mth potential gradient auxiliary function for a multipole of order n.
   */
  public double anm(int n, int m) {
    if (m == 0) {
      return an0[n];
    }
    var ret = 0.0;
    var coef = anmc[m];
    for (int i = 1; i <= m; i++) {
      ret += coef[i - 1] * fn[i] * anm(n + 1, m - i);
    }
    return ret;
  }

  /**
   * Compute the derivative with respect to a Born radius of the mth potential gradient auxiliary
   * function for a multipole of order n.
   *
   * @param n Multipole order.
   * @param m Mth potential gradient auxiliary function.
   * @return Returns the derivative with respect to a Born radius of the mth potential gradient
   *     auxiliary function for a multipole of order n.
   */
  public double bnm(int n, int m) {
    if (m == 0) {
      // return bn(0) * an0(n + 1);
      return bn[0] * an0[n + 1];
    }
    var ret = 0.0;
    var coef = anmc[m];
    for (int i = 1; i <= m; i++) {
      ret += coef[i - 1] * bn[i] * anm(n + 1, m - i);
      ret += coef[i - 1] * fn[i] * bnm(n + 1, m - i);
    }
    return ret;
  }

  /**
   * Returns nth value of the function f, which are chain rule terms from differentiating zeroth
   * order auxiliary functions (an0) with respect to x, y or z.
   *
   * @param n Multipole order.
   * @return Returns the nth value of the function f.
   */
  private double fn(int n) {
    switch (n) {
      case 0:
        return f;
      case 1:
        return 1.0 - expTerm / gc;
      default:
        var gcAiAj = gc * ai * aj;
        var f2 = 2.0 * expTerm / (gc * gcAiAj);
        var fr = -2.0 / gcAiAj;
        return pow(fr, n - 2) * f2;
    }
  }

  /**
   * Returns nth value of the function b, which are chain rule terms from differentiating zeroth
   * order auxiliary functions (an0) with respect to Ai or Aj.
   *
   * @param n Multipole order.
   * @return Returns the nth value of the function f.
   */
  protected double bn(int n) {
    var gcAiAj = gc * ai * aj;
    var ratio = -r2 / gcAiAj;
    switch (n) {
      case 0:
        return 0.5 * expTerm * (1.0 - ratio);
      case 1:
        return -r2 * expTerm / (gcAiAj * gcAiAj);
      default:
        var b2 = 2.0 * expTerm / (gcAiAj * gcAiAj) * (-ratio - 1.0);
        var br = 2.0 / (gcAiAj * ai * aj);
        var f2 = 2.0 / (gc * gcAiAj) * expTerm;
        var fr = -2.0 / (gcAiAj);
        return (n - 2) * pow(fr, n - 3) * br * f2 + pow(fr, n - 2) * b2;
    }
  }

  /**
   * Return coefficients needed when taking derivatives of auxiliary functions.
   *
   * @param n Multipole order.
   * @return Returns coefficients needed when taking derivatives of auxiliary functions.
   */
  public static double[] anmc(int n) {
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
