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
public class GKTensorQI extends CoulombTensorQI {

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
   * The GK tensor can be constructed for monopoles (GB), dipoles or quadrupoles.
   */
  protected final GK_MULTIPOLE_ORDER multipoleOrder;

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
   * The "mode" for the tensor (either POTENTIAL or BORN).
   */
  private GK_TENSOR_MODE mode = GK_TENSOR_MODE.POTENTIAL;

  /**
   * The "mode" for the tensor (either POTENTIAL or BORN).
   */
  public enum GK_TENSOR_MODE {POTENTIAL, BORN}

  /**
   * The GK tensor can be constructed for monopoles (GB), dipoles or quadrupoles.
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

  /**
   * @param multipoleOrder The multipole order.
   * @param order The number of derivatives to complete.
   * @param gc Generalized Kirkwood constant.
   * @param Eh Homogeneous dielectric constant.
   * @param Es Solvent dielectric constant.
   */
  public GKTensorQI(GK_MULTIPOLE_ORDER multipoleOrder, int order, double gc, double Eh, double Es) {
    super(order);
    this.multipoleOrder = multipoleOrder;
    this.order = order;
    this.gc = gc;
    this.Eh = Eh;
    this.Es = Es;

    // Load the dielectric function
    c = cn(multipoleOrder.getOrder(), Eh, Es);

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
   * GK Permanent multipole energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return the GK permanent multipole energy.
   */
  @Override
  public double multipoleEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    switch (multipoleOrder.getOrder()) {
      default:
      case 0:
        chargeIPotentialAtK(mI, 2);
        double eK = multipoleEnergy(mK);
        chargeKPotentialAtI(mK, 2);
        double eI = multipoleEnergy(mI);
        return 0.5 * (eK + eI);
      case 1:
        dipoleIPotentialAtK(mI, 2);
        eK = multipoleEnergy(mK);
        dipoleKPotentialAtI(mK, 2);
        eI = multipoleEnergy(mI);
        return 0.5 * (eK + eI);
      case 2:
        quadrupoleIPotentialAtK(mI, 2);
        eK = multipoleEnergy(mK);
        quadrupoleKPotentialAtI(mK, 2);
        eI = multipoleEnergy(mI);
        return 0.5 * (eK + eI);
    }
  }

  /**
   * GK Permanent multipole energy and gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  @Override
  public double multipoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {
    switch (multipoleOrder) {
      default:
      case MONOPOLE:
        return monopoleEnergyAndGradient(mI, mK, Gi, Gk, Ti, Tk);
      case DIPOLE:
        return dipoleEnergyAndGradient(mI, mK, Gi, Gk, Ti, Tk);
      case QUADRUPOLE:
        return quadrupoleEnergyAndGradient(mI, mK, Gi, Gk, Ti, Tk);
    }
  }

  /**
   * Permanent multipole energy and gradient using the GK monopole tensor.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  protected double monopoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {

    // Compute the potential due to a multipole component at site I.
    chargeIPotentialAtK(mI, 3);
    double eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Compute the potential due to a multipole component at site K.
    chargeKPotentialAtI(mK, 3);
    double eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    Gi[0] = 0.5 * (Gi[0] - Gk[0]);
    Gi[1] = 0.5 * (Gi[1] - Gk[1]);
    Gi[2] = 0.5 * (Gi[2] - Gk[2]);
    Gk[0] = -Gi[0];
    Gk[1] = -Gi[1];
    Gk[2] = -Gi[2];

    Ti[0] = 0.5 * Ti[0];
    Ti[1] = 0.5 * Ti[1];
    Ti[2] = 0.5 * Ti[2];
    Tk[0] = 0.5 * Tk[0];
    Tk[1] = 0.5 * Tk[1];
    Tk[2] = 0.5 * Tk[2];

    return 0.5 * (eK + eI);
  }

  /**
   * Permanent multipole energy and gradient using the GK dipole tensor.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  protected double dipoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {

    // Compute the potential due to a multipole component at site I.
    dipoleIPotentialAtK(mI, 3);
    double eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Need the torque on site I pole due to site K multipole.
    // Only torque on the site I dipole.
    multipoleKPotentialAtI(mK, 1);
    dipoleTorque(mI, Ti);

    // Compute the potential due to a multipole component at site K.
    dipoleKPotentialAtI(mK, 3);
    double eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    // Need the torque on site K pole due to site I multipole.
    // Only torque on the site K dipole.
    multipoleIPotentialAtK(mI, 1);
    dipoleTorque(mK, Tk);

    Gi[0] = 0.5 * (Gi[0] - Gk[0]);
    Gi[1] = 0.5 * (Gi[1] - Gk[1]);
    Gi[2] = 0.5 * (Gi[2] - Gk[2]);
    Gk[0] = -Gi[0];
    Gk[1] = -Gi[1];
    Gk[2] = -Gi[2];

    Ti[0] = 0.5 * Ti[0];
    Ti[1] = 0.5 * Ti[1];
    Ti[2] = 0.5 * Ti[2];
    Tk[0] = 0.5 * Tk[0];
    Tk[1] = 0.5 * Tk[1];
    Tk[2] = 0.5 * Tk[2];

    return 0.5 * (eK + eI);
  }

  /**
   * Permanent multipole energy and gradient using the GK quadrupole tensor.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole GK energy.
   */
  protected double quadrupoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
      double[] Gi, double[] Gk, double[] Ti, double[] Tk) {

    // Compute the potential due to a multipole component at site I.
    quadrupoleIPotentialAtK(mI, 3);
    double eK = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    multipoleTorque(mK, Tk);

    // Need the torque on site I pole due to site K multipole.
    // Only torque on the site I quadrupole.
    multipoleKPotentialAtI(mK, 2);
    quadrupoleTorque(mI, Ti);

    // Compute the potential due to a multipole component at site K.
    quadrupoleKPotentialAtI(mK, 3);
    double eI = multipoleEnergy(mI);
    multipoleGradient(mI, Gi);
    multipoleTorque(mI, Ti);

    // Need the torque on site K pole due to site I multipole.
    // Only torque on the site K quadrupole.
    multipoleIPotentialAtK(mI, 2);
    quadrupoleTorque(mK, Tk);

    Gi[0] = 0.5 * (Gi[0] - Gk[0]);
    Gi[1] = 0.5 * (Gi[1] - Gk[1]);
    Gi[2] = 0.5 * (Gi[2] - Gk[2]);
    Gk[0] = -Gi[0];
    Gk[1] = -Gi[1];
    Gk[2] = -Gi[2];

    Ti[0] = 0.5 * Ti[0];
    Ti[1] = 0.5 * Ti[1];
    Ti[2] = 0.5 * Ti[2];
    Tk[0] = 0.5 * Tk[0];
    Tk[1] = 0.5 * Tk[1];
    Tk[2] = 0.5 * Tk[2];

    return 0.5 * (eK + eI);
  }

  /**
   * GK Permanent multipole Born grad.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public double multipoleEnergyBornGrad(PolarizableMultipole mI, PolarizableMultipole mK) {
    GK_TENSOR_MODE currentMode = mode;
    setMode(GK_TENSOR_MODE.BORN);
    generateTensor();
    double db = multipoleEnergy(mI, mK);
    setMode(currentMode);
    return db;
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param scaleEnergy This is ignored, since masking/scaling is not applied to GK
   *     interactions.
   * @return a double.
   */
  @Override
  public double polarizationEnergy(PolarizableMultipole mI, PolarizableMultipole mK,
      double scaleEnergy) {
    return polarizationEnergy(mI, mK);
  }

  /**
   * GK Polarization Energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public double polarizationEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    switch (multipoleOrder) {
      default:
      case MONOPOLE:
        // Find the GK charge potential of site I at site K.
        chargeIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent charge I.
        double eK = polarizationEnergy(mK);
        // Find the GK charge potential of site K at site I.
        chargeKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent charge K.
        double eI = polarizationEnergy(mI);
        return 0.5 * (eK + eI);
      case DIPOLE:
        // Find the GK dipole potential of site I at site K.
        dipoleIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent dipole I.
        eK = polarizationEnergy(mK);
        // Find the GK induced dipole potential of site I at site K.
        inducedIPotentialAtK(mI);
        // Energy of permanent multipole K in the field of induced dipole I.
        eK += 0.5 * multipoleEnergy(mK);
        // Find the GK dipole potential of site K at site I.
        dipoleKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent dipole K.
        eI = polarizationEnergy(mI);
        // Find the GK induced dipole potential of site K at site I.
        inducedKPotentialAtI(mK);
        // Energy of permanent multipole I in the field of induced dipole K.
        eI += 0.5 * multipoleEnergy(mI);
        return 0.5 * (eK + eI);
      case QUADRUPOLE:
        // Find the GK quadrupole potential of site I at site K.
        quadrupoleIPotentialAtK(mI, 1);
        // Energy of induced dipole K in the field of permanent quadrupole I.
        eK = polarizationEnergy(mK);
        // Find the GK quadrupole potential of site K at site I.
        quadrupoleKPotentialAtI(mK, 1);
        // Energy of induced dipole I in the field of permanent quadrupole K.
        eI = polarizationEnergy(mI);
        return 0.5 * (eK + eI);
    }
  }

  /**
   * Set the separation vector.
   *
   * @param r Separation vector.
   * @param ai Born radius for Atom i.
   * @param aj Born radius for Atom j.
   */
  public void setR(double[] r, double ai, double aj) {
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
  public void setR(double dx, double dy, double dz, double ai, double aj) {
    setR(dx, dy, dz);
    this.ai = ai;
    this.aj = aj;
    expTerm = exp(-r2 / (gc * ai * aj));
    f = sqrt(r2 + ai * aj * expTerm);
  }

  /**
   * Set the "mode" for the tensor (either POTENTIAL or BORN).
   *
   * @param mode The mode for tensor generation.
   */
  protected void setMode(GK_TENSOR_MODE mode) {
    this.mode = mode;
  }

  /**
   * Generate source terms for the Kirkwood version of the Challacombe et al. recursion.
   */
  protected void source(double[] work) {
    int multipoleOrder = this.multipoleOrder.getOrder();
    if (mode == GK_TENSOR_MODE.POTENTIAL) {
      // Prepare the GK Potential tensor.
      for (int n = 0; n <= order; n++) {
        an0[n] = an0(n);
        fn[n] = fn(n);
      }
      for (int n = 0; n <= order; n++) {
        if (n < multipoleOrder) {
          work[n] = 0.0;
        } else {
          work[n] = anm(multipoleOrder, n - multipoleOrder);
        }
      }
    } else {
      // Prepare the GK Born-chain rule tensor.
      for (int n = 0; n <= order; n++) {
        an0[n] = an0(n);
        fn[n] = fn(n);
        bn[n] = bn(n);
      }
      // Only up to order - 1.
      for (int n = 0; n <= order - 1; n++) {
        if (n < multipoleOrder) {
          work[n] = 0.0;
        } else {
          work[n] = bnm(multipoleOrder, n - multipoleOrder);
        }
      }
    }
  }

  /**
   * Compute the potential auxiliary function for a multipole of order n.
   *
   * @param n Multipole order.
   * @return The potential auxiliary function for a multipole of order n.
   */
  protected double an0(int n) {
    return kirkwoodSource[n] / pow(f, 2 * n + 1);
  }

  /**
   * Compute the mth potential gradient auxiliary function for a multipole of order n.
   *
   * @param n Multipole order.
   * @param m Mth potential gradient auxiliary function.
   * @return Returns the mth potential gradient auxiliary function for a multipole of order n.
   */
  protected double anm(int n, int m) {
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
  protected double bnm(int n, int m) {
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
  protected double fn(int n) {
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
  protected static double[] anmc(int n) {
    return GKTensorGlobal.anmc(n);
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
    return GKTensorGlobal.cn(n, Eh, Es);
  }

}