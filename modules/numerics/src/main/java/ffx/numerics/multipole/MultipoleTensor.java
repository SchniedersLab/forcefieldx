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
package ffx.numerics.multipole;

import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.ScalarMath.doubleFactorial;
import static ffx.numerics.multipole.MultipoleUtilities.lmn;
import static ffx.numerics.multipole.MultipoleUtilities.loadTensor;
import static ffx.numerics.multipole.MultipoleUtilities.storePotential;
import static ffx.numerics.multipole.MultipoleUtilities.storePotentialNeg;
import static ffx.numerics.multipole.MultipoleUtilities.term;
import static java.lang.Math.fma;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The abstract MultipoleTensor is extended by classes that compute derivatives of 1/|<b>r</b>| via
 * recursion to arbitrary order using Cartesian multipoles in either a global frame or a
 * quasi-internal frame. <br> This class serves as the abstract parent to both and defines all shared
 * logic. Non-abstract methods are declared final to disallow unnecessary overrides.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class MultipoleTensor {

  /**
   * Logger for the MultipoleTensor class.
   */
  private static final Logger logger = Logger.getLogger(MultipoleTensor.class.getName());

  /**
   * Order of the tensor recursion (5th is needed for AMOEBA forces).
   */
  protected final int order;

  /**
   * Order plus 1.
   */
  protected final int o1;

  /**
   * Order plus one.
   */
  protected final int il;

  /**
   * im = (Order plus one)^2.
   */
  protected final int im;

  /**
   * in = (Order plus one)^3.
   */
  protected final int in;

  /**
   * Size = (order + 1) * (order + 2) * (order + 3) / 6;
   */
  protected final int size;

  /**
   * The OPERATOR in use.
   */
  protected Operator operator;

  /**
   * These are the "source" terms for the recursion for the Coulomb operator (1/R).
   */
  protected final double[] coulombSource;

  /**
   * Store the work array to avoid memory consumption. Note that rather than use an array for
   * intermediate values, a 4D matrix was tried. It was approximately 50% slower than the linear work
   * array.
   */
  protected final double[] work;

  /**
   * Store the auxiliary tensor memory to avoid memory consumption.
   */
  protected final double[] T000;

  /**
   * The coordinate system in use (global or QI).
   */
  protected final CoordinateSystem coordinates;

  /**
   * Separation distance.
   */
  protected double R;

  /**
   * Separation distance squared.
   */
  protected double r2;

  /**
   * Xk - Xi.
   */
  protected double x;

  /**
   * Yk - Yi.
   */
  protected double y;

  /**
   * Zk - Zi.
   */
  protected double z;

  private double dEdZ = 0.0;
  private double d2EdZ2 = 0.0;

  // Components of the potential, field and field gradient.
  double E000; // Potential
  // l + m + n = 1 (3)   4
  double E100; // d/dX
  double E010; // d/dY
  double E001; // d/dz
  // l + m + n = 2 (6)  10
  double E200; // d^2/dXdX
  double E020; // d^2/dYdY
  double E002; // d^2/dZdZ
  double E110; // d^2/dXdY
  double E101; // d^2/dXdZ
  double E011; // d^2/dYdZ
  // l + m + n = 3 (10) 20
  double E300; // d^3/dXdXdX
  double E030; // d^3/dYdYdY
  double E003; // d^3/dZdZdZ
  double E210; // d^3/dXdXdY
  double E201; // d^3/dXdXdZ
  double E120; // d^3/dXdYdY
  double E021; // d^3/dYdYdZ
  double E102; // d^3/dXdZdZ
  double E012; // d^3/dYdZdZ
  double E111; // d^3/dXdYdZ

  // Cartesian tensor elements (for 1/R, erfc(Beta*R)/R or Thole damping.
  double R000;
  // l + m + n = 1 (3)   4
  double R100;
  double R010;
  double R001;
  // l + m + n = 2 (6)  10
  double R200;
  double R020;
  double R002;
  double R110;
  double R101;
  double R011;
  // l + m + n = 3 (10) 20
  double R300;
  double R030;
  double R003;
  double R210;
  double R201;
  double R120;
  double R021;
  double R102;
  double R012;
  double R111;
  // l + m + n = 4 (15) 35
  double R400;
  double R040;
  double R004;
  double R310;
  double R301;
  double R130;
  double R031;
  double R103;
  double R013;
  double R220;
  double R202;
  double R022;
  double R211;
  double R121;
  double R112;
  // l + m + n = 5 (21) 56
  double R500;
  double R050;
  double R005;
  double R410;
  double R401;
  double R140;
  double R041;
  double R104;
  double R014;
  double R320;
  double R302;
  double R230;
  double R032;
  double R203;
  double R023;
  double R311;
  double R131;
  double R113;
  double R221;
  double R212;
  double R122;
  // l + m + n = 6 (28) 84
  double R006;
  double R402;
  double R042;
  double R204;
  double R024;
  double R222;
  double R600;
  double R060;
  double R510;
  double R501;
  double R150;
  double R051;
  double R105;
  double R015;
  double R420;
  double R240;
  double R411;
  double R141;
  double R114;
  double R330;
  double R303;
  double R033;
  double R321;
  double R231;
  double R213;
  double R312;
  double R132;
  double R123;

  /**
   * Constructor for MultipoleTensor.
   *
   * @param order       The order of the tensor.
   * @param coordinates a {@link CoordinateSystem} object.
   */
  public MultipoleTensor(CoordinateSystem coordinates, int order) {
    assert (order > 0);
    o1 = order + 1;
    il = o1;
    im = il * o1;
    in = im * o1;
    size = (order + 1) * (order + 2) * (order + 3) / 6;
    work = new double[in * o1];

    this.order = order;
    this.coordinates = coordinates;
    this.operator = Operator.COULOMB;

    // Auxiliary terms for Coulomb and Thole Screening.
    coulombSource = new double[o1];
    for (short n = 0; n <= order; n++) {
      /*
       Math.pow(-1.0, j) returns positive for all j, with -1.0 as the //
       argument rather than -1. This is a bug?
       Challacombe Eq. 21, first two factors.
      */
      coulombSource[n] = pow(-1, n) * doubleFactorial(2 * n - 1);
    }

    T000 = new double[order + 1];
    // l + m + n = 0 (1)
    t000 = MultipoleUtilities.ti(0, 0, 0, order);
    // l + m + n = 1 (3)   4
    t100 = MultipoleUtilities.ti(1, 0, 0, order);
    t010 = MultipoleUtilities.ti(0, 1, 0, order);
    t001 = MultipoleUtilities.ti(0, 0, 1, order);
    // l + m + n = 2 (6)  10
    t200 = MultipoleUtilities.ti(2, 0, 0, order);
    t020 = MultipoleUtilities.ti(0, 2, 0, order);
    t002 = MultipoleUtilities.ti(0, 0, 2, order);
    t110 = MultipoleUtilities.ti(1, 1, 0, order);
    t101 = MultipoleUtilities.ti(1, 0, 1, order);
    t011 = MultipoleUtilities.ti(0, 1, 1, order);
    // l + m + n = 3 (10) 20
    t300 = MultipoleUtilities.ti(3, 0, 0, order);
    t030 = MultipoleUtilities.ti(0, 3, 0, order);
    t003 = MultipoleUtilities.ti(0, 0, 3, order);
    t210 = MultipoleUtilities.ti(2, 1, 0, order);
    t201 = MultipoleUtilities.ti(2, 0, 1, order);
    t120 = MultipoleUtilities.ti(1, 2, 0, order);
    t021 = MultipoleUtilities.ti(0, 2, 1, order);
    t102 = MultipoleUtilities.ti(1, 0, 2, order);
    t012 = MultipoleUtilities.ti(0, 1, 2, order);
    t111 = MultipoleUtilities.ti(1, 1, 1, order);
    // l + m + n = 4 (15) 35
    t400 = MultipoleUtilities.ti(4, 0, 0, order);
    t040 = MultipoleUtilities.ti(0, 4, 0, order);
    t004 = MultipoleUtilities.ti(0, 0, 4, order);
    t310 = MultipoleUtilities.ti(3, 1, 0, order);
    t301 = MultipoleUtilities.ti(3, 0, 1, order);
    t130 = MultipoleUtilities.ti(1, 3, 0, order);
    t031 = MultipoleUtilities.ti(0, 3, 1, order);
    t103 = MultipoleUtilities.ti(1, 0, 3, order);
    t013 = MultipoleUtilities.ti(0, 1, 3, order);
    t220 = MultipoleUtilities.ti(2, 2, 0, order);
    t202 = MultipoleUtilities.ti(2, 0, 2, order);
    t022 = MultipoleUtilities.ti(0, 2, 2, order);
    t211 = MultipoleUtilities.ti(2, 1, 1, order);
    t121 = MultipoleUtilities.ti(1, 2, 1, order);
    t112 = MultipoleUtilities.ti(1, 1, 2, order);
    // l + m + n = 5 (21) 56
    t500 = MultipoleUtilities.ti(5, 0, 0, order);
    t050 = MultipoleUtilities.ti(0, 5, 0, order);
    t005 = MultipoleUtilities.ti(0, 0, 5, order);
    t410 = MultipoleUtilities.ti(4, 1, 0, order);
    t401 = MultipoleUtilities.ti(4, 0, 1, order);
    t140 = MultipoleUtilities.ti(1, 4, 0, order);
    t041 = MultipoleUtilities.ti(0, 4, 1, order);
    t104 = MultipoleUtilities.ti(1, 0, 4, order);
    t014 = MultipoleUtilities.ti(0, 1, 4, order);
    t320 = MultipoleUtilities.ti(3, 2, 0, order);
    t302 = MultipoleUtilities.ti(3, 0, 2, order);
    t230 = MultipoleUtilities.ti(2, 3, 0, order);
    t032 = MultipoleUtilities.ti(0, 3, 2, order);
    t203 = MultipoleUtilities.ti(2, 0, 3, order);
    t023 = MultipoleUtilities.ti(0, 2, 3, order);
    t311 = MultipoleUtilities.ti(3, 1, 1, order);
    t131 = MultipoleUtilities.ti(1, 3, 1, order);
    t113 = MultipoleUtilities.ti(1, 1, 3, order);
    t221 = MultipoleUtilities.ti(2, 2, 1, order);
    t212 = MultipoleUtilities.ti(2, 1, 2, order);
    t122 = MultipoleUtilities.ti(1, 2, 2, order);
    // l + m + n = 6 (28) 84
    t600 = MultipoleUtilities.ti(6, 0, 0, order);
    t060 = MultipoleUtilities.ti(0, 6, 0, order);
    t006 = MultipoleUtilities.ti(0, 0, 6, order);
    t510 = MultipoleUtilities.ti(5, 1, 0, order);
    t501 = MultipoleUtilities.ti(5, 0, 1, order);
    t150 = MultipoleUtilities.ti(1, 5, 0, order);
    t051 = MultipoleUtilities.ti(0, 5, 1, order);
    t105 = MultipoleUtilities.ti(1, 0, 5, order);
    t015 = MultipoleUtilities.ti(0, 1, 5, order);
    t420 = MultipoleUtilities.ti(4, 2, 0, order);
    t402 = MultipoleUtilities.ti(4, 0, 2, order);
    t240 = MultipoleUtilities.ti(2, 4, 0, order);
    t042 = MultipoleUtilities.ti(0, 4, 2, order);
    t204 = MultipoleUtilities.ti(2, 0, 4, order);
    t024 = MultipoleUtilities.ti(0, 2, 4, order);
    t411 = MultipoleUtilities.ti(4, 1, 1, order);
    t141 = MultipoleUtilities.ti(1, 4, 1, order);
    t114 = MultipoleUtilities.ti(1, 1, 4, order);
    t330 = MultipoleUtilities.ti(3, 3, 0, order);
    t303 = MultipoleUtilities.ti(3, 0, 3, order);
    t033 = MultipoleUtilities.ti(0, 3, 3, order);
    t321 = MultipoleUtilities.ti(3, 2, 1, order);
    t231 = MultipoleUtilities.ti(2, 3, 1, order);
    t213 = MultipoleUtilities.ti(2, 1, 3, order);
    t312 = MultipoleUtilities.ti(3, 1, 2, order);
    t132 = MultipoleUtilities.ti(1, 3, 2, order);
    t123 = MultipoleUtilities.ti(1, 2, 3, order);
    t222 = MultipoleUtilities.ti(2, 2, 2, order);
  }

  /**
   * getd2EdZ2.
   *
   * @return a double.
   */
  public double getd2EdZ2() {
    return d2EdZ2;
  }

  /**
   * getdEdZ.
   *
   * @return a double.
   */
  public double getdEdZ() {
    return dEdZ;
  }

  /**
   * log.
   *
   * @param tensor an array of double values.
   */
  public void log(double[] tensor) {
    log(this.operator, this.order, tensor);
  }

  /**
   * Permanent multipole energy.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return a double.
   */
  public double multipoleEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    multipoleIPotentialAtK(mI, 2);
    return multipoleEnergy(mK);
  }

  /**
   * Permanent multipole energy and gradient.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @param Gi Coordinate gradient at site I.
   * @param Gk Coordinate gradient at site K.
   * @param Ti Torque at site I.
   * @param Tk Torque at site K.
   * @return the permanent multipole energy.
   */
  public double multipoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
                                           double[] Gi, double[] Gk, double[] Ti, double[] Tk) {
    multipoleIPotentialAtK(mI, 3);
    double energy = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    Gi[0] = -Gk[0];
    Gi[1] = -Gk[1];
    Gi[2] = -Gk[2];

    // Torques
    multipoleTorque(mK, Tk);
    multipoleKPotentialAtI(mK, 2);
    multipoleTorque(mI, Ti);

    // dEdZ = -Gi[2];
    // if (order >= 6) {
    //  multipoleIdZ2(mI);
    //  d2EdZ2 = dotMultipole(mK);
    // }

    return energy;
  }

  /**
   * Polarization Energy.
   *
   * @param mI          PolarizableMultipole at site I.
   * @param mK          PolarizableMultipole at site K.
   * @param scaleEnergy a double.
   * @return a double.
   */
  public double polarizationEnergy(PolarizableMultipole mI, PolarizableMultipole mK,
                                   double scaleEnergy) {
    // Incorporate core charges into multipole potential
    if (mI.Z != 0 && mK.Z != 0 && operator == Operator.THOLE_DIRECT_FIELD) {
      mI.q += mI.Z;
      mK.q += mK.Z;
    }

    // Find the permanent multipole potential and derivatives at k.
    multipoleIPotentialAtK(mI, 1);
    // Energy of induced dipole k in the field of permanent multipole i.
    double eK = polarizationEnergy(mK);
    // Find the permanent multipole potential and derivatives at site i.
    multipoleKPotentialAtI(mK, 1);
    // Energy of induced dipole i in the field of permanent multipole k.
    double eI = polarizationEnergy(mI);

    if (mI.Z != 0 && mK.Z != 0 && operator == Operator.THOLE_DIRECT_FIELD) {
      mI.q -= mI.Z;
      mK.q -= mK.Z;
    }

    return scaleEnergy * (eI + eK);
  }

  /**
   * Polarization Energy and Gradient.
   *
   * @param mI            PolarizableMultipole at site I.
   * @param mK            PolarizableMultipole at site K.
   * @param inductionMask a double.
   * @param energyMask    a double.
   * @param mutualMask    a double.
   * @param Gi            an array of double values.
   * @param Ti            an array of double values.
   * @param Tk            an array of double values.
   * @return a double.
   */
  public double polarizationEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
                                              double inductionMask, double energyMask, double mutualMask,
                                              double[] Gi, double[] Ti, double[] Tk) {

    if (mI.Z != 0 && mK.Z != 0) {
      mI.q += mI.Z;
      mK.q += mK.Z;
    }

    // Add the induction and energy masks to create an "averaged" induced dipole (sx, sy, sz).
    mI.applyMasks(inductionMask, energyMask);
    mK.applyMasks(inductionMask, energyMask);

    // Find the permanent multipole potential and derivatives at site k.
    multipoleIPotentialAtK(mI, 2);
    // Energy of induced dipole k in the field of multipole i.
    // The field E_x = -E100.
    double eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = -(mK.sx * E200 + mK.sy * E110 + mK.sz * E101);
    Gi[1] = -(mK.sx * E110 + mK.sy * E020 + mK.sz * E011);
    Gi[2] = -(mK.sx * E101 + mK.sy * E011 + mK.sz * E002);

    // Find the permanent multipole potential and derivatives at site i.
    multipoleKPotentialAtI(mK, 2);
    // Energy of induced dipole i in the field of multipole k.
    double eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] += (mI.sx * E200 + mI.sy * E110 + mI.sz * E101);
    Gi[1] += (mI.sx * E110 + mI.sy * E020 + mI.sz * E011);
    Gi[2] += (mI.sx * E101 + mI.sy * E011 + mI.sz * E002);

    // Total polarization energy.
    double energy = energyMask * (eI + eK);

    // Get the induced-induced portion of the force (Ud . dC/dX . Up).
    // This contribution does not exist for direct polarization (mutualMask == 0.0).
    if (mutualMask != 0.0) {
      // Find the potential and its derivatives at k due to induced dipole i.
      dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
      Gi[0] -= 0.5 * mutualMask * (mK.px * E200 + mK.py * E110 + mK.pz * E101);
      Gi[1] -= 0.5 * mutualMask * (mK.px * E110 + mK.py * E020 + mK.pz * E011);
      Gi[2] -= 0.5 * mutualMask * (mK.px * E101 + mK.py * E011 + mK.pz * E002);

      // Find the potential and its derivatives at i due to induced dipole k.
      dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
      Gi[0] += 0.5 * mutualMask * (mI.px * E200 + mI.py * E110 + mI.pz * E101);
      Gi[1] += 0.5 * mutualMask * (mI.px * E110 + mI.py * E020 + mI.pz * E011);
      Gi[2] += 0.5 * mutualMask * (mI.px * E101 + mI.py * E011 + mI.pz * E002);
    }

    // Find the potential and its derivatives at K due to the averaged induced dipole at site i.
    dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 2);
    multipoleTorque(mK, Tk);

    // Find the potential and its derivatives at I due to the averaged induced dipole at site k.
    dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 2);
    multipoleTorque(mI, Ti);

    if (mI.Z != 0 && mK.Z != 0) {
      mI.q -= mI.Z;
      mK.q -= mK.Z;
    }

    return energy;
  }

  /**
   * Permanent Multipole + Polarization Energy.
   *
   * @param mI               PolarizableMultipole at site I.
   * @param mK               PolarizableMultipole at site K.
   * @param scaleEnergy      a double.
   * @param energyComponents The permanent and induced energy components.
   * @return a double.
   */
  public double totalEnergy(PolarizableMultipole mI, PolarizableMultipole mK, double scaleEnergy,
                            double[] energyComponents) {
    multipoleIPotentialAtK(mI, 2);
    double permanentEnergy = multipoleEnergy(mK);

    // Energy of multipole k in the field of induced dipole i.
    // The field E_x = -E100.
    double eK = -(mK.ux * E100 + mK.uy * E010 + mK.uz * E001);

    // Find the permanent multipole potential and derivatives at site i.
    multipoleKPotentialAtI(mK, 2);
    // Energy of multipole i in the field of induced dipole k.
    double eI = -(mI.ux * E100 + mI.uy * E010 + mI.uz * E001);

    double inducedEnergy = -0.5 * scaleEnergy * (eI + eK);

    // Story the energy components.
    energyComponents[0] = permanentEnergy;
    energyComponents[1] = inducedEnergy;

    return permanentEnergy + inducedEnergy;
  }

  /**
   * Set the separation vector.
   *
   * @param r The separation vector.
   */
  public final void setR(double[] r) {
    setR(r[0], r[1], r[2]);
  }

  /**
   * Set the separation vector.
   *
   * @param dx Separation along the X-axis.
   * @param dy Separation along the Y-axis.
   * @param dz Separation along the Z-axis.
   */
  public abstract void setR(double dx, double dy, double dz);

  /**
   * For the MultipoleTensorTest class and testing.
   *
   * @param r an array of double values.
   */
  public final void generateTensor(double[] r) {
    setR(r);
    generateTensor();
  }

  /**
   * Generate the tensor using hard-coded methods or via recursion.
   */
  public void generateTensor() {
    switch (order) {
      case 1 -> order1();
      case 2 -> order2();
      case 3 -> order3();
      case 4 -> order4();
      case 5 -> order5();
      case 6 -> order6();
      default -> {
        double[] r = {x, y, z};
        recursion(r, work);
      }
    }
  }

  /**
   * Generate source terms for the Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  protected abstract void source(double[] T000);

  /**
   * The index is based on the idea of filling tetrahedron.
   * <p>
   * 1/r has an index of 0.
   * <br>
   * derivatives of x are first; indices from 1..o for d/dx..(d/dx)^o
   * <br>
   * derivatives of x and y are second; base triangle of size (o+1)(o+2)/2
   * <br>
   * derivatives of x, y and z are last; total size (o+1)*(o+2)*(o+3)/6
   * <br>
   * <p>
   * This function is useful to set up masking constants:
   * <br>
   * static int Tlmn = ti(l,m,n,order)
   * <br>
   * For example the (d/dy)^2 (1/R) storage location:
   * <br>
   * static int T020 = ti(0,2,0,order)
   *
   * @param dx int The number of d/dx operations.
   * @param dy int The number of d/dy operations.
   * @param dz int The number of d/dz operations.
   * @return int in the range (0..binomial(order + 3, 3) - 1)
   */
  protected final int ti(int dx, int dy, int dz) {
    return MultipoleUtilities.ti(dx, dy, dz, order);
  }

  /**
   * Contract multipole moments with their respective electrostatic potential derivatives.
   *
   * @param mI PolarizableMultipole at site I.
   * @param T  array of electrostatic potential and partial derivatives
   * @param l  apply (d/dx)^l to the potential
   * @param m  apply (d/dy)^l to the potential
   * @param n  apply (d/dz)^l to the potential
   * @return the contracted interaction.
   */
  protected double contractMultipoleI(PolarizableMultipole mI, double[] T, int l, int m, int n) {
    double total = 0.0;
    total += mI.q * T[ti(l, m, n)];
    total -= mI.dx * T[ti(l + 1, m, n)];
    total -= mI.dy * T[ti(l, m + 1, n)];
    total -= mI.dz * T[ti(l, m, n + 1)];
    total += mI.qxx * T[ti(l + 2, m, n)];
    total += mI.qyy * T[ti(l, m + 2, n)];
    total += mI.qzz * T[ti(l, m, n + 2)];
    total += mI.qxy * T[ti(l + 1, m + 1, n)];
    total += mI.qxz * T[ti(l + 1, m, n + 1)];
    total += mI.qyz * T[ti(l, m + 1, n + 1)];
    return total;
  }

  /**
   * Collect the field at R due to Q multipole moments at the origin (site I).
   *
   * @param mI PolarizableMultipole at site I.
   * @param T  Electrostatic potential and partial derivatives
   * @param l  apply (d/dx)^l to the potential
   * @param m  apply (d/dy)^l to the potential
   * @param n  apply (d/dz)^l to the potential
   */
  protected final void potentialMultipoleI(PolarizableMultipole mI, double[] T, int l, int m,
                                           int n) {
    E000 = contractMultipoleI(mI, T, l, m, n);
    E100 = contractMultipoleI(mI, T, l + 1, m, n);
    E010 = contractMultipoleI(mI, T, l, m + 1, n);
    E001 = contractMultipoleI(mI, T, l, m, n + 1);
    E200 = contractMultipoleI(mI, T, l + 2, m, n);
    E020 = contractMultipoleI(mI, T, l, m + 2, n);
    E002 = contractMultipoleI(mI, T, l, m, n + 2);
    E110 = contractMultipoleI(mI, T, l + 1, n + 1, m);
    E101 = contractMultipoleI(mI, T, l + 1, m, n + 1);
    E011 = contractMultipoleI(mI, T, l, m + 1, n + 1);
  }

  /**
   * <p>This function is a driver to collect elements of the Cartesian multipole tensor. Collecting
   * all tensors scales slightly better than O(order^4).
   *
   * <p>For a multipole expansion truncated at quadrupole order, for example, up to order 5 is
   * needed for energy gradients. The number of terms this requires is binomial(5 + 3, 3) or 8! / (5!
   * * 3!), which is 56.
   *
   * <p>The packing of the tensor elements for order = 1<br>
   * tensor[0] = 1/|r| <br> tensor[1] = -x/|r|^3 <br> tensor[2] = -y/|r|^3 <br> tensor[3] = -z/|r|^3
   * <br>
   *
   * @param r      The separation vector.
   * @param tensor The tensor elements.
   * @return The tensor recursion code.
   * @since 1.0
   */
  protected abstract String codeTensorRecursion(final double[] r, final double[] tensor);

  /**
   * Contract multipole moments with their respective electrostatic potential derivatives.
   *
   * @param mI PolarizableMultipole at site I.
   * @param T  array of electrostatic potential and partial derivatives
   * @param l  apply (d/dx)^l to the potential
   * @param m  apply (d/dy)^l to the potential
   * @param n  apply (d/dz)^l to the potential
   * @param sb the code will be appended to the StringBuilder.
   * @return the contracted interaction.
   */
  private double codeContractMultipoleI(PolarizableMultipole mI, double[] T, int l, int m, int n,
                                        StringBuilder sb) {
    double total = 0.0;
    String name = term(l, m, n);
    sb.append(format("double %s = 0.0;\n", name));
    StringBuilder sb1 = new StringBuilder();
    double term = mI.q * T[ti(l, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mI.q, R%s, %s);\n", name, lmn(l, m, n), name));
    }
    term = mI.dx * T[ti(l + 1, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mI.dx, -R%s, %s);\n", name, lmn(l + 1, m, n), name));
    }
    term = mI.dy * T[ti(l, m + 1, n)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mI.dy, -R%s, %s);\n", name, lmn(l, m + 1, n), name));
    }
    term = mI.dz * T[ti(l, m, n + 1)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mI.dz, -R%s, %s);\n", name, lmn(l, m, n + 1), name));
    }
    StringBuilder traceSB = new StringBuilder();
    double trace = 0.0;
    term = mI.qxx * T[ti(l + 2, m, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(format("%s = fma(mI.qxx, R%s, %s);\n", name, lmn(l + 2, m, n), name));
    }
    term = mI.qyy * T[ti(l, m + 2, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(format("%s = fma(mI.qyy, R%s, %s);\n", name, lmn(l, m + 2, n), name));
    }
    term = mI.qzz * T[ti(l, m, n + 2)];
    if (term != 0) {
      trace += term;
      traceSB.append(format("%s = fma(mI.qzz, R%s, %s);\n", name, lmn(l, m, n + 2), name));
    }
    total += trace;
    if (total != 0) {
      sb.append(sb1);
      if (trace != 0) {
        sb.append(traceSB);
      }
    }
    term = mI.qxy * T[ti(l + 1, m + 1, n)];
    if (term != 0) {
      total += term;
      sb.append(format("%s = fma(mI.qxy, R%s, %s);\n", name, lmn(l + 1, m + 1, n), name));
    }
    term = mI.qxz * T[ti(l + 1, m, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(format("%s = fma(mI.qxz, R%s, %s);\n", name, lmn(l + 1, m, n + 1), name));
    }
    term = mI.qyz * T[ti(l, m + 1, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(format("%s = fma(mI.qyz, R%s, %s);\n", name, lmn(l, m + 1, n + 1), name));
    }
    return total;
  }


  /**
   * Contract multipole moments with their respective electrostatic potential derivatives
   * using SIMD instructions.
   *
   * @param mI PolarizableMultipole at site I.
   * @param T  array of electrostatic potential and partial derivatives
   * @param l  apply (d/dx)^l to the potential
   * @param m  apply (d/dy)^l to the potential
   * @param n  apply (d/dz)^l to the potential
   * @param sb the code will be appended to the StringBuilder.
   * @return the contracted interaction.
   */
  private double codeContractMultipoleISIMD(PolarizableMultipole mI, double[] T, int l, int m, int n,
                                            StringBuilder sb, HashMap<Integer, String> tensorMap) {
    double total = 0.0;
    String name = term(l, m, n);
    StringBuilder sb1 = new StringBuilder();
    double term = mI.q * T[ti(l, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l, m, n, tensorMap));
      sb1.append(format("\tDoubleVector %s = q.mul(t%s);\n", name, lmn(l, m, n)));
    } else {
      sb.append(format("\tDoubleVector %s = zero;\n", name));
    }
    term = mI.dx * T[ti(l + 1, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l + 1, m, n, tensorMap));
      sb1.append(format("\t%s = dx.fma(t%s.neg(), %s);\n", name, lmn(l + 1, m, n), name));
    }
    term = mI.dy * T[ti(l, m + 1, n)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l, m + 1, n, tensorMap));
      sb1.append(format("\t%s = dy.fma(t%s.neg(), %s);\n", name, lmn(l, m + 1, n), name));
    }
    term = mI.dz * T[ti(l, m, n + 1)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l, m, n + 1, tensorMap));
      sb1.append(format("\t%s = dz.fma(t%s.neg(), %s);\n", name, lmn(l, m, n + 1), name));
    }
    StringBuilder traceSB = new StringBuilder();
    double trace = 0.0;
    term = mI.qxx * T[ti(l + 2, m, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(loadTensor(l + 2, m, n, tensorMap));
      traceSB.append(format("\t%s = qxx.fma(t%s, %s);\n", name, lmn(l + 2, m, n), name));
    }
    term = mI.qyy * T[ti(l, m + 2, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(loadTensor(l, m + 2, n, tensorMap));
      traceSB.append(format("\t%s = qyy.fma(t%s, %s);\n", name, lmn(l, m + 2, n), name));
    }
    term = mI.qzz * T[ti(l, m, n + 2)];
    if (term != 0) {
      trace += term;
      traceSB.append(loadTensor(l, m, n + 2, tensorMap));
      traceSB.append(format("\t%s = qzz.fma(t%s, %s);\n", name, lmn(l, m, n + 2), name));
    }
    total += trace;
    if (total != 0) {
      sb.append(sb1);
      if (trace != 0) {
        sb.append(traceSB);
      }
    }
    term = mI.qxy * T[ti(l + 1, m + 1, n)];
    if (term != 0) {
      total += term;
      sb.append(loadTensor(l + 1, m + 1, n, tensorMap));
      sb.append(format("\t%s = qxy.fma(t%s, %s);\n", name, lmn(l + 1, m + 1, n), name));
    }
    term = mI.qxz * T[ti(l + 1, m, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(loadTensor(l + 1, m, n + 1, tensorMap));
      sb.append(format("\t%s = qxz.fma(t%s, %s);\n", name, lmn(l + 1, m, n + 1), name));
    }
    term = mI.qyz * T[ti(l, m + 1, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(loadTensor(l, m + 1, n + 1, tensorMap));
      sb.append(format("\t%s = qyz.fma(t%s, %s);\n", name, lmn(l, m + 1, n + 1), name));
    }
    return total;
  }

  /**
   * Contract multipole moments with their respective electrostatic potential derivatives.
   *
   * @param mK PolarizableMultipole at site K.
   * @param T  array of electrostatic potential and partial derivatives
   * @param l  apply (d/dx)^l to the potential
   * @param m  apply (d/dy)^l to the potential
   * @param n  apply (d/dz)^l to the potential
   * @param sb the code will be appended to the StringBuilder.
   * @return the contracted interaction.
   */
  private double codeContractMultipoleK(PolarizableMultipole mK, double[] T, int l, int m, int n,
                                        StringBuilder sb) {
    double total = 0.0;
    String name = term(l, m, n);
    sb.append(format("double %s = 0.0;\n", name));
    StringBuilder sb1 = new StringBuilder();
    double term = mK.q * T[ti(l, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mK.q, R%s, %s);\n", name, lmn(l, m, n), name));
    }
    term = mK.dx * T[ti(l + 1, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mK.dx, R%s, %s);\n", name, lmn(l + 1, m, n), name));
    }
    term = mK.dy * T[ti(l, m + 1, n)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mK.dy, R%s, %s);\n", name, lmn(l, m + 1, n), name));
    }
    term = mK.dz * T[ti(l, m, n + 1)];
    if (term != 0) {
      total += term;
      sb1.append(format("%s = fma(mK.dz, R%s, %s);\n", name, lmn(l, m, n + 1), name));
    }
    StringBuilder traceSB = new StringBuilder();
    double trace = 0.0;
    term = mK.qxx * T[ti(l + 2, m, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(format("%s = fma(mK.qxx, R%s, %s);\n", name, lmn(l + 2, m, n), name));
    }
    term = mK.qyy * T[ti(l, m + 2, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(format("%s = fma(mK.qyy, R%s, %s);\n", name, lmn(l, m + 2, n), name));
    }
    term = mK.qzz * T[ti(l, m, n + 2)];
    if (term != 0) {
      trace += term;
      traceSB.append(format("%s = fma(mK.qzz, R%s, %s);\n", name, lmn(l, m, n + 2), name));
    }
    total += trace;
    if (total != 0) {
      sb.append(sb1);
      if (trace != 0) {
        sb.append(traceSB);
      }
    }
    term = mK.qxy * T[ti(l + 1, m + 1, n)];
    if (term != 0) {
      total += term;
      sb.append(format("%s = fma(mK.qxy, R%s, %s);\n", name, lmn(l + 1, m + 1, n), name));
    }
    term = mK.qxz * T[ti(l + 1, m, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(format("%s = fma(mK.qxz, R%s, %s);\n", name, lmn(l + 1, m, n + 1), name));
    }
    term = mK.qyz * T[ti(l, m + 1, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(format("%s = fma(mK.qyz, R%s, %s);\n", name, lmn(l, m + 1, n + 1), name));
    }
    return total;
  }

  /**
   * Contract multipole moments with their respective electrostatic potential derivatives
   * using SIMD instructions.
   *
   * @param mK PolarizableMultipole at site K.
   * @param T  array of electrostatic potential and partial derivatives
   * @param l  apply (d/dx)^l to the potential
   * @param m  apply (d/dy)^l to the potential
   * @param n  apply (d/dz)^l to the potential
   * @param sb the code will be appended to the StringBuilder.
   * @return the contracted interaction.
   */
  private double codeContractMultipoleKSIMD(PolarizableMultipole mK, double[] T, int l, int m, int n,
                                            StringBuilder sb, HashMap<Integer, String> tensorHash) {
    double total = 0.0;
    String name = term(l, m, n);
    StringBuilder sb1 = new StringBuilder();
    double term = mK.q * T[ti(l, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l, m, n, tensorHash));
      sb1.append(format("\tDoubleVector %s = q.mul(t%s);\n", name, lmn(l, m, n)));
    } else {
      sb.append(format("\tDoubleVector %s = zero;\n", name));
    }
    term = mK.dx * T[ti(l + 1, m, n)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l + 1, m, n, tensorHash));
      sb1.append(format("\t%s = dx.fma(t%s, %s);\n", name, lmn(l + 1, m, n), name));
    }
    term = mK.dy * T[ti(l, m + 1, n)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l, m + 1, n, tensorHash));
      sb1.append(format("\t%s = dy.fma(t%s, %s);\n", name, lmn(l, m + 1, n), name));
    }
    term = mK.dz * T[ti(l, m, n + 1)];
    if (term != 0) {
      total += term;
      sb1.append(loadTensor(l, m, n + 1, tensorHash));
      sb1.append(format("\t%s = dz.fma(t%s, %s);\n", name, lmn(l, m, n + 1), name));
    }
    StringBuilder traceSB = new StringBuilder();
    double trace = 0.0;
    term = mK.qxx * T[ti(l + 2, m, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(loadTensor(l + 2, m, n, tensorHash));
      traceSB.append(format("\t%s = qxx.fma(t%s, %s);\n", name, lmn(l + 2, m, n), name));
    }
    term = mK.qyy * T[ti(l, m + 2, n)];
    if (term != 0) {
      trace += term;
      traceSB.append(loadTensor(l, m + 2, n, tensorHash));
      traceSB.append(format("\t%s = qyy.fma(t%s, %s);\n", name, lmn(l, m + 2, n), name));
    }
    term = mK.qzz * T[ti(l, m, n + 2)];
    if (term != 0) {
      trace += term;
      traceSB.append(loadTensor(l, m, n + 2, tensorHash));
      traceSB.append(format("\t%s = qzz.fma(t%s, %s);\n", name, lmn(l, m, n + 2), name));
    }
    total += trace;
    if (total != 0) {
      sb.append(sb1);
      if (trace != 0) {
        sb.append(traceSB);
      }
    }
    term = mK.qxy * T[ti(l + 1, m + 1, n)];
    if (term != 0) {
      total += term;
      sb.append(loadTensor(l + 1, m + 1, n, tensorHash));
      sb.append(format("\t%s = qxy.fma(t%s, %s);\n", name, lmn(l + 1, m + 1, n), name));
    }
    term = mK.qxz * T[ti(l + 1, m, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(loadTensor(l + 1, m, n + 1, tensorHash));
      sb.append(format("\t%s = qxz.fma(t%s, %s);\n", name, lmn(l + 1, m, n + 1), name));
    }
    term = mK.qyz * T[ti(l, m + 1, n + 1)];
    if (term != 0) {
      total += term;
      sb.append(loadTensor(l, m + 1, n + 1, tensorHash));
      sb.append(format("\t%s = qyz.fma(t%s, %s);\n", name, lmn(l, m + 1, n + 1), name));
    }
    return total;
  }

  /**
   * Collect the potential its partial derivatives at K due to multipole moments at the origin.
   *
   * @param mI PolarizableMultipole at site I.
   * @param T  Electrostatic potential and partial derivatives.
   * @param l  apply (d/dx)^l to the potential.
   * @param m  apply (d/dy)^l to the potential.
   * @param n  apply (d/dz)^l to the potential.
   * @param sb Append the code to the StringBuilder.
   */
  protected void codePotentialMultipoleI(PolarizableMultipole mI, double[] T, int l, int m, int n, StringBuilder sb) {
    E000 = codeContractMultipoleI(mI, T, l, m, n, sb);
    if (E000 != 0) {
      sb.append(format("E000 = %s;\n", term(l, m, n)));
    }
    // Order 1
    E100 = codeContractMultipoleI(mI, T, l + 1, m, n, sb);
    if (E100 != 0) {
      sb.append(format("E100 = %s;\n", term(l + 1, m, n)));
    }
    E010 = codeContractMultipoleI(mI, T, l, m + 1, n, sb);
    if (E100 != 0) {
      sb.append(format("E010 = %s;\n", term(l, m + 1, n)));
    }
    E001 = codeContractMultipoleI(mI, T, l, m, n + 1, sb);
    if (E001 != 0) {
      sb.append(format("E001 = %s;\n", term(l, m, n + 1)));
    }
    // Order 2
    E200 = codeContractMultipoleI(mI, T, l + 2, m, n, sb);
    if (E200 != 0) {
      sb.append(format("E200 = %s;\n", term(l + 2, m, n)));
    }
    E020 = codeContractMultipoleI(mI, T, l, m + 2, n, sb);
    if (E020 != 0) {
      sb.append(format("E020 = %s;\n", term(l, m + 2, n)));
    }
    E002 = codeContractMultipoleI(mI, T, l, m, n + 2, sb);
    if (E002 != 0) {
      sb.append(format("E002 = %s;\n", term(l, m, n + 2)));
    }
    E110 = codeContractMultipoleI(mI, T, l + 1, m + 1, n, sb);
    if (E110 != 0) {
      sb.append(format("E110 = %s;\n", term(l + 1, m + 1, n)));
    }
    E101 = codeContractMultipoleI(mI, T, l + 1, m, n + 1, sb);
    if (E101 != 0) {
      sb.append(format("E101 = %s;\n", term(l + 1, m, n + 1)));
    }
    E011 = codeContractMultipoleI(mI, T, l, m + 1, n + 1, sb);
    if (E011 != 0) {
      sb.append(format("E011 = %s;\n", term(l, m + 1, n + 1)));
    }
    // Order 3
    E300 = codeContractMultipoleI(mI, T, l + 3, m, n, sb);
    if (E300 != 0) {
      sb.append(format("E300 = %s;\n", term(l + 3, m, n)));
    }
    E030 = codeContractMultipoleI(mI, T, l, m + 3, n, sb);
    if (E030 != 0) {
      sb.append(format("E030 = %s;\n", term(l, m + 3, n)));
    }
    E003 = codeContractMultipoleI(mI, T, l, m, n + 3, sb);
    if (E003 != 0) {
      sb.append(format("E003 = %s;\n", term(l, m, n + 3)));
    }
    E210 = codeContractMultipoleI(mI, T, l + 2, m + 1, n, sb);
    if (E210 != 0) {
      sb.append(format("E210 = %s;\n", term(l + 2, m + 1, n)));
    }
    E201 = codeContractMultipoleI(mI, T, l + 2, m, n + 1, sb);
    if (E201 != 0) {
      sb.append(format("E201 = %s;\n", term(l + 2, m, n + 1)));
    }
    E120 = codeContractMultipoleI(mI, T, l + 1, m + 2, n, sb);
    if (E120 != 0) {
      sb.append(format("E120 = %s;\n", term(l + 1, m + 2, n)));
    }
    E021 = codeContractMultipoleI(mI, T, l, m + 2, n + 1, sb);
    if (E021 != 0) {
      sb.append(format("E021 = %s;\n", term(l, m + 2, n + 1)));
    }
    E102 = codeContractMultipoleI(mI, T, l + 1, m, n + 2, sb);
    if (E102 != 0) {
      sb.append(format("E102 = %s;\n", term(l + 1, m, n + 2)));
    }
    E012 = codeContractMultipoleI(mI, T, l, m + 1, n + 2, sb);
    if (E012 != 0) {
      sb.append(format("E012 = %s;\n", term(l, m + 1, n + 2)));
    }
    E111 = codeContractMultipoleI(mI, T, l + 1, m + 1, n + 1, sb);
    if (E111 != 0) {
      sb.append(format("E111 = %s;\n", term(l + 1, m + 1, n + 1)));
    }
  }

  /**
   * Collect the potential its partial derivatives at K due to multipole moments at the origin
   * using SIMD instructions.
   *
   * @param mI PolarizableMultipole at site I.
   * @param T  Electrostatic potential and partial derivatives.
   * @param l  apply (d/dx)^l to the potential.
   * @param m  apply (d/dy)^l to the potential.
   * @param n  apply (d/dz)^l to the potential.
   * @param sb Append the code to the StringBuilder.
   */
  protected void codePotentialMultipoleISIMD(PolarizableMultipole mI, double[] T, int l, int m, int n, StringBuilder sb) {
    String to = "e";
    HashMap<Integer, String> tensorHash = new HashMap<>();
    sb.append("\n// Order 0\n");
    E000 = codeContractMultipoleISIMD(mI, T, l, m, n, sb, tensorHash);
    if (E000 != 0) {
      sb.append(storePotential(to, l, m, n));
    }

    // Order 1
    sb.append("\n// Order 1\n");
    E100 = codeContractMultipoleISIMD(mI, T, l + 1, m, n, sb, tensorHash);
    if (E100 != 0) {
      sb.append(storePotential(to, l + 1, m, n));
    }
    E010 = codeContractMultipoleISIMD(mI, T, l, m + 1, n, sb, tensorHash);
    if (E100 != 0) {
      sb.append(storePotential(to, l, m + 1, n));
    }
    E001 = codeContractMultipoleISIMD(mI, T, l, m, n + 1, sb, tensorHash);
    if (E001 != 0) {
      sb.append(storePotential(to, l, m, n + 1));
    }

    // Order 2
    sb.append("\n// Order 2\n");
    E200 = codeContractMultipoleISIMD(mI, T, l + 2, m, n, sb, tensorHash);
    if (E200 != 0) {
      sb.append(storePotential(to, l + 2, m, n));
    }
    E020 = codeContractMultipoleISIMD(mI, T, l, m + 2, n, sb, tensorHash);
    if (E020 != 0) {
      sb.append(storePotential(to, l, m + 2, n));
    }
    E002 = codeContractMultipoleISIMD(mI, T, l, m, n + 2, sb, tensorHash);
    if (E002 != 0) {
      sb.append(storePotential(to, l, m, n + 2));
    }
    E110 = codeContractMultipoleISIMD(mI, T, l + 1, m + 1, n, sb, tensorHash);
    if (E110 != 0) {
      sb.append(storePotential(to, l + 1, m + 1, n));
    }
    E101 = codeContractMultipoleISIMD(mI, T, l + 1, m, n + 1, sb, tensorHash);
    if (E101 != 0) {
      sb.append(storePotential(to, l + 1, m, n + 1));
    }
    E011 = codeContractMultipoleISIMD(mI, T, l, m + 1, n + 1, sb, tensorHash);
    if (E011 != 0) {
      sb.append(storePotential(to, l, m + 1, n + 1));
    }

    // Order 3
    sb.append("\n// Order 3\n");
    E300 = codeContractMultipoleISIMD(mI, T, l + 3, m, n, sb, tensorHash);
    if (E300 != 0) {
      sb.append(storePotential(to, l + 3, m, n));
    }
    E030 = codeContractMultipoleISIMD(mI, T, l, m + 3, n, sb, tensorHash);
    if (E030 != 0) {
      sb.append(storePotential(to, l, m + 3, n));
    }
    E003 = codeContractMultipoleISIMD(mI, T, l, m, n + 3, sb, tensorHash);
    if (E003 != 0) {
      sb.append(storePotential(to, l, m, n + 3));
    }
    E210 = codeContractMultipoleISIMD(mI, T, l + 2, m + 1, n, sb, tensorHash);
    if (E210 != 0) {
      sb.append(storePotential(to, l + 2, m + 1, n));
    }
    E201 = codeContractMultipoleISIMD(mI, T, l + 2, m, n + 1, sb, tensorHash);
    if (E201 != 0) {
      sb.append(storePotential(to, l + 2, m, n + 1));
    }
    E120 = codeContractMultipoleISIMD(mI, T, l + 1, m + 2, n, sb, tensorHash);
    if (E120 != 0) {
      sb.append(storePotential(to, l + 1, m + 2, n));
    }
    E021 = codeContractMultipoleISIMD(mI, T, l, m + 2, n + 1, sb, tensorHash);
    if (E021 != 0) {
      sb.append(storePotential(to, l, m + 2, n + 1));
    }
    E102 = codeContractMultipoleISIMD(mI, T, l + 1, m, n + 2, sb, tensorHash);
    if (E102 != 0) {
      sb.append(storePotential(to, l + 1, m, n + 2));
    }
    E012 = codeContractMultipoleISIMD(mI, T, l, m + 1, n + 2, sb, tensorHash);
    if (E012 != 0) {
      sb.append(storePotential(to, l, m + 1, n + 2));
    }
    E111 = codeContractMultipoleISIMD(mI, T, l + 1, m + 1, n + 1, sb, tensorHash);
    if (E111 != 0) {
      sb.append(storePotential(to, l + 1, m + 1, n + 1));
    }
  }

  /**
   * Collect the potential its partial derivatives at the origin due to multipole moments at site K.
   *
   * @param mK PolarizableMultipole at site I.
   * @param T  Electrostatic potential and partial derivatives.
   * @param l  apply (d/dx)^l to the potential.
   * @param m  apply (d/dy)^l to the potential.
   * @param n  apply (d/dz)^l to the potential.
   * @param sb Append the code to the StringBuilder.
   */
  protected void codePotentialMultipoleK(PolarizableMultipole mK, double[] T, int l, int m, int n, StringBuilder sb) {
    E000 = codeContractMultipoleK(mK, T, l, m, n, sb);
    if (E000 != 0) {
      sb.append(format("E000 = %s;\n", term(l, m, n)));
    }
    // Order 1 (need a minus sign)
    E100 = codeContractMultipoleK(mK, T, l + 1, m, n, sb);
    if (E100 != 0) {
      sb.append(format("E100 = -%s;\n", term(l + 1, m, n)));
    }
    E010 = codeContractMultipoleK(mK, T, l, m + 1, n, sb);
    if (E100 != 0) {
      sb.append(format("E010 = -%s;\n", term(l, m + 1, n)));
    }
    E001 = codeContractMultipoleK(mK, T, l, m, n + 1, sb);
    if (E001 != 0) {
      sb.append(format("E001 = -%s;\n", term(l, m, n + 1)));
    }
    // Order 2
    E200 = codeContractMultipoleK(mK, T, l + 2, m, n, sb);
    if (E200 != 0) {
      sb.append(format("E200 = %s;\n", term(l + 2, m, n)));
    }
    E020 = codeContractMultipoleK(mK, T, l, m + 2, n, sb);
    if (E020 != 0) {
      sb.append(format("E020 = %s;\n", term(l, m + 2, n)));
    }
    E002 = codeContractMultipoleK(mK, T, l, m, n + 2, sb);
    if (E002 != 0) {
      sb.append(format("E002 = %s;\n", term(l, m, n + 2)));
    }
    E110 = codeContractMultipoleK(mK, T, l + 1, m + 1, n, sb);
    if (E110 != 0) {
      sb.append(format("E110 = %s;\n", term(l + 1, m + 1, n)));
    }
    E101 = codeContractMultipoleK(mK, T, l + 1, m, n + 1, sb);
    if (E101 != 0) {
      sb.append(format("E101 = %s;\n", term(l + 1, m, n + 1)));
    }
    E011 = codeContractMultipoleK(mK, T, l, m + 1, n + 1, sb);
    if (E011 != 0) {
      sb.append(format("E011 = %s;\n", term(l, m + 1, n + 1)));
    }
    // Order 3 (need a minus sign)
    E300 = codeContractMultipoleK(mK, T, l + 3, m, n, sb);
    if (E300 != 0) {
      sb.append(format("E300 = -%s;\n", term(l + 3, m, n)));
    }
    E030 = codeContractMultipoleK(mK, T, l, m + 3, n, sb);
    if (E030 != 0) {
      sb.append(format("E030 = -%s;\n", term(l, m + 3, n)));
    }
    E003 = codeContractMultipoleK(mK, T, l, m, n + 3, sb);
    if (E003 != 0) {
      sb.append(format("E003 = -%s;\n", term(l, m, n + 3)));
    }
    E210 = codeContractMultipoleK(mK, T, l + 2, m + 1, n, sb);
    if (E210 != 0) {
      sb.append(format("E210 = -%s;\n", term(l + 2, m + 1, n)));
    }
    E201 = codeContractMultipoleK(mK, T, l + 2, m, n + 1, sb);
    if (E201 != 0) {
      sb.append(format("E201 = -%s;\n", term(l + 2, m, n + 1)));
    }
    E120 = codeContractMultipoleK(mK, T, l + 1, m + 2, n, sb);
    if (E120 != 0) {
      sb.append(format("E120 = -%s;\n", term(l + 1, m + 2, n)));
    }
    E021 = codeContractMultipoleK(mK, T, l, m + 2, n + 1, sb);
    if (E021 != 0) {
      sb.append(format("E021 = -%s;\n", term(l, m + 2, n + 1)));
    }
    E102 = codeContractMultipoleK(mK, T, l + 1, m, n + 2, sb);
    if (E102 != 0) {
      sb.append(format("E102 = -%s;\n", term(l + 1, m, n + 2)));
    }
    E012 = codeContractMultipoleK(mK, T, l, m + 1, n + 2, sb);
    if (E012 != 0) {
      sb.append(format("E012 = -%s;\n", term(l, m + 1, n + 2)));
    }
    E111 = codeContractMultipoleK(mK, T, l + 1, m + 1, n + 1, sb);
    if (E111 != 0) {
      sb.append(format("E111 = -%s;\n", term(l + 1, m + 1, n + 1)));
    }
  }

  /**
   * Collect the potential its partial derivatives at the origin due to multipole moments at site K
   * using SIMD instructions.
   *
   * @param mK PolarizableMultipole at site I.
   * @param T  Electrostatic potential and partial derivatives.
   * @param l  apply (d/dx)^l to the potential.
   * @param m  apply (d/dy)^l to the potential.
   * @param n  apply (d/dz)^l to the potential.
   * @param sb Append the code to the StringBuilder.
   */
  protected void codePotentialMultipoleKSIMD(PolarizableMultipole mK, double[] T, int l, int m, int n, StringBuilder sb) {
    String to = "e";

    // Order 0.
    sb.append("\n// Order 0\n");
    // Store tensor names to avoid reloading them.
    HashMap<Integer, String> tensorHash = new HashMap<>();
    E000 = codeContractMultipoleKSIMD(mK, T, l, m, n, sb, tensorHash);
    if (E000 != 0) {
      sb.append(storePotential(to, l, m, n));
    }

    // Order 1 (need a minus sign)
    sb.append("\n// Order 1\n");
    E100 = codeContractMultipoleKSIMD(mK, T, l + 1, m, n, sb, tensorHash);
    if (E100 != 0) {
      sb.append(storePotentialNeg(to, l + 1, m, n));
    }
    E010 = codeContractMultipoleKSIMD(mK, T, l, m + 1, n, sb, tensorHash);
    if (E010 != 0) {
      sb.append(storePotentialNeg(to, l, m + 1, n));
    }
    E001 = codeContractMultipoleKSIMD(mK, T, l, m, n + 1, sb, tensorHash);
    if (E001 != 0) {
      sb.append(storePotentialNeg(to, l, m, n + 1));
    }

    // Order 2
    sb.append("\n// Order 2\n");
    E200 = codeContractMultipoleKSIMD(mK, T, l + 2, m, n, sb, tensorHash);
    if (E200 != 0) {
      sb.append(storePotential(to, l + 2, m, n));
    }
    E020 = codeContractMultipoleKSIMD(mK, T, l, m + 2, n, sb, tensorHash);
    if (E020 != 0) {
      sb.append(storePotential(to, l, m + 2, n));
    }
    E002 = codeContractMultipoleKSIMD(mK, T, l, m, n + 2, sb, tensorHash);
    if (E002 != 0) {
      sb.append(storePotential(to, l, m, n + 2));
    }
    E110 = codeContractMultipoleKSIMD(mK, T, l + 1, m + 1, n, sb, tensorHash);
    if (E110 != 0) {
      sb.append(storePotential(to, l + 1, m + 1, n));
    }
    E101 = codeContractMultipoleKSIMD(mK, T, l + 1, m, n + 1, sb, tensorHash);
    if (E101 != 0) {
      sb.append(storePotential(to, l + 1, m, n + 1));
    }
    E011 = codeContractMultipoleKSIMD(mK, T, l, m + 1, n + 1, sb, tensorHash);
    if (E011 != 0) {
      sb.append(storePotential(to, l, m + 1, n + 1));
    }

    // Order 3 (need a minus sign)
    sb.append("\n// Order 3\n");
    E300 = codeContractMultipoleKSIMD(mK, T, l + 3, m, n, sb, tensorHash);
    if (E300 != 0) {
      sb.append(storePotentialNeg(to, l + 3, m, n));
    }
    E030 = codeContractMultipoleKSIMD(mK, T, l, m + 3, n, sb, tensorHash);
    if (E030 != 0) {
      sb.append(storePotentialNeg(to, l, m + 3, n));
    }
    E003 = codeContractMultipoleKSIMD(mK, T, l, m, n + 3, sb, tensorHash);
    if (E003 != 0) {
      sb.append(storePotentialNeg(to, l, m, n + 3));
    }
    E210 = codeContractMultipoleKSIMD(mK, T, l + 2, m + 1, n, sb, tensorHash);
    if (E210 != 0) {
      sb.append(storePotentialNeg(to, l + 2, m + 1, n));
    }
    E201 = codeContractMultipoleKSIMD(mK, T, l + 2, m, n + 1, sb, tensorHash);
    if (E201 != 0) {
      sb.append(storePotentialNeg(to, l + 2, m, n + 1));
    }
    E120 = codeContractMultipoleKSIMD(mK, T, l + 1, m + 2, n, sb, tensorHash);
    if (E120 != 0) {
      sb.append(storePotentialNeg(to, l + 1, m + 2, n));
    }
    E021 = codeContractMultipoleKSIMD(mK, T, l, m + 2, n + 1, sb, tensorHash);
    if (E021 != 0) {
      sb.append(storePotentialNeg(to, l, m + 2, n + 1));
    }
    E102 = codeContractMultipoleKSIMD(mK, T, l + 1, m, n + 2, sb, tensorHash);
    if (E102 != 0) {
      sb.append(storePotentialNeg(to, l + 1, m, n + 2));
    }
    E012 = codeContractMultipoleKSIMD(mK, T, l, m + 1, n + 2, sb, tensorHash);
    if (E012 != 0) {
      sb.append(storePotentialNeg(to, l, m + 1, n + 2));
    }
    E111 = codeContractMultipoleKSIMD(mK, T, l + 1, m + 1, n + 1, sb, tensorHash);
    if (E111 != 0) {
      sb.append(storePotentialNeg(to, l + 1, m + 1, n + 1));
    }
  }

  /**
   * Contract a multipole with the potential and its derivatives.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @return The permanent multipole energy.
   */
  protected final double multipoleEnergy(PolarizableMultipole m) {
    double total = m.q * E000;
    total = fma(m.dx, E100, total);
    total = fma(m.dy, E010, total);
    total = fma(m.dz, E001, total);
    total = fma(m.qxx, E200, total);
    total = fma(m.qyy, E020, total);
    total = fma(m.qzz, E002, total);
    total = fma(m.qxy, E110, total);
    total = fma(m.qxz, E101, total);
    total = fma(m.qyz, E011, total);
    return total;
  }

  /**
   * Compute the permanent multipole gradient.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @param g The atomic gradient.
   */
  protected final void multipoleGradient(PolarizableMultipole m, double[] g) {
    // dEnergy/dY
    double total = m.q * E100;
    total = fma(m.dx, E200, total);
    total = fma(m.dy, E110, total);
    total = fma(m.dz, E101, total);
    total = fma(m.qxx, E300, total);
    total = fma(m.qyy, E120, total);
    total = fma(m.qzz, E102, total);
    total = fma(m.qxy, E210, total);
    total = fma(m.qxz, E201, total);
    total = fma(m.qyz, E111, total);
    g[0] = total;

    // dEnergy/dY
    total = m.q * E010;
    total = fma(m.dx, E110, total);
    total = fma(m.dy, E020, total);
    total = fma(m.dz, E011, total);
    total = fma(m.qxx, E210, total);
    total = fma(m.qyy, E030, total);
    total = fma(m.qzz, E012, total);
    total = fma(m.qxy, E120, total);
    total = fma(m.qxz, E111, total);
    total = fma(m.qyz, E021, total);
    g[1] = total;

    // dEnergy/dZ
    total = m.q * E001;
    total = fma(m.dx, E101, total);
    total = fma(m.dy, E011, total);
    total = fma(m.dz, E002, total);
    total = fma(m.qxx, E201, total);
    total = fma(m.qyy, E021, total);
    total = fma(m.qzz, E003, total);
    total = fma(m.qxy, E111, total);
    total = fma(m.qxz, E102, total);
    total = fma(m.qyz, E012, total);
    g[2] = total;
  }

  /**
   * Compute the torque on a permanent multipole.
   *
   * @param m      PolarizableMultipole at the site of the potential.
   * @param torque an array of double values.
   */
  protected final void multipoleTorque(PolarizableMultipole m, double[] torque) {
    // Torque on the permanent dipole due to the field.
    double dx = m.dy * E001 - m.dz * E010;
    double dy = m.dz * E100 - m.dx * E001;
    double dz = m.dx * E010 - m.dy * E100;

    // Torque on the permanent quadrupole due to the gradient of the field.
    double qx = m.qxy * E101 + 2.0 * m.qyy * E011 + m.qyz * E002
        - (m.qxz * E110 + m.qyz * E020 + 2.0 * m.qzz * E011);
    double qy = m.qxz * E200 + m.qyz * E110 + 2.0 * m.qzz * E101
        - (2.0 * m.qxx * E101 + m.qxy * E011 + m.qxz * E002);
    double qz = 2.0 * m.qxx * E110 + m.qxy * E020 + m.qxz * E011
        - (m.qxy * E200 + 2.0 * m.qyy * E110 + m.qyz * E101);

    // The field along X is -E001, so we need a negative sign.
    torque[0] -= (dx + qx);
    torque[1] -= (dy + qy);
    torque[2] -= (dz + qz);
  }

  /**
   * Compute the torque on a permanent dipole.
   *
   * @param m      PolarizableMultipole at the site of the potential.
   * @param torque an array of double values.
   */
  protected final void dipoleTorque(PolarizableMultipole m, double[] torque) {
    // Torque on the permanent dipole due to the field.
    double dx = m.dy * E001 - m.dz * E010;
    double dy = m.dz * E100 - m.dx * E001;
    double dz = m.dx * E010 - m.dy * E100;

    // The field along X is -E001, so we need a negative sign.
    torque[0] -= dx;
    torque[1] -= dy;
    torque[2] -= dz;
  }

  /**
   * Compute the torque on a permanent quadrupole.
   *
   * @param m      PolarizableMultipole at the site of the potential.
   * @param torque an array of double values.
   */
  protected final void quadrupoleTorque(PolarizableMultipole m, double[] torque) {
    // Torque on the permanent quadrupole due to the gradient of the field.
    double qx = m.qxy * E101 + 2.0 * m.qyy * E011 + m.qyz * E002
        - (m.qxz * E110 + m.qyz * E020 + 2.0 * m.qzz * E011);
    double qy = m.qxz * E200 + m.qyz * E110 + 2.0 * m.qzz * E101
        - (2.0 * m.qxx * E101 + m.qxy * E011 + m.qxz * E002);
    double qz = 2.0 * m.qxx * E110 + m.qxy * E020 + m.qxz * E011
        - (m.qxy * E200 + 2.0 * m.qyy * E110 + m.qyz * E101);

    // The field along X is -E001, so we need a negative sign.
    torque[0] -= qx;
    torque[1] -= qy;
    torque[2] -= qz;
  }

  /**
   * Contract an induced dipole with the potential and its derivatives.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @return The polarization energy.
   */
  protected final double polarizationEnergy(PolarizableMultipole m) {
    // E = -1/2 * u.E
    // No negative sign because the field E = [-E100, -E010, -E001].
    return 0.5 * (m.ux * E100 + m.uy * E010 + m.uz * E001);
  }

  /**
   * Contract an induced dipole with the potential and its derivatives.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @return The polarization energy.
   */
  protected final double polarizationEnergyS(PolarizableMultipole m) {
    // E = -1/2 * u.E
    // No negative sign because the field E = [-E100, -E010, -E001].
    return 0.5 * (m.sx * E100 + m.sy * E010 + m.sz * E001);
  }

  /**
   * Load the tensor components.
   *
   * @param T an array of double values.
   */
  @SuppressWarnings("fallthrough")
  protected final void getTensor(double[] T) {
    switch (order) {
      default:
      case 5:
        // l + m + n = 5 (21) 56
        T[t500] = R500;
        T[t050] = R050;
        T[t005] = R005;
        T[t410] = R410;
        T[t401] = R401;
        T[t140] = R140;
        T[t041] = R041;
        T[t104] = R104;
        T[t014] = R014;
        T[t320] = R320;
        T[t302] = R302;
        T[t230] = R230;
        T[t032] = R032;
        T[t203] = R203;
        T[t023] = R023;
        T[t311] = R311;
        T[t131] = R131;
        T[t113] = R113;
        T[t221] = R221;
        T[t212] = R212;
        T[t122] = R122;
        // Fall through to 4th order.
      case 4:
        // l + m + n = 4 (15) 35
        T[t400] = R400;
        T[t040] = R040;
        T[t004] = R004;
        T[t310] = R310;
        T[t301] = R301;
        T[t130] = R130;
        T[t031] = R031;
        T[t103] = R103;
        T[t013] = R013;
        T[t220] = R220;
        T[t202] = R202;
        T[t022] = R022;
        T[t211] = R211;
        T[t121] = R121;
        T[t112] = R112;
        // Fall through to 3rd order.
      case 3:
        // l + m + n = 3 (10) 20
        T[t300] = R300;
        T[t030] = R030;
        T[t003] = R003;
        T[t210] = R210;
        T[t201] = R201;
        T[t120] = R120;
        T[t021] = R021;
        T[t102] = R102;
        T[t012] = R012;
        T[t111] = R111;
        // Fall through to 2nd order.
      case 2:
        // l + m + n = 2 (6)  10
        T[t200] = R200;
        T[t020] = R020;
        T[t002] = R002;
        T[t110] = R110;
        T[t101] = R101;
        T[t011] = R011;
        // Fall through to 1st order.
      case 1:
        // l + m + n = 1 (3)   4
        T[t100] = R100;
        T[t010] = R010;
        T[t001] = R001;
        // Fall through to the potential.
      case 0:
        // l + m + n = 0 (1)
        T[t000] = R000;
    }
  }

  /**
   * Set the tensor components.
   *
   * @param T an array of double values.
   */
  @SuppressWarnings("fallthrough")
  protected final void setTensor(double[] T) {
    switch (order) {
      case 5:
        // l + m + n = 5 (21) 56
        R500 = T[t500];
        R050 = T[t050];
        R005 = T[t005];
        R410 = T[t410];
        R401 = T[t401];
        R140 = T[t140];
        R041 = T[t041];
        R104 = T[t104];
        R014 = T[t014];
        R320 = T[t320];
        R302 = T[t302];
        R230 = T[t230];
        R032 = T[t032];
        R203 = T[t203];
        R023 = T[t023];
        R311 = T[t311];
        R131 = T[t131];
        R113 = T[t113];
        R221 = T[t221];
        R212 = T[t212];
        R122 = T[t122];
        // Fall through to 4th order.
      case 4:
        // l + m + n = 4 (15) 35
        R400 = T[t400];
        R040 = T[t040];
        R004 = T[t004];
        R310 = T[t310];
        R301 = T[t301];
        R130 = T[t130];
        R031 = T[t031];
        R103 = T[t103];
        R013 = T[t013];
        R220 = T[t220];
        R202 = T[t202];
        R022 = T[t022];
        R211 = T[t211];
        R121 = T[t121];
        R112 = T[t112];
        // Fall through to 3rd order.
      case 3:
        // l + m + n = 3 (10) 20
        R300 = T[t300];
        R030 = T[t030];
        R003 = T[t003];
        R210 = T[t210];
        R201 = T[t201];
        R120 = T[t120];
        R021 = T[t021];
        R102 = T[t102];
        R012 = T[t012];
        R111 = T[t111];
        // Fall through to 2nd order.
      case 2:
        // l + m + n = 2 (6)  10
        R200 = T[t200];
        R020 = T[t020];
        R002 = T[t002];
        R110 = T[t110];
        R101 = T[t101];
        R011 = T[t011];
        // Fall through to 1st order.
      case 1:
        // l + m + n = 1 (3)   4
        R100 = T[t100];
        R010 = T[t010];
        R001 = T[t001];
        // Fall through to the potential.
      case 0:
        // l + m + n = 0 (1)
        R000 = T[t000];
    }
  }

  /**
   * This method is a driver to collect elements of the Cartesian multipole tensor given the
   * recursion relationships implemented by the method "Tlmnj", which can be called directly to get a
   * single tensor element. It does not store intermediate values of the recursion, causing it to
   * scale O(order^8). For order = 5, this approach is a factor of 10 slower than recursion.
   *
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   */
  protected abstract void noStorageRecursion(double[] tensor);

  /**
   * This method is a driver to collect elements of the Cartesian multipole tensor given the
   * recursion relationships implemented by the method "Tlmnj", which can be called directly to get a
   * single tensor element. It does not store intermediate values of the recursion, causing it to
   * scale O(order^8). For order = 5, this approach is a factor of 10 slower than recursion.
   *
   * @param r      double[] vector between two sites.
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   */
  protected abstract void noStorageRecursion(double[] r, double[] tensor);

  /**
   * This routine implements the recurrence relations for computation of any Cartesian multipole
   * tensor in ~O(L^8) time, where L is the total order l + m + n, given the auxiliary elements
   * T0000. <br> It implements the recursion relationships in brute force fashion, without saving
   * intermediate values. This is useful for finding a single tensor, rather than all binomial(L + 3,
   * 3). <br> The specific recursion equations (41-43) and set of auxiliary tensor elements from
   * equation (40) can be found in Challacombe et al.
   *
   * @param l    int The number of (d/dx) operations.
   * @param m    int The number of (d/dy) operations.
   * @param n    int The number of (d/dz) operations.
   * @param j    int j = 0 is the Tlmn tensor, j .GT. 0 is an intermediate.
   * @param r    double[] The {x,y,z} coordinates.
   * @param T000 double[] Initial auxiliary tensor elements from Eq. (40).
   * @return double The requested Tensor element (intermediate if j .GT. 0).
   * @since 1.0
   */
  protected abstract double Tlmnj(final int l, final int m, final int n, final int j,
                                  final double[] r, final double[] T000);

  /**
   * This method is a driver to collect elements of the Cartesian multipole tensor using recursion
   * relationships and storing intermediate values. It scales approximately O(order^4).
   *
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   */
  protected abstract void recursion(final double[] tensor);

  /**
   * This method is a driver to collect elements of the Cartesian multipole tensor using recursion
   * relationships and storing intermediate values. It scales approximately O(order^4).
   *
   * @param r      double[] vector between two sites.
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   */
  protected abstract void recursion(final double[] r, final double[] tensor);

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 1st order.
   */
  protected abstract void order1();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 2nd order.
   */
  protected abstract void order2();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 3rd order.
   */
  protected abstract void order3();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 4th order.
   */
  protected abstract void order4();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 5th order, which is needed for
   * quadrupole-quadrupole forces.
   */
  protected abstract void order5();

  /**
   * Hard coded computation of the Cartesian multipole tensors up to 6th order, which is needed for
   * quadrupole-quadrupole forces and orthogonal space sampling.
   */
  protected abstract void order6();

  /**
   * Compute the field components due to multipole I at site K.
   *
   * @param mI    PolarizableMultipole at site I.
   * @param order Compute derivatives of the potential up to this order. Order 0: Electrostatic
   *              potential (E000) Order 1: First derivatives: d/dX is E100, d/dY is E010, d/dZ is E001. Order
   *              2: Second derivatives: d2/dXdX is E200, d2/dXdY is E110 (needed for quadrupole energy) Order
   *              3: 3rd derivatives: (needed for quadrupole forces).
   */
  protected abstract void multipoleIPotentialAtK(PolarizableMultipole mI, int order);

  /**
   * Compute the field components due to charge I at site K.
   *
   * @param mI    PolarizableMultipole at site I.
   * @param order Compute derivatives of the potential up to this order. Order 0: Electrostatic
   *              potential (E000) Order 1: First derivatives: d/dX is E100, d/dY is E010, d/dZ is E001. Order
   *              2: Second derivatives: d2/dXdX is E200, d2/dXdY is E110 (needed for quadrupole energy) Order
   *              3: 3rd derivatives: (needed for quadrupole forces).
   */
  protected abstract void chargeIPotentialAtK(PolarizableMultipole mI, int order);

  /**
   * Compute the induced dipole field components due to site I at site K.
   *
   * @param uxi   X-dipole component.
   * @param uyi   Y-dipole component.
   * @param uzi   Z-dipole component.
   * @param order Potential order.
   */
  protected abstract void dipoleIPotentialAtK(double uxi, double uyi, double uzi, int order);

  /**
   * Compute the field components due to quadrupole I at site K.
   *
   * @param mI    PolarizableMultipole at site I.
   * @param order Compute derivatives of the potential up to this order. Order 0: Electrostatic
   *              potential (E000) Order 1: First derivatives: d/dX is E100, d/dY is E010, d/dZ is E001. Order
   *              2: Second derivatives: d2/dXdX is E200, d2/dXdY is E110 (needed for quadrupole energy) Order
   *              3: 3rd derivatives: (needed for quadrupole forces).
   */
  protected abstract void quadrupoleIPotentialAtK(PolarizableMultipole mI, int order);

  /**
   * Compute the field components due to multipole K at site I.
   *
   * @param mK    PolarizableMultipole at site K.
   * @param order Compute derivatives of the potential up to this order. Order 0: Electrostatic
   *              potential (E000) Order 1: First derivatives: d/dX is E100, d/dY is E010, d/dZ is E001. Order
   *              2: Second derivatives: d2/dXdX is E200, d2/dXdY is E110 (needed for quadrupole energy) Order
   *              3: 3rd derivatives: (needed for quadrupole forces).
   */
  protected abstract void multipoleKPotentialAtI(PolarizableMultipole mK, int order);

  /**
   * Compute the field components due to multipole K at site I.
   *
   * @param mK    PolarizableMultipole at site K.
   * @param order Compute derivatives of the potential up to this order. Order 0: Electrostatic
   *              potential (E000) Order 1: First derivatives: d/dX is E100, d/dY is E010, d/dZ is E001. Order
   *              2: Second derivatives: d2/dXdX is E200, d2/dXdY is E110 (needed for quadrupole energy) Order
   *              3: 3rd derivatives: (needed for quadrupole forces).
   */
  protected abstract void chargeKPotentialAtI(PolarizableMultipole mK, int order);

  /**
   * Compute the induced dipole field components due to site K at site I.
   *
   * @param uxk   X-dipole component.
   * @param uyk   Y-dipole component.
   * @param uzk   Z-dipole component.
   * @param order Potential order.
   */
  protected abstract void dipoleKPotentialAtI(double uxk, double uyk, double uzk, int order);

  /**
   * Compute the field components due to multipole K at site I.
   *
   * @param mK    PolarizableMultipole at site K.
   * @param order Compute derivatives of the potential up to this order. Order 0: Electrostatic
   *              potential (E000) Order 1: First derivatives: d/dX is E100, d/dY is E010, d/dZ is E001. Order
   *              2: Second derivatives: d2/dXdX is E200, d2/dXdY is E110 (needed for quadrupole energy) Order
   *              3: 3rd derivatives: (needed for quadrupole forces).
   */
  protected abstract void quadrupoleKPotentialAtI(PolarizableMultipole mK, int order);

  /**
   * Log the tensors.
   *
   * @param operator The OPERATOR to use.
   * @param order    The tensor order.
   * @param tensor   The tensor array.
   */
  private static void log(Operator operator, int order, double[] tensor) {
    final int o1 = order + 1;
    StringBuilder sb = new StringBuilder();

    sb.append(format("\n %s Operator to order %d:", operator, order));
    sb.append(format("\n%5s %4s %4s %4s %12s\n", "Index", "d/dx", "d/dy", "d/dz", "Tensor"));
    sb.append(format("%5d %4d %4d %4d %12.8f\n", 0, 0, 0, 0, tensor[0]));
    int count = 1;
    // Print (d/dx)^l for l = 1..order (m = 0, n = 0)
    for (int l = 1; l <= order; l++) {
      double value = tensor[MultipoleUtilities.ti(l, 0, 0, order)];
      if (value != 0.0) {
        sb.append(format("%5d %4d %4d %4d %12.8f\n", MultipoleUtilities.ti(l, 0, 0, order), l, 0, 0, value));
        count++;
      }
    }
    // Print (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
    for (int l = 0; l <= o1; l++) {
      for (int m = 1; m <= order - l; m++) {
        double value = tensor[MultipoleUtilities.ti(l, m, 0, order)];
        if (value != 0.0) {
          sb.append(format("%5d %4d %4d %4d %12.8f\n", MultipoleUtilities.ti(l, m, 0, order), l, m, 0, value));
          count++;
        }
      }
    }
    // Print (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
    for (int l = 0; l <= o1; l++) {
      for (int m = 0; m <= o1 - l; m++) {
        for (int n = 1; n <= order - l - m; n++) {
          double value = tensor[MultipoleUtilities.ti(l, m, n, order)];
          if (value != 0.0) {
            sb.append(format("%5d %4d %4d %4d %12.8f\n", MultipoleUtilities.ti(l, m, n, order), l, m, n, value));
            count++;
          }
        }
      }
    }
    sb.append(format("\n Total number of active tensors: %d\n", count));
    logger.log(Level.INFO, sb.toString());
  }

  // l + m + n = 0 (1)
  /**
   * No derivatives.
   */
  protected final int t000;
  // l + m + n = 1 (3)   4
  /**
   * First derivative with respect to x.
   */
  protected final int t100;
  /**
   * First derivative with respect to y.
   */
  protected final int t010;
  /**
   * First derivative with respect to z.
   */
  protected final int t001;
  // l + m + n = 2 (6)  10
  /**
   * Second derivative with respect to x.
   */
  protected final int t200;
  /**
   * Second derivative with respect to y.
   */
  protected final int t020;
  /**
   * Second derivative with respect to z.
   */
  protected final int t002;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t110;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t101;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t011;
  // l + m + n = 3 (10) 20
  /**
   * Third derivative with respect to x.
   */
  protected final int t300;
  /**
   * Third derivative with respect to y.
   */
  protected final int t030;
  /**
   * Third derivative with respect to z.
   */
  protected final int t003;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t210;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t201;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t120;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t021;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t102;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t012;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t111;
  // l + m + n = 4 (15) 35
  /**
   * Fourth derivative with respect to x.
   */
  protected final int t400;
  /**
   * Fourth derivative with respect to y.
   */
  protected final int t040;
  /**
   * Fourth derivative with respect to z.
   */
  protected final int t004;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t310;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t301;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t130;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t031;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t103;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t013;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t220;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t202;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t022;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t211;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t121;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t112;
  // l + m + n = 5 (21) 56
  /**
   * Fifth derivative with respect to x.
   */
  protected final int t500;
  /**
   * Fifth derivative with respect to y.
   */
  protected final int t050;
  /**
   * Fifth derivative with respect to z.
   */
  protected final int t005;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t410;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t401;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t140;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t041;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t104;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t014;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t320;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t302;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t230;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t032;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t203;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t023;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t311;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t131;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t113;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t221;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t212;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t122;

  // l + m + n = 6 (28) 84
  /**
   * Sixth derivative with respect to x.
   */
  protected final int t600;
  /**
   * Sixth derivative with respect to y.
   */
  protected final int t060;
  /**
   * Sixth derivative with respect to z.
   */
  protected final int t006;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t510;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t501;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t150;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t051;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t105;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t015;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t420;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t402;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t240;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t042;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t204;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t024;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t411;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t141;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t114;
  /**
   * Derivatives with respect to x and y.
   */
  protected final int t330;
  /**
   * Derivatives with respect to x and z.
   */
  protected final int t303;
  /**
   * Derivatives with respect to y and z.
   */
  protected final int t033;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t321;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t231;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t213;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t312;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t132;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t123;
  /**
   * Derivatives with respect to x, y and z.
   */
  protected final int t222;
}
