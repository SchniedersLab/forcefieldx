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

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.ScalarMath.doubleFactorial;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * The abstract MultipoleTensorSIMD is extended by classes that compute derivatives of 1/|<b>r</b>| via
 * recursion to arbitrary order using Cartesian multipoles in either a global frame or a
 * quasi-internal frame.
 * <br>
 * This class serves as the abstract parent to both and defines all shared logic. Non-abstract methods
 * are declared final to disallow unnecessary overrides.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public abstract class MultipoleTensorSIMD {

  /**
   * Logger for the MultipoleTensor class.
   */
  private static final Logger logger = Logger.getLogger(MultipoleTensorSIMD.class.getName());

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
   * The coordinate system in use (global or QI).
   */
  protected final CoordinateSystem coordinates;

  /**
   * Separation distance.
   */
  protected DoubleVector R;

  /**
   * Separation distance squared.
   */
  protected DoubleVector r2;

  /**
   * Xk - Xi.
   */
  protected DoubleVector x;

  /**
   * Yk - Yi.
   */
  protected DoubleVector y;

  /**
   * Zk - Zi.
   */
  protected DoubleVector z;

  /**
   * A work array with auxiliary terms for the recursion.
   */
  protected final DoubleVector[] work;

  // Cartesian tensor elements (for 1/R, erfc(Beta*R)/R or DoubleVectorhole damping.
  protected DoubleVector R000;
  // l + m + n = 1 (3)   4
  protected DoubleVector R100;
  protected DoubleVector R010;
  protected DoubleVector R001;
  // l + m + n = 2 (6)  10
  protected DoubleVector R200;
  protected DoubleVector R020;
  protected DoubleVector R002;
  protected DoubleVector R110;
  protected DoubleVector R101;
  protected DoubleVector R011;
  // l + m + n = 3 (10) 20
  protected DoubleVector R300;
  protected DoubleVector R030;
  protected DoubleVector R003;
  protected DoubleVector R210;
  protected DoubleVector R201;
  protected DoubleVector R120;
  protected DoubleVector R021;
  protected DoubleVector R102;
  protected DoubleVector R012;
  protected DoubleVector R111;
  // l + m + n = 4 (15) 35
  protected DoubleVector R400;
  protected DoubleVector R040;
  protected DoubleVector R004;
  protected DoubleVector R310;
  protected DoubleVector R301;
  protected DoubleVector R130;
  protected DoubleVector R031;
  protected DoubleVector R103;
  protected DoubleVector R013;
  protected DoubleVector R220;
  protected DoubleVector R202;
  protected DoubleVector R022;
  protected DoubleVector R211;
  protected DoubleVector R121;
  protected DoubleVector R112;
  // l + m + n = 5 (21) 56
  protected DoubleVector R500;
  protected DoubleVector R050;
  protected DoubleVector R005;
  protected DoubleVector R410;
  protected DoubleVector R401;
  protected DoubleVector R140;
  protected DoubleVector R041;
  protected DoubleVector R104;
  protected DoubleVector R014;
  protected DoubleVector R320;
  protected DoubleVector R302;
  protected DoubleVector R230;
  protected DoubleVector R032;
  protected DoubleVector R203;
  protected DoubleVector R023;
  protected DoubleVector R311;
  protected DoubleVector R131;
  protected DoubleVector R113;
  protected DoubleVector R221;
  protected DoubleVector R212;
  protected DoubleVector R122;
  // l + m + n = 6 (28) 84
  protected DoubleVector R006;
  protected DoubleVector R402;
  protected DoubleVector R042;
  protected DoubleVector R204;
  protected DoubleVector R024;
  protected DoubleVector R222;
  protected DoubleVector R600;
  protected DoubleVector R060;
  protected DoubleVector R510;
  protected DoubleVector R501;
  protected DoubleVector R150;
  protected DoubleVector R051;
  protected DoubleVector R105;
  protected DoubleVector R015;
  protected DoubleVector R420;
  protected DoubleVector R240;
  protected DoubleVector R411;
  protected DoubleVector R141;
  protected DoubleVector R114;
  protected DoubleVector R330;
  protected DoubleVector R303;
  protected DoubleVector R033;
  protected DoubleVector R321;
  protected DoubleVector R231;
  protected DoubleVector R213;
  protected DoubleVector R312;
  protected DoubleVector R132;
  protected DoubleVector R123;

  // Components of the potential, field and field gradient.
  protected DoubleVector E000; // Potential
  // l + m + n = 1 (3)   4
  protected DoubleVector E100; // d/dX
  protected DoubleVector E010; // d/dY
  protected DoubleVector E001; // d/dz
  // l + m + n = 2 (6)  10
  protected DoubleVector E200; // d^2/dXdX
  protected DoubleVector E020; // d^2/dYdY
  protected DoubleVector E002; // d^2/dZdZ
  protected DoubleVector E110; // d^2/dXdY
  protected DoubleVector E101; // d^2/dXdZ
  protected DoubleVector E011; // d^2/dYdZ
  // l + m + n = 3 (10) 20
  protected DoubleVector E300; // d^3/dXdXdX
  protected DoubleVector E030; // d^3/dYdYdY
  protected DoubleVector E003; // d^3/dZdZdZ
  protected DoubleVector E210; // d^3/dXdXdY
  protected DoubleVector E201; // d^3/dXdXdZ
  protected DoubleVector E120; // d^3/dXdYdY
  protected DoubleVector E021; // d^3/dYdYdZ
  protected DoubleVector E102; // d^3/dXdZdZ
  protected DoubleVector E012; // d^3/dYdZdZ
  protected DoubleVector E111; // d^3/dXdYdZ

  /**
   * Constructor for MultipoleTensor.
   *
   * @param order       The order of the tensor.
   * @param coordinates a {@link CoordinateSystem} object.
   */
  public MultipoleTensorSIMD(CoordinateSystem coordinates, int order) {
    assert (order > 0);
    o1 = order + 1;
    il = o1;
    im = il * o1;
    in = im * o1;
    size = (order + 1) * (order + 2) * (order + 3) / 6;

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

    work = new DoubleVector[in * o1];

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
   * Set the separation vector.
   *
   * @param r The separation vector.
   */
  public final void setR(DoubleVector[] r) {
    setR(r[0], r[1], r[2]);
  }

  /**
   * Set the separation vector.
   *
   * @param dx Separation along the X-axis.
   * @param dy Separation along the Y-axis.
   * @param dz Separation along the Z-axis.
   */
  public abstract void setR(DoubleVector dx, DoubleVector dy, DoubleVector dz);

  /**
   * Generate the tensor using hard-coded methods.
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
        if (logger.isLoggable(Level.WARNING)) {
          logger.severe("Order " + order + " not supported.");
        }
      }
    }
  }

  /**
   * Return the source terms.
   *
   * @return A DoubleVector array of source terms.
   */
  public DoubleVector[] getSource() {
    source(work);
    return work;
  }

  /**
   * Contract a multipole with the potential and its derivatives.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @return The permanent multipole energy.
   */
  protected final DoubleVector multipoleEnergy(PolarizableMultipoleSIMD m) {
    DoubleVector total = m.q.mul(E000);
    total = m.dx.fma(E100, total);
    total = m.dy.fma(E010, total);
    total = m.dz.fma(E001, total);
    total = m.qxx.fma(E200, total);
    total = m.qyy.fma(E020, total);
    total = m.qzz.fma(E002, total);
    total = m.qxy.fma(E110, total);
    total = m.qxz.fma(E101, total);
    total = m.qyz.fma(E011, total);
    return total;
  }

  /**
   * Contract a multipole with the potential and its derivatives.
   *
   * @param mI PolarizableMultipole at site I.
   * @param mK PolarizableMultipole at site K.
   * @return The permanent multipole energy.
   */
  public DoubleVector multipoleEnergy(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    multipoleIPotentialAtK(mI, 2);
    return multipoleEnergy(mK);
  }

  /**
   * Compute the permanent multipole gradient.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @param g The atomic gradient.
   */
  protected final void multipoleGradient(PolarizableMultipoleSIMD m, DoubleVector[] g) {
    // dEnergy/dX
    DoubleVector total = m.q.mul(E100);
    total = m.dx.fma(E200, total);
    total = m.dy.fma(E110, total);
    total = m.dz.fma(E101, total);
    total = m.qxx.fma(E300, total);
    total = m.qyy.fma(E120, total);
    total = m.qzz.fma(E102, total);
    total = m.qxy.fma(E210, total);
    total = m.qxz.fma(E201, total);
    total = m.qyz.fma(E111, total);
    g[0] = total;

    // dEnergy/dY
    total = m.q.mul(E010);
    total = m.dx.fma(E110, total);
    total = m.dy.fma(E020, total);
    total = m.dz.fma(E011, total);
    total = m.qxx.fma(E210, total);
    total = m.qyy.fma(E030, total);
    total = m.qzz.fma(E012, total);
    total = m.qxy.fma(E120, total);
    total = m.qxz.fma(E111, total);
    total = m.qyz.fma(E021, total);
    g[1] = total;

    // dEnergy/dZ
    total = m.q.mul(E001);
    total = m.dx.fma(E101, total);
    total = m.dy.fma(E011, total);
    total = m.dz.fma(E002, total);
    total = m.qxx.fma(E201, total);
    total = m.qyy.fma(E021, total);
    total = m.qzz.fma(E003, total);
    total = m.qxy.fma(E111, total);
    total = m.qxz.fma(E102, total);
    total = m.qyz.fma(E012, total);
    g[2] = total;
  }

  /**
   * Compute the torque on a permanent multipole.
   *
   * @param m      PolarizableMultipole at the site of the potential.
   * @param torque an array of double values.
   */
  protected final void multipoleTorque(PolarizableMultipoleSIMD m, DoubleVector[] torque) {
    // Torque on the permanent dipole due to the field.
    DoubleVector dx = m.dy.mul(E001).sub(m.dz.mul(E010));
    DoubleVector dy = m.dz.mul(E100).sub(m.dx.mul(E001));
    DoubleVector dz = m.dx.mul(E010).sub(m.dy.mul(E100));

    // Torque on the permanent quadrupole due to the gradient of the field.
    DoubleVector qx = m.qxy.mul(E101).add(m.qyy.mul(E011).mul(2.0)).add(m.qyz.mul(E002))
        .sub(m.qxz.mul(E110).add(m.qyz.mul(E020)).add(m.qzz.mul(E011).mul(2.0)));
    DoubleVector qy = m.qxz.mul(E200).add(m.qyz.mul(E110)).add(m.qzz.mul(E101).mul(2.0))
        .sub(m.qxx.mul(E101).mul(2.0).add(m.qxy.mul(E011)).add(m.qxz.mul(E002)));
    DoubleVector qz = m.qxx.mul(E110).mul(2.0).add(m.qxy.mul(E020)).add(m.qxz.mul(E011))
        .sub(m.qxy.mul(E200).add(m.qyy.mul(E110).mul(2.0)).add(m.qyz.mul(E101)));

    // The field along X is -E001, so we need a negative sign.
    torque[0] = torque[0].sub(dx.add(qx));
    torque[1] = torque[1].sub(dy.add(qy));
    torque[2] = torque[2].sub(dz.add(qz));
  }

  /**
   * Compute the torque on a permanent dipole.
   *
   * @param m      PolarizableMultipole at the site of the potential.
   * @param torque an array of double values.
   */
  protected final void dipoleTorque(PolarizableMultipoleSIMD m, DoubleVector[] torque) {
    // Torque on the permanent dipole due to the field.
    DoubleVector dx = m.dy.mul(E001).sub(m.dz.mul(E010));
    DoubleVector dy = m.dz.mul(E100).sub(m.dx.mul(E001));
    DoubleVector dz = m.dx.mul(E010).sub(m.dy.mul(E100));

    // The field along X is -E001, so we need a negative sign.
    torque[0] = torque[0].sub(dx);
    torque[1] = torque[1].sub(dy);
    torque[2] = torque[2].sub(dz);
  }

  /**
   * Compute the torque on a permanent quadrupole.
   *
   * @param m      PolarizableMultipole at the site of the potential.
   * @param torque an array of double values.
   */
  protected final void quadrupoleTorque(PolarizableMultipoleSIMD m, DoubleVector[] torque) {
    // Torque on the permanent quadrupole due to the gradient of the field.
    DoubleVector qx = m.qxy.mul(E101).add(m.qyy.mul(E011).mul(2.0)).add(m.qyz.mul(E002))
        .sub(m.qxz.mul(E110).add(m.qyz.mul(E020)).add(m.qzz.mul(E011).mul(2.0)));
    DoubleVector qy = m.qxz.mul(E200).add(m.qyz.mul(E110)).add(m.qzz.mul(E101).mul(2.0))
        .sub(m.qxx.mul(E101).mul(2.0).add(m.qxy.mul(E011)).add(m.qxz.mul(E002)));
    DoubleVector qz = m.qxx.mul(E110).mul(2.0).add(m.qxy.mul(E020)).add(m.qxz.mul(E011))
        .sub(m.qxy.mul(E200).add(m.qyy.mul(E110).mul(2.0)).add(m.qyz.mul(E101)));

    // The field along X is -E001, so we need a negative sign.
    torque[0] = torque[0].sub(qx);
    torque[1] = torque[1].sub(qy);
    torque[2] = torque[2].sub(qz);
  }

  /**
   * Contract an induced dipole with the potential and its derivatives.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @return The polarization energy.
   */
  protected final DoubleVector polarizationEnergy(PolarizableMultipoleSIMD m) {
    // E = -1/2 * u.E
    // No negative sign because the field E = [-E100, -E010, -E001].
    return (m.ux.mul(E100).add(m.uy.mul(E010)).add(m.uz.mul(E001))).mul(.5);
  }

  /**
   * Contract an induced dipole with the potential and its derivatives.
   *
   * @param m PolarizableMultipole at the site of the potential.
   * @return The polarization energy.
   */
  protected final DoubleVector polarizationEnergyS(PolarizableMultipoleSIMD m) {
    // E = -1/2 * u.E
    // No negative sign because the field E = [-E100, -E010, -E001].
    return (m.sx.mul(E100).add(m.sy.mul(E010)).add(m.sz.mul(E001))).mul(.5);
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
  public DoubleVector polarizationEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                    DoubleVector inductionMask, DoubleVector energyMask, DoubleVector mutualMask,
                                                    DoubleVector[] Gi, DoubleVector[] Ti, DoubleVector[] Tk) {

    // Add the induction and energy masks to create an "averaged" induced dipole (sx, sy, sz).
    mI.applyMasks(inductionMask, energyMask);
    mK.applyMasks(inductionMask, energyMask);

    // Find the permanent multipole potential and derivatives at site k.
    multipoleIPotentialAtK(mI, 2);
    // Energy of induced dipole k in the field of multipole i.
    // The field E_x = -E100.
    DoubleVector eK = polarizationEnergy(mK);
    // Derivative with respect to moving atom k.
    Gi[0] = mK.sx.mul(E200).add(mK.sy.mul(E110)).add(mK.sz.mul(E101)).neg();
    Gi[1] = mK.sx.mul(E110).add(mK.sy.mul(E020)).add(mK.sz.mul(E011)).neg();
    Gi[2] = mK.sx.mul(E101).add(mK.sy.mul(E011)).add(mK.sz.mul(E002)).neg();

    // Find the permanent multipole potential and derivatives at site i.
    multipoleKPotentialAtI(mK, 2);
    // Energy of induced dipole i in the field of multipole k.
    DoubleVector eI = polarizationEnergy(mI);
    // Derivative with respect to moving atom i.
    Gi[0] = Gi[0].add(mI.sx.mul(E200).add(mI.sy.mul(E110)).add(mI.sz.mul(E101)));
    Gi[1] = Gi[1].add(mI.sx.mul(E110).add(mI.sy.mul(E020)).add(mI.sz.mul(E011)));
    Gi[2] = Gi[2].add(mI.sx.mul(E101).add(mI.sy.mul(E011)).add(mI.sz.mul(E002)));

    // Total polarization energy.
    DoubleVector energy = eI.add(eK).mul(energyMask);

    // Get the induced-induced portion of the force (Ud . dC/dX . Up).
    // This contribution does not exist for direct polarization (mutualMask == 0.0).
    // For SIMD code, we are unable to hide this in an if statement, calculation is still correct with it
    // Find the potential and its derivatives at k due to induced dipole i.
    dipoleIPotentialAtK(mI.ux, mI.uy, mI.uz, 2);
    Gi[0] = Gi[0].sub(mK.px.mul(E200).add(mK.py.mul(E110)).add(mK.pz.mul(E101)).mul(0.5).mul(mutualMask));
    Gi[1] = Gi[1].sub(mK.px.mul(E110).add(mK.py.mul(E020)).add(mK.pz.mul(E011)).mul(0.5).mul(mutualMask));
    Gi[2] = Gi[2].sub(mK.px.mul(E101).add(mK.py.mul(E011)).add(mK.pz.mul(E002)).mul(0.5).mul(mutualMask));

    // Find the potential and its derivatives at i due to induced dipole k.
    dipoleKPotentialAtI(mK.ux, mK.uy, mK.uz, 2);
    Gi[0] = Gi[0].add(mI.px.mul(E200).add(mI.py.mul(E110)).add(mI.pz.mul(E101)).mul(0.5).mul(mutualMask));
    Gi[1] = Gi[1].add(mI.px.mul(E110).add(mI.py.mul(E020)).add(mI.pz.mul(E011)).mul(0.5).mul(mutualMask));
    Gi[2] = Gi[2].add(mI.px.mul(E101).add(mI.py.mul(E011)).add(mI.pz.mul(E002)).mul(0.5).mul(mutualMask));

    // Find the potential and its derivatives at K due to the averaged induced dipole at site i.
    dipoleIPotentialAtK(mI.sx, mI.sy, mI.sz, 2);
    multipoleTorque(mK, Tk);

    // Find the potential and its derivatives at I due to the averaged induced dipole at site k.
    dipoleKPotentialAtI(mK.sx, mK.sy, mK.sz, 2);
    multipoleTorque(mI, Ti);

    return energy;
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
  public DoubleVector multipoleEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                 DoubleVector[] Gi, DoubleVector[] Gk, DoubleVector[] Ti, DoubleVector[] Tk) {
    multipoleIPotentialAtK(mI, 3);
    DoubleVector energy = multipoleEnergy(mK);
    multipoleGradient(mK, Gk);
    Gi[0] = Gi[0].sub(Gk[0]);
    Gi[1] = Gi[1].sub(Gk[1]);
    Gi[2] = Gi[2].sub(Gk[2]);

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
   * Generate source terms for the Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  protected abstract void source(DoubleVector[] T000);

  /**
   * Hard coded tensor computation up to 1st order. This code is auto-generated for both Global and QI frames.
   */
  protected abstract void order1();

  /**
   * Hard coded tensor computation up to 2nd order. This code is auto-generated for both Global and QI frames.
   */
  protected abstract void order2();

  /**
   * Hard coded tensor computation up to 3rd order. This code is auto-generated for both Global and QI frames.
   */
  protected abstract void order3();

  /**
   * Hard coded tensor computation up to 4th order. This code is auto-generated for both Global and QI frames.
   */
  protected abstract void order4();

  /**
   * Hard coded tensor computation up to 5th order. This code is auto-generated for both Global and QI frames.
   * <br>
   * The 5th order recursion is needed for quadrupole-quadrupole forces.
   */
  protected abstract void order5();

  /**
   * Hard coded tensor computation up to 6th order. This code is auto-generated for both Global and QI frames.
   * <br>
   * This is needed for quadrupole-quadrupole forces and orthogonal space sampling.
   */
  protected abstract void order6();

  /**
   * Compute the field components due to site I multipole at site K.
   *
   * @param mI    PolarizableMultipoleSIMD at site I.
   * @param order Potential order.
   */
  @SuppressWarnings("fallthrough")
  protected abstract void multipoleIPotentialAtK(PolarizableMultipoleSIMD mI, int order);

  /**
   * Compute the field components due to site K multipole at site I.
   *
   * @param mK    PolarizableMultipoleSIMD at site I.
   * @param order Potential order.
   */
  @SuppressWarnings("fallthrough")
  protected abstract void multipoleKPotentialAtI(PolarizableMultipoleSIMD mK, int order);

  /**
   * Compute the field components due to site K charge at site I.
   *
   * @param mK    PolarizableMultipoleSIMD at site K.
   * @param order Potential order.
   */
  protected abstract void chargeKPotentialAtI(PolarizableMultipoleSIMD mK, int order);

  /**
   * Compute the induced dipole field components due to site K at site I.
   *
   * @param uxk   X-dipole component.
   * @param uyk   Y-dipole component.
   * @param uzk   Z-dipole component.
   * @param order Potential order.
   */
  protected abstract void dipoleKPotentialAtI(DoubleVector uxk, DoubleVector uyk, DoubleVector uzk, int order);

  /**
   * Compute the field components due to site K quadrupole at site I.
   *
   * @param mK    MultipoleTensorSIMD at site K.
   * @param order Potential order.
   */
  protected abstract void quadrupoleKPotentialAtI(PolarizableMultipoleSIMD mK, int order);

  /**
   * Compute the field components due to site I charge at site K.
   *
   * @param mI    PolarizableMultipoleSIMD at site I.
   * @param order Potential order.
   */
  protected abstract void chargeIPotentialAtK(PolarizableMultipoleSIMD mI, int order);

  /**
   * Compute the induced dipole field components due to site I at site K.
   *
   * @param uxi   X-dipole component.
   * @param uyi   Y-dipole component.
   * @param uzi   Z-dipole component.
   * @param order Potential order.
   */
  protected abstract void dipoleIPotentialAtK(DoubleVector uxi, DoubleVector uyi, DoubleVector uzi, int order);

  /**
   * Compute the field components due to site I quadrupole at site K.
   *
   * @param mI    MultipoleTensorSIMD at site I.
   * @param order Potential order.
   */
  protected abstract void quadrupoleIPotentialAtK(PolarizableMultipoleSIMD mI, int order);

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
