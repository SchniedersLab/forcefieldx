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

import static java.lang.Math.fma;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The CoulombTensorQI class computes derivatives of 1/|<b>r</b>| via recursion to arbitrary order
 * for Cartesian multipoles in a quasi-internal frame.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 *     Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 *     computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 *     Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @since 1.0
 */
public class CoulombTensorQI extends MultipoleTensor {

  public CoulombTensorQI(int order) {
    super(COORDINATES.QI, order);
    operator = OPERATOR.COULOMB;
  }

  /** {@inheritDoc} */
  @Override
  public void setR(double dx, double dy, double dz) {
    x = 0.0;
    y = 0.0;
    r2 = dx * dx + dy * dy + dz * dz;
    z = sqrt(r2);
    R = z;
  }

  /**
   * Generate source terms for the Coulomb Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  protected void source(double[] T000) {
    // Challacombe et al. Equation 21, last factor.
    // == (1/r) * (1/r^3) * (1/r^5) * (1/r^7) * ...
    double ir = 1.0 / R;
    double ir2 = ir * ir;
    for (int n = 0; n < o1; n++) {
      T000[n] = coulombSource[n] * ir;
      ir *= ir2;
    }
  }

  /** {@inheritDoc} */
  @Override
  protected void order1() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    R000 = term0000;
    R001 = z * term0001;
  }

  /** {@inheritDoc} */
  @Override
  protected void order2() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    R000 = term0000;
    R200 = term0001;
    R020 = term0001;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = fma(z, term0011, term0001);
  }

  /** {@inheritDoc} */
  @Override
  protected void order3() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    R000 = term0000;
    R200 = term0001;
    double term2001 = term0002;
    R020 = term0001;
    double term0201 = term0002;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = fma(z, term0011, term0001);
    double term0012 = z * term0003;
    double term0021 = fma(z, term0012, term0002);
    R003 = fma(z, term0021, 2 * term0011);
    R021 = z * term0201;
    R201 = z * term2001;
  }

  /** {@inheritDoc} */
  @Override
  protected void order4() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    double term0004 = work[4];
    R000 = term0000;
    R200 = term0001;
    double term2001 = term0002;
    double term2002 = term0003;
    R400 = 3 * term2001;
    R020 = term0001;
    double term0201 = term0002;
    double term0202 = term0003;
    R040 = 3 * term0201;
    R220 = term2001;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = fma(z, term0011, term0001);
    double term0012 = z * term0003;
    double term0021 = fma(z, term0012, term0002);
    R003 = fma(z, term0021, 2 * term0011);
    double term0013 = z * term0004;
    double term0022 = fma(z, term0013, term0003);
    double term0031 = fma(z, term0022, 2 * term0012);
    R004 = fma(z, term0031, 3 * term0021);
    R021 = z * term0201;
    double term0211 = z * term0202;
    R022 = fma(z, term0211, term0201);
    R201 = z * term2001;
    double term2011 = z * term2002;
    R202 = fma(z, term2011, term2001);
  }

  /** {@inheritDoc} */
  @Override
  protected void order5() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    double term0004 = work[4];
    double term0005 = work[5];
    R000 = term0000;
    R200 = term0001;
    double term2001 = term0002;
    double term2002 = term0003;
    R400 = 3 * term2001;
    double term2003 = term0004;
    double term4001 = 3 * term2002;
    R020 = term0001;
    double term0201 = term0002;
    double term0202 = term0003;
    R040 = 3 * term0201;
    double term0203 = term0004;
    double term0401 = 3 * term0202;
    R220 = term2001;
    double term2201 = term2002;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = fma(z, term0011, term0001);
    double term0012 = z * term0003;
    double term0021 = fma(z, term0012, term0002);
    R003 = fma(z, term0021, 2 * term0011);
    double term0013 = z * term0004;
    double term0022 = fma(z, term0013, term0003);
    double term0031 = fma(z, term0022, 2 * term0012);
    R004 = fma(z, term0031, 3 * term0021);
    double term0014 = z * term0005;
    double term0023 = fma(z, term0014, term0004);
    double term0032 = fma(z, term0023, 2 * term0013);
    double term0041 = fma(z, term0032, 3 * term0022);
    R005 = fma(z, term0041, 4 * term0031);
    R021 = z * term0201;
    double term0211 = z * term0202;
    R022 = fma(z, term0211, term0201);
    double term0212 = z * term0203;
    double term0221 = fma(z, term0212, term0202);
    R023 = fma(z, term0221, 2 * term0211);
    R041 = z * term0401;
    R201 = z * term2001;
    double term2011 = z * term2002;
    R202 = fma(z, term2011, term2001);
    double term2012 = z * term2003;
    double term2021 = fma(z, term2012, term2002);
    R203 = fma(z, term2021, 2 * term2011);
    R221 = z * term2201;
    R401 = z * term4001;
  }

  /** {@inheritDoc} */
  @Override
  protected void order6() {
    source(work);
    double term0000 = work[0];
    double term0001 = work[1];
    double term0002 = work[2];
    double term0003 = work[3];
    double term0004 = work[4];
    double term0005 = work[5];
    double term0006 = work[6];
    R000 = term0000;
    R200 = term0001;
    double term2001 = term0002;
    double term2002 = term0003;
    R400 = 3 * term2001;
    double term2003 = term0004;
    double term4001 = 3 * term2002;
    double term2004 = term0005;
    double term4002 = 3 * term2003;
    R600 = 5 * term4001;
    R020 = term0001;
    double term0201 = term0002;
    double term0202 = term0003;
    R040 = 3 * term0201;
    double term0203 = term0004;
    double term0401 = 3 * term0202;
    double term0204 = term0005;
    double term0402 = 3 * term0203;
    R060 = 5 * term0401;
    R220 = term2001;
    double term2201 = term2002;
    double term2202 = term2003;
    R240 = 3 * term2201;
    R420 = term4001;
    R001 = z * term0001;
    double term0011 = z * term0002;
    R002 = fma(z, term0011, term0001);
    double term0012 = z * term0003;
    double term0021 = fma(z, term0012, term0002);
    R003 = fma(z, term0021, 2 * term0011);
    double term0013 = z * term0004;
    double term0022 = fma(z, term0013, term0003);
    double term0031 = fma(z, term0022, 2 * term0012);
    R004 = fma(z, term0031, 3 * term0021);
    double term0014 = z * term0005;
    double term0023 = fma(z, term0014, term0004);
    double term0032 = fma(z, term0023, 2 * term0013);
    double term0041 = fma(z, term0032, 3 * term0022);
    R005 = fma(z, term0041, 4 * term0031);
    double term0015 = z * term0006;
    double term0024 = fma(z, term0015, term0005);
    double term0033 = fma(z, term0024, 2 * term0014);
    double term0042 = fma(z, term0033, 3 * term0023);
    double term0051 = fma(z, term0042, 4 * term0032);
    R006 = fma(z, term0051, 5 * term0041);
    R021 = z * term0201;
    double term0211 = z * term0202;
    R022 = fma(z, term0211, term0201);
    double term0212 = z * term0203;
    double term0221 = fma(z, term0212, term0202);
    R023 = fma(z, term0221, 2 * term0211);
    double term0213 = z * term0204;
    double term0222 = fma(z, term0213, term0203);
    double term0231 = fma(z, term0222, 2 * term0212);
    R024 = fma(z, term0231, 3 * term0221);
    R041 = z * term0401;
    double term0411 = z * term0402;
    R042 = fma(z, term0411, term0401);
    R201 = z * term2001;
    double term2011 = z * term2002;
    R202 = fma(z, term2011, term2001);
    double term2012 = z * term2003;
    double term2021 = fma(z, term2012, term2002);
    R203 = fma(z, term2021, 2 * term2011);
    double term2013 = z * term2004;
    double term2022 = fma(z, term2013, term2003);
    double term2031 = fma(z, term2022, 2 * term2012);
    R204 = fma(z, term2031, 3 * term2021);
    R221 = z * term2201;
    double term2211 = z * term2202;
    R222 = fma(z, term2211, term2201);
    R401 = z * term4001;
    double term4011 = z * term4002;
    R402 = fma(z, term4011, term4001);
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void multipoleIPotentialAtK(PolarizableMultipole mI, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3.
        double term300 = 0.0;
        term300 = fma(mI.dx, -R400, term300);
        term300 = fma(mI.qxz, R401, term300);
        E300 = term300;
        double term030 = 0.0;
        term030 = fma(mI.dy, -R040, term030);
        term030 = fma(mI.qyz, R041, term030);
        E030 = term030;
        double term003 = 0.0;
        term003 = fma(mI.q, R003, term003);
        term003 = fma(mI.dz, -R004, term003);
        term003 = fma(mI.qxx, R203, term003);
        term003 = fma(mI.qyy, R023, term003);
        term003 = fma(mI.qzz, R005, term003);
        E003 = term003;
        double term210 = 0.0;
        term210 = fma(mI.dy, -R220, term210);
        term210 = fma(mI.qyz, R221, term210);
        E210 = term210;
        double term201 = 0.0;
        term201 = fma(mI.q, R201, term201);
        term201 = fma(mI.dz, -R202, term201);
        term201 = fma(mI.qxx, R401, term201);
        term201 = fma(mI.qyy, R221, term201);
        term201 = fma(mI.qzz, R203, term201);
        E201 = term201;
        double term120 = 0.0;
        term120 = fma(mI.dx, -R220, term120);
        term120 = fma(mI.qxz, R221, term120);
        E120 = term120;
        double term021 = 0.0;
        term021 = fma(mI.q, R021, term021);
        term021 = fma(mI.dz, -R022, term021);
        term021 = fma(mI.qxx, R221, term021);
        term021 = fma(mI.qyy, R041, term021);
        term021 = fma(mI.qzz, R023, term021);
        E021 = term021;
        double term102 = 0.0;
        term102 = fma(mI.dx, -R202, term102);
        term102 = fma(mI.qxz, R203, term102);
        E102 = term102;
        double term012 = 0.0;
        term012 = fma(mI.dy, -R022, term012);
        term012 = fma(mI.qyz, R023, term012);
        E012 = term012;
        double term111 = 0.0;
        term111 = fma(mI.qxy, R221, term111);
        E111 = term111;
        // Fall through to 2nd order.
      case 2:
        // Order 2.
        double term200 = 0.0;
        term200 = fma(mI.q, R200, term200);
        term200 = fma(mI.dz, -R201, term200);
        term200 = fma(mI.qxx, R400, term200);
        term200 = fma(mI.qyy, R220, term200);
        term200 = fma(mI.qzz, R202, term200);
        E200 = term200;
        double term020 = 0.0;
        term020 = fma(mI.q, R020, term020);
        term020 = fma(mI.dz, -R021, term020);
        term020 = fma(mI.qxx, R220, term020);
        term020 = fma(mI.qyy, R040, term020);
        term020 = fma(mI.qzz, R022, term020);
        E020 = term020;
        double term002 = 0.0;
        term002 = fma(mI.q, R002, term002);
        term002 = fma(mI.dz, -R003, term002);
        term002 = fma(mI.qxx, R202, term002);
        term002 = fma(mI.qyy, R022, term002);
        term002 = fma(mI.qzz, R004, term002);
        E002 = term002;
        double term110 = 0.0;
        term110 = fma(mI.qxy, R220, term110);
        E110 = term110;
        double term101 = 0.0;
        term101 = fma(mI.dx, -R201, term101);
        term101 = fma(mI.qxz, R202, term101);
        E101 = term101;
        double term011 = 0.0;
        term011 = fma(mI.dy, -R021, term011);
        term011 = fma(mI.qyz, R022, term011);
        E011 = term011;
        // Fall through to 1st order.
      case 1:
        // Order 1.
        double term100 = 0.0;
        term100 = fma(mI.dx, -R200, term100);
        term100 = fma(mI.qxz, R201, term100);
        E100 = term100;
        double term010 = 0.0;
        term010 = fma(mI.dy, -R020, term010);
        term010 = fma(mI.qyz, R021, term010);
        E010 = term010;
        double term001 = 0.0;
        term001 = fma(mI.q, R001, term001);
        term001 = fma(mI.dz, -R002, term001);
        term001 = fma(mI.qxx, R201, term001);
        term001 = fma(mI.qyy, R021, term001);
        term001 = fma(mI.qzz, R003, term001);
        E001 = term001;
        // Fall through to the potential.
      case 0:
        double term000 = 0.0;
        term000 = fma(mI.q, R000, term000);
        term000 = fma(mI.dz, -R001, term000);
        term000 = fma(mI.qxx, R200, term000);
        term000 = fma(mI.qyy, R020, term000);
        term000 = fma(mI.qzz, R002, term000);
        E000 = term000;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void chargeIPotentialAtK(PolarizableMultipole mI, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3.
        E300 = 0.0;
        E030 = 0.0;
        E003 = mI.q * R003;
        E210 = 0.0;
        E201 = mI.q * R201;
        E120 = 0.0;
        E021 = mI.q * R021;
        E102 = 0.0;
        E012 = 0.0;
        E111 = 0.0;
        // Fall through to 2nd order.
      case 2:
        // Order 2.
        E200 = mI.q * R200;
        E020 = mI.q * R020;
        E002 = mI.q * R002;
        E110 = 0.0;
        E101 = 0.0;
        E011 = 0.0;
        // Fall through to 1st order.
      case 1:
        // Order 1.
        E100 = 0.0;
        E010 = 0.0;
        E001 = mI.q * R001;
        // Fall through to the potential.
      case 0:
        E000 = mI.q * R000;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void dipoleIPotentialAtK(double uxi, double uyi, double uzi, int order) {
    switch (order) {
      case 3:
      default:
        // Order 3
        E300 = -uxi * R400;
        E030 = -uyi * R040;
        E003 = -uzi * R004;
        E210 = -uyi * R220;
        E201 = -uzi * R202;
        E120 = -uxi * R220;
        E021 = -uzi * R022;
        E102 = -uxi * R202;
        E012 = -uyi * R022;
        E111 = 0.0;
        // Fall through to 2nd order.
      case 2:
        // Order 2.
        E200 = -uzi * R201;
        E020 = -uzi * R021;
        E002 = -uzi * R003;
        E110 = 0.0;
        E101 = -uxi * R201;
        E011 = -uyi * R021;
        // Fall through to 1st order.
      case 1:
        // Order 1.
        E100 = -uxi * R200;
        E010 = -uyi * R020;
        E001 = -uzi * R002;
        // Fall through to the potential.
      case 0:
        E000 = -uzi * R001;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void quadrupoleIPotentialAtK(PolarizableMultipole mI, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3.
        E300 = mI.qxz * R401;
        E030 = mI.qyz * R041;
        double term003 = 0.0;
        term003 = fma(mI.qxx, R203, term003);
        term003 = fma(mI.qyy, R023, term003);
        term003 = fma(mI.qzz, R005, term003);
        E003 = term003;
        E210 = mI.qyz * R221;
        double term201 = 0.0;
        term201 = fma(mI.qxx, R401, term201);
        term201 = fma(mI.qyy, R221, term201);
        term201 = fma(mI.qzz, R203, term201);
        E201 = term201;
        E120 = mI.qxz * R221;
        double term021 = 0.0;
        term021 = fma(mI.qxx, R221, term021);
        term021 = fma(mI.qyy, R041, term021);
        term021 = fma(mI.qzz, R023, term021);
        E021 = term021;
        E102 = mI.qxz * R203;
        E012 = mI.qyz * R023;
        E111 = mI.qxy * R221;
        // Fall through to 2nd order.
      case 2:
        // Order 2.
        double term200 = 0.0;
        term200 = fma(mI.qxx, R400, term200);
        term200 = fma(mI.qyy, R220, term200);
        term200 = fma(mI.qzz, R202, term200);
        E200 = term200;
        double term020 = 0.0;
        term020 = fma(mI.qxx, R220, term020);
        term020 = fma(mI.qyy, R040, term020);
        term020 = fma(mI.qzz, R022, term020);
        E020 = term020;
        double term002 = 0.0;
        term002 = fma(mI.qxx, R202, term002);
        term002 = fma(mI.qyy, R022, term002);
        term002 = fma(mI.qzz, R004, term002);
        E002 = term002;
        E110 = mI.qxy * R220;
        E101 = mI.qxz * R202;
        E011 = mI.qyz * R022;
        // Fall through to 1st order.
      case 1:
        // Order 1.
        E100 = mI.qxz * R201;
        E010 = mI.qyz * R021;
        double term001 = 0.0;
        term001 = fma(mI.qxx, R201, term001);
        term001 = fma(mI.qyy, R021, term001);
        term001 = fma(mI.qzz, R003, term001);
        E001 = term001;
        // Fall through to the potential.
      case 0:
        double term000 = 0.0;
        term000 = fma(mI.qxx, R200, term000);
        term000 = fma(mI.qyy, R020, term000);
        term000 = fma(mI.qzz, R002, term000);
        E000 = term000;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void multipoleKPotentialAtI(PolarizableMultipole mK, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        double term300 = 0.0;
        term300 = fma(mK.dx, R400, term300);
        term300 = fma(mK.qxz, R401, term300);
        E300 = -term300;
        double term030 = 0.0;
        term030 = fma(mK.dy, R040, term030);
        term030 = fma(mK.qyz, R041, term030);
        E030 = -term030;
        double term003 = 0.0;
        term003 = fma(mK.q, R003, term003);
        term003 = fma(mK.dz, R004, term003);
        term003 = fma(mK.qxx, R203, term003);
        term003 = fma(mK.qyy, R023, term003);
        term003 = fma(mK.qzz, R005, term003);
        E003 = -term003;
        double term210 = 0.0;
        term210 = fma(mK.dy, R220, term210);
        term210 = fma(mK.qyz, R221, term210);
        E210 = -term210;
        double term201 = 0.0;
        term201 = fma(mK.q, R201, term201);
        term201 = fma(mK.dz, R202, term201);
        term201 = fma(mK.qxx, R401, term201);
        term201 = fma(mK.qyy, R221, term201);
        term201 = fma(mK.qzz, R203, term201);
        E201 = -term201;
        double term120 = 0.0;
        term120 = fma(mK.dx, R220, term120);
        term120 = fma(mK.qxz, R221, term120);
        E120 = -term120;
        double term021 = 0.0;
        term021 = fma(mK.q, R021, term021);
        term021 = fma(mK.dz, R022, term021);
        term021 = fma(mK.qxx, R221, term021);
        term021 = fma(mK.qyy, R041, term021);
        term021 = fma(mK.qzz, R023, term021);
        E021 = -term021;
        double term102 = 0.0;
        term102 = fma(mK.dx, R202, term102);
        term102 = fma(mK.qxz, R203, term102);
        E102 = -term102;
        double term012 = 0.0;
        term012 = fma(mK.dy, R022, term012);
        term012 = fma(mK.qyz, R023, term012);
        E012 = -term012;
        double term111 = 0.0;
        term111 = fma(mK.qxy, R221, term111);
        E111 = -term111;
        // Fall through to 2nd order.
      case 2:
        // Order 2
        double term200 = 0.0;
        term200 = fma(mK.q, R200, term200);
        term200 = fma(mK.dz, R201, term200);
        term200 = fma(mK.qxx, R400, term200);
        term200 = fma(mK.qyy, R220, term200);
        term200 = fma(mK.qzz, R202, term200);
        E200 = term200;
        double term020 = 0.0;
        term020 = fma(mK.q, R020, term020);
        term020 = fma(mK.dz, R021, term020);
        term020 = fma(mK.qxx, R220, term020);
        term020 = fma(mK.qyy, R040, term020);
        term020 = fma(mK.qzz, R022, term020);
        E020 = term020;
        double term002 = 0.0;
        term002 = fma(mK.q, R002, term002);
        term002 = fma(mK.dz, R003, term002);
        term002 = fma(mK.qxx, R202, term002);
        term002 = fma(mK.qyy, R022, term002);
        term002 = fma(mK.qzz, R004, term002);
        E002 = term002;
        double term110 = 0.0;
        term110 = fma(mK.qxy, R220, term110);
        E110 = term110;
        double term101 = 0.0;
        term101 = fma(mK.dx, R201, term101);
        term101 = fma(mK.qxz, R202, term101);
        E101 = term101;
        double term011 = 0.0;
        term011 = fma(mK.dy, R021, term011);
        term011 = fma(mK.qyz, R022, term011);
        E011 = term011;
        // Fall through to 1st order.
      case 1:
        // Order 1
        double term100 = 0.0;
        term100 = fma(mK.dx, R200, term100);
        term100 = fma(mK.qxz, R201, term100);
        E100 = -term100;
        double term010 = 0.0;
        term010 = fma(mK.dy, R020, term010);
        term010 = fma(mK.qyz, R021, term010);
        E010 = -term010;
        double term001 = 0.0;
        term001 = fma(mK.q, R001, term001);
        term001 = fma(mK.dz, R002, term001);
        term001 = fma(mK.qxx, R201, term001);
        term001 = fma(mK.qyy, R021, term001);
        term001 = fma(mK.qzz, R003, term001);
        E001 = -term001;
      case 0:
        double term000 = 0.0;
        term000 = fma(mK.q, R000, term000);
        term000 = fma(mK.dz, R001, term000);
        term000 = fma(mK.qxx, R200, term000);
        term000 = fma(mK.qyy, R020, term000);
        term000 = fma(mK.qzz, R002, term000);
        E000 = term000;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void chargeKPotentialAtI(PolarizableMultipole mK, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        E300 = 0.0;
        E030 = 0.0;
        E003 = -mK.q * R003;
        E210 = 0.0;
        E201 = -mK.q * R201;
        E120 = 0.0;
        E021 = -mK.q * R021;
        E102 = 0.0;
        E012 = 0.0;
        E111 = 0.0;
        // Fall through to 2nd order.
      case 2:
        // Order 2
        E200 = mK.q * R200;
        E020 = mK.q * R020;
        E002 = mK.q * R002;
        E110 = 0.0;
        E101 = 0.0;
        E011 = 0.0;
        // Fall through to 1st order.
      case 1:
        // Order 1
        E100 = 0.0;
        E010 = 0.0;
        E001 = -mK.q * R001;
      case 0:
        E000 = mK.q * R000;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void dipoleKPotentialAtI(double uxk, double uyk, double uzk, int order) {
    switch (order) {
      case 3:
      default:
        // Order 3
        E300 = -uxk * R400;
        E030 = -uyk * R040;
        E003 = -uzk * R004;
        E210 = -uyk * R220;
        E201 = -uzk * R202;
        E120 = -uxk * R220;
        E021 = -uzk * R022;
        E102 = -uxk * R202;
        E012 = -uyk * R022;
        E111 = 0.0;
        // Fall through to 2nd order.
      case 2:
        E200 = uzk * R201;
        E020 = uzk * R021;
        E002 = uzk * R003;
        E110 = 0.0;
        E101 = uxk * R201;
        E011 = uyk * R021;
      case 1:
        E100 = -uxk * R200;
        E010 = -uyk * R020;
        E001 = -uzk * R002;
      case 0:
        E000 = uzk * R001;
    }
  }

  /** {@inheritDoc} */
  @SuppressWarnings("fallthrough")
  @Override
  protected void quadrupoleKPotentialAtI(PolarizableMultipole mK, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        E300 = -mK.qxz * R401;
        E030 = -mK.qyz * R041;
        double term003 = 0.0;
        term003 = fma(mK.qxx, R203, term003);
        term003 = fma(mK.qyy, R023, term003);
        term003 = fma(mK.qzz, R005, term003);
        E003 = -term003;
        E210 = -mK.qyz * R221;
        double term201 = 0.0;
        term201 = fma(mK.qxx, R401, term201);
        term201 = fma(mK.qyy, R221, term201);
        term201 = fma(mK.qzz, R203, term201);
        E201 = -term201;
        E120 = -mK.qxz * R221;
        double term021 = 0.0;
        term021 = fma(mK.qxx, R221, term021);
        term021 = fma(mK.qyy, R041, term021);
        term021 = fma(mK.qzz, R023, term021);
        E021 = -term021;
        E102 = -mK.qxz * R203;
        E012 = -mK.qyz * R023;
        E111 = -mK.qxy * R221;
        // Fall through to 2nd order.
      case 2:
        // Order 2
        double term200 = 0.0;
        term200 = fma(mK.qxx, R400, term200);
        term200 = fma(mK.qyy, R220, term200);
        term200 = fma(mK.qzz, R202, term200);
        E200 = term200;
        double term020 = 0.0;
        term020 = fma(mK.qxx, R220, term020);
        term020 = fma(mK.qyy, R040, term020);
        term020 = fma(mK.qzz, R022, term020);
        E020 = term020;
        double term002 = 0.0;
        term002 = fma(mK.qxx, R202, term002);
        term002 = fma(mK.qyy, R022, term002);
        term002 = fma(mK.qzz, R004, term002);
        E002 = term002;
        E110 = mK.qxy * R220;
        E101 = mK.qxz * R202;
        E011 = mK.qyz * R022;
        // Fall through to 1st order.
      case 1:
        // Order 1
        E100 = -mK.qxz * R201;
        E010 = -mK.qyz * R021;
        double term001 = 0.0;
        term001 = fma(mK.qxx, R201, term001);
        term001 = fma(mK.qyy, R021, term001);
        term001 = fma(mK.qzz, R003, term001);
        E001 = -term001;
      case 0:
        double term000 = 0.0;
        term000 = fma(mK.qxx, R200, term000);
        term000 = fma(mK.qyy, R020, term000);
        term000 = fma(mK.qzz, R002, term000);
        E000 = term000;
    }
  }

  /** {@inheritDoc} */
  @Override
  protected double Tlmnj(
      final int l, final int m, final int n, final int j, final double[] r, final double[] T000) {
    double z = r[2];
    assert (r[0] == 0.0 && r[1] == 0.0);

    if (m == 0 && n == 0) {
      if (l > 1) {
        return (l - 1) * Tlmnj(l - 2, 0, 0, j + 1, r, T000);
      } else if (l == 1) { // l == 1, d/dx is done.
        return 0.0;
      } else { // l = m = n = 0. Recursion is done.
        return T000[j];
      }
    } else if (n == 0) { // m >= 1
      if (m > 1) {
        return (m - 1) * Tlmnj(l, m - 2, 0, j + 1, r, T000);
      }
      return 0.0;
    } else { // n >= 1
      if (n > 1) {
        return z * Tlmnj(l, m, n - 1, j + 1, r, T000)
            + (n - 1) * Tlmnj(l, m, n - 2, j + 1, r, T000);
      }
      return z * Tlmnj(l, m, 0, j + 1, r, T000);
    }
  }

  /**
   * This method is a driver to collect elements of the Cartesian multipole tensor given the
   * recursion relationships implemented by the method "Tlmnj", which can be called directly to get a
   * single tensor element. It does not store intermediate values of the recursion, causing it to
   * scale O(order^8). For order = 5, this approach is a factor of 10 slower than recursion that
   * stores intermediates.
   *
   * @param r double[] vector between two sites. r[0] and r[1] must equal 0.0.
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   */
  @Override
  protected void noStorageRecursion(double[] r, double[] tensor) {
    setR(r);
    noStorageRecursion(tensor);
  }

  /**
   * This method is a driver to collect elements of the Cartesian multipole tensor given the
   * recursion relationships implemented by the method "Tlmnj", which can be called directly to get a
   * single tensor element. It does not store intermediate values of the recursion, causing it to
   * scale O(order^8). For order = 5, this approach is a factor of 10 slower than recursion that
   * stores intermediates.
   *
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   */
  @Override
  protected void noStorageRecursion(double[] tensor) {
    assert (x == 0.0 && y == 0.0);
    double[] r = {x, y, z};
    source(T000);
    // 1/r
    tensor[0] = T000[0];
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    for (int l = 1; l <= order; l++) {
      tensor[ti(l, 0, 0)] = Tlmnj(l, 0, 0, 0, r, T000);
    }
    // Find (d/dx)^l * (d/dy)^m for l + m = 1..order (m >= 1, n = 0)
    for (int l = 0; l <= o1; l++) {
      for (int m = 1; m <= order - l; m++) {
        tensor[ti(l, m, 0)] = Tlmnj(l, m, 0, 0, r, T000);
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
    for (int l = 0; l <= o1; l++) {
      for (int m = 0; m <= o1 - l; m++) {
        for (int n = 1; n <= order - l - m; n++) {
          tensor[ti(l, m, n)] = Tlmnj(l, m, n, 0, r, T000);
        }
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  protected void recursion(final double[] r, final double[] tensor) {
    setR(r);
    recursion(tensor);
  }

  /** {@inheritDoc} */
  @Override
  protected void recursion(final double[] tensor) {
    assert (x == 0.0 && y == 0.0);
    source(work);
    tensor[0] = work[0];
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    // Any (d/dx) term can be formed as
    // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
    // All intermediate terms are indexed as l*il + m*im + n*in + j;
    // Store the l=1 tensor T100 (d/dx)
    // Starting the loop at l=2 avoids an if statement.
    double current;
    double previous = work[1];
    tensor[ti(1, 0, 0)] = 0.0;
    for (int l = 2; l < o1; l++) {
      // Initial condition for the inner loop is formation of T100(l-1).
      // Starting the inner loop at a=1 avoids an if statement.
      // T100(l-1) = 0.0 * T000(l)
      current = 0.0;
      int iw = il + l - 1;
      work[iw] = current;
      for (int a = 1; a < l - 1; a++) {
        // T200(l-2) = 0.0 * T100(l-1) + (2 - 1) * T000(l-1)
        // T300(l-3) = 0.0 * T200(l-2) + (3 - 1) * T100(l-2)
        // ...
        // T(l-1)001 = 0.0 * T(l-2)002 + (l - 2) * T(l-3)002
        current = a * work[iw - il];
        iw += il - 1;
        work[iw] = current;
      }
      // Store the Tl00 tensor (d/dx)^l
      // Tl00 = 0.0 * [[ T(l-1)001 ]] + (l - 1) * T(l-2)001
      tensor[ti(l, 0, 0)] = (l - 1) * previous;
      previous = current;
    }
    // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
    // Any (d/dy) term can be formed as:
    // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
    for (int l = 0; l < order; l++) {
      // Store the m=1 tensor (d/dx)^l *(d/dy)
      // Tl10 = y * Tl001
      previous = work[l * il + 1];
      tensor[ti(l, 1, 0)] = 0.0;
      for (int m = 2; m + l < o1; m++) {
        // Tl10(m-1) = y * Tl00m;
        int iw = l * il + m;
        current = 0.0;
        iw += im - 1;
        work[iw] = current;
        for (int a = 1; a < m - 1; a++) {
          // Tl20(m-2) = 0.0 * Tl10(m-1) + (2 - 1) * T100(m-1)
          // Tl30(m-3) = 0.0 * Tl20(m-2) + (3 - 1) * Tl10(m-2)
          // ...
          // Tl(m-1)01 = 0.0 * Tl(m-2)02 + (m - 2) * Tl(m-3)02
          current = a * work[iw - im];
          iw += im - 1;
          work[iw] = current;
        }
        // Store the tensor (d/dx)^l * (d/dy)^m
        // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
        tensor[ti(l, m, 0)] = (m - 1) * previous;
        previous = current;
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
    // Any (d/dz) term can be formed as:
    // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m + l < order; m++) {
        // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
        // Tlmn = z * Tlm01
        final int lm = m + l;
        final int lilmim = l * il + m * im;
        previous = work[lilmim + 1];
        tensor[ti(l, m, 1)] = z * previous;
        for (int n = 2; lm + n < o1; n++) {
          // Tlm1(n-1) = z * Tlm0n;
          int iw = lilmim + n;
          current = z * work[iw];
          iw += in - 1;
          work[iw] = current;
          final int n1 = n - 1;
          for (int a = 1; a < n1; a++) {
            // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
            // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
            // ...
            // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
            current = z * current + a * work[iw - in];
            iw += in - 1;
            work[iw] = current;
          }
          // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
          // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
          tensor[ti(l, m, n)] = z * current + n1 * previous;
          previous = current;
        }
      }
    }
  }

  /**
   * This function is a driver to collect elements of the Cartesian multipole tensor. Collecting all
   * tensors scales slightly better than O(order^4).
   *
   * <p>For a multipole expansion truncated at quadrupole order, for example, up to order 5 is
   * needed for energy gradients. The number of terms this requires is binomial(5 + 3, 3) or 8! / (5!
   * * 3!), which is 56.
   *
   * <p>The packing of the tensor elements for order = 1<br>
   * tensor[0] = 1/|r| <br> tensor[1] = -x/|r|^3 <br> tensor[2] = -y/|r|^3 <br> tensor[3] = -z/|r|^3
   * <br>
   *
   * <p>
   *
   * @param r double[] vector between two sites.
   * @param tensor double[] length must be at least binomial(order + 3, 3).
   * @return Java code for the tensor recursion.
   * @since 1.0
   */
  @Override
  protected String codeTensorRecursion(final double[] r, final double[] tensor) {
    setR(r);
    assert (x == 0.0 && y == 0.0);
    source(work);
    StringBuilder sb = new StringBuilder();
    tensor[0] = work[0];
    if (work[0] > 0) {
      sb.append(format("%s = %s;\n", rlmn(0, 0, 0), term(0, 0, 0, 0)));
    }
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    // Any (d/dx) term can be formed as
    // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
    // All intermediate terms are indexed as l*il + m*im + n*in + j;
    // Store the l=1 tensor T100 (d/dx)
    // Starting the loop at l=2 avoids an if statement.
    double current;
    double previous = work[1];
    tensor[ti(1, 0, 0)] = 0.0;
    for (int l = 2; l < o1; l++) {
      // Initial condition for the inner loop is formation of T100(l-1).
      // Starting the inner loop at a=1 avoids an if statement.
      // T100(l-1) = 0.0 * T000(l)
      current = 0.0;
      int iw = il + (l - 1);
      work[iw] = current;
      // sb.append(format("double %s = 0.0;\n", term(1, 0, 0, l - 1)));
      for (int a = 1; a < l - 1; a++) {
        // T200(l-2) = 0.0 * T100(l-1) + (2 - 1) * T000(l-1)
        // T300(l-3) = 0.0 * T200(l-2) + (3 - 1) * T100(l-2)
        // ...
        // T(l-1)001 = 0.0 * T(l-2)002 + (l - 2) * T(l-3)002
        // iw = (a - 1) * il + (l - a)
        current = a * work[iw - il];
        iw += il - 1;
        // iw = (a + 1) * il + (l - a - 1)
        work[iw] = current;
        if (current != 0) {
          if (a > 2) {
            sb.append(
                format(
                    "double %s = %d * %s;\n",
                    term(a + 1, 0, 0, l - a - 1), a, term(a - 1, 0, 0, l - a)));
          } else {
            sb.append(
                format(
                    "double %s = %s;\n", term(a + 1, 0, 0, l - a - 1), term(a - 1, 0, 0, l - a)));
          }
        }
      }
      // Store the Tl00 tensor (d/dx)^l
      // Tl00 = 0.0 * T(l-1)001 + (l - 1) * T(l-2)001
      int index = ti(l, 0, 0);
      tensor[index] = (l - 1) * previous;
      previous = current;
      if (tensor[index] != 0) {
        if (l > 2) {
          sb.append(format("%s = %d * %s;\n", rlmn(l, 0, 0), (l - 1), term(l - 2, 0, 0, 1)));
        } else {
          sb.append(format("%s = %s;\n", rlmn(l, 0, 0), term(l - 2, 0, 0, 1)));
        }
      }
    }
    // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
    // Any (d/dy) term can be formed as:
    // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
    for (int l = 0; l < order; l++) {
      // Store the m=1 tensor (d/dx)^l *(d/dy)
      // Tl10 = y * Tl001
      previous = work[l * il + 1];
      tensor[ti(l, 1, 0)] = 0.0;
      for (int m = 2; m + l < o1; m++) {
        // Tl10(m-1) = y * Tl00m;
        int iw = l * il + m;
        current = 0.0;
        iw += im - 1;
        work[iw] = current;
        for (int a = 1; a < m - 1; a++) {
          // Tl20(m-2) = 0.0 * Tl10(m-1) + (2 - 1) * Tl00(m-1)
          // Tl30(m-3) = 0.0 * Tl20(m-2) + (3 - 1) * Tl10(m-2)
          // ...
          // Tl(m-1)01 = 0.0 * Tl(m-2)02 + (m - 2) * Tl(m-3)02
          current = a * work[iw - im];
          iw += im - 1;
          work[iw] = current;
          if (current != 0) {
            if (a > 1) {
              sb.append(
                  format(
                      "double %s = %d * %s;\n",
                      term(l, a + 1, 0, m - a - 1), a, term(l, a - 1, 0, m - a)));
            } else {
              sb.append(
                  format(
                      "double %s = %s;\n", term(l, a + 1, 0, m - a - 1), term(l, a - 1, 0, m - a)));
            }
          }
        }
        // Store the tensor (d/dx)^l * (d/dy)^m
        // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
        int index = ti(l, m, 0);
        tensor[index] = (m - 1) * previous;
        previous = current;
        if (tensor[index] != 0) {
          if (m > 2) {
            sb.append(format("%s = %d * %s;\n", rlmn(l, m, 0), (m - 1), term(l, m - 2, 0, 1)));
          } else {
            sb.append(format("%s = %s;\n", rlmn(l, m, 0), term(l, m - 2, 0, 1)));
          }
        }
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
    // Any (d/dz) term can be formed as:
    // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m + l < order; m++) {
        // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
        // Tlmn = z * Tlm01
        final int lm = m + l;
        final int lilmim = l * il + m * im;
        previous = work[lilmim + 1];
        int index = ti(l, m, 1);
        tensor[index] = z * previous;
        if (tensor[index] != 0) {
          sb.append(format("%s = z * %s;\n", rlmn(l, m, 1), term(l, m, 0, 1)));
        }
        for (int n = 2; lm + n < o1; n++) {
          // Tlm1(n-1) = z * Tlm0n;
          int iw = lilmim + n;
          current = z * work[iw];
          iw += in - 1;
          work[iw] = current;
          if (current != 0) {
            sb.append(format("double %s = z * %s;\n", term(l, m, 1, n - 1), term(l, m, 0, n)));
          }
          final int n1 = n - 1;
          for (int a = 1; a < n1; a++) {
            // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
            // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
            // ...
            // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
            current = z * current + a * work[iw - in];
            iw += in - 1;
            work[iw] = current;
            if (current != 0) {
              if (a > 1) {
                sb.append(
                    format(
                        "double %s = fma(z, %s, %d * %s);\n",
                        term(l, m, a + 1, n - a - 1),
                        term(l, m, a, n - a),
                        a,
                        term(l, m, a - 1, n - a)));
              } else {
                sb.append(
                    format(
                        "double %s = fma(z, %s, %s);\n",
                        term(l, m, a + 1, n - a - 1),
                        term(l, m, a, n - a),
                        term(l, m, a - 1, n - a)));
              }
            }
          }
          // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
          // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
          index = ti(l, m, n);
          tensor[index] = z * current + n1 * previous;
          previous = current;
          if (tensor[index] != 0) {
            if (n > 2) {
              sb.append(
                  format(
                      "%s = fma(z, %s, %d * %s);\n",
                      rlmn(l, m, n), term(l, m, n - 1, 1), (n - 1), term(l, m, n - 2, 1)));
            } else {
              sb.append(
                  format(
                      "%s = fma(z, %s, %s);\n",
                      rlmn(l, m, n), term(l, m, n - 1, 1), term(l, m, n - 2, 1)));
            }
          }
        }
      }
    }
    return sb.toString();
  }

}
