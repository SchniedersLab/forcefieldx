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

/**
 * The CoulombTensorQISIMD class computes derivatives of 1/|<b>r</b>| via recursion to arbitrary order
 * for Cartesian multipoles in a quasi-internal frame using SIMD instructions.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 * Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 * computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 * Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @since 1.0
 */
public class CoulombTensorQISIMD extends MultipoleTensorSIMD {

  private static final DoubleVector zero = DoubleVector.zero(DoubleVector.SPECIES_PREFERRED);

  /**
   * Create a new CoulombTensorQI object.
   *
   * @param order The tensor order.
   */
  public CoulombTensorQISIMD(int order) {
    super(CoordinateSystem.QI, order);
    operator = Operator.COULOMB;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setR(DoubleVector dx, DoubleVector dy, DoubleVector dz) {
    // QI separation along x is zero.
    x = DoubleVector.broadcast(dz.species(), 0.0);
    // QI separation along y is zero.
    y = DoubleVector.broadcast(dz.species(), 0.0);
    // Compute the separation distance squared.
    DoubleVector dx2 = dx.mul(dx);
    DoubleVector dy2 = dy.mul(dy);
    DoubleVector dz2 = dz.mul(dz);
    r2 = dx2.add(dy2).add(dz2);
    // QI separation vector is along z.
    z = r2.sqrt();
    R = z;
  }

  /**
   * Generate source terms for the Coulomb Challacombe et al. recursion.
   *
   * @param T000 Location to store the source terms.
   */
  protected void source(DoubleVector[] T000) {
    // Challacombe et al. Equation 21, last factor.
    // == (1/r) * (1/r^3) * (1/r^5) * (1/r^7) * ...
    DoubleVector ONE = DoubleVector.broadcast(R.species(), 1.0);
    DoubleVector ir = ONE.div(R);
    DoubleVector ir2 = ir.mul(ir);
    for (int n = 0; n < o1; n++) {
      T000[n] = ir.mul(coulombSource[n]);
      ir = ir.mul(ir2);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  protected void order1() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    R000 = term0000;
    R001 = z.mul(term0001);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  protected void order2() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    R000 = term0000;
    R200 = term0001;
    R020 = term0001;
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  protected void order3() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    DoubleVector term0003 = work[3];
    R000 = term0000;
    R200 = term0001;
    R020 = term0001;
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
    DoubleVector term0012 = z.mul(term0003);
    DoubleVector term0021 = z.fma(term0012, term0002);
    R003 = z.fma(term0021, term0011.mul(2));
    R021 = z.mul(term0002);
    R201 = z.mul(term0002);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  protected void order4() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    DoubleVector term0003 = work[3];
    DoubleVector term0004 = work[4];
    R000 = term0000;
    R200 = term0001;
    R400 = term0002.mul(3);
    R020 = term0001;
    R040 = term0002.mul(3);
    R220 = term0002;
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
    DoubleVector term0012 = z.mul(term0003);
    DoubleVector term0021 = z.fma(term0012, term0002);
    R003 = z.fma(term0021, term0011.mul(2));
    DoubleVector term0013 = z.mul(term0004);
    DoubleVector term0022 = z.fma(term0013, term0003);
    DoubleVector term0031 = z.fma(term0022, term0012.mul(2));
    R004 = z.fma(term0031, term0021.mul(3));
    R021 = z.mul(term0002);
    DoubleVector term0211 = z.mul(term0003);
    R022 = z.fma(term0211, term0002);
    R201 = z.mul(term0002);
    DoubleVector term2011 = z.mul(term0003);
    R202 = z.fma(term2011, term0002);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  protected void order5() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    DoubleVector term0003 = work[3];
    DoubleVector term0004 = work[4];
    DoubleVector term0005 = work[5];
    R000 = term0000;
    R200 = term0001;
    R400 = term0002.mul(3);
    DoubleVector term4001 = term0003.mul(3);
    R020 = term0001;
    R040 = term0002.mul(3);
    DoubleVector term0401 = term0003.mul(3);
    R220 = term0002;
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
    DoubleVector term0012 = z.mul(term0003);
    DoubleVector term0021 = z.fma(term0012, term0002);
    R003 = z.fma(term0021, term0011.mul(2));
    DoubleVector term0013 = z.mul(term0004);
    DoubleVector term0022 = z.fma(term0013, term0003);
    DoubleVector term0031 = z.fma(term0022, term0012.mul(2));
    R004 = z.fma(term0031, term0021.mul(3));
    DoubleVector term0014 = z.mul(term0005);
    DoubleVector term0023 = z.fma(term0014, term0004);
    DoubleVector term0032 = z.fma(term0023, term0013.mul(2));
    DoubleVector term0041 = z.fma(term0032, term0022.mul(3));
    R005 = z.fma(term0041, term0031.mul(4));
    R021 = z.mul(term0002);
    DoubleVector term0211 = z.mul(term0003);
    R022 = z.fma(term0211, term0002);
    DoubleVector term0212 = z.mul(term0004);
    DoubleVector term0221 = z.fma(term0212, term0003);
    R023 = z.fma(term0221, term0211.mul(2));
    R041 = z.mul(term0401);
    R201 = z.mul(term0002);
    DoubleVector term2011 = z.mul(term0003);
    R202 = z.fma(term2011, term0002);
    DoubleVector term2012 = z.mul(term0004);
    DoubleVector term2021 = z.fma(term2012, term0003);
    R203 = z.fma(term2021, term2011.mul(2));
    R221 = z.mul(term0003);
    R401 = z.mul(term4001);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  protected void order6() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    DoubleVector term0003 = work[3];
    DoubleVector term0004 = work[4];
    DoubleVector term0005 = work[5];
    DoubleVector term0006 = work[6];
    R000 = term0000;
    R200 = term0001;
    R400 = term0002.mul(3);
    DoubleVector term4001 = term0003.mul(3);
    DoubleVector term4002 = term0004.mul(3);
    R600 = term4001.mul(5);
    R020 = term0001;
    R040 = term0002.mul(3);
    DoubleVector term0401 = term0003.mul(3);
    DoubleVector term0402 = term0004.mul(3);
    R060 = term0401.mul(5);
    R220 = term0002;
    R240 = term0003.mul(3);
    R420 = term4001;
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
    DoubleVector term0012 = z.mul(term0003);
    DoubleVector term0021 = z.fma(term0012, term0002);
    R003 = z.fma(term0021, term0011.mul(2));
    DoubleVector term0013 = z.mul(term0004);
    DoubleVector term0022 = z.fma(term0013, term0003);
    DoubleVector term0031 = z.fma(term0022, term0012.mul(2));
    R004 = z.fma(term0031, term0021.mul(3));
    DoubleVector term0014 = z.mul(term0005);
    DoubleVector term0023 = z.fma(term0014, term0004);
    DoubleVector term0032 = z.fma(term0023, term0013.mul(2));
    DoubleVector term0041 = z.fma(term0032, term0022.mul(3));
    R005 = z.fma(term0041, term0031.mul(4));
    DoubleVector term0015 = z.mul(term0006);
    DoubleVector term0024 = z.fma(term0015, term0005);
    DoubleVector term0033 = z.fma(term0024, term0014.mul(2));
    DoubleVector term0042 = z.fma(term0033, term0023.mul(3));
    DoubleVector term0051 = z.fma(term0042, term0032.mul(4));
    R006 = z.fma(term0051, term0041.mul(5));
    R021 = z.mul(term0002);
    DoubleVector term0211 = z.mul(term0003);
    R022 = z.fma(term0211, term0002);
    DoubleVector term0212 = z.mul(term0004);
    DoubleVector term0221 = z.fma(term0212, term0003);
    R023 = z.fma(term0221, term0211.mul(2));
    DoubleVector term0213 = z.mul(term0005);
    DoubleVector term0222 = z.fma(term0213, term0004);
    DoubleVector term0231 = z.fma(term0222, term0212.mul(2));
    R024 = z.fma(term0231, term0221.mul(3));
    R041 = z.mul(term0401);
    DoubleVector term0411 = z.mul(term0402);
    R042 = z.fma(term0411, term0401);
    R201 = z.mul(term0002);
    DoubleVector term2011 = z.mul(term0003);
    R202 = z.fma(term2011, term0002);
    DoubleVector term2012 = z.mul(term0004);
    DoubleVector term2021 = z.fma(term2012, term0003);
    R203 = z.fma(term2021, term2011.mul(2));
    DoubleVector term2013 = z.mul(term0005);
    DoubleVector term2022 = z.fma(term2013, term0004);
    DoubleVector term2031 = z.fma(term2022, term2012.mul(2));
    R204 = z.fma(term2031, term2021.mul(3));
    R221 = z.mul(term0003);
    DoubleVector term2211 = z.mul(term0004);
    R222 = z.fma(term2211, term0003);
    R401 = z.mul(term4001);
    DoubleVector term4011 = z.mul(term4002);
    R402 = z.fma(term4011, term4001);
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void multipoleKPotentialAtI(PolarizableMultipoleSIMD mK, int order) {
    switch (order) {
      default:
      case 3:
        DoubleVector term300 = mK.dx.mul(R400);
        term300 = mK.qxz.fma(R401, term300);
        E300 = term300.neg();
        DoubleVector term030 = mK.dy.mul(R040);
        term030 = mK.qyz.fma(R041, term030);
        E030 = term030.neg();
        DoubleVector term003 = mK.q.mul(R003);
        term003 = mK.dz.fma(R004, term003);
        term003 = mK.qxx.fma(R203, term003);
        term003 = mK.qyy.fma(R023, term003);
        term003 = mK.qzz.fma(R005, term003);
        E003 = term003.neg();
        DoubleVector term210 = mK.dy.mul(R220);
        term210 = mK.qyz.fma(R221, term210);
        E210 = term210.neg();
        DoubleVector term201 = mK.q.mul(R201);
        term201 = mK.dz.fma(R202, term201);
        term201 = mK.qxx.fma(R401, term201);
        term201 = mK.qyy.fma(R221, term201);
        term201 = mK.qzz.fma(R203, term201);
        E201 = term201.neg();
        DoubleVector term120 = mK.dx.mul(R220);
        term120 = mK.qxz.fma(R221, term120);
        E120 = term120.neg();
        DoubleVector term021 = mK.q.mul(R021);
        term021 = mK.dz.fma(R022, term021);
        term021 = mK.qxx.fma(R221, term021);
        term021 = mK.qyy.fma(R041, term021);
        term021 = mK.qzz.fma(R023, term021);
        E021 = term021.neg();
        DoubleVector term102 = mK.dx.mul(R202);
        term102 = mK.qxz.fma(R203, term102);
        E102 = term102.neg();
        DoubleVector term012 = mK.dy.mul(R022);
        term012 = mK.qyz.fma(R023, term012);
        E012 = term012.neg();
        DoubleVector term111 = mK.qxy.mul(R221);
        E111 = term111.neg();
      case 2:
        DoubleVector term200 = mK.q.mul(R200);
        term200 = mK.dz.fma(R201, term200);
        term200 = mK.qxx.fma(R400, term200);
        term200 = mK.qyy.fma(R220, term200);
        term200 = mK.qzz.fma(R202, term200);
        E200 = term200;
        DoubleVector term020 = mK.q.mul(R020);
        term020 = mK.dz.fma(R021, term020);
        term020 = mK.qxx.fma(R220, term020);
        term020 = mK.qyy.fma(R040, term020);
        term020 = mK.qzz.fma(R022, term020);
        E020 = term020;
        DoubleVector term002 = mK.q.mul(R002);
        term002 = mK.dz.fma(R003, term002);
        term002 = mK.qxx.fma(R202, term002);
        term002 = mK.qyy.fma(R022, term002);
        term002 = mK.qzz.fma(R004, term002);
        E002 = term002;
        E110 = mK.qxy.mul(R220);
        DoubleVector term101 = mK.dx.mul(R201);
        term101 = mK.qxz.fma(R202, term101);
        E101 = term101;
        DoubleVector term011 = mK.dy.mul(R021);
        term011 = mK.qyz.fma(R022, term011);
        E011 = term011;
      case 1:
        DoubleVector term100 = mK.dx.mul(R200);
        term100 = mK.qxz.fma(R201, term100);
        E100 = term100.neg();
        DoubleVector term010 = mK.dy.mul(R020);
        term010 = mK.qyz.fma(R021, term010);
        E010 = term010.neg();
        DoubleVector term001 = mK.q.mul(R001);
        term001 = mK.dz.fma(R002, term001);
        term001 = mK.qxx.fma(R201, term001);
        term001 = mK.qyy.fma(R021, term001);
        term001 = mK.qzz.fma(R003, term001);
        E001 = term001.neg();
      case 0:
        DoubleVector term000 = mK.q.mul(R000);
        term000 = mK.dz.fma(R001, term000);
        term000 = mK.qxx.fma(R200, term000);
        term000 = mK.qyy.fma(R020, term000);
        term000 = mK.qzz.fma(R002, term000);
        E000 = term000;
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void chargeKPotentialAtI(PolarizableMultipoleSIMD mK, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        E300 = zero;
        E030 = zero;
        E003 = mK.q.mul(R003).neg();
        E210 = zero;
        E201 = mK.q.mul(R201).neg();
        E120 = zero;
        E021 = mK.q.mul(R021).neg();
        E102 = zero;
        E012 = zero;
        E111 = zero;
        // Fall through to 2nd order.
      case 2:
        // Order 2
        E200 = mK.q.mul(R200);
        E020 = mK.q.mul(R020);
        E002 = mK.q.mul(R002);
        E110 = zero;
        E101 = zero;
        E011 = zero;
        // Fall through to 1st order.
      case 1:
        // Order 1
        E100 = zero;
        E010 = zero;
        E001 = mK.q.mul(R001).neg();
      case 0:
        E000 = mK.q.mul(R000);
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void dipoleKPotentialAtI(DoubleVector uxk, DoubleVector uyk, DoubleVector uzk, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        E300 = uxk.mul(R400).neg();
        E030 = uyk.mul(R040).neg();
        E003 = uzk.mul(R004).neg();
        E210 = uyk.mul(R220).neg();
        E201 = uzk.mul(R202).neg();
        E120 = uxk.mul(R220).neg();
        E021 = uzk.mul(R022).neg();
        E102 = uxk.mul(R202).neg();
        E012 = uyk.mul(R022).neg();
        E111 = zero;
        // Fall through to 2nd order.
      case 2:
        // Order 2
        E200 = uzk.mul(R201);
        E020 = uzk.mul(R021);
        E002 = uzk.mul(R003);
        E110 = zero;
        E101 = uxk.mul(R201);
        E011 = uyk.mul(R021);
        // Fall through to 1st order.
      case 1:
        // Order 1
        E100 = uxk.mul(R200).neg();
        E010 = uyk.mul(R020).neg();
        E001 = uzk.mul(R002).neg();
        // Fall through to the potential.
      case 0:
        E000 = uzk.mul(R001);
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void quadrupoleKPotentialAtI(PolarizableMultipoleSIMD mK, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        E300 = mK.qxz.mul(R401).neg();
        E030 = mK.qyz.mul(R041).neg();
        DoubleVector term003 = mK.qxx.mul(R203);
        term003 = mK.qyy.fma(R023, term003);
        term003 = mK.qzz.fma(R005, term003);
        E003 = term003;
        E210 = mK.qyz.mul(R221).neg();
        DoubleVector term201 = mK.qxx.mul(R401);
        term201 = mK.qyy.fma(R221, term201);
        term201 = mK.qzz.fma(R203, term201);
        E201 = term201.neg();
        E120 = mK.qxz.mul(R221).neg();
        DoubleVector term021 = mK.qxx.mul(R221);
        term021 = mK.qyy.fma(R041, term021);
        term021 = mK.qzz.fma(R023, term021);
        E021 = term021.neg();
        E102 = mK.qxz.mul(R203).neg();
        E012 = mK.qyz.mul(R023).neg();
        E111 = mK.qxy.mul(R221).neg();
        // Fall through to 2nd order.
      case 2:
        // Order 2
        DoubleVector term200 = mK.qxx.mul(R400);
        term200 = mK.qyy.fma(R220, term200);
        term200 = mK.qzz.fma(R202, term200);
        E200 = term200;
        DoubleVector term020 = mK.qxx.mul(R220);
        term020 = mK.qyy.fma(R040, term020);
        term020 = mK.qzz.fma(R022, term020);
        E020 = term020;
        DoubleVector term002 = mK.qxx.mul(R202);
        term002 = mK.qyy.fma(R022, term002);
        term002 = mK.qzz.fma(R004, term002);
        E002 = term002;
        E110 = mK.qxy.mul(R220);
        E101 = mK.qxz.mul(R202);
        E011 = mK.qyz.mul(R022);
        // Fall through to 1st order.
      case 1:
        // Order 1
        E100 = mK.qxz.mul(R201).neg();
        E010 = mK.qyz.mul(R021).neg();
        DoubleVector term001 = mK.qxx.mul(R201);
        term001 = mK.qyy.fma(R021, term001);
        term001 = mK.qzz.fma(R003, term001);
        E001 = term001.neg();
      case 0:
        DoubleVector term000 = mK.qxx.mul(R200);
        term000 = mK.qyy.fma(R020, term000);
        term000 = mK.qzz.fma(R002, term000);
        E000 = term000;
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void multipoleIPotentialAtK(PolarizableMultipoleSIMD mI, int order) {
    switch (order) {
      default:
      case 3:
        DoubleVector term300 = mI.dx.mul(R400.neg());
        term300 = mI.qxz.fma(R401, term300);
        E300 = term300;
        DoubleVector term030 = mI.dy.mul(R040.neg());
        term030 = mI.qyz.fma(R041, term030);
        E030 = term030;
        DoubleVector term003 = mI.q.mul(R003);
        term003 = mI.dz.fma(R004.neg(), term003);
        term003 = mI.qxx.fma(R203, term003);
        term003 = mI.qyy.fma(R023, term003);
        term003 = mI.qzz.fma(R005, term003);
        E003 = term003;
        DoubleVector term210 = mI.dy.mul(R220.neg());
        term210 = mI.qyz.fma(R221, term210);
        E210 = term210;
        DoubleVector term201 = mI.q.mul(R201);
        term201 = mI.dz.fma(R202.neg(), term201);
        term201 = mI.qxx.fma(R401, term201);
        term201 = mI.qyy.fma(R221, term201);
        term201 = mI.qzz.fma(R203, term201);
        E201 = term201;
        DoubleVector term120 = mI.dx.mul(R220.neg());
        term120 = mI.qxz.fma(R221, term120);
        E120 = term120;
        DoubleVector term021 = mI.q.mul(R021);
        term021 = mI.dz.fma(R022.neg(), term021);
        term021 = mI.qxx.fma(R221, term021);
        term021 = mI.qyy.fma(R041, term021);
        term021 = mI.qzz.fma(R023, term021);
        E021 = term021;
        DoubleVector term102 = mI.dx.mul(R202.neg());
        term102 = mI.qxz.fma(R203, term102);
        E102 = term102;
        DoubleVector term012 = mI.dy.mul(R022.neg());
        term012 = mI.qyz.fma(R023, term012);
        E012 = term012;
        E111 = mI.qxy.mul(R221);
      case 2:
        DoubleVector term200 = mI.q.mul(R200);
        term200 = mI.dz.fma(R201.neg(), term200);
        term200 = mI.qxx.fma(R400, term200);
        term200 = mI.qyy.fma(R220, term200);
        term200 = mI.qzz.fma(R202, term200);
        E200 = term200;
        DoubleVector term020 = mI.q.mul(R020);
        term020 = mI.dz.fma(R021.neg(), term020);
        term020 = mI.qxx.fma(R220, term020);
        term020 = mI.qyy.fma(R040, term020);
        term020 = mI.qzz.fma(R022, term020);
        E020 = term020;
        DoubleVector term002 = mI.q.mul(R002);
        term002 = mI.dz.fma(R003.neg(), term002);
        term002 = mI.qxx.fma(R202, term002);
        term002 = mI.qyy.fma(R022, term002);
        term002 = mI.qzz.fma(R004, term002);
        E002 = term002;
        E110 = mI.qxy.mul(R220);
        DoubleVector term101 = mI.dx.mul(R201.neg());
        term101 = mI.qxz.fma(R202, term101);
        E101 = term101;
        DoubleVector term011 = mI.dy.mul(R021.neg());
        term011 = mI.qyz.fma(R022, term011);
        E011 = term011;
      case 1:
        DoubleVector term100 = mI.dx.mul(R200.neg());
        term100 = mI.qxz.fma(R201, term100);
        E100 = term100;
        DoubleVector term010 = mI.dy.mul(R020.neg());
        term010 = mI.qyz.fma(R021, term010);
        E010 = term010;
        DoubleVector term001 = mI.q.mul(R001);
        term001 = mI.dz.fma(R002.neg(), term001);
        term001 = mI.qxx.fma(R201, term001);
        term001 = mI.qyy.fma(R021, term001);
        term001 = mI.qzz.fma(R003, term001);
        E001 = term001;
      case 0:
        DoubleVector term000 = mI.q.mul(R000);
        term000 = mI.dz.fma(R001.neg(), term000);
        term000 = mI.qxx.fma(R200, term000);
        term000 = mI.qyy.fma(R020, term000);
        term000 = mI.qzz.fma(R002, term000);
        E000 = term000;
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void chargeIPotentialAtK(PolarizableMultipoleSIMD mI, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3.
        E300 = zero;
        E030 = zero;
        E003 = mI.q.mul(R003);
        E210 = zero;
        E201 = mI.q.mul(R201);
        E120 = zero;
        E021 = mI.q.mul(R021);
        E102 = zero;
        E012 = zero;
        E111 = zero;
        // Fall through to 2nd order.
      case 2:
        // Order 2.
        E200 = mI.q.mul(R200);
        E020 = mI.q.mul(R020);
        E002 = mI.q.mul(R002);
        E110 = zero;
        E101 = zero;
        E011 = zero;
        // Fall through to 1st order.
      case 1:
        // Order 1.
        E100 = zero;
        E010 = zero;
        E001 = mI.q.mul(R001);
        // Fall through to the potential.
      case 0:
        E000 = mI.q.mul(R000);
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void dipoleIPotentialAtK(DoubleVector uxi, DoubleVector uyi, DoubleVector uzi, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3
        E300 = uxi.mul(R400).neg();
        E030 = uyi.mul(R040).neg();
        E003 = uzi.mul(R004).neg();
        E210 = uyi.mul(R220).neg();
        E201 = uzi.mul(R202).neg();
        E120 = uxi.mul(R220).neg();
        E021 = uzi.mul(R022).neg();
        E102 = uxi.mul(R202).neg();
        E012 = uyi.mul(R022).neg();
        E111 = zero;
        // Fall through to 2nd order.
      case 2:
        E200 = uzi.mul(R201.neg());
        E020 = uzi.mul(R021.neg());
        E002 = uzi.mul(R003.neg());
        E110 = zero;
        E101 = uxi.mul(R201.neg());
        E011 = uyi.mul(R021.neg());
        // Fall through to 1st order.
      case 1:
        E100 = uxi.mul(R200.neg());
        E010 = uyi.mul(R020.neg());
        E001 = uzi.mul(R002.neg());
        // Fall through to the potential.
      case 0:
        E000 = uzi.mul(R001.neg());
    }
  }

  /**
   * {@inheritDoc}
   */
  @SuppressWarnings("fallthrough")
  @Override
  protected void quadrupoleIPotentialAtK(PolarizableMultipoleSIMD mI, int order) {
    switch (order) {
      default:
      case 3:
        // Order 3.
        E300 = mI.qxz.mul(R401);
        E030 = mI.qyz.mul(R041);
        DoubleVector term003 = mI.qxx.mul(R203);
        term003 = mI.qyy.fma(R023, term003);
        term003 = mI.qzz.fma(R005, term003);
        E003 = term003;
        E210 = mI.qyz.mul(R221);
        DoubleVector term201 = mI.qxx.mul(R401);
        term201 = mI.qyy.fma(R221, term201);
        term201 = mI.qzz.fma(R203, term201);
        E201 = term201;
        E120 = mI.qxz.mul(R221);
        DoubleVector term021 = mI.qxx.mul(R221);
        term021 = mI.qyy.fma(R041, term021);
        term021 = mI.qzz.fma(R023, term021);
        E021 = term021;
        E102 = mI.qxz.mul(R203);
        E012 = mI.qyz.mul(R023);
        E111 = mI.qxy.mul(R221);
        // Fall through to 2nd order.
      case 2:
        // Order 2.
        DoubleVector term200 = mI.qxx.mul(R400);
        term200 = mI.qyy.fma(R220, term200);
        term200 = mI.qzz.fma(R202, term200);
        E200 = term200;
        DoubleVector term020 = mI.qxx.mul(R220);
        term020 = mI.qyy.fma(R040, term020);
        term020 = mI.qzz.fma(R022, term020);
        E020 = term020;
        DoubleVector term002 = mI.qxx.mul(R202);
        term002 = mI.qyy.fma(R022, term002);
        term002 = mI.qzz.fma(R004, term002);
        E002 = term002;
        E110 = mI.qxy.mul(R220);
        E101 = mI.qxz.mul(R202);
        E011 = mI.qyz.mul(R022);
        // Fall through to 1st order.
      case 1:
        // Order 1.
        E100 = mI.qxz.mul(R201);
        E010 = mI.qyz.mul(R021);
        DoubleVector term001 = mI.qxx.mul(R201);
        term001 = mI.qyy.fma(R021, term001);
        term001 = mI.qzz.fma(R003, term001);
        E001 = term001;
        // Fall through to the potential.
      case 0:
        DoubleVector term000 = mI.qxx.mul(R200);
        term000 = mI.qyy.fma(R020, term000);
        term000 = mI.qzz.fma(R002, term000);
        E000 = term000;
    }
  }
}
