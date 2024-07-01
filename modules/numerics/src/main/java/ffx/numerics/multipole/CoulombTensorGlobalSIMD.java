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

import static ffx.numerics.multipole.CoordinateSystem.GLOBAL;

/**
 * The CoulombTensorGlobal class computes derivatives of 1/|<b>r</b>| via recursion to arbitrary
 * order for Cartesian multipoles in the global frame using SIMD instructions.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://doi.org/10.1142/9789812830364_0002" target="_blank"> Matt Challacombe, Eric
 * Schwegler and Jan Almlof, Modern developments in Hartree-Fock theory: Fast methods for
 * computing the Coulomb matrix. Computational Chemistry: Review of Current Trends. pp. 53-107,
 * Ed. J. Leczszynski, World Scientifc, 1996. </a>
 * @since 1.0
 */
public class CoulombTensorGlobalSIMD extends MultipoleTensorSIMD {

  /**
   * Constructor for CoulombTensorGlobalSIMD.
   *
   * @param order The order of the tensor.
   */
  public CoulombTensorGlobalSIMD(int order) {
    super(GLOBAL, order);
    operator = Operator.COULOMB;
  }

  @Override
  public void setR(DoubleVector dx, DoubleVector dy, DoubleVector dz) {
    x = dx;
    y = dy;
    z = dz;
    DoubleVector x2 = x.mul(x);
    DoubleVector y2 = y.mul(y);
    DoubleVector z2 = z.mul(z);
    r2 = x2.add(y2).add(z2);
    R = r2.sqrt();
  }

  @Override
  protected void source(DoubleVector[] T000) {
    DoubleVector ONE = DoubleVector.broadcast(R.species(), 1.0);
    DoubleVector ir = ONE.div(R);
    DoubleVector ir2 = ir.mul(ir);
    for (int n = 0; n < o1; n++) {
      T000[n] = ir.mul(coulombSource[n]);
      ir = ir.mul(ir2);
    }
  }

  @Override
  protected void order1() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    R000 = term0000;
    R100 = x.mul(term0001);
    R010 = y.mul(term0001);
    R001 = z.mul(term0001);
  }

  @Override
  protected void order2() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    R000 = term0000;
    R100 = x.mul(term0001);
    DoubleVector term1001 = x.mul(term0002);
    R200 = x.fma(term1001, term0001);
    R010 = y.mul(term0001);
    DoubleVector term0101 = y.mul(term0002);
    R020 = y.fma(term0101, term0001);
    R110 = y.mul(term1001);
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
    R011 = z.mul(term0101);
    R101 = z.mul(term1001);
  }

  @Override
  protected void order3() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    DoubleVector term0003 = work[3];
    R000 = term0000;
    R100 = x.mul(term0001);
    DoubleVector term1001 = x.mul(term0002);
    R200 = x.fma(term1001, term0001);
    DoubleVector term1002 = x.mul(term0003);
    DoubleVector term2001 = x.fma(term1002, term0002);
    R300 = x.fma(term2001, term1001.mul(2));
    R010 = y.mul(term0001);
    DoubleVector term0101 = y.mul(term0002);
    R020 = y.fma(term0101, term0001);
    DoubleVector term0102 = y.mul(term0003);
    DoubleVector term0201 = y.fma(term0102, term0002);
    R030 = y.fma(term0201, term0101.mul(2));
    R110 = y.mul(term1001);
    DoubleVector term1101 = y.mul(term1002);
    R120 = y.fma(term1101, term1001);
    R210 = y.mul(term2001);
    R001 = z.mul(term0001);
    DoubleVector term0011 = z.mul(term0002);
    R002 = z.fma(term0011, term0001);
    DoubleVector term0012 = z.mul(term0003);
    DoubleVector term0021 = z.fma(term0012, term0002);
    R003 = z.fma(term0021, term0011.mul(2));
    R011 = z.mul(term0101);
    DoubleVector term0111 = z.mul(term0102);
    R012 = z.fma(term0111, term0101);
    R021 = z.mul(term0201);
    R101 = z.mul(term1001);
    DoubleVector term1011 = z.mul(term1002);
    R102 = z.fma(term1011, term1001);
    R111 = z.mul(term1101);
    R201 = z.mul(term2001);
  }

  @Override
  protected void order4() {
    source(work);
    DoubleVector term0000 = work[0];
    DoubleVector term0001 = work[1];
    DoubleVector term0002 = work[2];
    DoubleVector term0003 = work[3];
    DoubleVector term0004 = work[4];
    R000 = term0000;
    R100 = x.mul(term0001);
    DoubleVector term1001 = x.mul(term0002);
    R200 = x.fma(term1001, term0001);
    DoubleVector term1002 = x.mul(term0003);
    DoubleVector term2001 = x.fma(term1002, term0002);
    R300 = x.fma(term2001, term1001.mul(2));
    DoubleVector term1003 = x.mul(term0004);
    DoubleVector term2002 = x.fma(term1003, term0003);
    DoubleVector term3001 = x.fma(term2002, term1002.mul(2));
    R400 = x.fma(term3001, term2001.mul(3));
    R010 = y.mul(term0001);
    DoubleVector term0101 = y.mul(term0002);
    R020 = y.fma(term0101, term0001);
    DoubleVector term0102 = y.mul(term0003);
    DoubleVector term0201 = y.fma(term0102, term0002);
    R030 = y.fma(term0201, term0101.mul(2));
    DoubleVector term0103 = y.mul(term0004);
    DoubleVector term0202 = y.fma(term0103, term0003);
    DoubleVector term0301 = y.fma(term0202, term0102.mul(2));
    R040 = y.fma(term0301, term0201.mul(3));
    R110 = y.mul(term1001);
    DoubleVector term1101 = y.mul(term1002);
    R120 = y.fma(term1101, term1001);
    DoubleVector term1102 = y.mul(term1003);
    DoubleVector term1201 = y.fma(term1102, term1002);
    R130 = y.fma(term1201, term1101.mul(2));
    R210 = y.mul(term2001);
    DoubleVector term2101 = y.mul(term2002);
    R220 = y.fma(term2101, term2001);
    R310 = y.mul(term3001);
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
    R011 = z.mul(term0101);
    DoubleVector term0111 = z.mul(term0102);
    R012 = z.fma(term0111, term0101);
    DoubleVector term0112 = z.mul(term0103);
    DoubleVector term0121 = z.fma(term0112, term0102);
    R013 = z.fma(term0121, term0111.mul(2));
    R021 = z.mul(term0201);
    DoubleVector term0211 = z.mul(term0202);
    R022 = z.fma(term0211, term0201);
    R031 = z.mul(term0301);
    R101 = z.mul(term1001);
    DoubleVector term1011 = z.mul(term1002);
    R102 = z.fma(term1011, term1001);
    DoubleVector term1012 = z.mul(term1003);
    DoubleVector term1021 = z.fma(term1012, term1002);
    R103 = z.fma(term1021, term1011.mul(2));
    R111 = z.mul(term1101);
    DoubleVector term1111 = z.mul(term1102);
    R112 = z.fma(term1111, term1101);
    R121 = z.mul(term1201);
    R201 = z.mul(term2001);
    DoubleVector term2011 = z.mul(term2002);
    R202 = z.fma(term2011, term2001);
    R211 = z.mul(term2101);
    R301 = z.mul(term3001);
  }

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
    R100 = x.mul(term0001);
    DoubleVector term1001 = x.mul(term0002);
    R200 = x.fma(term1001, term0001);
    DoubleVector term1002 = x.mul(term0003);
    DoubleVector term2001 = x.fma(term1002, term0002);
    R300 = x.fma(term2001, term1001.mul(2));
    DoubleVector term1003 = x.mul(term0004);
    DoubleVector term2002 = x.fma(term1003, term0003);
    DoubleVector term3001 = x.fma(term2002, term1002.mul(2));
    R400 = x.fma(term3001, term2001.mul(3));
    DoubleVector term1004 = x.mul(term0005);
    DoubleVector term2003 = x.fma(term1004, term0004);
    DoubleVector term3002 = x.fma(term2003, term1003.mul(2));
    DoubleVector term4001 = x.fma(term3002, term2002.mul(3));
    R500 = x.fma(term4001, term3001.mul(4));
    R010 = y.mul(term0001);
    DoubleVector term0101 = y.mul(term0002);
    R020 = y.fma(term0101, term0001);
    DoubleVector term0102 = y.mul(term0003);
    DoubleVector term0201 = y.fma(term0102, term0002);
    R030 = y.fma(term0201, term0101.mul(2));
    DoubleVector term0103 = y.mul(term0004);
    DoubleVector term0202 = y.fma(term0103, term0003);
    DoubleVector term0301 = y.fma(term0202, term0102.mul(2));
    R040 = y.fma(term0301, term0201.mul(3));
    DoubleVector term0104 = y.mul(term0005);
    DoubleVector term0203 = y.fma(term0104, term0004);
    DoubleVector term0302 = y.fma(term0203, term0103.mul(2));
    DoubleVector term0401 = y.fma(term0302, term0202.mul(3));
    R050 = y.fma(term0401, term0301.mul(4));
    R110 = y.mul(term1001);
    DoubleVector term1101 = y.mul(term1002);
    R120 = y.fma(term1101, term1001);
    DoubleVector term1102 = y.mul(term1003);
    DoubleVector term1201 = y.fma(term1102, term1002);
    R130 = y.fma(term1201, term1101.mul(2));
    DoubleVector term1103 = y.mul(term1004);
    DoubleVector term1202 = y.fma(term1103, term1003);
    DoubleVector term1301 = y.fma(term1202, term1102.mul(2));
    R140 = y.fma(term1301, term1201.mul(3));
    R210 = y.mul(term2001);
    DoubleVector term2101 = y.mul(term2002);
    R220 = y.fma(term2101, term2001);
    DoubleVector term2102 = y.mul(term2003);
    DoubleVector term2201 = y.fma(term2102, term2002);
    R230 = y.fma(term2201, term2101.mul(2));
    R310 = y.mul(term3001);
    DoubleVector term3101 = y.mul(term3002);
    R320 = y.fma(term3101, term3001);
    R410 = y.mul(term4001);
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
    R011 = z.mul(term0101);
    DoubleVector term0111 = z.mul(term0102);
    R012 = z.fma(term0111, term0101);
    DoubleVector term0112 = z.mul(term0103);
    DoubleVector term0121 = z.fma(term0112, term0102);
    R013 = z.fma(term0121, term0111.mul(2));
    DoubleVector term0113 = z.mul(term0104);
    DoubleVector term0122 = z.fma(term0113, term0103);
    DoubleVector term0131 = z.fma(term0122, term0112.mul(2));
    R014 = z.fma(term0131, term0121.mul(3));
    R021 = z.mul(term0201);
    DoubleVector term0211 = z.mul(term0202);
    R022 = z.fma(term0211, term0201);
    DoubleVector term0212 = z.mul(term0203);
    DoubleVector term0221 = z.fma(term0212, term0202);
    R023 = z.fma(term0221, term0211.mul(2));
    R031 = z.mul(term0301);
    DoubleVector term0311 = z.mul(term0302);
    R032 = z.fma(term0311, term0301);
    R041 = z.mul(term0401);
    R101 = z.mul(term1001);
    DoubleVector term1011 = z.mul(term1002);
    R102 = z.fma(term1011, term1001);
    DoubleVector term1012 = z.mul(term1003);
    DoubleVector term1021 = z.fma(term1012, term1002);
    R103 = z.fma(term1021, term1011.mul(2));
    DoubleVector term1013 = z.mul(term1004);
    DoubleVector term1022 = z.fma(term1013, term1003);
    DoubleVector term1031 = z.fma(term1022, term1012.mul(2));
    R104 = z.fma(term1031, term1021.mul(3));
    R111 = z.mul(term1101);
    DoubleVector term1111 = z.mul(term1102);
    R112 = z.fma(term1111, term1101);
    DoubleVector term1112 = z.mul(term1103);
    DoubleVector term1121 = z.fma(term1112, term1102);
    R113 = z.fma(term1121, term1111.mul(2));
    R121 = z.mul(term1201);
    DoubleVector term1211 = z.mul(term1202);
    R122 = z.fma(term1211, term1201);
    R131 = z.mul(term1301);
    R201 = z.mul(term2001);
    DoubleVector term2011 = z.mul(term2002);
    R202 = z.fma(term2011, term2001);
    DoubleVector term2012 = z.mul(term2003);
    DoubleVector term2021 = z.fma(term2012, term2002);
    R203 = z.fma(term2021, term2011.mul(2));
    R211 = z.mul(term2101);
    DoubleVector term2111 = z.mul(term2102);
    R212 = z.fma(term2111, term2101);
    R221 = z.mul(term2201);
    R301 = z.mul(term3001);
    DoubleVector term3011 = z.mul(term3002);
    R302 = z.fma(term3011, term3001);
    R311 = z.mul(term3101);
    R401 = z.mul(term4001);
  }

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
    R100 = x.mul(term0001);
    DoubleVector term1001 = x.mul(term0002);
    R200 = x.fma(term1001, term0001);
    DoubleVector term1002 = x.mul(term0003);
    DoubleVector term2001 = x.fma(term1002, term0002);
    R300 = x.fma(term2001, term1001.mul(2));
    DoubleVector term1003 = x.mul(term0004);
    DoubleVector term2002 = x.fma(term1003, term0003);
    DoubleVector term3001 = x.fma(term2002, term1002.mul(2));
    R400 = x.fma(term3001, term2001.mul(3));
    DoubleVector term1004 = x.mul(term0005);
    DoubleVector term2003 = x.fma(term1004, term0004);
    DoubleVector term3002 = x.fma(term2003, term1003.mul(2));
    DoubleVector term4001 = x.fma(term3002, term2002.mul(3));
    R500 = x.fma(term4001, term3001.mul(4));
    DoubleVector term1005 = x.mul(term0006);
    DoubleVector term2004 = x.fma(term1005, term0005);
    DoubleVector term3003 = x.fma(term2004, term1004.mul(2));
    DoubleVector term4002 = x.fma(term3003, term2003.mul(3));
    DoubleVector term5001 = x.fma(term4002, term3002.mul(4));
    R600 = x.fma(term5001, term4001.mul(5));
    R010 = y.mul(term0001);
    DoubleVector term0101 = y.mul(term0002);
    R020 = y.fma(term0101, term0001);
    DoubleVector term0102 = y.mul(term0003);
    DoubleVector term0201 = y.fma(term0102, term0002);
    R030 = y.fma(term0201, term0101.mul(2));
    DoubleVector term0103 = y.mul(term0004);
    DoubleVector term0202 = y.fma(term0103, term0003);
    DoubleVector term0301 = y.fma(term0202, term0102.mul(2));
    R040 = y.fma(term0301, term0201.mul(3));
    DoubleVector term0104 = y.mul(term0005);
    DoubleVector term0203 = y.fma(term0104, term0004);
    DoubleVector term0302 = y.fma(term0203, term0103.mul(2));
    DoubleVector term0401 = y.fma(term0302, term0202.mul(3));
    R050 = y.fma(term0401, term0301.mul(4));
    DoubleVector term0105 = y.mul(term0006);
    DoubleVector term0204 = y.fma(term0105, term0005);
    DoubleVector term0303 = y.fma(term0204, term0104.mul(2));
    DoubleVector term0402 = y.fma(term0303, term0203.mul(3));
    DoubleVector term0501 = y.fma(term0402, term0302.mul(4));
    R060 = y.fma(term0501, term0401.mul(5));
    R110 = y.mul(term1001);
    DoubleVector term1101 = y.mul(term1002);
    R120 = y.fma(term1101, term1001);
    DoubleVector term1102 = y.mul(term1003);
    DoubleVector term1201 = y.fma(term1102, term1002);
    R130 = y.fma(term1201, term1101.mul(2));
    DoubleVector term1103 = y.mul(term1004);
    DoubleVector term1202 = y.fma(term1103, term1003);
    DoubleVector term1301 = y.fma(term1202, term1102.mul(2));
    R140 = y.fma(term1301, term1201.mul(3));
    DoubleVector term1104 = y.mul(term1005);
    DoubleVector term1203 = y.fma(term1104, term1004);
    DoubleVector term1302 = y.fma(term1203, term1103.mul(2));
    DoubleVector term1401 = y.fma(term1302, term1202.mul(3));
    R150 = y.fma(term1401, term1301.mul(4));
    R210 = y.mul(term2001);
    DoubleVector term2101 = y.mul(term2002);
    R220 = y.fma(term2101, term2001);
    DoubleVector term2102 = y.mul(term2003);
    DoubleVector term2201 = y.fma(term2102, term2002);
    R230 = y.fma(term2201, term2101.mul(2));
    DoubleVector term2103 = y.mul(term2004);
    DoubleVector term2202 = y.fma(term2103, term2003);
    DoubleVector term2301 = y.fma(term2202, term2102.mul(2));
    R240 = y.fma(term2301, term2201.mul(3));
    R310 = y.mul(term3001);
    DoubleVector term3101 = y.mul(term3002);
    R320 = y.fma(term3101, term3001);
    DoubleVector term3102 = y.mul(term3003);
    DoubleVector term3201 = y.fma(term3102, term3002);
    R330 = y.fma(term3201, term3101.mul(2));
    R410 = y.mul(term4001);
    DoubleVector term4101 = y.mul(term4002);
    R420 = y.fma(term4101, term4001);
    R510 = y.mul(term5001);
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
    R011 = z.mul(term0101);
    DoubleVector term0111 = z.mul(term0102);
    R012 = z.fma(term0111, term0101);
    DoubleVector term0112 = z.mul(term0103);
    DoubleVector term0121 = z.fma(term0112, term0102);
    R013 = z.fma(term0121, term0111.mul(2));
    DoubleVector term0113 = z.mul(term0104);
    DoubleVector term0122 = z.fma(term0113, term0103);
    DoubleVector term0131 = z.fma(term0122, term0112.mul(2));
    R014 = z.fma(term0131, term0121.mul(3));
    DoubleVector term0114 = z.mul(term0105);
    DoubleVector term0123 = z.fma(term0114, term0104);
    DoubleVector term0132 = z.fma(term0123, term0113.mul(2));
    DoubleVector term0141 = z.fma(term0132, term0122.mul(3));
    R015 = z.fma(term0141, term0131.mul(4));
    R021 = z.mul(term0201);
    DoubleVector term0211 = z.mul(term0202);
    R022 = z.fma(term0211, term0201);
    DoubleVector term0212 = z.mul(term0203);
    DoubleVector term0221 = z.fma(term0212, term0202);
    R023 = z.fma(term0221, term0211.mul(2));
    DoubleVector term0213 = z.mul(term0204);
    DoubleVector term0222 = z.fma(term0213, term0203);
    DoubleVector term0231 = z.fma(term0222, term0212.mul(2));
    R024 = z.fma(term0231, term0221.mul(3));
    R031 = z.mul(term0301);
    DoubleVector term0311 = z.mul(term0302);
    R032 = z.fma(term0311, term0301);
    DoubleVector term0312 = z.mul(term0303);
    DoubleVector term0321 = z.fma(term0312, term0302);
    R033 = z.fma(term0321, term0311.mul(2));
    R041 = z.mul(term0401);
    DoubleVector term0411 = z.mul(term0402);
    R042 = z.fma(term0411, term0401);
    R051 = z.mul(term0501);
    R101 = z.mul(term1001);
    DoubleVector term1011 = z.mul(term1002);
    R102 = z.fma(term1011, term1001);
    DoubleVector term1012 = z.mul(term1003);
    DoubleVector term1021 = z.fma(term1012, term1002);
    R103 = z.fma(term1021, term1011.mul(2));
    DoubleVector term1013 = z.mul(term1004);
    DoubleVector term1022 = z.fma(term1013, term1003);
    DoubleVector term1031 = z.fma(term1022, term1012.mul(2));
    R104 = z.fma(term1031, term1021.mul(3));
    DoubleVector term1014 = z.mul(term1005);
    DoubleVector term1023 = z.fma(term1014, term1004);
    DoubleVector term1032 = z.fma(term1023, term1013.mul(2));
    DoubleVector term1041 = z.fma(term1032, term1022.mul(3));
    R105 = z.fma(term1041, term1031.mul(4));
    R111 = z.mul(term1101);
    DoubleVector term1111 = z.mul(term1102);
    R112 = z.fma(term1111, term1101);
    DoubleVector term1112 = z.mul(term1103);
    DoubleVector term1121 = z.fma(term1112, term1102);
    R113 = z.fma(term1121, term1111.mul(2));
    DoubleVector term1113 = z.mul(term1104);
    DoubleVector term1122 = z.fma(term1113, term1103);
    DoubleVector term1131 = z.fma(term1122, term1112.mul(2));
    R114 = z.fma(term1131, term1121.mul(3));
    R121 = z.mul(term1201);
    DoubleVector term1211 = z.mul(term1202);
    R122 = z.fma(term1211, term1201);
    DoubleVector term1212 = z.mul(term1203);
    DoubleVector term1221 = z.fma(term1212, term1202);
    R123 = z.fma(term1221, term1211.mul(2));
    R131 = z.mul(term1301);
    DoubleVector term1311 = z.mul(term1302);
    R132 = z.fma(term1311, term1301);
    R141 = z.mul(term1401);
    R201 = z.mul(term2001);
    DoubleVector term2011 = z.mul(term2002);
    R202 = z.fma(term2011, term2001);
    DoubleVector term2012 = z.mul(term2003);
    DoubleVector term2021 = z.fma(term2012, term2002);
    R203 = z.fma(term2021, term2011.mul(2));
    DoubleVector term2013 = z.mul(term2004);
    DoubleVector term2022 = z.fma(term2013, term2003);
    DoubleVector term2031 = z.fma(term2022, term2012.mul(2));
    R204 = z.fma(term2031, term2021.mul(3));
    R211 = z.mul(term2101);
    DoubleVector term2111 = z.mul(term2102);
    R212 = z.fma(term2111, term2101);
    DoubleVector term2112 = z.mul(term2103);
    DoubleVector term2121 = z.fma(term2112, term2102);
    R213 = z.fma(term2121, term2111.mul(2));
    R221 = z.mul(term2201);
    DoubleVector term2211 = z.mul(term2202);
    R222 = z.fma(term2211, term2201);
    R231 = z.mul(term2301);
    R301 = z.mul(term3001);
    DoubleVector term3011 = z.mul(term3002);
    R302 = z.fma(term3011, term3001);
    DoubleVector term3012 = z.mul(term3003);
    DoubleVector term3021 = z.fma(term3012, term3002);
    R303 = z.fma(term3021, term3011.mul(2));
    R311 = z.mul(term3101);
    DoubleVector term3111 = z.mul(term3102);
    R312 = z.fma(term3111, term3101);
    R321 = z.mul(term3201);
    R401 = z.mul(term4001);
    DoubleVector term4011 = z.mul(term4002);
    R402 = z.fma(term4011, term4001);
    R411 = z.mul(term4101);
    R501 = z.mul(term5001);
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
        // Order 3 (need a minus sign).
        DoubleVector term300 = mK.q.mul(R300);
        term300 = mK.dx.fma(R400, term300);
        term300 = mK.dy.fma(R310, term300);
        term300 = mK.dz.fma(R301, term300);
        term300 = mK.qxx.fma(R500, term300);
        term300 = mK.qyy.fma(R320, term300);
        term300 = mK.qzz.fma(R302, term300);
        term300 = mK.qxy.fma(R410, term300);
        term300 = mK.qxz.fma(R401, term300);
        term300 = mK.qyz.fma(R311, term300);
        E300 = term300.neg();
        DoubleVector term030 = mK.q.mul(R030);
        term030 = mK.dx.fma(R130, term030);
        term030 = mK.dy.fma(R040, term030);
        term030 = mK.dz.fma(R031, term030);
        term030 = mK.qxx.fma(R230, term030);
        term030 = mK.qyy.fma(R050, term030);
        term030 = mK.qzz.fma(R032, term030);
        term030 = mK.qxy.fma(R140, term030);
        term030 = mK.qxz.fma(R131, term030);
        term030 = mK.qyz.fma(R041, term030);
        E030 = term030.neg();
        DoubleVector term003 = mK.q.mul(R003);
        term003 = mK.dx.fma(R103, term003);
        term003 = mK.dy.fma(R013, term003);
        term003 = mK.dz.fma(R004, term003);
        term003 = mK.qxx.fma(R203, term003);
        term003 = mK.qyy.fma(R023, term003);
        term003 = mK.qzz.fma(R005, term003);
        term003 = mK.qxy.fma(R113, term003);
        term003 = mK.qxz.fma(R104, term003);
        term003 = mK.qyz.fma(R014, term003);
        E003 = term003.neg();
        DoubleVector term210 = mK.q.mul(R210);
        term210 = mK.dx.fma(R310, term210);
        term210 = mK.dy.fma(R220, term210);
        term210 = mK.dz.fma(R211, term210);
        term210 = mK.qxx.fma(R410, term210);
        term210 = mK.qyy.fma(R230, term210);
        term210 = mK.qzz.fma(R212, term210);
        term210 = mK.qxy.fma(R320, term210);
        term210 = mK.qxz.fma(R311, term210);
        term210 = mK.qyz.fma(R221, term210);
        E210 = term210.neg();
        DoubleVector term201 = mK.q.mul(R201);
        term201 = mK.dx.fma(R301, term201);
        term201 = mK.dy.fma(R211, term201);
        term201 = mK.dz.fma(R202, term201);
        term201 = mK.qxx.fma(R401, term201);
        term201 = mK.qyy.fma(R221, term201);
        term201 = mK.qzz.fma(R203, term201);
        term201 = mK.qxy.fma(R311, term201);
        term201 = mK.qxz.fma(R302, term201);
        term201 = mK.qyz.fma(R212, term201);
        E201 = term201.neg();
        DoubleVector term120 = mK.q.mul(R120);
        term120 = mK.dx.fma(R220, term120);
        term120 = mK.dy.fma(R130, term120);
        term120 = mK.dz.fma(R121, term120);
        term120 = mK.qxx.fma(R320, term120);
        term120 = mK.qyy.fma(R140, term120);
        term120 = mK.qzz.fma(R122, term120);
        term120 = mK.qxy.fma(R230, term120);
        term120 = mK.qxz.fma(R221, term120);
        term120 = mK.qyz.fma(R131, term120);
        E120 = term120.neg();
        DoubleVector term021 = mK.q.mul(R021);
        term021 = mK.dx.fma(R121, term021);
        term021 = mK.dy.fma(R031, term021);
        term021 = mK.dz.fma(R022, term021);
        term021 = mK.qxx.fma(R221, term021);
        term021 = mK.qyy.fma(R041, term021);
        term021 = mK.qzz.fma(R023, term021);
        term021 = mK.qxy.fma(R131, term021);
        term021 = mK.qxz.fma(R122, term021);
        term021 = mK.qyz.fma(R032, term021);
        E021 = term021.neg();
        DoubleVector term102 = mK.q.mul(R102);
        term102 = mK.dx.fma(R202, term102);
        term102 = mK.dy.fma(R112, term102);
        term102 = mK.dz.fma(R103, term102);
        term102 = mK.qxx.fma(R302, term102);
        term102 = mK.qyy.fma(R122, term102);
        term102 = mK.qzz.fma(R104, term102);
        term102 = mK.qxy.fma(R212, term102);
        term102 = mK.qxz.fma(R203, term102);
        term102 = mK.qyz.fma(R113, term102);
        E102 = term102.neg();
        DoubleVector term012 = mK.q.mul(R012);
        term012 = mK.dx.fma(R112, term012);
        term012 = mK.dy.fma(R022, term012);
        term012 = mK.dz.fma(R013, term012);
        term012 = mK.qxx.fma(R212, term012);
        term012 = mK.qyy.fma(R032, term012);
        term012 = mK.qzz.fma(R014, term012);
        term012 = mK.qxy.fma(R122, term012);
        term012 = mK.qxz.fma(R113, term012);
        term012 = mK.qyz.fma(R023, term012);
        E012 = term012.neg();
        DoubleVector term111 = mK.q.mul(R111);
        term111 = mK.dx.fma(R211, term111);
        term111 = mK.dy.fma(R121, term111);
        term111 = mK.dz.fma(R112, term111);
        term111 = mK.qxx.fma(R311, term111);
        term111 = mK.qyy.fma(R131, term111);
        term111 = mK.qzz.fma(R113, term111);
        term111 = mK.qxy.fma(R221, term111);
        term111 = mK.qxz.fma(R212, term111);
        term111 = mK.qyz.fma(R122, term111);
        E111 = term111.neg();
      case 2:
        // Order 2.
        DoubleVector term200 = mK.q.mul(R200);
        term200 = mK.dx.fma(R300, term200);
        term200 = mK.dy.fma(R210, term200);
        term200 = mK.dz.fma(R201, term200);
        term200 = mK.qxx.fma(R400, term200);
        term200 = mK.qyy.fma(R220, term200);
        term200 = mK.qzz.fma(R202, term200);
        term200 = mK.qxy.fma(R310, term200);
        term200 = mK.qxz.fma(R301, term200);
        term200 = mK.qyz.fma(R211, term200);
        E200 = term200;
        DoubleVector term020 = mK.q.mul(R020);
        term020 = mK.dx.fma(R120, term020);
        term020 = mK.dy.fma(R030, term020);
        term020 = mK.dz.fma(R021, term020);
        term020 = mK.qxx.fma(R220, term020);
        term020 = mK.qyy.fma(R040, term020);
        term020 = mK.qzz.fma(R022, term020);
        term020 = mK.qxy.fma(R130, term020);
        term020 = mK.qxz.fma(R121, term020);
        term020 = mK.qyz.fma(R031, term020);
        E020 = term020;
        DoubleVector term002 = mK.q.mul(R002);
        term002 = mK.dx.fma(R102, term002);
        term002 = mK.dy.fma(R012, term002);
        term002 = mK.dz.fma(R003, term002);
        term002 = mK.qxx.fma(R202, term002);
        term002 = mK.qyy.fma(R022, term002);
        term002 = mK.qzz.fma(R004, term002);
        term002 = mK.qxy.fma(R112, term002);
        term002 = mK.qxz.fma(R103, term002);
        term002 = mK.qyz.fma(R013, term002);
        E002 = term002;
        DoubleVector term110 = mK.q.mul(R110);
        term110 = mK.dx.fma(R210, term110);
        term110 = mK.dy.fma(R120, term110);
        term110 = mK.dz.fma(R111, term110);
        term110 = mK.qxx.fma(R310, term110);
        term110 = mK.qyy.fma(R130, term110);
        term110 = mK.qzz.fma(R112, term110);
        term110 = mK.qxy.fma(R220, term110);
        term110 = mK.qxz.fma(R211, term110);
        term110 = mK.qyz.fma(R121, term110);
        E110 = term110;
        DoubleVector term101 = mK.q.mul(R101);
        term101 = mK.dx.fma(R201, term101);
        term101 = mK.dy.fma(R111, term101);
        term101 = mK.dz.fma(R102, term101);
        term101 = mK.qxx.fma(R301, term101);
        term101 = mK.qyy.fma(R121, term101);
        term101 = mK.qzz.fma(R103, term101);
        term101 = mK.qxy.fma(R211, term101);
        term101 = mK.qxz.fma(R202, term101);
        term101 = mK.qyz.fma(R112, term101);
        E101 = term101;
        DoubleVector term011 = mK.q.mul(R011);
        term011 = mK.dx.fma(R111, term011);
        term011 = mK.dy.fma(R021, term011);
        term011 = mK.dz.fma(R012, term011);
        term011 = mK.qxx.fma(R211, term011);
        term011 = mK.qyy.fma(R031, term011);
        term011 = mK.qzz.fma(R013, term011);
        term011 = mK.qxy.fma(R121, term011);
        term011 = mK.qxz.fma(R112, term011);
        term011 = mK.qyz.fma(R022, term011);
        E011 = term011;
      case 1:
        // Order 1 (need a minus sign).
        DoubleVector term100 = mK.q.mul(R100);
        term100 = mK.dx.fma(R200, term100);
        term100 = mK.dy.fma(R110, term100);
        term100 = mK.dz.fma(R101, term100);
        term100 = mK.qxx.fma(R300, term100);
        term100 = mK.qyy.fma(R120, term100);
        term100 = mK.qzz.fma(R102, term100);
        term100 = mK.qxy.fma(R210, term100);
        term100 = mK.qxz.fma(R201, term100);
        term100 = mK.qyz.fma(R111, term100);
        E100 = term100.neg();
        DoubleVector term010 = mK.q.mul(R010);
        term010 = mK.dx.fma(R110, term010);
        term010 = mK.dy.fma(R020, term010);
        term010 = mK.dz.fma(R011, term010);
        term010 = mK.qxx.fma(R210, term010);
        term010 = mK.qyy.fma(R030, term010);
        term010 = mK.qzz.fma(R012, term010);
        term010 = mK.qxy.fma(R120, term010);
        term010 = mK.qxz.fma(R111, term010);
        term010 = mK.qyz.fma(R021, term010);
        E010 = term010.neg();
        DoubleVector term001 = mK.q.mul(R001);
        term001 = mK.dx.fma(R101, term001);
        term001 = mK.dy.fma(R011, term001);
        term001 = mK.dz.fma(R002, term001);
        term001 = mK.qxx.fma(R201, term001);
        term001 = mK.qyy.fma(R021, term001);
        term001 = mK.qzz.fma(R003, term001);
        term001 = mK.qxy.fma(R111, term001);
        term001 = mK.qxz.fma(R102, term001);
        term001 = mK.qyz.fma(R012, term001);
        E001 = term001.neg();
      case 0:
        DoubleVector term000 = mK.q.mul(R000);
        term000 = mK.dx.fma(R100, term000);
        term000 = mK.dy.fma(R010, term000);
        term000 = mK.dz.fma(R001, term000);
        term000 = mK.qxx.fma(R200, term000);
        term000 = mK.qyy.fma(R020, term000);
        term000 = mK.qzz.fma(R002, term000);
        term000 = mK.qxy.fma(R110, term000);
        term000 = mK.qxz.fma(R101, term000);
        term000 = mK.qyz.fma(R011, term000);
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
        // Order 3 (need a minus sign).
        E300 = mK.q.mul(R300).neg();
        E030 = mK.q.mul(R030).neg();
        E003 = mK.q.mul(R003).neg();
        E210 = mK.q.mul(R210).neg();
        E201 = mK.q.mul(R201).neg();
        E120 = mK.q.mul(R120).neg();
        E021 = mK.q.mul(R021).neg();
        E102 = mK.q.mul(R102).neg();
        E012 = mK.q.mul(R012).neg();
        E111 = mK.q.mul(R111).neg();
      case 2:
        // Order 2.
        E200 = mK.q.mul(R200);
        E020 = mK.q.mul(R020);
        E002 = mK.q.mul(R002);
        E110 = mK.q.mul(R110);
        E101 = mK.q.mul(R101);
        E011 = mK.q.mul(R011);
      case 1:
        // Order 1 (need a minus sign).
        E100 = mK.q.mul(R100).neg();
        E010 = mK.q.mul(R010).neg();
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
        DoubleVector term300 = uxk.mul(R400);
        term300 = uyk.fma(R310, term300);
        term300 = uzk.fma(R301, term300);
        E300 = term300.neg();
        DoubleVector term030 = uxk.mul(R130);
        term030 = uyk.fma(R040, term030);
        term030 = uzk.fma(R031, term030);
        E030 = term030.neg();
        DoubleVector term003 = uxk.mul(R103);
        term003 = uyk.fma(R013, term003);
        term003 = uzk.fma(R004, term003);
        E003 = term003.neg();
        DoubleVector term210 = uxk.mul(R310);
        term210 = uyk.fma(R220, term210);
        term210 = uzk.fma(R211, term210);
        E210 = term210.neg();
        DoubleVector term201 = uxk.mul(R301);
        term201 = uyk.fma(R211, term201);
        term201 = uzk.fma(R202, term201);
        E201 = term201.neg();
        DoubleVector term120 = uxk.mul(R220);
        term120 = uyk.fma(R130, term120);
        term120 = uzk.fma(R121, term120);
        E120 = term120.neg();
        DoubleVector term021 = uxk.mul(R121);
        term021 = uyk.fma(R031, term021);
        term021 = uzk.fma(R022, term021);
        E021 = term021.neg();
        DoubleVector term102 = uxk.mul(R202);
        term102 = uyk.fma(R112, term102);
        term102 = uzk.fma(R103, term102);
        E102 = term102.neg();
        DoubleVector term012 = uxk.mul(R112);
        term012 = uyk.fma(R022, term012);
        term012 = uzk.fma(R013, term012);
        E012 = term012.neg();
        DoubleVector term111 = uxk.mul(R211);
        term111 = uyk.fma(R121, term111);
        term111 = uzk.fma(R112, term111);
        E111 = term111.neg();
        // Foll through to 2nd order.
      case 2:
        // Order 2.
        DoubleVector term200 = uxk.mul(R300);
        term200 = uyk.fma(R210, term200);
        term200 = uzk.fma(R201, term200);
        E200 = term200;
        DoubleVector term020 = uxk.mul(R120);
        term020 = uyk.fma(R030, term020);
        term020 = uzk.fma(R021, term020);
        E020 = term020;
        DoubleVector term002 = uxk.mul(R102);
        term002 = uyk.fma(R012, term002);
        term002 = uzk.fma(R003, term002);
        E002 = term002;
        DoubleVector term110 = uxk.mul(R210);
        term110 = uyk.fma(R120, term110);
        term110 = uzk.fma(R111, term110);
        E110 = term110;
        DoubleVector term101 = uxk.mul(R201);
        term101 = uyk.fma(R111, term101);
        term101 = uzk.fma(R102, term101);
        E101 = term101;
        DoubleVector term011 = uxk.mul(R111);
        term011 = uyk.fma(R021, term011);
        term011 = uzk.fma(R012, term011);
        E011 = term011;
      case 1:
        // Order 1 (need a minus sign).
        DoubleVector term100 = uxk.mul(R200);
        term100 = uyk.fma(R110, term100);
        term100 = uzk.fma(R101, term100);
        E100 = term100.neg();
        DoubleVector term010 = uxk.mul(R110);
        term010 = uyk.fma(R020, term010);
        term010 = uzk.fma(R011, term010);
        E010 = term010.neg();
        DoubleVector term001 = uxk.mul(R101);
        term001 = uyk.fma(R011, term001);
        term001 = uzk.fma(R002, term001);
        E001 = term001.neg();
      case 0:
        DoubleVector term000 = uxk.mul(R100);
        term000 = uyk.fma(R010, term000);
        term000 = uzk.fma(R001, term000);
        E000 = term000;
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
        // Order 3 (need a minus sign).
        DoubleVector term300 = mK.qxx.mul(R500);
        term300 = mK.qyy.fma(R320, term300);
        term300 = mK.qzz.fma(R302, term300);
        term300 = mK.qxy.fma(R410, term300);
        term300 = mK.qxz.fma(R401, term300);
        term300 = mK.qyz.fma(R311, term300);
        E300 = term300.neg();
        DoubleVector term030 = mK.qxx.mul(R230);
        term030 = mK.qyy.fma(R050, term030);
        term030 = mK.qzz.fma(R032, term030);
        term030 = mK.qxy.fma(R140, term030);
        term030 = mK.qxz.fma(R131, term030);
        term030 = mK.qyz.fma(R041, term030);
        E030 = term030.neg();
        DoubleVector term003 = mK.qxx.mul(R203);
        term003 = mK.qyy.fma(R023, term003);
        term003 = mK.qzz.fma(R005, term003);
        term003 = mK.qxy.fma(R113, term003);
        term003 = mK.qxz.fma(R104, term003);
        term003 = mK.qyz.fma(R014, term003);
        E003 = term003.neg();
        DoubleVector term210 = mK.qxx.mul(R410);
        term210 = mK.qyy.fma(R230, term210);
        term210 = mK.qzz.fma(R212, term210);
        term210 = mK.qxy.fma(R320, term210);
        term210 = mK.qxz.fma(R311, term210);
        term210 = mK.qyz.fma(R221, term210);
        E210 = term210.neg();
        DoubleVector term201 = mK.qxx.mul(R401);
        term201 = mK.qyy.fma(R221, term201);
        term201 = mK.qzz.fma(R203, term201);
        term201 = mK.qxy.fma(R311, term201);
        term201 = mK.qxz.fma(R302, term201);
        term201 = mK.qyz.fma(R212, term201);
        E201 = term201.neg();
        DoubleVector term120 = mK.qxx.mul(R320);
        term120 = mK.qyy.fma(R140, term120);
        term120 = mK.qzz.fma(R122, term120);
        term120 = mK.qxy.fma(R230, term120);
        term120 = mK.qxz.fma(R221, term120);
        term120 = mK.qyz.fma(R131, term120);
        E120 = term120.neg();
        DoubleVector term021 = mK.qxx.mul(R221);
        term021 = mK.qyy.fma(R041, term021);
        term021 = mK.qzz.fma(R023, term021);
        term021 = mK.qxy.fma(R131, term021);
        term021 = mK.qxz.fma(R122, term021);
        term021 = mK.qyz.fma(R032, term021);
        E021 = term021.neg();
        DoubleVector term102 = mK.qxx.mul(R302);
        term102 = mK.qyy.fma(R122, term102);
        term102 = mK.qzz.fma(R104, term102);
        term102 = mK.qxy.fma(R212, term102);
        term102 = mK.qxz.fma(R203, term102);
        term102 = mK.qyz.fma(R113, term102);
        E102 = term102.neg();
        DoubleVector term012 = mK.qxx.mul(R212);
        term012 = mK.qyy.fma(R032, term012);
        term012 = mK.qzz.fma(R014, term012);
        term012 = mK.qxy.fma(R122, term012);
        term012 = mK.qxz.fma(R113, term012);
        term012 = mK.qyz.fma(R023, term012);
        E012 = term012.neg();
        DoubleVector term111 = mK.qxx.mul(R311);
        term111 = mK.qyy.fma(R131, term111);
        term111 = mK.qzz.fma(R113, term111);
        term111 = mK.qxy.fma(R221, term111);
        term111 = mK.qxz.fma(R212, term111);
        term111 = mK.qyz.fma(R122, term111);
        E111 = term111.neg();
      case 2:
        // Order 2.
        DoubleVector term200 = mK.qxx.mul(R400);
        term200 = mK.qyy.fma(R220, term200);
        term200 = mK.qzz.fma(R202, term200);
        term200 = mK.qxy.fma(R310, term200);
        term200 = mK.qxz.fma(R301, term200);
        term200 = mK.qyz.fma(R211, term200);
        E200 = term200;
        DoubleVector term020 = mK.qxx.mul(R220);
        term020 = mK.qyy.fma(R040, term020);
        term020 = mK.qzz.fma(R022, term020);
        term020 = mK.qxy.fma(R130, term020);
        term020 = mK.qxz.fma(R121, term020);
        term020 = mK.qyz.fma(R031, term020);
        E020 = term020;
        DoubleVector term002 = mK.qxx.mul(R202);
        term002 = mK.qyy.fma(R022, term002);
        term002 = mK.qzz.fma(R004, term002);
        term002 = mK.qxy.fma(R112, term002);
        term002 = mK.qxz.fma(R103, term002);
        term002 = mK.qyz.fma(R013, term002);
        E002 = term002;
        DoubleVector term110 = mK.qxx.mul(R310);
        term110 = mK.qyy.fma(R130, term110);
        term110 = mK.qzz.fma(R112, term110);
        term110 = mK.qxy.fma(R220, term110);
        term110 = mK.qxz.fma(R211, term110);
        term110 = mK.qyz.fma(R121, term110);
        E110 = term110;
        DoubleVector term101 = mK.qxx.mul(R301);
        term101 = mK.qyy.fma(R121, term101);
        term101 = mK.qzz.fma(R103, term101);
        term101 = mK.qxy.fma(R211, term101);
        term101 = mK.qxz.fma(R202, term101);
        term101 = mK.qyz.fma(R112, term101);
        E101 = term101;
        DoubleVector term011 = mK.qxx.mul(R211);
        term011 = mK.qyy.fma(R031, term011);
        term011 = mK.qzz.fma(R013, term011);
        term011 = mK.qxy.fma(R121, term011);
        term011 = mK.qxz.fma(R112, term011);
        term011 = mK.qyz.fma(R022, term011);
        E011 = term011;
      case 1:
        // Order 1 (need a minus sign).
        DoubleVector term100 = mK.qxx.mul(R300);
        term100 = mK.qyy.fma(R120, term100);
        term100 = mK.qzz.fma(R102, term100);
        term100 = mK.qxy.fma(R210, term100);
        term100 = mK.qxz.fma(R201, term100);
        term100 = mK.qyz.fma(R111, term100);
        E100 = term100.neg();
        DoubleVector term010 = mK.qxx.mul(R210);
        term010 = mK.qyy.fma(R030, term010);
        term010 = mK.qzz.fma(R012, term010);
        term010 = mK.qxy.fma(R120, term010);
        term010 = mK.qxz.fma(R111, term010);
        term010 = mK.qyz.fma(R021, term010);
        E010 = term010.neg();
        DoubleVector term001 = mK.qxx.mul(R201);
        term001 = mK.qyy.fma(R021, term001);
        term001 = mK.qzz.fma(R003, term001);
        term001 = mK.qxy.fma(R111, term001);
        term001 = mK.qxz.fma(R102, term001);
        term001 = mK.qyz.fma(R012, term001);
        E001 = term001.neg();
      case 0:
        DoubleVector term000 = mK.qxx.mul(R200);
        term000 = mK.qyy.fma(R020, term000);
        term000 = mK.qzz.fma(R002, term000);
        term000 = mK.qxy.fma(R110, term000);
        term000 = mK.qxz.fma(R101, term000);
        term000 = mK.qyz.fma(R011, term000);
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
        DoubleVector term300 = mI.q.mul(R300);
        term300 = mI.dx.fma(R400.neg(), term300);
        term300 = mI.dy.fma(R310.neg(), term300);
        term300 = mI.dz.fma(R301.neg(), term300);
        term300 = mI.qxx.fma(R500, term300);
        term300 = mI.qyy.fma(R320, term300);
        term300 = mI.qzz.fma(R302, term300);
        term300 = mI.qxy.fma(R410, term300);
        term300 = mI.qxz.fma(R401, term300);
        term300 = mI.qyz.fma(R311, term300);
        E300 = term300;
        DoubleVector term030 = mI.q.mul(R030);
        term030 = mI.dx.fma(R130.neg(), term030);
        term030 = mI.dy.fma(R040.neg(), term030);
        term030 = mI.dz.fma(R031.neg(), term030);
        term030 = mI.qxx.fma(R230, term030);
        term030 = mI.qyy.fma(R050, term030);
        term030 = mI.qzz.fma(R032, term030);
        term030 = mI.qxy.fma(R140, term030);
        term030 = mI.qxz.fma(R131, term030);
        term030 = mI.qyz.fma(R041, term030);
        E030 = term030;
        DoubleVector term003 = mI.q.mul(R003);
        term003 = mI.dx.fma(R103.neg(), term003);
        term003 = mI.dy.fma(R013.neg(), term003);
        term003 = mI.dz.fma(R004.neg(), term003);
        term003 = mI.qxx.fma(R203, term003);
        term003 = mI.qyy.fma(R023, term003);
        term003 = mI.qzz.fma(R005, term003);
        term003 = mI.qxy.fma(R113, term003);
        term003 = mI.qxz.fma(R104, term003);
        term003 = mI.qyz.fma(R014, term003);
        E003 = term003;
        DoubleVector term210 = mI.q.mul(R210);
        term210 = mI.dx.fma(R310.neg(), term210);
        term210 = mI.dy.fma(R220.neg(), term210);
        term210 = mI.dz.fma(R211.neg(), term210);
        term210 = mI.qxx.fma(R410, term210);
        term210 = mI.qyy.fma(R230, term210);
        term210 = mI.qzz.fma(R212, term210);
        term210 = mI.qxy.fma(R320, term210);
        term210 = mI.qxz.fma(R311, term210);
        term210 = mI.qyz.fma(R221, term210);
        E210 = term210;
        DoubleVector term201 = mI.q.mul(R201);
        term201 = mI.dx.fma(R301.neg(), term201);
        term201 = mI.dy.fma(R211.neg(), term201);
        term201 = mI.dz.fma(R202.neg(), term201);
        term201 = mI.qxx.fma(R401, term201);
        term201 = mI.qyy.fma(R221, term201);
        term201 = mI.qzz.fma(R203, term201);
        term201 = mI.qxy.fma(R311, term201);
        term201 = mI.qxz.fma(R302, term201);
        term201 = mI.qyz.fma(R212, term201);
        E201 = term201;
        DoubleVector term120 = mI.q.mul(R120);
        term120 = mI.dx.fma(R220.neg(), term120);
        term120 = mI.dy.fma(R130.neg(), term120);
        term120 = mI.dz.fma(R121.neg(), term120);
        term120 = mI.qxx.fma(R320, term120);
        term120 = mI.qyy.fma(R140, term120);
        term120 = mI.qzz.fma(R122, term120);
        term120 = mI.qxy.fma(R230, term120);
        term120 = mI.qxz.fma(R221, term120);
        term120 = mI.qyz.fma(R131, term120);
        E120 = term120;
        DoubleVector term021 = mI.q.mul(R021);
        term021 = mI.dx.fma(R121.neg(), term021);
        term021 = mI.dy.fma(R031.neg(), term021);
        term021 = mI.dz.fma(R022.neg(), term021);
        term021 = mI.qxx.fma(R221, term021);
        term021 = mI.qyy.fma(R041, term021);
        term021 = mI.qzz.fma(R023, term021);
        term021 = mI.qxy.fma(R131, term021);
        term021 = mI.qxz.fma(R122, term021);
        term021 = mI.qyz.fma(R032, term021);
        E021 = term021;
        DoubleVector term102 = mI.q.mul(R102);
        term102 = mI.dx.fma(R202.neg(), term102);
        term102 = mI.dy.fma(R112.neg(), term102);
        term102 = mI.dz.fma(R103.neg(), term102);
        term102 = mI.qxx.fma(R302, term102);
        term102 = mI.qyy.fma(R122, term102);
        term102 = mI.qzz.fma(R104, term102);
        term102 = mI.qxy.fma(R212, term102);
        term102 = mI.qxz.fma(R203, term102);
        term102 = mI.qyz.fma(R113, term102);
        E102 = term102;
        DoubleVector term012 = mI.q.mul(R012);
        term012 = mI.dx.fma(R112.neg(), term012);
        term012 = mI.dy.fma(R022.neg(), term012);
        term012 = mI.dz.fma(R013.neg(), term012);
        term012 = mI.qxx.fma(R212, term012);
        term012 = mI.qyy.fma(R032, term012);
        term012 = mI.qzz.fma(R014, term012);
        term012 = mI.qxy.fma(R122, term012);
        term012 = mI.qxz.fma(R113, term012);
        term012 = mI.qyz.fma(R023, term012);
        E012 = term012;
        DoubleVector term111 = mI.q.mul(R111);
        term111 = mI.dx.fma(R211.neg(), term111);
        term111 = mI.dy.fma(R121.neg(), term111);
        term111 = mI.dz.fma(R112.neg(), term111);
        term111 = mI.qxx.fma(R311, term111);
        term111 = mI.qyy.fma(R131, term111);
        term111 = mI.qzz.fma(R113, term111);
        term111 = mI.qxy.fma(R221, term111);
        term111 = mI.qxz.fma(R212, term111);
        term111 = mI.qyz.fma(R122, term111);
        E111 = term111;
      case 2:
        DoubleVector term200 = mI.q.mul(R200);
        term200 = mI.dx.fma(R300.neg(), term200);
        term200 = mI.dy.fma(R210.neg(), term200);
        term200 = mI.dz.fma(R201.neg(), term200);
        term200 = mI.qxx.fma(R400, term200);
        term200 = mI.qyy.fma(R220, term200);
        term200 = mI.qzz.fma(R202, term200);
        term200 = mI.qxy.fma(R310, term200);
        term200 = mI.qxz.fma(R301, term200);
        term200 = mI.qyz.fma(R211, term200);
        E200 = term200;
        DoubleVector term020 = mI.q.mul(R020);
        term020 = mI.dx.fma(R120.neg(), term020);
        term020 = mI.dy.fma(R030.neg(), term020);
        term020 = mI.dz.fma(R021.neg(), term020);
        term020 = mI.qxx.fma(R220, term020);
        term020 = mI.qyy.fma(R040, term020);
        term020 = mI.qzz.fma(R022, term020);
        term020 = mI.qxy.fma(R130, term020);
        term020 = mI.qxz.fma(R121, term020);
        term020 = mI.qyz.fma(R031, term020);
        E020 = term020;
        DoubleVector term002 = mI.q.mul(R002);
        term002 = mI.dx.fma(R102.neg(), term002);
        term002 = mI.dy.fma(R012.neg(), term002);
        term002 = mI.dz.fma(R003.neg(), term002);
        term002 = mI.qxx.fma(R202, term002);
        term002 = mI.qyy.fma(R022, term002);
        term002 = mI.qzz.fma(R004, term002);
        term002 = mI.qxy.fma(R112, term002);
        term002 = mI.qxz.fma(R103, term002);
        term002 = mI.qyz.fma(R013, term002);
        E002 = term002;
        DoubleVector term110 = mI.q.mul(R110);
        term110 = mI.dx.fma(R210.neg(), term110);
        term110 = mI.dy.fma(R120.neg(), term110);
        term110 = mI.dz.fma(R111.neg(), term110);
        term110 = mI.qxx.fma(R310, term110);
        term110 = mI.qyy.fma(R130, term110);
        term110 = mI.qzz.fma(R112, term110);
        term110 = mI.qxy.fma(R220, term110);
        term110 = mI.qxz.fma(R211, term110);
        term110 = mI.qyz.fma(R121, term110);
        E110 = term110;
        DoubleVector term101 = mI.q.mul(R101);
        term101 = mI.dx.fma(R201.neg(), term101);
        term101 = mI.dy.fma(R111.neg(), term101);
        term101 = mI.dz.fma(R102.neg(), term101);
        term101 = mI.qxx.fma(R301, term101);
        term101 = mI.qyy.fma(R121, term101);
        term101 = mI.qzz.fma(R103, term101);
        term101 = mI.qxy.fma(R211, term101);
        term101 = mI.qxz.fma(R202, term101);
        term101 = mI.qyz.fma(R112, term101);
        E101 = term101;
        DoubleVector term011 = mI.q.mul(R011);
        term011 = mI.dx.fma(R111.neg(), term011);
        term011 = mI.dy.fma(R021.neg(), term011);
        term011 = mI.dz.fma(R012.neg(), term011);
        term011 = mI.qxx.fma(R211, term011);
        term011 = mI.qyy.fma(R031, term011);
        term011 = mI.qzz.fma(R013, term011);
        term011 = mI.qxy.fma(R121, term011);
        term011 = mI.qxz.fma(R112, term011);
        term011 = mI.qyz.fma(R022, term011);
        E011 = term011;
      case 1:
        DoubleVector term100 = mI.q.mul(R100);
        term100 = mI.dx.fma(R200.neg(), term100);
        term100 = mI.dy.fma(R110.neg(), term100);
        term100 = mI.dz.fma(R101.neg(), term100);
        term100 = mI.qxx.fma(R300, term100);
        term100 = mI.qyy.fma(R120, term100);
        term100 = mI.qzz.fma(R102, term100);
        term100 = mI.qxy.fma(R210, term100);
        term100 = mI.qxz.fma(R201, term100);
        term100 = mI.qyz.fma(R111, term100);
        E100 = term100;
        DoubleVector term010 = mI.q.mul(R010);
        term010 = mI.dx.fma(R110.neg(), term010);
        term010 = mI.dy.fma(R020.neg(), term010);
        term010 = mI.dz.fma(R011.neg(), term010);
        term010 = mI.qxx.fma(R210, term010);
        term010 = mI.qyy.fma(R030, term010);
        term010 = mI.qzz.fma(R012, term010);
        term010 = mI.qxy.fma(R120, term010);
        term010 = mI.qxz.fma(R111, term010);
        term010 = mI.qyz.fma(R021, term010);
        E010 = term010;
        DoubleVector term001 = mI.q.mul(R001);
        term001 = mI.dx.fma(R101.neg(), term001);
        term001 = mI.dy.fma(R011.neg(), term001);
        term001 = mI.dz.fma(R002.neg(), term001);
        term001 = mI.qxx.fma(R201, term001);
        term001 = mI.qyy.fma(R021, term001);
        term001 = mI.qzz.fma(R003, term001);
        term001 = mI.qxy.fma(R111, term001);
        term001 = mI.qxz.fma(R102, term001);
        term001 = mI.qyz.fma(R012, term001);
        E001 = term001;
      case 0:
        DoubleVector term000 = mI.q.mul(R000);
        term000 = mI.dx.fma(R100.neg(), term000);
        term000 = mI.dy.fma(R010.neg(), term000);
        term000 = mI.dz.fma(R001.neg(), term000);
        term000 = mI.qxx.fma(R200, term000);
        term000 = mI.qyy.fma(R020, term000);
        term000 = mI.qzz.fma(R002, term000);
        term000 = mI.qxy.fma(R110, term000);
        term000 = mI.qxz.fma(R101, term000);
        term000 = mI.qyz.fma(R011, term000);
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
        E300 = mI.q.mul(R300);
        E030 = mI.q.mul(R030);
        E003 = mI.q.mul(R003);
        E210 = mI.q.mul(R210);
        E201 = mI.q.mul(R201);
        E120 = mI.q.mul(R120);
        E021 = mI.q.mul(R021);
        E102 = mI.q.mul(R102);
        E012 = mI.q.mul(R012);
        E111 = mI.q.mul(R111);
      case 2:
        E200 = mI.q.mul(R200);
        E020 = mI.q.mul(R020);
        E002 = mI.q.mul(R002);
        E110 = mI.q.mul(R110);
        E101 = mI.q.mul(R101);
        E011 = mI.q.mul(R011);
      case 1:
        E100 = mI.q.mul(R100);
        E010 = mI.q.mul(R010);
        E001 = mI.q.mul(R001);
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
        DoubleVector term300 = uxi.mul(R400.neg());
        term300 = uyi.fma(R310.neg(), term300);
        term300 = uzi.fma(R301.neg(), term300);
        E300 = term300;
        DoubleVector term030 = uxi.mul(R130.neg());
        term030 = uyi.fma(R040.neg(), term030);
        term030 = uzi.fma(R031.neg(), term030);
        E030 = term030;
        DoubleVector term003 = uxi.mul(R103.neg());
        term003 = uyi.fma(R013.neg(), term003);
        term003 = uzi.fma(R004.neg(), term003);
        E003 = term003;
        DoubleVector term210 = uxi.mul(R310.neg());
        term210 = uyi.fma(R220.neg(), term210);
        term210 = uzi.fma(R211.neg(), term210);
        E210 = term210;
        DoubleVector term201 = uxi.mul(R301.neg());
        term201 = uyi.fma(R211.neg(), term201);
        term201 = uzi.fma(R202.neg(), term201);
        E201 = term201;
        DoubleVector term120 = uxi.mul(R220.neg());
        term120 = uyi.fma(R130.neg(), term120);
        term120 = uzi.fma(R121.neg(), term120);
        E120 = term120;
        DoubleVector term021 = uxi.mul(R121.neg());
        term021 = uyi.fma(R031.neg(), term021);
        term021 = uzi.fma(R022.neg(), term021);
        E021 = term021;
        DoubleVector term102 = uxi.mul(R202.neg());
        term102 = uyi.fma(R112.neg(), term102);
        term102 = uzi.fma(R103.neg(), term102);
        E102 = term102;
        DoubleVector term012 = uxi.mul(R112.neg());
        term012 = uyi.fma(R022.neg(), term012);
        term012 = uzi.fma(R013.neg(), term012);
        E012 = term012;
        DoubleVector term111 = uxi.mul(R211.neg());
        term111 = uyi.fma(R121.neg(), term111);
        term111 = uzi.fma(R112.neg(), term111);
        E111 = term111;
        // Fall through to 2nd order.
      case 2:
        // Order 2
        DoubleVector term200 = uxi.mul(R300.neg());
        term200 = uyi.fma(R210.neg(), term200);
        term200 = uzi.fma(R201.neg(), term200);
        E200 = term200;
        DoubleVector term020 = uxi.mul(R120.neg());
        term020 = uyi.fma(R030.neg(), term020);
        term020 = uzi.fma(R021.neg(), term020);
        E020 = term020;
        DoubleVector term002 = uxi.mul(R102.neg());
        term002 = uyi.fma(R012.neg(), term002);
        term002 = uzi.fma(R003.neg(), term002);
        E002 = term002;
        DoubleVector term110 = uxi.mul(R210.neg());
        term110 = uyi.fma(R120.neg(), term110);
        term110 = uzi.fma(R111.neg(), term110);
        E110 = term110;
        DoubleVector term101 = uxi.mul(R201.neg());
        term101 = uyi.fma(R111.neg(), term101);
        term101 = uzi.fma(R102.neg(), term101);
        E101 = term101;
        DoubleVector term011 = uxi.mul(R111.neg());
        term011 = uyi.fma(R021.neg(), term011);
        term011 = uzi.fma(R012.neg(), term011);
        E011 = term011;
        // Fall through to 1st order.
      case 1:
        // Order 1
        // This is d/dX of equation 3.1.3 in the Stone book.
        DoubleVector term100 = uxi.mul(R200.neg());
        term100 = uyi.fma(R110.neg(), term100);
        term100 = uzi.fma(R101.neg(), term100);
        E100 = term100;
        // This is d/dY of equation 3.1.3 in the Stone book.
        DoubleVector term010 = uxi.mul(R110.neg());
        term010 = uyi.fma(R020.neg(), term010);
        term010 = uzi.fma(R011.neg(), term010);
        E010 = term010;
        DoubleVector term001 = uxi.mul(R101.neg());
        term001 = uyi.fma(R011.neg(), term001);
        term001 = uzi.fma(R002.neg(), term001);
        E001 = term001;
        // Fall through to the potential.
      case 0:
        // This is equation 3.1.3 in the Stone book.
        DoubleVector term000 = uxi.mul(R100.neg());
        term000 = uyi.fma(R010.neg(), term000);
        term000 = uzi.fma(R001.neg(), term000);
        E000 = term000;
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
        DoubleVector term300 = mI.qxx.mul(R500);
        term300 = mI.qyy.fma(R320, term300);
        term300 = mI.qzz.fma(R302, term300);
        term300 = mI.qxy.fma(R410, term300);
        term300 = mI.qxz.fma(R401, term300);
        term300 = mI.qyz.fma(R311, term300);
        E300 = term300;
        DoubleVector term030 = mI.qxx.mul(R230);
        term030 = mI.qyy.fma(R050, term030);
        term030 = mI.qzz.fma(R032, term030);
        term030 = mI.qxy.fma(R140, term030);
        term030 = mI.qxz.fma(R131, term030);
        term030 = mI.qyz.fma(R041, term030);
        E030 = term030;
        DoubleVector term003 = mI.qxx.mul(R203);
        term003 = mI.qyy.fma(R023, term003);
        term003 = mI.qzz.fma(R005, term003);
        term003 = mI.qxy.fma(R113, term003);
        term003 = mI.qxz.fma(R104, term003);
        term003 = mI.qyz.fma(R014, term003);
        E003 = term003;
        DoubleVector term210 = mI.qxx.mul(R410);
        term210 = mI.qyy.fma(R230, term210);
        term210 = mI.qzz.fma(R212, term210);
        term210 = mI.qxy.fma(R320, term210);
        term210 = mI.qxz.fma(R311, term210);
        term210 = mI.qyz.fma(R221, term210);
        E210 = term210;
        DoubleVector term201 = mI.qxx.mul(R401);
        term201 = mI.qyy.fma(R221, term201);
        term201 = mI.qzz.fma(R203, term201);
        term201 = mI.qxy.fma(R311, term201);
        term201 = mI.qxz.fma(R302, term201);
        term201 = mI.qyz.fma(R212, term201);
        E201 = term201;
        DoubleVector term120 = mI.qxx.mul(R320);
        term120 = mI.qyy.fma(R140, term120);
        term120 = mI.qzz.fma(R122, term120);
        term120 = mI.qxy.fma(R230, term120);
        term120 = mI.qxz.fma(R221, term120);
        term120 = mI.qyz.fma(R131, term120);
        E120 = term120;
        DoubleVector term021 = mI.qxx.mul(R221);
        term021 = mI.qyy.fma(R041, term021);
        term021 = mI.qzz.fma(R023, term021);
        term021 = mI.qxy.fma(R131, term021);
        term021 = mI.qxz.fma(R122, term021);
        term021 = mI.qyz.fma(R032, term021);
        E021 = term021;
        DoubleVector term102 = mI.qxx.mul(R302);
        term102 = mI.qyy.fma(R122, term102);
        term102 = mI.qzz.fma(R104, term102);
        term102 = mI.qxy.fma(R212, term102);
        term102 = mI.qxz.fma(R203, term102);
        term102 = mI.qyz.fma(R113, term102);
        E102 = term102;
        DoubleVector term012 = mI.qxx.mul(R212);
        term012 = mI.qyy.fma(R032, term012);
        term012 = mI.qzz.fma(R014, term012);
        term012 = mI.qxy.fma(R122, term012);
        term012 = mI.qxz.fma(R113, term012);
        term012 = mI.qyz.fma(R023, term012);
        E012 = term012;
        DoubleVector term111 = mI.qxx.mul(R311);
        term111 = mI.qyy.fma(R131, term111);
        term111 = mI.qzz.fma(R113, term111);
        term111 = mI.qxy.fma(R221, term111);
        term111 = mI.qxz.fma(R212, term111);
        term111 = mI.qyz.fma(R122, term111);
        E111 = term111;
      case 2:
        DoubleVector term200 = mI.qxx.mul(R400);
        term200 = mI.qyy.fma(R220, term200);
        term200 = mI.qzz.fma(R202, term200);
        term200 = mI.qxy.fma(R310, term200);
        term200 = mI.qxz.fma(R301, term200);
        term200 = mI.qyz.fma(R211, term200);
        E200 = term200;
        DoubleVector term020 = mI.qxx.mul(R220);
        term020 = mI.qyy.fma(R040, term020);
        term020 = mI.qzz.fma(R022, term020);
        term020 = mI.qxy.fma(R130, term020);
        term020 = mI.qxz.fma(R121, term020);
        term020 = mI.qyz.fma(R031, term020);
        E020 = term020;
        DoubleVector term002 = mI.qxx.mul(R202);
        term002 = mI.qyy.fma(R022, term002);
        term002 = mI.qzz.fma(R004, term002);
        term002 = mI.qxy.fma(R112, term002);
        term002 = mI.qxz.fma(R103, term002);
        term002 = mI.qyz.fma(R013, term002);
        E002 = term002;
        DoubleVector term110 = mI.qxx.mul(R310);
        term110 = mI.qyy.fma(R130, term110);
        term110 = mI.qzz.fma(R112, term110);
        term110 = mI.qxy.fma(R220, term110);
        term110 = mI.qxz.fma(R211, term110);
        term110 = mI.qyz.fma(R121, term110);
        E110 = term110;
        DoubleVector term101 = mI.qxx.mul(R301);
        term101 = mI.qyy.fma(R121, term101);
        term101 = mI.qzz.fma(R103, term101);
        term101 = mI.qxy.fma(R211, term101);
        term101 = mI.qxz.fma(R202, term101);
        term101 = mI.qyz.fma(R112, term101);
        E101 = term101;
        DoubleVector term011 = mI.qxx.mul(R211);
        term011 = mI.qyy.fma(R031, term011);
        term011 = mI.qzz.fma(R013, term011);
        term011 = mI.qxy.fma(R121, term011);
        term011 = mI.qxz.fma(R112, term011);
        term011 = mI.qyz.fma(R022, term011);
        E011 = term011;
      case 1:
        DoubleVector term100 = mI.qxx.mul(R300);
        term100 = mI.qyy.fma(R120, term100);
        term100 = mI.qzz.fma(R102, term100);
        term100 = mI.qxy.fma(R210, term100);
        term100 = mI.qxz.fma(R201, term100);
        term100 = mI.qyz.fma(R111, term100);
        E100 = term100;
        DoubleVector term010 = mI.qxx.mul(R210);
        term010 = mI.qyy.fma(R030, term010);
        term010 = mI.qzz.fma(R012, term010);
        term010 = mI.qxy.fma(R120, term010);
        term010 = mI.qxz.fma(R111, term010);
        term010 = mI.qyz.fma(R021, term010);
        E010 = term010;
        DoubleVector term001 = mI.qxx.mul(R201);
        term001 = mI.qyy.fma(R021, term001);
        term001 = mI.qzz.fma(R003, term001);
        term001 = mI.qxy.fma(R111, term001);
        term001 = mI.qxz.fma(R102, term001);
        term001 = mI.qyz.fma(R012, term001);
        E001 = term001;
      case 0:
        DoubleVector term000 = mI.qxx.mul(R200);
        term000 = mI.qyy.fma(R020, term000);
        term000 = mI.qzz.fma(R002, term000);
        term000 = mI.qxy.fma(R110, term000);
        term000 = mI.qxz.fma(R101, term000);
        term000 = mI.qyz.fma(R011, term000);
        E000 = term000;
    }
  }
}
