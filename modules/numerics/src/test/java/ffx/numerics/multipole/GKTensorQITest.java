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

import static ffx.numerics.math.DoubleMath.length;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GKTensorQITest {

  private final double tolerance = 1.0e-9;
  private final double fdTolerance = 1.0e-6;

  private static final double gc = 2.455;
  private static final double Eh = 1.0;
  private static final double Es = 78.3;
  private static final double Ai = 2.0;
  private static final double Aj = 2.0;
  private static final double[] r = {1.0, 1.0, 1.0};

  @Test
  public void tensorAuxiliaryTest() {

    double x = r[0];
    double y = r[1];
    double z = r[2];

    int order = 3;
    int multipoleOrder = 0;
    double[] work = new double[order + 1];

    GKTensorQI tensorQI = new GKTensorQI(multipoleOrder, order, gc, Eh, Es);
    tensorQI.setR(r[0], r[1], r[2], Ai, Aj);
    tensorQI.source(work);

    // Test the "bn" method.
    double bn0 = tensorQI.bn(0);
    double bn1 = tensorQI.bn(1);
    double bn2 = tensorQI.bn(2);
    double bn3 = tensorQI.bn(3);
    double bn4 = tensorQI.bn(4);
    double mapleBN0 = 0.4809168788;
    double mapleBN1 = -0.02292037702;
    double mapleBN2 = -0.01061215012;
    double mapleBN3 = 0.005273401519;
    double mapleBN4 = -0.001707834828;
    assertEquals(" bn0", mapleBN0, bn0, tolerance);
    assertEquals(" bn1", mapleBN1, bn1, tolerance);
    assertEquals(" bn2", mapleBN2, bn2, tolerance);
    assertEquals(" bn3", mapleBN3, bn3, tolerance);
    assertEquals(" bn4", mapleBN4, bn4, tolerance);

    // Monopole potential and derivatives.
    double A00 = tensorQI.anm(0, 0);
    double A01 = tensorQI.anm(0, 1);
    double A02 = tensorQI.anm(0, 2);
    double A03 = tensorQI.anm(0, 3);
    final double mapleA00 = -0.4048255705;
    final double mapleA01 = 0.04764329318;
    final double mapleA02 = -0.01266056821;
    final double mapleA03 = 0.004644003643;
    assertEquals(" A00", mapleA00, A00, tolerance);
    assertEquals(" A01", mapleA01, A01, tolerance);
    assertEquals(" A02", mapleA02, A02, tolerance);
    assertEquals(" A03", mapleA03, A03, tolerance);
    double B00 = tensorQI.bnm(0, 0);
    double B01 = tensorQI.bnm(0, 1);
    double B02 = tensorQI.bnm(0, 2);
    // Monpole potential Born chain-rule derivatives.
    final double mapleB00 = 0.03273696138;
    final double mapleB01 = -0.01311852185;
    final double mapleB02 = 0.006171353772;
    assertEquals(" B00", mapleB00, B00, tolerance);
    assertEquals(" B01", mapleB01, B01, tolerance);
    assertEquals(" B02", mapleB02, B02, tolerance);

    // Dipole potential and derivatives.
    order = 4;
    multipoleOrder = 1;
    work = new double[order + 1];
    tensorQI = new GKTensorQI(multipoleOrder, order, gc, Eh, Es);
    tensorQI.setR(r[0], r[1], r[2], Ai, Aj);
    tensorQI.source(work);
    double A10 = tensorQI.anm(1, 0);
    double A11 = tensorQI.anm(1, 1);
    double A12 = tensorQI.anm(1, 2);
    double A13 = tensorQI.anm(1, 3);
    final double mapleA10 = 0.06764004545;
    final double mapleA11 = -0.02388135594;
    final double mapleA12 = 0.01196727051;
    final double mapleA13 = -0.007470574812;
    assertEquals(" A10", mapleA10, A10, tolerance);
    assertEquals(" A11", mapleA11, A11, tolerance);
    assertEquals(" A12", mapleA12, A12, tolerance);
    assertEquals(" A13", mapleA13, A13, tolerance);
    // Dipole potential Born chain-rule derivatives.
    double B10 = tensorQI.bnm(1, 0);
    double B11 = tensorQI.bnm(1, 1);
    double B12 = tensorQI.bnm(1, 2);
    final double mapleB10 = -0.01640950855;
    final double mapleB11 = 0.01043812102;
    final double mapleB12 = -0.007669896149;
    assertEquals(" B10", mapleB10, B10, tolerance);
    assertEquals(" B11", mapleB11, B11, tolerance);
    assertEquals(" B12", mapleB12, B12, tolerance);

    // Quadrupole potential and derivatives.
    order = 5;
    multipoleOrder = 2;
    work = new double[order + 1];
    tensorQI = new GKTensorQI(multipoleOrder, order, gc, Eh, Es);
    tensorQI.setR(r[0], r[1], r[2], Ai, Aj);
    tensorQI.source(work);
    double A20 = tensorQI.anm(2, 0);
    double A21 = tensorQI.anm(2, 1);
    double A22 = tensorQI.anm(2, 2);
    double A23 = tensorQI.anm(2, 3);
    final double mapleA20 = -0.03404928262;
    final double mapleA21 = 0.02003603617;
    final double mapleA22 = -0.01475634878;
    final double mapleA23 = 0.01280244361;
    assertEquals(" A20", mapleA20, A20, tolerance);
    assertEquals(" A21", mapleA21, A21, tolerance);
    assertEquals(" A22", mapleA22, A22, tolerance);
    assertEquals(" A23", mapleA23, A23, tolerance);
    // Quadrupole potential Born chain-rule derivatives.
    double B20 = tensorQI.bnm(2, 0);
    double B21 = tensorQI.bnm(2, 1);
    double B22 = tensorQI.bnm(2, 2);
    final double mapleB20 = 0.01376728808;
    final double mapleB21 = -0.01199790086;
    final double mapleB22 = 0.01179997601;
    assertEquals(" B20", mapleB20, B20, tolerance);
    assertEquals(" B21", mapleB21, B21, tolerance);
    assertEquals(" B22", mapleB22, B22, tolerance);
  }

  @Test
  public void chargeTensorTest() {
    double[] r = {0.7, 0.8, 0.9};
    int order = 3;
    int multipoleOrder = 0;
    double[] work = new double[order + 1];

    // Monopole potential and derivatives.
    GKTensorQI gkMonopoleTensor = new GKTensorQI(multipoleOrder, order, gc, Eh, Es);
    gkMonopoleTensor.setR(r[0], r[1], r[2], Ai, Aj);
    gkMonopoleTensor.source(work);

    double A00 = gkMonopoleTensor.anm(0, 0);
    double A01 = gkMonopoleTensor.anm(0, 1);
    double A02 = gkMonopoleTensor.anm(0, 2);
    double A03 = gkMonopoleTensor.anm(0, 3);
    gkMonopoleTensor.generateTensor();

    double z = length(r);
    assertEquals(" R000", A00, gkMonopoleTensor.R000, tolerance);
    assertEquals(" R001", z * A01, gkMonopoleTensor.R001, tolerance);
    assertEquals(" R002", z * z * A02 + A01, gkMonopoleTensor.R002, tolerance);
    assertEquals(" R003", z * z * z * A03 + 3.0 * z * A02, gkMonopoleTensor.R003, tolerance);
  }

  @Test
  public void dipoleTensorTest() {
    double[] r = {0.7, 0.8, 0.9};
    int order = 4;
    int multipoleOrder = 1;
    double[] work = new double[order + 1];

    // Dipole potential and derivatives.
    GKTensorQI gkDipoleTensor = new GKTensorQI(multipoleOrder, order, gc, Eh, Es);
    gkDipoleTensor.setR(r[0], r[1], r[2], Ai, Aj);
    gkDipoleTensor.source(work);

    double A10 = gkDipoleTensor.anm(1, 0);
    double A11 = gkDipoleTensor.anm(1, 1);
    double A12 = gkDipoleTensor.anm(1, 2);
    double A13 = gkDipoleTensor.anm(1, 3);
    gkDipoleTensor.generateTensor();

    // No charge potential.
    assertEquals(" R000", 0.0, gkDipoleTensor.R000, tolerance);

    double z = length(r);
    assertEquals(" R001", z * A10, gkDipoleTensor.R001, tolerance);
    assertEquals(" R002", z * z * A11 + A10, gkDipoleTensor.R002, tolerance);
    assertEquals(" R003", z * z * z * A12 + 3.0 * z * A11, gkDipoleTensor.R003, tolerance);
    assertEquals(" R004", z * z * z * z * A13 + 6.0 * z * z * A12 + 3.0 * A11, gkDipoleTensor.R004,
        tolerance);
  }

  @Test
  public void quadrupoleTensorTest() {
    double[] r = {0.7, 0.8, 0.9};
    int order = 5;
    int multipoleOrder = 2;
    double[] work = new double[order + 1];

    // Quadrupole potential and derivatives.
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(multipoleOrder, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(r[0], r[1], r[2], Ai, Aj);
    gkQuadrupoleTensor.source(work);

    double A20 = gkQuadrupoleTensor.anm(2, 0);
    double A21 = gkQuadrupoleTensor.anm(2, 1);
    double A22 = gkQuadrupoleTensor.anm(2, 2);
    double A23 = gkQuadrupoleTensor.anm(2, 3);
    gkQuadrupoleTensor.generateTensor();

    // No charge potential.
    assertEquals(" R000", 0.0, gkQuadrupoleTensor.R000, tolerance);
    // No dipole potential.
    assertEquals(" R100", 0.0, gkQuadrupoleTensor.R100, tolerance);

    double z = length(r);
    double x = 0.0;
    double y = 0.0;

    // Trace.
    assertEquals(" R200", x * x * A20, gkQuadrupoleTensor.R200, tolerance);
    assertEquals(" R020", y * y * A20, gkQuadrupoleTensor.R020, tolerance);
    assertEquals(" R002", z * z * A20, gkQuadrupoleTensor.R002, tolerance);

    // Gradient of the trace with respect to z (note that the x*A20 term sums to 2*x*A20 over the trace).
    assertEquals(" R003", z * z * z * A21 + 3.0 * z * A20, gkQuadrupoleTensor.R003, tolerance);
    assertEquals(" R201", z * A20, gkQuadrupoleTensor.R201, tolerance);
    assertEquals(" R021", z * A20, gkQuadrupoleTensor.R021, tolerance);

    // Higher order Qzz terms.
    assertEquals(" R004", z * z * z * z * A22 + 6.0 * z * z * A21 + 3.0 * A20,
        gkQuadrupoleTensor.R004, tolerance);
    assertEquals(" R005", z * z * z * z * z * A23 + 10.0 * z * z * z * A22 + 15.0 * z * A21,
        gkQuadrupoleTensor.R005, tolerance);

    // Qxy
    assertEquals(" R110", 0.0, gkQuadrupoleTensor.R110, tolerance);
    assertEquals(" R111", 0.0, gkQuadrupoleTensor.R111, tolerance);
    assertEquals(" R310", 0.0, gkQuadrupoleTensor.R310, tolerance);
    assertEquals(" R220", A20, gkQuadrupoleTensor.R220, tolerance);
    assertEquals(" R211", 0.0, gkQuadrupoleTensor.R211, tolerance);
    assertEquals(" R401", 3.0 * z * A21, gkQuadrupoleTensor.R401, tolerance);
    assertEquals(" R023", z * z * z * A22 + 3.0 * z * A21, gkQuadrupoleTensor.R023, tolerance);
    assertEquals(" R311", 0.0, gkQuadrupoleTensor.R311, tolerance);
  }

  @Test
  public void chargeFiniteDifferenceTest() {
    int order = 6;
    double[] r = {0.7, 0.8, 0.9};

    r[2] = length(r);
    r[0] = 0.0;
    r[1] = 0.0;

    GKTensorQI gkTensorQI = new GKTensorQI(0, order, gc, Eh, Es);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.noStorageRecursion(tensor);
    double[] tensorsPz = new double[tensorCount];
    double[] tensorsNz = new double[tensorCount];
    double delta = 1.0e-5;
    double delta2 = delta * 2;

    r[2] += delta;
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.noStorageRecursion(tensorsPz);
    r[2] -= delta2;
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.noStorageRecursion(tensorsNz);
    r[2] += delta;

    tensorFiniteDifference(gkTensorQI, delta2, order, tensor, tensorsPz, tensorsNz);

    // Order(L^4) recursion.
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.recursion(tensor);
    tensorFiniteDifference(gkTensorQI, delta2, order, tensor, tensorsPz, tensorsNz);

    // Machine generated code.
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.generateTensor();
    gkTensorQI.getTensor(tensor);
    tensorFiniteDifference(gkTensorQI, delta2, order, tensor, tensorsPz, tensorsNz);
  }

  @Test
  public void dipoleFiniteDifferenceTest() {
    int order = 6;
    int mulitpoleOrder = 1;
    double[] r = {0.7, 0.8, 0.9};
    r[2] = length(r);
    r[0] = 0.0;
    r[1] = 0.0;

    GKTensorQI gkTensorQI = new GKTensorQI(mulitpoleOrder, order, gc, Eh, Es);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.noStorageRecursion(tensor);
    double[] tensorsPz = new double[tensorCount];
    double[] tensorsNz = new double[tensorCount];
    double delta = 1.0e-5;
    double delta2 = delta * 2;

    r[2] += delta;
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.noStorageRecursion(tensorsPz);
    r[2] -= delta2;
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.noStorageRecursion(tensorsNz);
    r[2] += delta;

    tensorFiniteDifference(gkTensorQI, delta2, order, tensor, tensorsPz, tensorsNz);

    // Order(L^4) recursion.
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.recursion(tensor);
    tensorFiniteDifference(gkTensorQI, delta2, order, tensor, tensorsPz, tensorsNz);

    // Machine generated code.
    gkTensorQI.setR(r, Ai, Aj);
    gkTensorQI.generateTensor();
    gkTensorQI.getTensor(tensor);
    tensorFiniteDifference(gkTensorQI, delta2, order, tensor, tensorsPz, tensorsNz);
  }

  private void tensorFiniteDifference(GKTensorQI gkTensorQI, double delta2, int order,
      double[] tensor, double[] tensorsPz, double[] tensorsNz) {

    int start = gkTensorQI.multipoleOrder;

    String info = "QK QI";

    // Test the partial derivatives for all tensor components.
    for (int l = start; l < order; l++) {
      // Test Z derivative
      double expect = tensor[gkTensorQI.ti(l, 0, 1)];
      double actual =
          (tensorsPz[gkTensorQI.ti(l, 0, 0)] - tensorsNz[gkTensorQI.ti(l, 0, 0)]) / delta2;
      assertEquals(info + " (d/dx): " + l, expect, actual, fdTolerance);
    }
    for (int l = 0; l < order; l++) {
      for (int m = 1; m < order - l; m++) {
        // Test Z derivative
        double expect = tensor[gkTensorQI.ti(l, m, 1)];
        double actual =
            (tensorsPz[gkTensorQI.ti(l, m, 0)] - tensorsNz[gkTensorQI.ti(l, m, 0)])
                / delta2;
        assertEquals(info + " (d/dx): " + l + " (d/dy): " + m, expect, actual, fdTolerance);
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m < order - l; m++) {
        for (int n = 1; n < order - l - m; n++) {
          // Test Z derivative
          double expect = tensor[gkTensorQI.ti(l, m, n + 1)];
          double actual =
              (tensorsPz[gkTensorQI.ti(l, m, n)] - tensorsNz[gkTensorQI.ti(l, m, n)])
                  / delta2;
          assertEquals(info + " (d/dx): " + l + " (d/dy): " + m + " (d/dz):" + n, expect, actual,
              fdTolerance);
        }
      }
    }
  }
}
