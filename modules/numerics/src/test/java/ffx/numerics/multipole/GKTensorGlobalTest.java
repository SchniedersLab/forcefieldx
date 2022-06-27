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

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class GKTensorGlobalTest {

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

    GKTensorGlobal GKTensorGlobal =
        new GKTensorGlobal(multipoleOrder, order, gc, Eh, Es);
    GKTensorGlobal.setR(r[0], r[1], r[2], Ai, Aj);
    GKTensorGlobal.source(work);

    // Test the "bn" method.
    double bn0 = GKTensorGlobal.bn(0);
    double bn1 = GKTensorGlobal.bn(1);
    double bn2 = GKTensorGlobal.bn(2);
    double bn3 = GKTensorGlobal.bn(3);
    double bn4 = GKTensorGlobal.bn(4);
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
    double A00 = GKTensorGlobal.anm(0, 0);
    double A01 = GKTensorGlobal.anm(0, 1);
    double A02 = GKTensorGlobal.anm(0, 2);
    double A03 = GKTensorGlobal.anm(0, 3);
    final double mapleA00 = -0.4048255705;
    final double mapleA01 = 0.04764329318;
    final double mapleA02 = -0.01266056821;
    final double mapleA03 = 0.004644003643;
    assertEquals(" A00", mapleA00, A00, tolerance);
    assertEquals(" A01", mapleA01, A01, tolerance);
    assertEquals(" A02", mapleA02, A02, tolerance);
    assertEquals(" A03", mapleA03, A03, tolerance);
    double B00 = GKTensorGlobal.bnm(0, 0);
    double B01 = GKTensorGlobal.bnm(0, 1);
    double B02 = GKTensorGlobal.bnm(0, 2);
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
    GKTensorGlobal = new GKTensorGlobal(multipoleOrder, order, gc, Eh, Es);
    GKTensorGlobal.setR(r[0], r[1], r[2], Ai, Aj);
    GKTensorGlobal.source(work);
    double A10 = GKTensorGlobal.anm(1, 0);
    double A11 = GKTensorGlobal.anm(1, 1);
    double A12 = GKTensorGlobal.anm(1, 2);
    double A13 = GKTensorGlobal.anm(1, 3);
    final double mapleA10 = 0.06764004545;
    final double mapleA11 = -0.02388135594;
    final double mapleA12 = 0.01196727051;
    final double mapleA13 = -0.007470574812;
    assertEquals(" A10", mapleA10, A10, tolerance);
    assertEquals(" A11", mapleA11, A11, tolerance);
    assertEquals(" A12", mapleA12, A12, tolerance);
    assertEquals(" A13", mapleA13, A13, tolerance);
    // Dipole potential Born chain-rule derivatives.
    double B10 = GKTensorGlobal.bnm(1, 0);
    double B11 = GKTensorGlobal.bnm(1, 1);
    double B12 = GKTensorGlobal.bnm(1, 2);
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
    GKTensorGlobal = new GKTensorGlobal(multipoleOrder, order, gc, Eh, Es);
    GKTensorGlobal.setR(r[0], r[1], r[2], Ai, Aj);
    GKTensorGlobal.source(work);
    double A20 = GKTensorGlobal.anm(2, 0);
    double A21 = GKTensorGlobal.anm(2, 1);
    double A22 = GKTensorGlobal.anm(2, 2);
    double A23 = GKTensorGlobal.anm(2, 3);
    final double mapleA20 = -0.03404928262;
    final double mapleA21 = 0.02003603617;
    final double mapleA22 = -0.01475634878;
    final double mapleA23 = 0.01280244361;
    assertEquals(" A20", mapleA20, A20, tolerance);
    assertEquals(" A21", mapleA21, A21, tolerance);
    assertEquals(" A22", mapleA22, A22, tolerance);
    assertEquals(" A23", mapleA23, A23, tolerance);
    // Quadrupole potential Born chain-rule derivatives.
    double B20 = GKTensorGlobal.bnm(2, 0);
    double B21 = GKTensorGlobal.bnm(2, 1);
    double B22 = GKTensorGlobal.bnm(2, 2);
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

    GKTensorGlobal gkMonopoleTensor =
        new GKTensorGlobal(multipoleOrder, order, gc, Eh, Es);
    gkMonopoleTensor.setR(r[0], r[1], r[2], Ai, Aj);
    gkMonopoleTensor.source(work);

    // Monopole potential and derivatives.
    double A00 = gkMonopoleTensor.anm(0, 0);
    double A01 = gkMonopoleTensor.anm(0, 1);
    double A02 = gkMonopoleTensor.anm(0, 2);
    double A03 = gkMonopoleTensor.anm(0, 3);
    gkMonopoleTensor.generateTensor();

    double x = r[0];
    assertEquals(" R000", A00, gkMonopoleTensor.R000, tolerance);
    assertEquals(" R100", x * A01, gkMonopoleTensor.R100, tolerance);
    assertEquals(" R200", x * x * A02 + A01, gkMonopoleTensor.R200, tolerance);
    assertEquals(" R300", x * x * x * A03 + 3.0 * x * A02, gkMonopoleTensor.R300, tolerance);
  }

  @Test
  public void dipoleTensorTest() {
    double[] r = {0.7, 0.8, 0.9};
    double x = r[0];
    double y = r[1];
    double z = r[2];

    // Dipole potential and derivatives.
    int order = 4;
    int multipoleOrder = 1;
    double[] work = new double[order + 1];
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(multipoleOrder, order, gc, Eh, Es);
    gkDipoleTensor.setR(r[0], r[1], r[2], Ai, Aj);
    gkDipoleTensor.source(work);
    double A10 = gkDipoleTensor.anm(1, 0);
    double A11 = gkDipoleTensor.anm(1, 1);
    double A12 = gkDipoleTensor.anm(1, 2);
    double A13 = gkDipoleTensor.anm(1, 3);
    gkDipoleTensor.generateTensor();

    // No charge potential.
    assertEquals(" R000", 0.0, gkDipoleTensor.R000, tolerance);

    assertEquals(" R100", x * A10, gkDipoleTensor.R100, tolerance);
    assertEquals(" R200", x * x * A11 + A10, gkDipoleTensor.R200, tolerance);
    assertEquals(" R300", x * x * x * A12 + 3.0 * x * A11, gkDipoleTensor.R300, tolerance);
    assertEquals(" R400", x * x * x * x * A13 + 6.0 * x * x * A12 + 3.0 * A11, gkDipoleTensor.R400,
        tolerance);
  }

  @Test
  public void quadrupoleTensorTest() {
    double[] r = {0.7, 0.8, 0.9};
    double x = r[0];
    double y = r[1];
    double z = r[2];

    // Quadrupole potential and derivatives.
    int order = 5;
    int multipoleOrder = 2;
    double[] work = new double[order + 1];
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(multipoleOrder, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(r[0], r[1], r[2], Ai, Aj);
    gkQuadrupoleTensor.source(work);
    double A20 = gkQuadrupoleTensor.anm(2, 0);
    double A21 = gkQuadrupoleTensor.anm(2, 1);
    double A22 = gkQuadrupoleTensor.anm(2, 2);
    double A23 = gkQuadrupoleTensor.anm(2, 3);

    // int tensorCount = MultipoleTensor.tensorCount(order);
    // double[] tensor = new double[tensorCount];
    // gkQuadrupoleTensor.noStorageRecursion(tensor);

    gkQuadrupoleTensor.generateTensor();

    // assertEquals(" R200", tensor[MultipoleTensor.ti(2, 0, 0, order)], gkQuadrupoleTensor.R200, tolerance);

    // No charge potential.
    assertEquals(" R000", 0.0, gkQuadrupoleTensor.R000, tolerance);
    // No dipole potential.
    assertEquals(" R100", 0.0, gkQuadrupoleTensor.R100, tolerance);

    // Potential for the quadrupole trace
    assertEquals(" R200", x * x * A20, gkQuadrupoleTensor.R200, tolerance);
    assertEquals(" R020", y * y * A20, gkQuadrupoleTensor.R020, tolerance);
    assertEquals(" R002", z * z * A20, gkQuadrupoleTensor.R002, tolerance);

    // Gradient of the trace with respect to X (note that the x*A20 term sums to 2*x*A20 over the trace).
    assertEquals(" R300", x * x * x * A21 + 3.0 * x * A20, gkQuadrupoleTensor.R300, tolerance);
    assertEquals(" R120", x * y * y * A21 + x * A20, gkQuadrupoleTensor.R120, tolerance);
    assertEquals(" R102", x * z * z * A21 + x * A20, gkQuadrupoleTensor.R102, tolerance);

    // Higher order Qxx terms.
    assertEquals(" R400", x * x * x * x * A22 + 6.0 * x * x * A21 + 3.0 * A20, gkQuadrupoleTensor.R400, tolerance);
    assertEquals(" R500", x * x * x * x * x * A23 + 10.0 * x * x * x * A22 + 15.0 * x * A21, gkQuadrupoleTensor.R500, tolerance);

    // Qxy
    assertEquals(" R110", x * y * A20, gkQuadrupoleTensor.R110, tolerance);
    assertEquals(" R111", x * y * z * A21, gkQuadrupoleTensor.R111, tolerance);
    assertEquals(" R310", x * x * x * y * A22 + 3.0 * x * y * A21, gkQuadrupoleTensor.R310, tolerance);
    assertEquals(" R220", x * x * y * y * A22 + x * x * A21 + y * y * A21 + A20, gkQuadrupoleTensor.R220, tolerance);
    assertEquals(" R211", x * x * y * z * A22 + y * z * A21, gkQuadrupoleTensor.R211, tolerance);
    assertEquals(" R410", x * x * x * x * y * A23 + 6.0 * x * x * y * A22 + 3.0 * y * A21, gkQuadrupoleTensor.R410, tolerance);
    assertEquals(" R320", x * x * x * y * y * A23 + x * x * x * A22 + 3.0 * x * y * y * A22 + 3.0 * x * A21, gkQuadrupoleTensor.R320, tolerance);
    assertEquals(" R311", x * x * x * y * z * A23 + 3.0 * x * y * z * A22, gkQuadrupoleTensor.R311, tolerance);
  }

  @Test
  public void chargeFiniteDifferenceTest() {
    int order = 6;
    double[] r = {0.7, 0.8, 0.9};
    double x = r[0];
    double y = r[1];
    double z = r[2];

    GKTensorGlobal gkTensorGlobal = new GKTensorGlobal(0, order, gc, Eh, Es);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
    gkTensorGlobal.setR(x, y, z, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensor);
    double[] tensorsPx = new double[tensorCount];
    double[] tensorsNx = new double[tensorCount];
    double[] tensorsPy = new double[tensorCount];
    double[] tensorsNy = new double[tensorCount];
    double[] tensorsPz = new double[tensorCount];
    double[] tensorsNz = new double[tensorCount];
    double delta = 1.0e-5;
    double delta2 = delta * 2;
    r[0] += delta;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsPx);
    r[0] -= delta2;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsNx);
    r[0] += delta;

    r[1] += delta;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsPy);
    r[1] -= delta2;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsNy);
    r[1] += delta;

    r[2] += delta;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsPz);
    r[2] -= delta2;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsNz);
    r[2] += delta;

    tensorFiniteDifference(gkTensorGlobal, delta2, order, tensor, tensorsPx, tensorsNx, tensorsPy,
        tensorsNy, tensorsPz, tensorsNz);

    // Order(L^4) recursion.
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.recursion(tensor);
    tensorFiniteDifference(gkTensorGlobal, delta2, order, tensor, tensorsPx, tensorsNx, tensorsPy,
        tensorsNy, tensorsPz, tensorsNz);

    // Machine generated code.
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.generateTensor();
    gkTensorGlobal.getTensor(tensor);
    tensorFiniteDifference(gkTensorGlobal, delta2, order, tensor, tensorsPx, tensorsNx, tensorsPy,
        tensorsNy, tensorsPz, tensorsNz);
  }

  @Test
  public void dipoleFiniteDifferenceTest() {
    int order = 6;
    int mulitpoleOrder = 1;
    double[] r = {0.7, 0.8, 0.9};
    double x = r[0];
    double y = r[1];
    double z = r[2];

    GKTensorGlobal gkTensorGlobal = new GKTensorGlobal(mulitpoleOrder, order, gc, Eh, Es);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
    gkTensorGlobal.setR(x, y, z, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensor);
    double[] tensorsPx = new double[tensorCount];
    double[] tensorsNx = new double[tensorCount];
    double[] tensorsPy = new double[tensorCount];
    double[] tensorsNy = new double[tensorCount];
    double[] tensorsPz = new double[tensorCount];
    double[] tensorsNz = new double[tensorCount];
    double delta = 1.0e-5;
    double delta2 = delta * 2;
    r[0] += delta;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsPx);
    r[0] -= delta2;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsNx);
    r[0] += delta;

    r[1] += delta;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsPy);
    r[1] -= delta2;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsNy);
    r[1] += delta;

    r[2] += delta;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsPz);
    r[2] -= delta2;
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.noStorageRecursion(tensorsNz);
    r[2] += delta;

    tensorFiniteDifference(gkTensorGlobal, delta2, order, tensor, tensorsPx, tensorsNx, tensorsPy,
        tensorsNy, tensorsPz, tensorsNz);

    // Order(L^4) recursion.
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.recursion(tensor);
    tensorFiniteDifference(gkTensorGlobal, delta2, order, tensor, tensorsPx, tensorsNx, tensorsPy,
        tensorsNy, tensorsPz, tensorsNz);

    // Machine generated code.
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.generateTensor();
    gkTensorGlobal.getTensor(tensor);
    tensorFiniteDifference(gkTensorGlobal, delta2, order, tensor, tensorsPx, tensorsNx, tensorsPy,
        tensorsNy, tensorsPz, tensorsNz);
  }

  private void tensorFiniteDifference(GKTensorGlobal multipoleTensor,
      double delta2, int order, double[] tensor,
      double[] tensorsPx, double[] tensorsNx,
      double[] tensorsPy, double[] tensorsNy,
      double[] tensorsPz, double[] tensorsNz) {

    int start = multipoleTensor.multipoleOrder;

    String info = "GK Global";
    // Test the partial derivatives for all tensor components.
    for (int l = start; l < order; l++) {
      // Test X derivative
      double expect = tensor[multipoleTensor.ti(l + 1, 0, 0)];
      double actual = (tensorsPx[multipoleTensor.ti(l, 0, 0)] - tensorsNx[multipoleTensor.ti(l, 0, 0)]) / delta2;
      assertEquals(info + " @ " + l, expect, actual, fdTolerance);
      // Test Y derivative
      expect = tensor[multipoleTensor.ti(l, 1, 0)];
      actual =
          (tensorsPy[multipoleTensor.ti(l, 0, 0)] - tensorsNy[multipoleTensor.ti(l, 0, 0)]) / delta2;
      assertEquals(info + " @ " + l, expect, actual, fdTolerance);
      // Test Z derivative
      expect = tensor[multipoleTensor.ti(l, 0, 1)];
      actual =
          (tensorsPz[multipoleTensor.ti(l, 0, 0)] - tensorsNz[multipoleTensor.ti(l, 0, 0)]) / delta2;
      assertEquals(info + " @ " + l, expect, actual, fdTolerance);
    }
    for (int l = 0; l < order; l++) {
      for (int m = 1; m < order - l; m++) {
        // Test X derivative
        double expect = tensor[multipoleTensor.ti(l + 1, m, 0)];
        double actual = (tensorsPx[multipoleTensor.ti(l, m, 0)] - tensorsNx[multipoleTensor.ti(l, m, 0)]) / delta2;
        assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
        // Test Y derivative
        expect = tensor[multipoleTensor.ti(l, m + 1, 0)];
        actual = (tensorsPy[multipoleTensor.ti(l, m, 0)] - tensorsNy[multipoleTensor.ti(l, m, 0)]) / delta2;
        assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
        // Test Z derivative
        expect = tensor[multipoleTensor.ti(l, m, 1)];
        actual = (tensorsPz[multipoleTensor.ti(l, m, 0)] - tensorsNz[multipoleTensor.ti(l, m, 0)]) / delta2;
        assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l + m + n = 1..o (n >= 1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m < order - l; m++) {
        for (int n = 1; n < order - l - m; n++) {
          // Test X derivative
          double expect = tensor[multipoleTensor.ti(l + 1, m, n)];
          double actual =
              (tensorsPx[multipoleTensor.ti(l, m, n)] - tensorsNx[multipoleTensor.ti(l, m, n)])
                  / delta2;
          assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
          // Test Y derivative
          expect = tensor[multipoleTensor.ti(l, m + 1, n)];
          actual =
              (tensorsPy[multipoleTensor.ti(l, m, n)] - tensorsNy[multipoleTensor.ti(l, m, n)])
                  / delta2;
          assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
          // Test Z derivative
          expect = tensor[multipoleTensor.ti(l, m, n + 1)];
          actual =
              (tensorsPz[multipoleTensor.ti(l, m, n)] - tensorsNz[multipoleTensor.ti(l, m, n)])
                  / delta2;
          assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
        }
      }
    }
  }
}
