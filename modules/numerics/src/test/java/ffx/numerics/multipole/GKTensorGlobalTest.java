package ffx.numerics.multipole;

import static org.junit.Assert.assertEquals;

import java.util.logging.Logger;
import org.junit.Test;

public class GKTensorGlobalTest {

  /** Logger for the MultipoleTensor class. */
  private static final Logger logger = Logger.getLogger(GKTensorGlobalTest.class.getName());

  private final double tolerance = 1.0e-9;

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
    assertEquals(" R400", x * x * x * x * A13 + 6.0 * x * x * A12 + 3.0 * A11, gkDipoleTensor.R400, tolerance);
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
    gkQuadrupoleTensor.generateTensor();

    // No charge potential.
    assertEquals(" R000", 0.0, gkQuadrupoleTensor.R000, tolerance);
    // No dipole potential.
    assertEquals(" R100", 0.0, gkQuadrupoleTensor.R100, tolerance);

    // Qxx (checks de-tracing).
    assertEquals(" R200", x * x * A20, gkQuadrupoleTensor.R200, tolerance);
    assertEquals(" R300", x * x * x * A21 + 2.0 * x * A20, gkQuadrupoleTensor.R300, tolerance);
    assertEquals(" R400", x * x * x * x * A22 + 5.0 * x * x * A21 + 2.0 * A20, gkQuadrupoleTensor.R400, tolerance);
    assertEquals(" R500", x * x * x * x * x * A23 + 9.0 * x * x * x * A22 + 12.0 * x * A21, gkQuadrupoleTensor.R500, tolerance);

    // Qxy
    assertEquals(" R110", x * y * A20, gkQuadrupoleTensor.R110, tolerance);
    assertEquals(" R210", x * x * y * A21 + y * A20, gkQuadrupoleTensor.R210, tolerance);
    assertEquals(" R111", x * y * z * A21, gkQuadrupoleTensor.R111, tolerance);
    assertEquals(" R310", x * x * x * y * A22 + 3.0 * x * y * A21, gkQuadrupoleTensor.R310, tolerance);
    assertEquals(" R220", x * x * y * y * A22 + x * x * A21 + y * y * A21 + A20, gkQuadrupoleTensor.R220, tolerance);
    assertEquals(" R211", x * x * y * z * A22 + y * z * A21, gkQuadrupoleTensor.R211, tolerance);
    assertEquals(" R410", x * x * x * x * y * A23 + 6.0 * x * x * y * A22 + 3.0 *y * A21, gkQuadrupoleTensor.R410, tolerance);
    assertEquals(" R320", x * x * x * y * y * A23 + x * x * x * A22 + 3.0 * x * y * y * A22 + 3.0 * x * A21, gkQuadrupoleTensor.R320, tolerance);
    assertEquals(" R311", x * x * x * y * z * A23 + 3.0 * x * y * z * A22, gkQuadrupoleTensor.R311, tolerance);
  }
}
