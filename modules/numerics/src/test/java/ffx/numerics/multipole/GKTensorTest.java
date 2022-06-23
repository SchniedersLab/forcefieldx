package ffx.numerics.multipole;

import static org.junit.Assert.assertEquals;

import java.util.logging.Logger;
import org.junit.Test;

public class GKTensorTest {

  /** Logger for the MultipoleTensor class. */
  private static final Logger logger = Logger.getLogger(GKTensorTest.class.getName());

  private final double tolerance = 1.0e-9;

  private static final double gc = 2.455;
  private static final double Eh = 1.0;
  private static final double Es = 78.3;
  private static final double Ai = 2.0;
  private static final double Aj = 2.0;
  private static final double[] r = {1.0, 1.0, 1.0};

  @Test
  public void tensorTest() {

    GeneralizedKirkwoodTensor generalizedKirkwoodTensor =
        new GeneralizedKirkwoodTensor(0, 3, gc, Eh, Es);
    generalizedKirkwoodTensor.setR(r[0], r[1], r[2], Ai, Aj);
    generalizedKirkwoodTensor.source();

    // Test the "bn" method.
    double bn0 = generalizedKirkwoodTensor.bn(0);
    double bn1 = generalizedKirkwoodTensor.bn(1);
    double bn2 = generalizedKirkwoodTensor.bn(2);
    double bn3 = generalizedKirkwoodTensor.bn(3);
    double bn4 = generalizedKirkwoodTensor.bn(4);
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
    double A00 = generalizedKirkwoodTensor.anm(0, 0);
    double A01 = generalizedKirkwoodTensor.anm(0, 1);
    double A02 = generalizedKirkwoodTensor.anm(0, 2);
    double A03 = generalizedKirkwoodTensor.anm(0, 3);
    final double mapleA00 = -0.4048255705;
    final double mapleA01 = 0.04764329318;
    final double mapleA02 = -0.01266056821;
    final double mapleA03 = 0.004644003643;
    assertEquals(" A00", mapleA00, A00, tolerance);
    assertEquals(" A01", mapleA01, A01, tolerance);
    assertEquals(" A02", mapleA02, A02, tolerance);
    assertEquals(" A03", mapleA03, A03, tolerance);
    double B00 = generalizedKirkwoodTensor.bnm(0, 0);
    double B01 = generalizedKirkwoodTensor.bnm(0, 1);
    double B02 = generalizedKirkwoodTensor.bnm(0, 2);
    // Monpole potential Born chain-rule derivatives.
    final double mapleB00 = 0.03273696138;
    final double mapleB01 = -0.01311852185;
    final double mapleB02 = 0.006171353772;
    assertEquals(" B00", mapleB00, B00, tolerance);
    assertEquals(" B01", mapleB01, B01, tolerance);
    assertEquals(" B02", mapleB02, B02, tolerance);

    // Dipole potential and derivatives.
    generalizedKirkwoodTensor = new GeneralizedKirkwoodTensor(1, 4, gc, Eh, Es);
    generalizedKirkwoodTensor.setR(r[0], r[1], r[2], Ai, Aj);
    generalizedKirkwoodTensor.source();
    double A10 = generalizedKirkwoodTensor.anm(1, 0);
    double A11 = generalizedKirkwoodTensor.anm(1, 1);
    double A12 = generalizedKirkwoodTensor.anm(1, 2);
    double A13 = generalizedKirkwoodTensor.anm(1, 3);
    final double mapleA10 = 0.06764004545;
    final double mapleA11 = -0.02388135594;
    final double mapleA12 = 0.01196727051;
    final double mapleA13 = -0.007470574812;
    assertEquals(" A10", mapleA10, A10, tolerance);
    assertEquals(" A11", mapleA11, A11, tolerance);
    assertEquals(" A12", mapleA12, A12, tolerance);
    assertEquals(" A13", mapleA13, A13, tolerance);
    // Dipole potential Born chain-rule derivatives.
    double B10 = generalizedKirkwoodTensor.bnm(1, 0);
    double B11 = generalizedKirkwoodTensor.bnm(1, 1);
    double B12 = generalizedKirkwoodTensor.bnm(1, 2);
    final double mapleB10 = -0.01640950855;
    final double mapleB11 = 0.01043812102;
    final double mapleB12 = -0.007669896149;
    assertEquals(" B10", mapleB10, B10, tolerance);
    assertEquals(" B11", mapleB11, B11, tolerance);
    assertEquals(" B12", mapleB12, B12, tolerance);

    // Quadrupole potential and derivatives.
    generalizedKirkwoodTensor = new GeneralizedKirkwoodTensor(2, 5, gc, Eh, Es);
    generalizedKirkwoodTensor.setR(r[0], r[1], r[2], Ai, Aj);
    generalizedKirkwoodTensor.source();
    double A20 = generalizedKirkwoodTensor.anm(2, 0);
    double A21 = generalizedKirkwoodTensor.anm(2, 1);
    double A22 = generalizedKirkwoodTensor.anm(2, 2);
    double A23 = generalizedKirkwoodTensor.anm(2, 3);
    final double mapleA20 = -0.03404928262;
    final double mapleA21 = 0.02003603617;
    final double mapleA22 = -0.01475634878;
    final double mapleA23 = 0.01280244361;
    assertEquals(" A20", mapleA20, A20, tolerance);
    assertEquals(" A21", mapleA21, A21, tolerance);
    assertEquals(" A22", mapleA22, A22, tolerance);
    assertEquals(" A23", mapleA23, A23, tolerance);
    // Quadrupole potential Born chain-rule derivatives.
    double B20 = generalizedKirkwoodTensor.bnm(2, 0);
    double B21 = generalizedKirkwoodTensor.bnm(2, 1);
    double B22 = generalizedKirkwoodTensor.bnm(2, 2);
    final double mapleB20 = 0.01376728808;
    final double mapleB21 = -0.01199790086;
    final double mapleB22 = 0.01179997601;
    assertEquals(" B20", mapleB20, B20, tolerance);
    assertEquals(" B21", mapleB21, B21, tolerance);
    assertEquals(" B22", mapleB22, B22, tolerance);
  }

}
