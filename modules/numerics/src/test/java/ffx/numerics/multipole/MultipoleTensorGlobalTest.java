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

import ffx.utilities.FFXTest;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import java.util.Arrays;
import java.util.Collection;
import java.util.logging.Logger;

import static ffx.numerics.multipole.MultipoleTensorTest.Qi;
import static ffx.numerics.multipole.MultipoleTensorTest.Qk;
import static ffx.numerics.multipole.MultipoleTensorTest.Ui;
import static ffx.numerics.multipole.MultipoleTensorTest.UiEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.Uk;
import static ffx.numerics.multipole.MultipoleTensorTest.UkEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueI;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueIEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueK;
import static ffx.numerics.multipole.MultipoleTensorTest.permTorqueKEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.permanentEnergy;
import static ffx.numerics.multipole.MultipoleTensorTest.permanentEnergyEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarGradICoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarGradIEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarGradIThole;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueICoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueIEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueIThole;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueKCoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueKEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarTorqueKThole;
import static ffx.numerics.multipole.MultipoleTensorTest.polarizationEnergyCoulomb;
import static ffx.numerics.multipole.MultipoleTensorTest.polarizationEnergyEwald;
import static ffx.numerics.multipole.MultipoleTensorTest.polarizationEnergyThole;
import static ffx.numerics.multipole.MultipoleTensorTest.scaleMutual;
import static ffx.numerics.multipole.MultipoleTensorTest.xyz;
import static java.lang.String.format;
import static org.junit.Assert.assertEquals;

/**
 * Parameterized Test of the MultipoleTensor class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class MultipoleTensorGlobalTest extends FFXTest {

  /**
   * Logger for the MultipoleTensor class.
   */
  private static final Logger logger = Logger.getLogger(MultipoleTensorGlobalTest.class.getName());

  private final double tolerance = 1.0e-13;
  private final double fdTolerance = 1.0e-6;
  private final double[] r = new double[3];
  private final double[] tensor;
  private final double[] noStorageTensor;
  private final double[] fastTensor;
  private final int order;
  private final int tensorCount;
  private final String info;
  private final Operator operator;
  private final double beta = 0.545;
  private final double thole = 0.39;
  private final double AiAk = 1.061104559485911;

  public MultipoleTensorGlobalTest(
      String info,
      int order,
      Operator operator) {
    this.info = info;
    this.order = order;
    r[0] = xyz[0];
    r[1] = xyz[1];
    r[2] = xyz[2];
    this.tensorCount = MultipoleUtilities.tensorCount(order);
    tensor = new double[tensorCount];
    noStorageTensor = new double[tensorCount];
    fastTensor = new double[tensorCount];
    this.operator = operator;
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][]{
            {
                "Order 5 Coulomb",
                5,      // Order
                Operator.COULOMB
            },
            {
                "Order 5 Screened Coulomb",
                5,
                Operator.SCREENED_COULOMB
            },
            {
                "Order 4 Thole Field",
                4,
                Operator.THOLE_FIELD
            }
        });
  }

  @Test
  public void codingGenerationTest() {
    // Only run code generation for the Coulomb operator.
    if (operator != Operator.COULOMB) {
      return;
    }

    // For code generation, make sure no mutipole terms are zero.
    double[] Qi = {0.11, 0.22, 0.33, 0.44, 0.11, 0.22, -0.33, 0.12, 0.13, 0.14};
    double[] r = {2.11, 2.12, 2.13};
    double[] tensor = new double[tensorCount];
    CoulombTensorGlobal multipoleTensor = new CoulombTensorGlobal(order);
    logger.info(format(" Writing Cartesian Order %d tensor recursion code:", order));
    String code = multipoleTensor.codeTensorRecursion(r, tensor);
    logger.info(format("\n%s", code));

    int filter = 5;
    logger.info(format(" Writing Cartesian Order %d vectorized tensor recursion code filtered at O=%d:", order, filter));
    code = multipoleTensor.codeVectorTensorRecursion(filter);
    logger.info(format("\n%s", code));

    PolarizableMultipole polarizableMultipole = new PolarizableMultipole(Qi, Ui, Ui);
    StringBuilder sb = new StringBuilder();
    logger.info(" Writing potential code due to multipole I:");
    multipoleTensor.codePotentialMultipoleI(polarizableMultipole, tensor, 0, 0, 0, sb);
    logger.info("\n\n" + sb);

    sb = new StringBuilder();
    logger.info(" Writing potential vectorized code due to multipole I:");
    multipoleTensor.codePotentialMultipoleISIMD(polarizableMultipole, tensor, 0, 0, 0, sb);
    logger.info("\n\n" + sb);

    sb = new StringBuilder();
    logger.info(" Writing potential code due to multipole K:");
    multipoleTensor.codePotentialMultipoleK(polarizableMultipole, tensor, 0, 0, 0, sb);
    logger.info("\n\n" + sb);

    sb = new StringBuilder();
    logger.info(" Writing potential vectorized code due to multipole K:");
    multipoleTensor.codePotentialMultipoleKSIMD(polarizableMultipole, tensor, 0, 0, 0, sb);
    logger.info("\n\n" + sb);
  }

  @Test
  public void permanentMultipoleEnergyAndGradTest() {
    // Thole damping is not used for permanent AMOEBA electrostatics.
    if (operator == Operator.THOLE_FIELD) {
      return;
    }

    MultipoleTensor multipoleTensor = new CoulombTensorGlobal(order);
    if (operator == Operator.SCREENED_COULOMB) {
      multipoleTensor = new EwaldTensorGlobal(order, beta);
    }

    double delta = 1.0e-5;
    double delta2 = 2.0 * 1.0e-5;
    double[] Fi = new double[3];
    double[] Fk = new double[3];
    double[] Ti = new double[3];
    double[] Tk = new double[3];

    PolarizableMultipole mI = new PolarizableMultipole(Qi, Ui, Ui);
    PolarizableMultipole mK = new PolarizableMultipole(Qk, Uk, Uk);
    multipoleTensor.generateTensor(r);
    double e = multipoleTensor.multipoleEnergyAndGradient(mI, mK, Fi, Fk, Ti, Tk);

    if (operator == Operator.COULOMB) {
      assertEquals(info + " Cartesian Multipole Energy", permanentEnergy, e, tolerance);
      assertEquals(info + " Cartesian Multipole Torque I", permTorqueI[1], Ti[1], tolerance);
      assertEquals(info + " Cartesian Multipole Torque K", permTorqueK[1], Tk[1], tolerance);
    } else if (operator == Operator.SCREENED_COULOMB) {
      assertEquals(info + " Cartesian Ewald Permanent Energy", permanentEnergyEwald, e, tolerance);
      assertEquals(info + " Cartesian Ewald Multipole Torque I", permTorqueIEwald[1], Ti[1],
          tolerance);
      assertEquals(info + " Cartesian Ewald Multipole Torque K", permTorqueKEwald[1], Tk[1],
          tolerance);
    }

    // Analytic gradient.
    double aX = Fk[0];
    double aY = Fk[1];
    double aZ = Fk[2];

    r[0] += delta;
    multipoleTensor.generateTensor(r);
    double posX = multipoleTensor.multipoleEnergy(mI, mK);
    r[0] -= delta2;
    multipoleTensor.generateTensor(r);
    double negX = multipoleTensor.multipoleEnergy(mI, mK);
    r[0] += delta;
    r[1] += delta;
    multipoleTensor.generateTensor(r);
    double posY = multipoleTensor.multipoleEnergy(mI, mK);
    r[1] -= delta2;
    multipoleTensor.generateTensor(r);
    double negY = multipoleTensor.multipoleEnergy(mI, mK);
    r[1] += delta;
    r[2] += delta;
    multipoleTensor.generateTensor(r);
    double posZ = multipoleTensor.multipoleEnergy(mI, mK);
    r[2] -= delta2;
    multipoleTensor.generateTensor(r);
    double negZ = multipoleTensor.multipoleEnergy(mI, mK);
    r[2] += delta;

    double expect = aX;
    double actual = (posX - negX) / delta2;
    assertEquals(info + " Force X", expect, actual, fdTolerance);
    expect = aY;
    actual = (posY - negY) / delta2;
    assertEquals(info + " Force Y", expect, actual, fdTolerance);
    expect = aZ;
    actual = (posZ - negZ) / delta2;
    assertEquals(info + " Force Z", expect, actual, fdTolerance);
  }

  @Test
  public void polarizationEnergyAndGradTest() {

    MultipoleTensor multipoleTensor = new CoulombTensorGlobal(order);

    if (operator == Operator.THOLE_FIELD) {
      multipoleTensor = new TholeTensorGlobal(order, thole, AiAk);
    } else if (operator == Operator.SCREENED_COULOMB) {
      multipoleTensor = new EwaldTensorGlobal(order, beta);
    }

    double delta = 1.0e-5;
    double delta2 = 2.0 * 1.0e-5;
    double[] Gi = new double[3];
    double[] Ti = new double[3];
    double[] Tk = new double[3];

    PolarizableMultipole mI = new PolarizableMultipole(Qi, Ui, Ui);
    PolarizableMultipole mK = new PolarizableMultipole(Qk, Uk, Uk);
    if (operator == Operator.SCREENED_COULOMB) {
      mI.setInducedDipole(UiEwald, UiEwald);
      mK.setInducedDipole(UkEwald, UkEwald);
    }

    multipoleTensor.generateTensor(r);
    double e = multipoleTensor.polarizationEnergyAndGradient(mI, mK, 1.0, 1.0, scaleMutual, Gi, Ti,
        Tk);

    // Analytic gradient.
    double aX = Gi[0];
    double aY = Gi[1];
    double aZ = Gi[2];

    double energy = polarizationEnergyCoulomb;
    double[] gradI = polarGradICoulomb;
    double[] torqueI = polarTorqueICoulomb;
    double[] torqueK = polarTorqueKCoulomb;
    if (operator == Operator.THOLE_FIELD) {
      energy = polarizationEnergyThole;
      gradI = polarGradIThole;
      torqueI = polarTorqueIThole;
      torqueK = polarTorqueKThole;
    } else if (operator == Operator.SCREENED_COULOMB) {
      energy = polarizationEnergyEwald;
      gradI = polarGradIEwald;
      torqueI = polarTorqueIEwald;
      torqueK = polarTorqueKEwald;
    }
    assertEquals(info + " Polarization Energy", energy, e, tolerance);
    assertEquals(info + " Polarization GradX", gradI[0], aX, tolerance);
    assertEquals(info + " Polarization GradY", gradI[1], aY, tolerance);
    assertEquals(info + " Polarization GradZ", gradI[2], aZ, tolerance);
    assertEquals(info + " Polarization TorqueX I", torqueI[0], Ti[0], tolerance);
    assertEquals(info + " Polarization TorqueY I", torqueI[1], Ti[1], tolerance);
    assertEquals(info + " Polarization TorqueZ I", torqueI[2], Ti[2], tolerance);
    assertEquals(info + " Polarization TorqueX K", torqueK[0], Tk[0], tolerance);
    assertEquals(info + " Polarization TorqueY K", torqueK[1], Tk[1], tolerance);
    assertEquals(info + " Polarization TorqueZ K", torqueK[2], Tk[2], tolerance);

    // Finite-Difference test only works for direct polarization.
    if (scaleMutual == 0.0) {
      // Compute the gradient for atom K.
      r[0] += delta;
      multipoleTensor.generateTensor(r);
      double posX = multipoleTensor.polarizationEnergy(mI, mK, 1.0);
      r[0] -= delta2;
      multipoleTensor.generateTensor(r);
      double negX = multipoleTensor.polarizationEnergy(mI, mK, 1.0);
      r[0] += delta;

      r[1] += delta;
      multipoleTensor.generateTensor(r);
      double posY = multipoleTensor.polarizationEnergy(mI, mK, 1.0);
      r[1] -= delta2;
      multipoleTensor.generateTensor(r);
      double negY = multipoleTensor.polarizationEnergy(mI, mK, 1.0);
      r[1] += delta;
      r[2] += delta;
      multipoleTensor.generateTensor(r);
      double posZ = multipoleTensor.polarizationEnergy(mI, mK, 1.0);
      r[2] -= delta2;
      multipoleTensor.generateTensor(r);
      double negZ = multipoleTensor.polarizationEnergy(mI, mK, 1.0);
      r[2] += delta;

      // Note that the factor of 2 below accounts for the induced dipoles not being updated (i.e. the analytic gradient is for
      // induced dipoles, but the finite-difference is effectively for permanent dipoles due to lack of re-induction).
      double expect = -aX;
      double actual = 2.0 * (posX - negX) / delta2;
      assertEquals(info + " Grad X", expect, actual, fdTolerance);
      expect = -aY;
      actual = 2.0 * (posY - negY) / delta2;
      assertEquals(info + " Grad Y", expect, actual, fdTolerance);
      expect = -aZ;
      actual = 2.0 * (posZ - negZ) / delta2;
      assertEquals(info + " Grad Z", expect, actual, fdTolerance);
    }
  }

  @Test
  public void finiteDifferenceTest() {

    MultipoleTensor multipoleTensor = new CoulombTensorGlobal(order);
    if (operator == Operator.THOLE_FIELD) {
      multipoleTensor = new TholeTensorGlobal(order, thole, AiAk);
    } else if (operator == Operator.SCREENED_COULOMB) {
      multipoleTensor = new EwaldTensorGlobal(order, beta);
    }

    multipoleTensor.recursion(r, tensor);

    double[] tensorsPx = new double[tensorCount];
    double[] tensorsNx = new double[tensorCount];
    double[] tensorsPy = new double[tensorCount];
    double[] tensorsNy = new double[tensorCount];
    double[] tensorsPz = new double[tensorCount];
    double[] tensorsNz = new double[tensorCount];
    double delta = 1.0e-5;
    double delta2 = delta * 2;
    r[0] += delta;
    multipoleTensor.noStorageRecursion(r, tensorsPx);
    r[0] -= delta2;
    multipoleTensor.noStorageRecursion(r, tensorsNx);
    r[0] += delta;

    r[1] += delta;
    multipoleTensor.noStorageRecursion(r, tensorsPy);
    r[1] -= delta2;
    multipoleTensor.noStorageRecursion(r, tensorsNy);
    r[1] += delta;

    r[2] += delta;
    multipoleTensor.noStorageRecursion(r, tensorsPz);
    r[2] -= delta2;
    multipoleTensor.noStorageRecursion(r, tensorsNz);
    r[2] += delta;

    tensorFiniteDifference(
        multipoleTensor, delta2, tensorsPx, tensorsNx, tensorsPy, tensorsNy, tensorsPz, tensorsNz);

    multipoleTensor.recursion(r, tensor);
    tensorFiniteDifference(
        multipoleTensor, delta2, tensorsPx, tensorsNx, tensorsPy, tensorsNy, tensorsPz, tensorsNz);

    multipoleTensor.generateTensor();
    multipoleTensor.getTensor(tensor);
    tensorFiniteDifference(
        multipoleTensor, delta2, tensorsPx, tensorsNx, tensorsPy, tensorsNy, tensorsPz, tensorsNz);
  }

  /**
   * Test of recursion and noStorageRecursion methods, of class MultipoleTensor.
   *
   * @since 1.0
   */
  @Test
  public void multipoleTensorTest() {
    MultipoleTensor multipoleTensor = new CoulombTensorGlobal(order);
    if (operator == Operator.THOLE_FIELD) {
      multipoleTensor = new TholeTensorGlobal(order, thole, AiAk);
    } else if (operator == Operator.SCREENED_COULOMB) {
      multipoleTensor = new EwaldTensorGlobal(order, beta);
    }

    // Check Cartesian Tensors in the Global frame.
    multipoleTensor.noStorageRecursion(r, noStorageTensor);
    multipoleTensor.recursion(r, tensor);
    multipoleTensor.generateTensor(r);
    multipoleTensor.getTensor(fastTensor);

    for (int i = 0; i < tensorCount; i++) {
      double expect = noStorageTensor[i];
      double actual = tensor[i];
      assertEquals(info + " @ " + i, expect, actual, tolerance);
      if (order >= 2 || order <= 6) {
        expect = noStorageTensor[i];
        actual = fastTensor[i];
        assertEquals(info + " @ " + i, expect, actual, tolerance);
      }
    }
  }

  private void tensorFiniteDifference(MultipoleTensor multipoleTensor,
                                      double delta2,
                                      double[] tensorsPx, double[] tensorsNx,
                                      double[] tensorsPy, double[] tensorsNy,
                                      double[] tensorsPz, double[] tensorsNz) {

    int start = 0;

    // We do not calculate the zeroth term for Thole damping.
    if (operator == Operator.THOLE_FIELD) {
      start = 1;
    }

    // Test the partial derivatives for all tensor components.
    for (int l = start; l < order; l++) {
      // Test X derivative
      double expect = tensor[multipoleTensor.ti(l + 1, 0, 0)];
      double actual =
          (tensorsPx[multipoleTensor.ti(l, 0, 0)] - tensorsNx[multipoleTensor.ti(l, 0, 0)]) / delta2;
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
        double actual =
            (tensorsPx[multipoleTensor.ti(l, m, 0)] - tensorsNx[multipoleTensor.ti(l, m, 0)])
                / delta2;
        assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
        // Test Y derivative
        expect = tensor[multipoleTensor.ti(l, m + 1, 0)];
        actual =
            (tensorsPy[multipoleTensor.ti(l, m, 0)] - tensorsNy[multipoleTensor.ti(l, m, 0)])
                / delta2;
        assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
        // Test Z derivative
        expect = tensor[multipoleTensor.ti(l, m, 1)];
        actual =
            (tensorsPz[multipoleTensor.ti(l, m, 0)] - tensorsNz[multipoleTensor.ti(l, m, 0)])
                / delta2;
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
