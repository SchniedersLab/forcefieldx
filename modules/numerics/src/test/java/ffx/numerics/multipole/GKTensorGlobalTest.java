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

import static ffx.numerics.multipole.GKTensorGlobal.GK_MULTIPOLE_ORDER.DIPOLE;
import static ffx.numerics.multipole.GKTensorGlobal.GK_MULTIPOLE_ORDER.MONOPOLE;
import static ffx.numerics.multipole.GKTensorGlobal.GK_MULTIPOLE_ORDER.QUADRUPOLE;
import static java.lang.String.format;
import static org.junit.Assert.assertEquals;

import ffx.numerics.multipole.GKTensorGlobal.GK_TENSOR_MODE;
import org.junit.Test;

/**
 * Test the GK tensor evaluated in the global coordinate frame.
 * <p>
 * There is no quadrupole finite-difference test because the trace of the quadrupole potential is
 * neglected; thus the derivatives of the quadrupole potential are correct when summed over the
 * trace, but not on a per-element basis.
 */
public class GKTensorGlobalTest {

  private final double tolerance = 1.0e-9;
  private final double fdTolerance = 1.0e-6;
  private static final double gc = 2.455;
  private static final double Eh = 1.0;
  private static final double Es = 78.3;
  private static final double Ai = 2.0;
  private static final double Aj = 2.0;
  private static final double[] r = {0.7, 0.8, 0.9};

  public static final double[] rWater = {2.973385290000000, 0.000000000000000, 0.035464520000000};
  public static final double bornI = 1.403320706062232;
  public static final double bornK = 1.403320706062232;

  public static final double watPermEnergy = -0.086259327778797;
  public static final double[] multI = {-0.51966, 0.06979198988239577, 0.0, 0.0289581620819011,
      0.024871041044109393, -0.1170771231287098, 0.09220608208460039,
      0.0, -0.03374891685535346, 0.0};
  public static final double[] multK = {-0.51966, 0.05872406108747119, 0.0, 0.047549780780788455,
      0.048623413357888695, -0.1170771231287098, 0.06845370977082109,
      0.0, -0.04662811558421081, 0.0};
  public static final double[] permTorqueI = {0.0, -0.000497809865297, 0.0};
  public static final double[] permTorqueK = {0.0, 0.003535773811802, 0.0};
  public static final double[] permGradI = {-0.025569107173925, 0.0, 0.000716747957596};
  public static final double permGradBorn = 0.002598186313550;

  // Mutual polarization.
  public static final double watTotalEnergy = -0.086278549962107;
  public static final double watPolEnergy = watTotalEnergy - watPermEnergy;
  public static final double[] uI = {0.07723465799892439, 0.0, 0.027914344367154516};
  public static final double[] uK = {0.08743555236690667, 0.0, 0.05830830340255934};
  public static final double[] polGradI = {0.000254342542095, 0.0, -0.000047781304691};
  public static final double[] polTorqueI = {0.0, -0.000050312084142, 0.0};
  public static final double[] polTorqueK = {0.0, -0.000480106909619, 0.0};
  public static final double polGradBorn = -0.000262931043709;

  // Direct polarization.
  public static final double watTotalDirect = -0.086133659894249;
  public static final double watPolDirect = watTotalDirect - watPermEnergy;
  public static final double[] uIDirect = {0.05851012079839759, 0.0, 0.02312407977363764};
  public static final double[] uKDirect = {0.0603673253364741, 0.0, 0.04226502224741869};
  public static final double[] polGradIDirect = {0.000161932495081, 0.0, 0.000062539680741};
  public static final double[] polTorqueIDirect = {0.0, -0.000036645740123, 0.0};
  public static final double[] polTorqueKDirect = {0.0, -0.000368952557395, 0.0};
  public static final double polGradBornDirect = -0.000125792583010;

  @Test
  public void permanentEnergyTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);

    // Apply the GK monopole tensor.
    int order = 2; // Charge - Quadrupole Energy.
    GKTensorGlobal gkMonopoleTensor = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.multipoleEnergy(mI, mK);

    // Apply the GK dipole tensor.
    order = 3; // Dipole - Quadrupole Energy.
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.multipoleEnergy(mI, mK);

    // Apply the GK quadrupole tensor.
    order = 4; // Quadrupole - Quadrupole Energy.
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.multipoleEnergy(mI, mK);

    assertEquals("GK Permanent Energy", watPermEnergy, c + d + q, tolerance);
  }

  @Test
  public void permanentEnergyAndGradientTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);

    double[] gradI = new double[3];
    double[] gradK = new double[3];
    double[] torqueI = new double[3];
    double[] torqueK = new double[3];

    double[] tempGradI = new double[3];
    double[] tempGradK = new double[3];
    double[] tempTorqueI = new double[3];
    double[] tempTorqueK = new double[3];

    // Apply the GK monopole tensor.
    int order = 3; // Charge - Quadrupole Energy Gradient.
    GKTensorGlobal gkMonopoleTensor = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double e = gkMonopoleTensor.multipoleEnergyAndGradient(mI, mK, gradI, gradK, torqueI, torqueK);
    double db = gkMonopoleTensor.multipoleEnergyBornGrad(mI, mK);

    // Apply the GK dipole tensor.
    order = 4; // Dipole - Quadrupole Energy Gradient.
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    e += gkDipoleTensor.multipoleEnergyAndGradient(mI, mK, tempGradI, tempGradK, tempTorqueI, tempTorqueK);
    db += gkDipoleTensor.multipoleEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      gradK[i] += tempGradK[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
      tempTorqueI[i] = 0.0;
      tempTorqueK[i] = 0.0;
    }

    // Apply the GK quadrupole tensor.
    order = 5; // Quadrupole - Quadrupole Energy Gradient.
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    e += gkQuadrupoleTensor.multipoleEnergyAndGradient(mI, mK, tempGradI, tempGradK, tempTorqueI, tempTorqueK);
    db += gkQuadrupoleTensor.multipoleEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      gradK[i] += tempGradK[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
    }

    assertEquals("GK Permanent Energy", watPermEnergy, e, tolerance);
    assertEquals("GK Permanent Grad X", permGradI[0], gradI[0], tolerance);
    assertEquals("GK Permanent Grad Y", permGradI[1], gradI[1], tolerance);
    assertEquals("GK Permanent Grad Z", permGradI[2], gradI[2], tolerance);
    assertEquals("GK Permanent Torque I X", permTorqueI[0], torqueI[0], tolerance);
    assertEquals("GK Permanent Torque I Y", permTorqueI[1], torqueI[1], tolerance);
    assertEquals("GK Permanent Torque I Z", permTorqueI[2], torqueI[2], tolerance);
    assertEquals("GK Permanent Torque K X", permTorqueK[0], torqueK[0], tolerance);
    assertEquals("GK Permanent Torque K Y", permTorqueK[1], torqueK[1], tolerance);
    assertEquals("GK Permanent Torque K Z", permTorqueK[2], torqueK[2], tolerance);
    assertEquals("GK Born Grad I", permGradBorn, db * bornI, tolerance);
  }

  @Test
  public void polarizationEnergyTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);

    // Apply the GK monopole tensor.
    int order = 1; // Induced Dipole - Monopole Energy.
    GKTensorGlobal gkMonopoleTensor = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK dipole tensor.
    order = 3; // Induced Dipole - Quadrupole Energy.
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK quadrupole tensor.
    order = 3; // Induced Dipole - Quadrupole Energy.
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.polarizationEnergy(mI, mK);

    assertEquals("GK Polarization Energy", watPolEnergy, c + d + q, tolerance);
  }

  @Test
  public void polarizationEnergyDirectTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uIDirect, uIDirect);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uKDirect, uKDirect);

    // Apply the GK monopole tensor.
    int order = 1; // Induced Dipole - Monopole Energy.
    GKTensorGlobal gkMonopoleTensor = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK dipole tensor.
    order = 3; // Induced Dipole - Quadrupole Energy.
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK quadrupole tensor.
    order = 3; // Induced Dipole - Quadrupole Energy.
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.polarizationEnergy(mI, mK);

    assertEquals("GK Direct Polarization Energy", watPolDirect, c + d + q, tolerance);
  }

  @Test
  public void polarizationEnergyAndGradientTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);

    mI.applyMasks(1.0, 1.0);
    mK.applyMasks(1.0, 1.0);

    double[] gradI = new double[3];
    double[] torqueI = new double[3];
    double[] torqueK = new double[3];
    double[] tempGradI = new double[3];
    double[] tempTorqueI = new double[3];
    double[] tempTorqueK = new double[3];

    // Apply the GK monopole tensor.
    int order = 2; // Induced Dipole - Monopole Energy Gradient.
    GKTensorGlobal gkMonopoleTensor = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 1.0, gradI, torqueI, torqueK);
    double db = gkMonopoleTensor.polarizationEnergyBornGrad(mI, mK);

    // Apply the GK dipole tensor.
    order = 4; // Induced Dipole - Quadrupole Energy Gradient.
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 1.0, tempGradI, tempTorqueI, tempTorqueK);
    db += gkDipoleTensor.polarizationEnergyBornGrad(mI, mK);
    // Add the mutual polarization contribution to Born chain-rule term.
    db += gkDipoleTensor.mutualPolarizationEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
      tempGradI[i] = 0.0;
      tempTorqueI[i] = 0.0;
      tempTorqueK[i] = 0.0;
    }

    // Apply the GK quadrupole tensor.
    order = 4; // Induced Dipole - Quadrupole Energy Gradient.
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 1.0, tempGradI, tempTorqueI, tempTorqueK);
    db += gkQuadrupoleTensor.polarizationEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
    }

    assertEquals("GK Polarization Energy", watPolEnergy, c + d + q, tolerance);
    assertEquals("GK Polarization Grad X", polGradI[0], gradI[0], tolerance);
    assertEquals("GK Polarization Grad Y", polGradI[1], gradI[1], tolerance);
    assertEquals("GK Polarization Grad Z", polGradI[2], gradI[2], tolerance);
    assertEquals("GK Polarization Torque I X", polTorqueI[0], torqueI[0], tolerance);
    assertEquals("GK Polarization Torque I Y", polTorqueI[1], torqueI[1], tolerance);
    assertEquals("GK Polarization Torque I Z", polTorqueI[2], torqueI[2], tolerance);
    assertEquals("GK Polarization Torque K X", polTorqueK[0], torqueK[0], tolerance);
    assertEquals("GK Polarization Torque K Y", polTorqueK[1], torqueK[1], tolerance);
    assertEquals("GK Polarization Torque K Z", polTorqueK[2], torqueK[2], tolerance);
    assertEquals("GK Born Grad I", polGradBorn, db * bornI, tolerance);
  }

  @Test
  public void polarizationEnergyAndGradientDirectTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uIDirect, uIDirect);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uKDirect, uKDirect);

    mI.applyMasks(1.0, 1.0);
    mK.applyMasks(1.0, 1.0);

    double[] gradI = new double[3];
    double[] torqueI = new double[3];
    double[] torqueK = new double[3];
    double[] tempGradI = new double[3];
    double[] tempTorqueI = new double[3];
    double[] tempTorqueK = new double[3];

    // Apply the GK monopole tensor.
    int order = 2; // Induced Dipole - Monopole Energy Gradient.
    GKTensorGlobal gkMonopoleTensor = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 0.0, gradI, torqueI, torqueK);
    double db = gkMonopoleTensor.polarizationEnergyBornGrad(mI, mK);

    // Apply the GK dipole tensor.
    order = 4; // Induced Dipole - Quadrupole Energy Gradient.
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 0.0, tempGradI, tempTorqueI, tempTorqueK);
    db += gkDipoleTensor.polarizationEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
      tempGradI[i] = 0.0;
      tempTorqueI[i] = 0.0;
      tempTorqueK[i] = 0.0;
    }

    // Apply the GK quadrupole tensor.
    order = 4; // Induced Dipole - Quadrupole Energy Gradient.
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 0.0, tempGradI, tempTorqueI, tempTorqueK);
    db += gkQuadrupoleTensor.polarizationEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
    }

    assertEquals("GK Polarization Energy", watPolDirect, c + d + q, tolerance);
    assertEquals("GK Polarization Grad X", polGradIDirect[0], gradI[0], tolerance);
    assertEquals("GK Polarization Grad Y", polGradIDirect[1], gradI[1], tolerance);
    assertEquals("GK Polarization Grad Z", polGradIDirect[2], gradI[2], tolerance);
    assertEquals("GK Polarization Torque I X", polTorqueIDirect[0], torqueI[0], tolerance);
    assertEquals("GK Polarization Torque I Y", polTorqueIDirect[1], torqueI[1], tolerance);
    assertEquals("GK Polarization Torque I Z", polTorqueIDirect[2], torqueI[2], tolerance);
    assertEquals("GK Polarization Torque K X", polTorqueKDirect[0], torqueK[0], tolerance);
    assertEquals("GK Polarization Torque K Y", polTorqueKDirect[1], torqueK[1], tolerance);
    assertEquals("GK Polarization Torque K Z", polTorqueKDirect[2], torqueK[2], tolerance);
    assertEquals("GK Born Grad I", polGradBornDirect,db * bornI, tolerance);
  }

  @Test
  public void tensorAuxiliaryTest() {
    int order = 3;
    double[] work = new double[order + 1];

    GKTensorGlobal gkTensorGlobal = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.source(work);

    // Test the "bn" method.
    double bn0 = gkTensorGlobal.bn(0);
    double bn1 = gkTensorGlobal.bn(1);
    double bn2 = gkTensorGlobal.bn(2);
    double bn3 = gkTensorGlobal.bn(3);
    double bn4 = gkTensorGlobal.bn(4);
    double mapleBN0 = 0.4914375691;
    double mapleBN1 = -0.01651130007;
    double mapleBN2 = -0.01365916860;
    double mapleBN3 = 0.006248702126;
    double mapleBN4 = -0.001978716128;
    assertEquals("bn0", mapleBN0, bn0, tolerance);
    assertEquals("bn1", mapleBN1, bn1, tolerance);
    assertEquals("bn2", mapleBN2, bn2, tolerance);
    assertEquals("bn3", mapleBN3, bn3, tolerance);
    assertEquals("bn4", mapleBN4, bn4, tolerance);

    // Monopole potential and derivatives.
    double A00 = gkTensorGlobal.anm(0, 0);
    double A01 = gkTensorGlobal.anm(0, 1);
    double A02 = gkTensorGlobal.anm(0, 2);
    double A03 = gkTensorGlobal.anm(0, 3);
    final double mapleA00 = -0.4319767286;
    final double mapleA01 = 0.05505753880;
    final double mapleA02 = -0.01542067105;
    final double mapleA03 = 0.005809287830;
    assertEquals("A00", mapleA00, A00, tolerance);
    assertEquals("A01", mapleA01, A01, tolerance);
    assertEquals("A02", mapleA02, A02, tolerance);
    assertEquals("A03", mapleA03, A03, tolerance);

    gkTensorGlobal.setMode(GK_TENSOR_MODE.BORN);
    gkTensorGlobal.source(work);
    double B00 = gkTensorGlobal.bnm(0, 0);
    double B01 = gkTensorGlobal.bnm(0, 1);
    double B02 = gkTensorGlobal.bnm(0, 2);
    // Monpole potential Born chain-rule derivatives.
    final double mapleB00 = 0.04064563792;
    final double mapleB01 = -0.01690706430;
    final double mapleB02 = 0.008229167114;
    assertEquals("B00", mapleB00, B00, tolerance);
    assertEquals("B01", mapleB01, B01, tolerance);
    assertEquals("B02", mapleB02, B02, tolerance);

    // Dipole potential and derivatives.
    order = 4;
    work = new double[order + 1];
    gkTensorGlobal = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.source(work);
    double A10 = gkTensorGlobal.anm(1, 0);
    double A11 = gkTensorGlobal.anm(1, 1);
    double A12 = gkTensorGlobal.anm(1, 2);
    double A13 = gkTensorGlobal.anm(1, 3);
    final double mapleA10 = 0.08218283800;
    final double mapleA11 = -0.03142380936;
    final double mapleA12 = 0.01681150454;
    final double mapleA13 = -0.01106715285;
    assertEquals("A10", mapleA10, A10, tolerance);
    assertEquals("A11", mapleA11, A11, tolerance);
    assertEquals("A12", mapleA12, A12, tolerance);
    assertEquals("A13", mapleA13, A13, tolerance);
    // Dipole potential Born chain-rule derivatives.
    gkTensorGlobal.setMode(GK_TENSOR_MODE.BORN);
    gkTensorGlobal.source(work);
    double B10 = gkTensorGlobal.bnm(1, 0);
    double B11 = gkTensorGlobal.bnm(1, 1);
    double B12 = gkTensorGlobal.bnm(1, 2);
    final double mapleB10 = -0.02319829046;
    final double mapleB11 = 0.01556309100;
    final double mapleB12 = -0.01202628190;
    assertEquals("B10", mapleB10, B10, tolerance);
    assertEquals("B11", mapleB11, B11, tolerance);
    assertEquals("B12", mapleB12, B12, tolerance);

    // Quadrupole potential and derivatives.
    order = 5;
    work = new double[order + 1];
    gkTensorGlobal = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkTensorGlobal.setR(r, Ai, Aj);
    gkTensorGlobal.source(work);
    double A20 = gkTensorGlobal.anm(2, 0);
    double A21 = gkTensorGlobal.anm(2, 1);
    double A22 = gkTensorGlobal.anm(2, 2);
    double A23 = gkTensorGlobal.anm(2, 3);
    final double mapleA20 = -0.04710532877;
    final double mapleA21 = 0.03001901831;
    final double mapleA22 = -0.02371209206;
    final double mapleA23 = 0.02187861148;
    assertEquals("A20", mapleA20, A20, tolerance);
    assertEquals("A21", mapleA21, A21, tolerance);
    assertEquals("A22", mapleA22, A22, tolerance);
    assertEquals("A23", mapleA23, A23, tolerance);

    // Quadrupole potential Born chain-rule derivatives.
    gkTensorGlobal.setMode(GK_TENSOR_MODE.BORN);
    gkTensorGlobal.source(work);
    double B20 = gkTensorGlobal.bnm(2, 0);
    double B21 = gkTensorGlobal.bnm(2, 1);
    double B22 = gkTensorGlobal.bnm(2, 2);
    final double mapleB20 = 0.02216121853;
    final double mapleB21 = -0.02051645869;
    final double mapleB22 = 0.02137054028;
    assertEquals("B20", mapleB20, B20, tolerance);
    assertEquals("B21", mapleB21, B21, tolerance);
    assertEquals("B22", mapleB22, B22, tolerance);
  }

  @Test
  public void chargeTensorTest() {
    int order = 3;
    double[] work = new double[order + 1];

    GKTensorGlobal gkMonopoleTensor =
        new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(r, Ai, Aj);
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

    // Check Born radii chain rule terms.
    gkMonopoleTensor.setMode(GK_TENSOR_MODE.BORN);
    gkMonopoleTensor.generateTensor();
    double B00 = gkMonopoleTensor.bnm(0, 0);
    double B01 = gkMonopoleTensor.bnm(0, 1);
    double B02 = gkMonopoleTensor.bnm(0, 2);
    assertEquals(" B000", B00, gkMonopoleTensor.R000, tolerance);
    assertEquals(" B100", x * B01, gkMonopoleTensor.R100, tolerance);
    assertEquals(" B200", x * x * B02 + B01, gkMonopoleTensor.R200, tolerance);
  }

  @Test
  public void dipoleTensorTest() {
    double x = r[0];

    // Dipole potential and derivatives.
    int order = 4;
    double[] work = new double[order + 1];
    GKTensorGlobal gkDipoleTensor = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(r, Ai, Aj);
    gkDipoleTensor.source(work);
    double A10 = gkDipoleTensor.anm(1, 0);
    double A11 = gkDipoleTensor.anm(1, 1);
    double A12 = gkDipoleTensor.anm(1, 2);
    double A13 = gkDipoleTensor.anm(1, 3);
    gkDipoleTensor.generateTensor();

    // No charge potential.
    assertEquals(" R000", 0.0, gkDipoleTensor.R000, tolerance);
    // Ux Dipole potential
    assertEquals(" R100", x * A10, gkDipoleTensor.R100, tolerance);
    // Ux Dipole potential gradient
    assertEquals(" R200", x * x * A11 + A10, gkDipoleTensor.R200, tolerance);
    // Ux Dipole 2nd potential gradient
    assertEquals(" R300", x * x * x * A12 + 3.0 * x * A11, gkDipoleTensor.R300, tolerance);
    // Ux Dipole 3rd potential gradient
    assertEquals(" R400", x * x * x * x * A13 + 6.0 * x * x * A12 + 3.0 * A11, gkDipoleTensor.R400,
        tolerance);

    // Check Born radii chain rule terms.
    gkDipoleTensor.setMode(GK_TENSOR_MODE.BORN);
    gkDipoleTensor.generateTensor();
    double B10 = gkDipoleTensor.bnm(1, 0);
    double B11 = gkDipoleTensor.bnm(1, 1);
    double B12 = gkDipoleTensor.bnm(1, 2);
    assertEquals(" B100", x * B10, gkDipoleTensor.R100, tolerance);
    assertEquals(" B200", x * x * B11 + B10, gkDipoleTensor.R200, tolerance);
    assertEquals(" B300", x * x * x * B12 + 3.0 * x * B11, gkDipoleTensor.R300, tolerance);
  }

  @Test
  public void quadrupoleTensorTest() {
    double x = r[0];
    double y = r[1];
    double z = r[2];

    // Quadrupole potential and derivatives.
    int order = 5;
    double[] work = new double[order + 1];
    GKTensorGlobal gkQuadrupoleTensor = new GKTensorGlobal(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(r, Ai, Aj);
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

    // Potential for the quadrupole trace
    assertEquals(" R200", x * x * A20, gkQuadrupoleTensor.R200, tolerance);
    assertEquals(" R020", y * y * A20, gkQuadrupoleTensor.R020, tolerance);
    assertEquals(" R002", z * z * A20, gkQuadrupoleTensor.R002, tolerance);

    // Partial derivative of the trace with respect to X (note that the x*A20 term sums to 2*x*A20 over the trace).
    assertEquals(" R300", x * x * x * A21 + 3.0 * x * A20, gkQuadrupoleTensor.R300, tolerance);
    assertEquals(" R120", x * y * y * A21 + x * A20, gkQuadrupoleTensor.R120, tolerance);
    assertEquals(" R102", x * z * z * A21 + x * A20, gkQuadrupoleTensor.R102, tolerance);

    // Higher order Qxx terms.
    assertEquals(" R400", x * x * x * x * A22 + 6.0 * x * x * A21 + 3.0 * A20,
        gkQuadrupoleTensor.R400, tolerance);
    assertEquals(" R500", x * x * x * x * x * A23 + 10.0 * x * x * x * A22 + 15.0 * x * A21,
        gkQuadrupoleTensor.R500, tolerance);

    // Qxy
    assertEquals(" R110", x * y * A20, gkQuadrupoleTensor.R110, tolerance);
    assertEquals(" R111", x * y * z * A21, gkQuadrupoleTensor.R111, tolerance);
    assertEquals(" R310", x * x * x * y * A22 + 3.0 * x * y * A21, gkQuadrupoleTensor.R310,
        tolerance);
    assertEquals(" R220", x * x * y * y * A22 + x * x * A21 + y * y * A21 + A20,
        gkQuadrupoleTensor.R220, tolerance);
    assertEquals(" R211", x * x * y * z * A22 + y * z * A21, gkQuadrupoleTensor.R211, tolerance);
    assertEquals(" R410", x * x * x * x * y * A23 + 6.0 * x * x * y * A22 + 3.0 * y * A21,
        gkQuadrupoleTensor.R410, tolerance);
    assertEquals(" R320",
        x * x * x * y * y * A23 + x * x * x * A22 + 3.0 * x * y * y * A22 + 3.0 * x * A21,
        gkQuadrupoleTensor.R320, tolerance);
    assertEquals(" R311", x * x * x * y * z * A23 + 3.0 * x * y * z * A22, gkQuadrupoleTensor.R311,
        tolerance);

    // Check Born radii chain rule terms.
    gkQuadrupoleTensor.setMode(GK_TENSOR_MODE.BORN);
    gkQuadrupoleTensor.generateTensor();
    double B20 = gkQuadrupoleTensor.bnm(2, 0);
    double B21 = gkQuadrupoleTensor.bnm(2, 1);
    double B22 = gkQuadrupoleTensor.bnm(2, 2);
    // No charge or dipole potential.
    assertEquals(" B000", 0.0, gkQuadrupoleTensor.R000, tolerance);
    assertEquals(" B100", 0.0, gkQuadrupoleTensor.R100, tolerance);
    // Born chain rule for the potential for the quadrupole trace
    assertEquals(" B200", x * x * B20, gkQuadrupoleTensor.R200, tolerance);
    assertEquals(" B020", y * y * B20, gkQuadrupoleTensor.R020, tolerance);
    assertEquals(" B002", z * z * B20, gkQuadrupoleTensor.R002, tolerance);
    // Born chain rule for the partial derivative of the trace with respect to X.
    assertEquals(" B300", x * x * x * B21 + 3.0 * x * B20, gkQuadrupoleTensor.R300, tolerance);
    assertEquals(" B120", x * y * y * B21 + x * B20, gkQuadrupoleTensor.R120, tolerance);
    assertEquals(" B102", x * z * z * B21 + x * B20, gkQuadrupoleTensor.R102, tolerance);
    // Higher order term for Qxx.
    assertEquals(" B400", x * x * x * x * B22 + 6.0 * x * x * B21 + 3.0 * B20,
        gkQuadrupoleTensor.R400, tolerance);
  }

  @Test
  public void chargeFiniteDifferenceTest() {
    int order = 6;

    GKTensorGlobal gkTensorGlobal = new GKTensorGlobal(MONOPOLE, order, gc, Eh, Es);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
    gkTensorGlobal.setR(r, Ai, Aj);
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
    GKTensorGlobal gkTensorGlobal = new GKTensorGlobal(DIPOLE, order, gc, Eh, Es);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
    gkTensorGlobal.setR(r, Ai, Aj);
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

    int start = multipoleTensor.multipoleOrder.getOrder();

    String info = "GK Global";
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
        actual = (tensorsPy[multipoleTensor.ti(l, m, 0)] - tensorsNy[multipoleTensor.ti(l, m, 0)])
            / delta2;
        assertEquals(info + " " + l + " " + m, expect, actual, fdTolerance);
        // Test Z derivative
        expect = tensor[multipoleTensor.ti(l, m, 1)];
        actual = (tensorsPz[multipoleTensor.ti(l, m, 0)] - tensorsNz[multipoleTensor.ti(l, m, 0)])
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
