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
import static ffx.numerics.multipole.GKTensorGlobalTest.bornI;
import static ffx.numerics.multipole.GKTensorGlobalTest.bornK;
import static ffx.numerics.multipole.GKTensorGlobalTest.multI;
import static ffx.numerics.multipole.GKTensorGlobalTest.multK;
import static ffx.numerics.multipole.GKTensorGlobalTest.permGradBorn;
import static ffx.numerics.multipole.GKTensorGlobalTest.permGradI;
import static ffx.numerics.multipole.GKTensorGlobalTest.permTorqueI;
import static ffx.numerics.multipole.GKTensorGlobalTest.permTorqueK;
import static ffx.numerics.multipole.GKTensorGlobalTest.polGradBorn;
import static ffx.numerics.multipole.GKTensorGlobalTest.polGradBornDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.polGradI;
import static ffx.numerics.multipole.GKTensorGlobalTest.polGradIDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.polTorqueI;
import static ffx.numerics.multipole.GKTensorGlobalTest.polTorqueIDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.polTorqueK;
import static ffx.numerics.multipole.GKTensorGlobalTest.polTorqueKDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.rWater;
import static ffx.numerics.multipole.GKTensorGlobalTest.uI;
import static ffx.numerics.multipole.GKTensorGlobalTest.uIDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.uK;
import static ffx.numerics.multipole.GKTensorGlobalTest.uKDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.watPermEnergy;
import static ffx.numerics.multipole.GKTensorGlobalTest.watPolDirect;
import static ffx.numerics.multipole.GKTensorGlobalTest.watPolEnergy;
import static ffx.numerics.multipole.GKTensorQI.GK_MULTIPOLE_ORDER.DIPOLE;
import static ffx.numerics.multipole.GKTensorQI.GK_MULTIPOLE_ORDER.MONOPOLE;
import static ffx.numerics.multipole.GKTensorQI.GK_MULTIPOLE_ORDER.QUADRUPOLE;
import static org.junit.Assert.assertEquals;

import ffx.numerics.multipole.GKTensorQI.GK_MULTIPOLE_ORDER;
import ffx.numerics.multipole.GKTensorQI.GK_TENSOR_MODE;
import org.junit.Test;

/**
 * Test the GK tensor evaluated in the quasi-internal coordinate frame.
 * <p>
 * There is no quadrupole finite-difference test because the trace of the quadrupole potential is
 * neglected; thus the derivatives of the quadrupole potential are correct when summed over the
 * trace, but not on a per-element basis.
 */
public class GKTensorQITest {

  private final double tolerance = 1.0e-9;
  private final double fdTolerance = 1.0e-6;

  private static final double gc = 2.455;
  private static final double Eh = 1.0;
  private static final double Es = 78.3;
  private static final double Ai = 2.0;
  private static final double Aj = 2.0;
  private static final double[] r = {0.7, 0.8, 0.9};

  @Test
  public void permanentEnergyTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);
    QIFrame qiFrame = new QIFrame(rWater);
    qiFrame.rotatePolarizableMultipole(mI);
    qiFrame.rotatePolarizableMultipole(mK);

    // Apply the GK monopole tensor.
    int order = 2;
    GKTensorQI monopoleGK = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    monopoleGK.setR(rWater, bornI, bornK);
    monopoleGK.generateTensor();
    double c = monopoleGK.multipoleEnergy(mI, mK);

    // Apply the GK dipole tensor.
    order = 3;
    GKTensorQI dipoleGK = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    dipoleGK.setR(rWater, bornI, bornK);
    dipoleGK.generateTensor();
    double d = dipoleGK.multipoleEnergy(mI, mK);

    // Apply the GK quadrupole tensor.
    order = 4;
    GKTensorQI quadrupoleGK = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
    quadrupoleGK.setR(rWater, bornI, bornK);
    quadrupoleGK.generateTensor();
    double q = quadrupoleGK.multipoleEnergy(mI, mK);

    assertEquals("GK Permanent Energy", watPermEnergy, c + d + q, tolerance);
  }

  @Test
  public void permanentEnergyAndGradientTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);
    QIFrame qiFrame = new QIFrame(rWater);
    qiFrame.rotatePolarizableMultipole(mI);
    qiFrame.rotatePolarizableMultipole(mK);

    double[] gradI = new double[3];
    double[] gradK = new double[3];
    double[] torqueI = new double[3];
    double[] torqueK = new double[3];

    double[] tempGradI = new double[3];
    double[] tempGradK = new double[3];
    double[] tempTorqueI = new double[3];
    double[] tempTorqueK = new double[3];

    // Apply the GK monopole tensor.
    int order = 3;
    GKTensorQI gkMonopoleTensor = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double e = gkMonopoleTensor.multipoleEnergyAndGradient(mI, mK, gradI, gradK, torqueI, torqueK);
    double db = gkMonopoleTensor.multipoleEnergyBornGrad(mI, mK);

    // Apply the GK dipole tensor.
    order = 4;
    GKTensorQI gkDipoleTensor = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    e += gkDipoleTensor.multipoleEnergyAndGradient(mI, mK, tempGradI, tempGradK, tempTorqueI,
        tempTorqueK);
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
    order = 5;
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(GK_MULTIPOLE_ORDER.QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    e += gkQuadrupoleTensor.multipoleEnergyAndGradient(mI, mK, tempGradI, tempGradK, tempTorqueI,
        tempTorqueK);
    db += gkQuadrupoleTensor.multipoleEnergyBornGrad(mI, mK);

    for (int i = 0; i < 3; i++) {
      gradI[i] += tempGradI[i];
      gradK[i] += tempGradK[i];
      torqueI[i] += tempTorqueI[i];
      torqueK[i] += tempTorqueK[i];
    }

    qiFrame.toGlobal(gradI);
    qiFrame.toGlobal(torqueI);
    qiFrame.toGlobal(torqueK);

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

    QIFrame qiFrame = new QIFrame(rWater);
    qiFrame.rotatePolarizableMultipole(mI);
    qiFrame.rotatePolarizableMultipole(mK);

    // Apply the GK monopole tensor.
    int order = 2;
    GKTensorQI gkMonopoleTensor = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK dipole tensor.
    order = 3;
    GKTensorQI gkDipoleTensor = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK quadrupole tensor.
    order = 4;
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.polarizationEnergy(mI, mK);

    assertEquals("GK Polarization Energy", watPolEnergy, c + d + q, tolerance);
  }

  @Test
  public void polarizationEnergyDirectTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uIDirect, uIDirect);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uKDirect, uKDirect);

    QIFrame qiFrame = new QIFrame(rWater);
    qiFrame.rotatePolarizableMultipole(mI);
    qiFrame.rotatePolarizableMultipole(mK);

    // Apply the GK monopole tensor.
    int order = 1; // Induced Dipole - Monopole Energy.
    GKTensorQI gkMonopoleTensor = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK dipole tensor.
    order = 3; // Induced Dipole - Quadrupole Energy.
    GKTensorQI gkDipoleTensor = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(rWater, bornI, bornK);
    gkDipoleTensor.generateTensor();
    double d = gkDipoleTensor.polarizationEnergy(mI, mK);

    // Apply the GK quadrupole tensor.
    order = 3; // Induced Dipole - Quadrupole Energy.
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(rWater, bornI, bornK);
    gkQuadrupoleTensor.generateTensor();
    double q = gkQuadrupoleTensor.polarizationEnergy(mI, mK);

    assertEquals("GK Direct Polarization Energy", watPolDirect, c + d + q, tolerance);
  }

  @Test
  public void polarizationEnergyAndGradientTest() {
    PolarizableMultipole mI = new PolarizableMultipole(multI, uI, uI);
    PolarizableMultipole mK = new PolarizableMultipole(multK, uK, uK);

    QIFrame qiFrame = new QIFrame(rWater);
    qiFrame.rotatePolarizableMultipole(mI);
    qiFrame.rotatePolarizableMultipole(mK);

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
    GKTensorQI gkMonopoleTensor = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 1.0, gradI, torqueI, torqueK);
    double db = gkMonopoleTensor.polarizationEnergyBornGrad(mI, mK);

    // Apply the GK dipole tensor.
    order = 4; // Induced Dipole - Quadrupole Energy Gradient.
    GKTensorQI gkDipoleTensor = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
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
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
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

    qiFrame.toGlobal(gradI);
    qiFrame.toGlobal(torqueI);
    qiFrame.toGlobal(torqueK);

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

    QIFrame qiFrame = new QIFrame(rWater);
    qiFrame.rotatePolarizableMultipole(mI);
    qiFrame.rotatePolarizableMultipole(mK);

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
    GKTensorQI gkMonopoleTensor = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(rWater, bornI, bornK);
    gkMonopoleTensor.generateTensor();
    double c = gkMonopoleTensor.polarizationEnergyAndGradient(mI, mK,
        1.0, 1.0, 0.0, gradI, torqueI, torqueK);
    double db = gkMonopoleTensor.polarizationEnergyBornGrad(mI, mK);

    // Apply the GK dipole tensor.
    order = 4; // Induced Dipole - Quadrupole Energy Gradient.
    GKTensorQI gkDipoleTensor = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
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
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
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

    qiFrame.toGlobal(gradI);
    qiFrame.toGlobal(torqueI);
    qiFrame.toGlobal(torqueK);

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
    GKTensorQI tensorQI = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    tensorQI.setR(r, Ai, Aj);
    double[] work = new double[order + 1];
    tensorQI.source(work);

    // Test the "bn" method.
    double bn0 = tensorQI.bn(0);
    double bn1 = tensorQI.bn(1);
    double bn2 = tensorQI.bn(2);
    double bn3 = tensorQI.bn(3);
    double bn4 = tensorQI.bn(4);
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
    double A00 = tensorQI.anm(0, 0);
    double A01 = tensorQI.anm(0, 1);
    double A02 = tensorQI.anm(0, 2);
    double A03 = tensorQI.anm(0, 3);
    final double mapleA00 = -0.4319767286;
    final double mapleA01 = 0.05505753880;
    final double mapleA02 = -0.01542067105;
    final double mapleA03 = 0.005809287830;
    assertEquals("A00", mapleA00, A00, tolerance);
    assertEquals("A01", mapleA01, A01, tolerance);
    assertEquals("A02", mapleA02, A02, tolerance);
    assertEquals("A03", mapleA03, A03, tolerance);
    // Monpole potential Born chain-rule derivatives.
    tensorQI.setMode(GK_TENSOR_MODE.BORN);
    tensorQI.source(work);
    double B00 = tensorQI.bnm(0, 0);
    double B01 = tensorQI.bnm(0, 1);
    double B02 = tensorQI.bnm(0, 2);
    final double mapleB00 = 0.04064563792;
    final double mapleB01 = -0.01690706430;
    final double mapleB02 = 0.008229167114;
    assertEquals("B00", mapleB00, B00, tolerance);
    assertEquals("B01", mapleB01, B01, tolerance);
    assertEquals("B02", mapleB02, B02, tolerance);

    // Dipole potential and derivatives.
    order = 4;
    work = new double[order + 1];
    tensorQI = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    tensorQI.setR(r, Ai, Aj);
    tensorQI.source(work);
    double A10 = tensorQI.anm(1, 0);
    double A11 = tensorQI.anm(1, 1);
    double A12 = tensorQI.anm(1, 2);
    double A13 = tensorQI.anm(1, 3);
    final double mapleA10 = 0.08218283800;
    final double mapleA11 = -0.03142380936;
    final double mapleA12 = 0.01681150454;
    final double mapleA13 = -0.01106715285;
    assertEquals("A10", mapleA10, A10, tolerance);
    assertEquals("A11", mapleA11, A11, tolerance);
    assertEquals("A12", mapleA12, A12, tolerance);
    assertEquals("A13", mapleA13, A13, tolerance);
    // Dipole potential Born chain-rule derivatives.
    tensorQI.setMode(GK_TENSOR_MODE.BORN);
    tensorQI.source(work);
    double B10 = tensorQI.bnm(1, 0);
    double B11 = tensorQI.bnm(1, 1);
    double B12 = tensorQI.bnm(1, 2);
    final double mapleB10 = -0.02319829046;
    final double mapleB11 = 0.01556309100;
    final double mapleB12 = -0.01202628190;
    assertEquals("B10", mapleB10, B10, tolerance);
    assertEquals("B11", mapleB11, B11, tolerance);
    assertEquals("B12", mapleB12, B12, tolerance);

    // Quadrupole potential and derivatives.
    order = 5;
    work = new double[order + 1];
    tensorQI = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
    tensorQI.setR(r, Ai, Aj);
    tensorQI.source(work);
    double A20 = tensorQI.anm(2, 0);
    double A21 = tensorQI.anm(2, 1);
    double A22 = tensorQI.anm(2, 2);
    double A23 = tensorQI.anm(2, 3);
    final double mapleA20 = -0.04710532877;
    final double mapleA21 = 0.03001901831;
    final double mapleA22 = -0.02371209206;
    final double mapleA23 = 0.02187861148;
    assertEquals("A20", mapleA20, A20, tolerance);
    assertEquals("A21", mapleA21, A21, tolerance);
    assertEquals("A22", mapleA22, A22, tolerance);
    assertEquals("A23", mapleA23, A23, tolerance);
    // Quadrupole potential Born chain-rule derivatives.
    tensorQI.setMode(GK_TENSOR_MODE.BORN);
    tensorQI.source(work);
    double B20 = tensorQI.bnm(2, 0);
    double B21 = tensorQI.bnm(2, 1);
    double B22 = tensorQI.bnm(2, 2);
    final double mapleB20 = 0.02216121853;
    final double mapleB21 = -0.02051645869;
    final double mapleB22 = 0.02137054028;
    assertEquals("B20", mapleB20, B20, tolerance);
    assertEquals("B21", mapleB21, B21, tolerance);
    assertEquals("B22", mapleB22, B22, tolerance);
  }

  @Test
  public void chargeTensorTest() {
    // Monopole potential and derivatives.
    int order = 3;
    GKTensorQI gkMonopoleTensor = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    gkMonopoleTensor.setR(r, Ai, Aj);
    double[] work = new double[order + 1];
    gkMonopoleTensor.source(work);

    double A00 = gkMonopoleTensor.anm(0, 0);
    double A01 = gkMonopoleTensor.anm(0, 1);
    double A02 = gkMonopoleTensor.anm(0, 2);
    double A03 = gkMonopoleTensor.anm(0, 3);
    gkMonopoleTensor.generateTensor();

    double z = length(r);
    assertEquals("R000", A00, gkMonopoleTensor.R000, tolerance);
    assertEquals("R001", z * A01, gkMonopoleTensor.R001, tolerance);
    assertEquals("R002", z * z * A02 + A01, gkMonopoleTensor.R002, tolerance);
    assertEquals("R003", z * z * z * A03 + 3.0 * z * A02, gkMonopoleTensor.R003, tolerance);

    // Check Born radii chain rule terms.
    gkMonopoleTensor.setMode(GK_TENSOR_MODE.BORN);
    gkMonopoleTensor.generateTensor();
    double B00 = gkMonopoleTensor.bnm(0, 0);
    double B01 = gkMonopoleTensor.bnm(0, 1);
    double B02 = gkMonopoleTensor.bnm(0, 2);
    assertEquals("B000", B00, gkMonopoleTensor.R000, tolerance);
    assertEquals("B001", z * B01, gkMonopoleTensor.R001, tolerance);
    assertEquals("B002", z * z * B02 + B01, gkMonopoleTensor.R002, tolerance);
  }

  @Test
  public void dipoleTensorTest() {
    // Dipole potential and derivatives.
    int order = 4;
    GKTensorQI gkDipoleTensor = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    gkDipoleTensor.setR(r, Ai, Aj);
    double[] work = new double[order + 1];
    gkDipoleTensor.source(work);

    double A10 = gkDipoleTensor.anm(1, 0);
    double A11 = gkDipoleTensor.anm(1, 1);
    double A12 = gkDipoleTensor.anm(1, 2);
    double A13 = gkDipoleTensor.anm(1, 3);
    gkDipoleTensor.generateTensor();

    // No charge potential.
    assertEquals("R000", 0.0, gkDipoleTensor.R000, tolerance);

    double z = length(r);
    assertEquals("R001", z * A10, gkDipoleTensor.R001, tolerance);
    assertEquals("R002", z * z * A11 + A10, gkDipoleTensor.R002, tolerance);
    assertEquals("R003", z * z * z * A12 + 3.0 * z * A11, gkDipoleTensor.R003, tolerance);
    assertEquals("R004", z * z * z * z * A13 + 6.0 * z * z * A12 + 3.0 * A11, gkDipoleTensor.R004,
        tolerance);

    // Check Born radii chain rule terms.
    gkDipoleTensor.setMode(GK_TENSOR_MODE.BORN);
    gkDipoleTensor.generateTensor();
    double B10 = gkDipoleTensor.bnm(1, 0);
    double B11 = gkDipoleTensor.bnm(1, 1);
    double B12 = gkDipoleTensor.bnm(1, 2);
    assertEquals("B001", z * B10, gkDipoleTensor.R001, tolerance);
    assertEquals("B002", z * z * B11 + B10, gkDipoleTensor.R002, tolerance);
    assertEquals("B003", z * z * z * B12 + 3.0 * z * B11, gkDipoleTensor.R003, tolerance);
  }

  @Test
  public void quadrupoleTensorTest() {
    // Quadrupole potential and derivatives.
    int order = 5;
    GKTensorQI gkQuadrupoleTensor = new GKTensorQI(QUADRUPOLE, order, gc, Eh, Es);
    gkQuadrupoleTensor.setR(r, Ai, Aj);
    double[] work = new double[order + 1];
    gkQuadrupoleTensor.source(work);

    double A20 = gkQuadrupoleTensor.anm(2, 0);
    double A21 = gkQuadrupoleTensor.anm(2, 1);
    double A22 = gkQuadrupoleTensor.anm(2, 2);
    double A23 = gkQuadrupoleTensor.anm(2, 3);
    gkQuadrupoleTensor.generateTensor();

    // No charge potential.
    assertEquals("R000", 0.0, gkQuadrupoleTensor.R000, tolerance);
    // No dipole potential.
    assertEquals("R100", 0.0, gkQuadrupoleTensor.R100, tolerance);

    double z = length(r);
    double x = 0.0;
    double y = 0.0;

    // Trace.
    assertEquals("R200", x * x * A20, gkQuadrupoleTensor.R200, tolerance);
    assertEquals("R020", y * y * A20, gkQuadrupoleTensor.R020, tolerance);
    assertEquals("R002", z * z * A20, gkQuadrupoleTensor.R002, tolerance);

    // Gradient of the trace with respect to z (note that the x*A20 term sums to 2*x*A20 over the trace).
    assertEquals("R003", z * z * z * A21 + 3.0 * z * A20, gkQuadrupoleTensor.R003, tolerance);
    assertEquals("R201", z * A20, gkQuadrupoleTensor.R201, tolerance);
    assertEquals("R021", z * A20, gkQuadrupoleTensor.R021, tolerance);

    // Higher order Qzz terms.
    assertEquals("R004", z * z * z * z * A22 + 6.0 * z * z * A21 + 3.0 * A20,
        gkQuadrupoleTensor.R004, tolerance);
    assertEquals("R005", z * z * z * z * z * A23 + 10.0 * z * z * z * A22 + 15.0 * z * A21,
        gkQuadrupoleTensor.R005, tolerance);

    // Qxy
    assertEquals("R110", 0.0, gkQuadrupoleTensor.R110, tolerance);
    assertEquals("R111", 0.0, gkQuadrupoleTensor.R111, tolerance);
    assertEquals("R310", 0.0, gkQuadrupoleTensor.R310, tolerance);
    assertEquals("R220", A20, gkQuadrupoleTensor.R220, tolerance);
    assertEquals("R211", 0.0, gkQuadrupoleTensor.R211, tolerance);
    assertEquals("R401", 3.0 * z * A21, gkQuadrupoleTensor.R401, tolerance);
    assertEquals("R023", z * z * z * A22 + 3.0 * z * A21, gkQuadrupoleTensor.R023, tolerance);
    assertEquals("R311", 0.0, gkQuadrupoleTensor.R311, tolerance);

    // Check Born radii chain rule terms.
    gkQuadrupoleTensor.setMode(GK_TENSOR_MODE.BORN);
    gkQuadrupoleTensor.generateTensor();
    double B20 = gkQuadrupoleTensor.bnm(2, 0);
    double B21 = gkQuadrupoleTensor.bnm(2, 1);
    double B22 = gkQuadrupoleTensor.bnm(2, 2);
    // No charge or dipole potential.
    assertEquals("B000", 0.0, gkQuadrupoleTensor.R000, tolerance);
    assertEquals("B100", 0.0, gkQuadrupoleTensor.R100, tolerance);
    // Born chain rule for the potential for the quadrupole trace
    assertEquals("B200", 0.0, gkQuadrupoleTensor.R200, tolerance);
    assertEquals("B020", 0.0, gkQuadrupoleTensor.R020, tolerance);
    assertEquals("B002", z * z * B20, gkQuadrupoleTensor.R002, tolerance);
    // Born chain rule for the partial derivative of the trace with respect to Z.
    assertEquals("B003", z * z * z * B21 + 3.0 * z * B20, gkQuadrupoleTensor.R003, tolerance);
    assertEquals("B201", z * B20, gkQuadrupoleTensor.R201, tolerance);
    assertEquals("B021", z * B20, gkQuadrupoleTensor.R021, tolerance);
    // Higher order term for Qzz.
    assertEquals("B004", z * z * z * z * B22 + 6.0 * z * z * B21 + 3.0 * B20,
        gkQuadrupoleTensor.R004, tolerance);
  }

  @Test
  public void chargeFiniteDifferenceTest() {
    int order = 6;
    GKTensorQI gkTensorQI = new GKTensorQI(MONOPOLE, order, gc, Eh, Es);
    double[] r = {0.7, 0.8, 0.9};
    r[2] = length(r);
    r[0] = 0.0;
    r[1] = 0.0;
    gkTensorQI.setR(r, Ai, Aj);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
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
    GKTensorQI gkTensorQI = new GKTensorQI(DIPOLE, order, gc, Eh, Es);
    double[] r = {0.7, 0.8, 0.9};
    r[2] = length(r);
    r[0] = 0.0;
    r[1] = 0.0;
    gkTensorQI.setR(r, Ai, Aj);
    int tensorCount = MultipoleTensor.tensorCount(order);
    double[] tensor = new double[tensorCount];
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

    int start = gkTensorQI.multipoleOrder.getOrder();

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
