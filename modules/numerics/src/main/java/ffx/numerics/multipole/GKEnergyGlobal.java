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

import static ffx.numerics.multipole.GKSource.GK_MULTIPOLE_ORDER.DIPOLE;
import static ffx.numerics.multipole.GKSource.GK_MULTIPOLE_ORDER.MONOPOLE;
import static ffx.numerics.multipole.GKSource.GK_MULTIPOLE_ORDER.QUADRUPOLE;
import static ffx.numerics.multipole.GKSource.GK_TENSOR_MODE.BORN;
import static ffx.numerics.multipole.GKSource.GK_TENSOR_MODE.POTENTIAL;
import static java.util.Arrays.fill;

/**
 * GKEnergyGlobal computes the generalized Kirkwood energy and forces in a global frame.
 */
public class GKEnergyGlobal {

  private final GKSource gkSource;
  private final GKTensorGlobal gkMonopole;
  private final GKTensorGlobal gkDipole;
  private final GKTensorGlobal gkQuadrupole;

  /**
   * Constructor for GKEnergyGlobal.
   *
   * @param gkc      The GK generalizing function constant.
   * @param epsilon  The solvent dielectric.
   * @param gradient If true, compute the gradient and torque.
   */
  public GKEnergyGlobal(double gkc, double epsilon, boolean gradient) {
    int monopoleOrder = 2;
    int dipoleOrder = 3;
    int quadrupoleOrder = 4;
    if (gradient) {
      monopoleOrder = 3;
      dipoleOrder = 4;
      quadrupoleOrder = 5;
    }
    gkSource = new GKSource(quadrupoleOrder, gkc);
    gkMonopole = new GKTensorGlobal(MONOPOLE, monopoleOrder, gkSource, 1.0, epsilon);
    gkDipole = new GKTensorGlobal(DIPOLE, dipoleOrder, gkSource, 1.0, epsilon);
    gkQuadrupole = new GKTensorGlobal(QUADRUPOLE, quadrupoleOrder, gkSource, 1.0, epsilon);
  }

  /**
   * Initialize the potential.
   *
   * @param r   The separation. vector.
   * @param r2  The squared separation.
   * @param rbi The Born radius of atom i.
   * @param rbk The Born radius of atom k.
   */
  public void initPotential(double[] r, double r2, double rbi, double rbk) {
    gkSource.generateSource(POTENTIAL, QUADRUPOLE, r2, rbi, rbk);
    gkMonopole.setR(r);
    gkDipole.setR(r);
    gkQuadrupole.setR(r);
    gkMonopole.generateTensor();
    gkDipole.generateTensor();
    gkQuadrupole.generateTensor();
  }

  /**
   * Initialize for computing Born chain-rule terms.
   *
   * @param r   The separation vector.
   * @param r2  The squared separation.
   * @param rbi The Born radius of atom i.
   * @param rbk The Born radius of atom k.
   */
  public void initBorn(double[] r, double r2, double rbi, double rbk) {
    gkSource.generateSource(BORN, QUADRUPOLE, r2, rbi, rbk);
    gkMonopole.setR(r);
    gkDipole.setR(r);
    gkQuadrupole.setR(r);
    gkMonopole.generateTensor();
    gkDipole.generateTensor();
    gkQuadrupole.generateTensor();
  }

  /**
   * Compute the multipole energy.
   *
   * @param mI The polarizable multipole of atom i.
   * @param mK The polarizable multipole of atom k.
   * @return The multipole energy.
   */
  public double multipoleEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    double em = gkMonopole.multipoleEnergy(mI, mK);
    double ed = gkDipole.multipoleEnergy(mI, mK);
    double eq = gkQuadrupole.multipoleEnergy(mI, mK);
    return em + ed + eq;
  }

  /**
   * Compute the polarization energy.
   *
   * @param mI The polarizable multipole of atom i.
   * @param mK The polarizable multipole of atom k.
   * @return The polarization energy.
   */
  public double polarizationEnergy(PolarizableMultipole mI, PolarizableMultipole mK) {
    double emp = gkMonopole.polarizationEnergy(mI, mK);
    double edp = gkDipole.polarizationEnergy(mI, mK);
    double eqp = gkQuadrupole.polarizationEnergy(mI, mK);
    return emp + edp + eqp;
  }

  /**
   * Compute the multipole energy and gradient.
   *
   * @param mI      The polarizable multipole of atom i.
   * @param mK      The polarizable multipole of atom k.
   * @param gradI   The gradient for atom i.
   * @param torqueI The torque on atom i.
   * @param torqueK The torque on atom k.
   * @return The multipole energy.
   */
  public double multipoleEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
                                           double[] gradI, double[] torqueI, double[] torqueK) {
    double[] gI = new double[3];
    double[] gK = new double[3];
    double[] tI = new double[3];
    double[] tK = new double[3];
    double em = gkMonopole.multipoleEnergyAndGradient(mI, mK, gI, gK, tI, tK);
    for (int j = 0; j < 3; j++) {
      gradI[j] = gI[j];
      torqueI[j] = tI[j];
      torqueK[j] = tK[j];
    }
    fill(gI, 0.0);
    fill(gK, 0.0);
    fill(tI, 0.0);
    fill(tK, 0.0);
    double ed = gkDipole.multipoleEnergyAndGradient(mI, mK, gI, gK, tI, tK);
    for (int j = 0; j < 3; j++) {
      gradI[j] += gI[j];
      torqueI[j] += tI[j];
      torqueK[j] += tK[j];
    }
    fill(gI, 0.0);
    fill(gK, 0.0);
    fill(tI, 0.0);
    fill(tK, 0.0);
    double eq = gkQuadrupole.multipoleEnergyAndGradient(mI, mK, gI, gK, tI, tK);
    for (int j = 0; j < 3; j++) {
      gradI[j] += gI[j];
      torqueI[j] += tI[j];
      torqueK[j] += tK[j];
    }

    return em + ed + eq;
  }

  /**
   * Compute the polarization energy and gradient.
   *
   * @param mI         The polarizable multipole of atom i.
   * @param mK         The polarizable multipole of atom k.
   * @param mutualMask The mutual polarization mask.
   * @param gradI      The gradient for atom i.
   * @param torqueI    The torque on atom i.
   * @param torqueK    The torque on atom k.
   * @return The polarization energy.
   */
  public double polarizationEnergyAndGradient(PolarizableMultipole mI, PolarizableMultipole mK,
                                              double mutualMask,
                                              double[] gradI, double[] torqueI, double[] torqueK) {
    double[] gI = new double[3];
    double[] tI = new double[3];
    double[] tK = new double[3];
    double emp = gkMonopole.polarizationEnergyAndGradient(mI, mK, 1.0, 1.0, mutualMask, gI, tI, tK);
    for (int j = 0; j < 3; j++) {
      gradI[j] = gI[j];
      torqueI[j] = tI[j];
      torqueK[j] = tK[j];
    }
    fill(gI, 0.0);
    fill(tI, 0.0);
    fill(tK, 0.0);
    double edp = gkDipole.polarizationEnergyAndGradient(mI, mK, 1.0, 1.0, mutualMask, gI, tI, tK);
    for (int j = 0; j < 3; j++) {
      gradI[j] += gI[j];
      torqueI[j] += tI[j];
      torqueK[j] += tK[j];
    }
    fill(gI, 0.0);
    fill(tI, 0.0);
    fill(tK, 0.0);
    double eqp = gkQuadrupole.polarizationEnergyAndGradient(mI, mK, 1.0, 1.0, mutualMask, gI, tI,
        tK);
    for (int j = 0; j < 3; j++) {
      gradI[j] += gI[j];
      torqueI[j] += tI[j];
      torqueK[j] += tK[j];
    }

    // Sum the GK polarization interaction energy.
    return emp + edp + eqp;
  }

  /**
   * Compute the Born chain-rule term for the multipole energy.
   *
   * @param mI The polarizable multipole of atom i.
   * @param mK The polarizable multipole of atom k.
   * @return The Born chain-rule term.
   */
  public double multipoleEnergyBornGrad(PolarizableMultipole mI, PolarizableMultipole mK) {
    double db = gkMonopole.multipoleEnergyBornGrad(mI, mK);
    db += gkDipole.multipoleEnergyBornGrad(mI, mK);
    db += gkQuadrupole.multipoleEnergyBornGrad(mI, mK);
    return db;
  }

  /**
   * Compute the Born chain-rule term for the polarization energy.
   *
   * @param mI     The polarizable multipole of atom i.
   * @param mK     The polarizable multipole of atom k.
   * @param mutual If true, compute the mutual polarization contribution.
   * @return The Born chain-rule term.
   */
  public double polarizationEnergyBornGrad(PolarizableMultipole mI, PolarizableMultipole mK,
                                           boolean mutual) {
    // Compute the GK polarization Born chain-rule term.
    double db = gkMonopole.polarizationEnergyBornGrad(mI, mK);
    db += gkDipole.polarizationEnergyBornGrad(mI, mK);
    db += gkQuadrupole.polarizationEnergyBornGrad(mI, mK);
    // Add the mutual polarization contribution to Born chain-rule term.
    if (mutual) {
      db += gkDipole.mutualPolarizationEnergyBornGrad(mI, mK);
    }
    return db;
  }

}
