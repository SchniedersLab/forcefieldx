// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import static ffx.numerics.multipole.GKMultipoleOrder.DIPOLE;
import static ffx.numerics.multipole.GKMultipoleOrder.MONOPOLE;
import static ffx.numerics.multipole.GKMultipoleOrder.QUADRUPOLE;
import static ffx.numerics.multipole.GKTensorMode.BORN;
import static ffx.numerics.multipole.GKTensorMode.POTENTIAL;

/**
 * GKEnergyGlobal computes the generalized Kirkwood energy and forces in a global frame.
 */
public class GKEnergyGlobalSIMD {

  private final GKSourceSIMD gkSource;
  private final GKTensorGlobalSIMD gkMonopole;
  private final GKTensorGlobalSIMD gkDipole;
  private final GKTensorGlobalSIMD gkQuadrupole;

  private final DoubleVector one = DoubleVector.zero(DoubleVector.SPECIES_PREFERRED).add(1.0);

  /**
   * Constructor for GKEnergyGlobal.
   *
   * @param gkc      The GK generalizing function constant.
   * @param epsilon  The solvent dielectric.
   * @param gradient If true, compute the gradient and torque.
   */
  public GKEnergyGlobalSIMD(double gkc, double epsilon, boolean gradient) {
    int monopoleOrder = 2;
    int dipoleOrder = 3;
    int quadrupoleOrder = 4;
    if (gradient) {
      monopoleOrder = 3;
      dipoleOrder = 4;
      quadrupoleOrder = 5;
    }
    gkSource = new GKSourceSIMD(quadrupoleOrder, gkc);
    gkMonopole = new GKTensorGlobalSIMD(MONOPOLE, monopoleOrder, gkSource, 1.0, epsilon);
    gkDipole = new GKTensorGlobalSIMD(DIPOLE, dipoleOrder, gkSource, 1.0, epsilon);
    gkQuadrupole = new GKTensorGlobalSIMD(QUADRUPOLE, quadrupoleOrder, gkSource, 1.0, epsilon);
  }

  /**
   * Initialize the potential.
   *
   * @param r   The separation. vector.
   * @param r2  The squared separation.
   * @param rbi The Born radius of atom i.
   * @param rbk The Born radius of atom k.
   */
  public void initPotential(DoubleVector[] r, DoubleVector r2, DoubleVector rbi, DoubleVector rbk) {
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
  public void initBorn(DoubleVector[] r, DoubleVector r2, DoubleVector rbi, DoubleVector rbk) {
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
  public DoubleVector multipoleEnergy(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    DoubleVector em = gkMonopole.multipoleEnergy(mI, mK);
    DoubleVector ed = gkDipole.multipoleEnergy(mI, mK);
    DoubleVector eq = gkQuadrupole.multipoleEnergy(mI, mK);
    return em.add(ed).add(eq);
  }

  /**
   * Compute the polarization energy.
   *
   * @param mI The polarizable multipole of atom i.
   * @param mK The polarizable multipole of atom k.
   * @return The polarization energy.
   */
  public DoubleVector polarizationEnergy(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    DoubleVector emp = gkMonopole.polarizationEnergy(mI, mK);
    DoubleVector edp = gkDipole.polarizationEnergy(mI, mK);
    DoubleVector eqp = gkQuadrupole.polarizationEnergy(mI, mK);
    return emp.add(edp).add(eqp);
  }

  /**
   * Compute the multipole energy and gradient.
   *
   * @param mI The polarizable multipole of atom i.
   * @param mK The polarizable multipole of atom k.
   * @param gI The gradient for atom i.
   * @param tI The torque on atom i.
   * @param tK The torque on atom k.
   * @return The multipole energy.
   */
  public DoubleVector multipoleEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK,
                                                 DoubleVector[] gI, DoubleVector[] tI, DoubleVector[] tK) {
    DoubleVector[] gK = new DoubleVector[3];
    DoubleVector em = gkMonopole.multipoleEnergyAndGradient(mI, mK, gI, gK, tI, tK);
    DoubleVector ed = gkDipole.multipoleEnergyAndGradient(mI, mK, gI, gK, tI, tK);
    DoubleVector eq = gkQuadrupole.multipoleEnergyAndGradient(mI, mK, gI, gK, tI, tK);
    return em.add(ed).add(eq);
  }

  /**
   * Compute the polarization energy and gradient.
   *
   * @param mI         The polarizable multipole of atom i.
   * @param mK         The polarizable multipole of atom k.
   * @param mutualMask The mutual polarization mask.
   * @param gI         The gradient for atom i.
   * @param tI         The torque on atom i.
   * @param tK         The torque on atom k.
   * @return The polarization energy.
   */
  public DoubleVector polarizationEnergyAndGradient(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK, DoubleVector mutualMask,
                                                    DoubleVector[] gI, DoubleVector[] tI, DoubleVector[] tK) {
    DoubleVector emp = gkMonopole.polarizationEnergyAndGradient(mI, mK, one, one, mutualMask, gI, tI, tK);
    DoubleVector edp = gkDipole.polarizationEnergyAndGradient(mI, mK, one, one, mutualMask, gI, tI, tK);
    DoubleVector eqp = gkQuadrupole.polarizationEnergyAndGradient(mI, mK, one, one, mutualMask, gI, tI, tK);
    // Sum the GK polarization interaction energy.
    return emp.add(edp).add(eqp);
  }

  /**
   * Compute the Born chain-rule term for the multipole energy.
   *
   * @param mI The polarizable multipole of atom i.
   * @param mK The polarizable multipole of atom k.
   * @return The Born chain-rule term.
   */
  public DoubleVector multipoleEnergyBornGrad(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK) {
    DoubleVector db = gkMonopole.multipoleEnergyBornGrad(mI, mK);
    db = db.add(gkDipole.multipoleEnergyBornGrad(mI, mK));
    db = db.add(gkQuadrupole.multipoleEnergyBornGrad(mI, mK));
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
  public DoubleVector polarizationEnergyBornGrad(PolarizableMultipoleSIMD mI, PolarizableMultipoleSIMD mK, boolean mutual) {
    // Compute the GK polarization Born chain-rule term.
    DoubleVector db = gkMonopole.polarizationEnergyBornGrad(mI, mK);
    db = db.add(gkDipole.polarizationEnergyBornGrad(mI, mK));
    db = db.add(gkQuadrupole.polarizationEnergyBornGrad(mI, mK));
    // Add the mutual polarization contribution to Born chain-rule term.
    if (mutual) {
      db = db.add(gkDipole.mutualPolarizationEnergyBornGrad(mI, mK));
    }
    return db;
  }

}
