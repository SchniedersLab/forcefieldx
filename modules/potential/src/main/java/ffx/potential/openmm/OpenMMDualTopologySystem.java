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
package ffx.potential.openmm;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;

import javax.annotation.Nullable;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static java.lang.String.format;

/**
 * Create and manage an OpenMM Dual-Topology System.
 *
 * <p>The definition of a System involves four elements:
 *
 * <p>The particles and constraints are defined directly by the System object, while forces are
 * defined by objects that extend the Force class. After creating a System, call addParticle() once
 * for each particle, addConstraint() for each constraint, and addForce() for each Force.
 *
 * <p>In addition, particles may be designated as "virtual sites". These are particles whose
 * positions are computed automatically based on the positions of other particles. To define a
 * virtual site, call setVirtualSite(), passing in a VirtualSite object that defines the rules for
 * computing its position.
 */
public class OpenMMDualTopologySystem extends OpenMMSystem {

  private static final Logger logger = Logger.getLogger(OpenMMDualTopologySystem.class.getName());

  /**
   * The OpenMMDualTopologyEnergy instance.
   */
  private final OpenMMDualTopologyEnergy openMMDualTopologyEnergy;

  /**
   * The OpenMMEnergy instance for the second topology.
   */
  protected OpenMMEnergy openMMEnergy2;
  /**
   * OpenMM Custom Bond Force for topology 2.
   */
  protected BondForce bondForce2 = null;
  /**
   * OpenMM AMOEBA van der Waals Force for topology 2.
   */
  private AmoebaVdwForce amoebaVDWForce2 = null;
  /**
   * OpenMM AMOEBA Multipole Force for topology 2.
   */
  private AmoebaMultipoleForce amoebaMultipoleForce2 = null;

  /**
   * OpenMMDualTopologyEnergy constructor.
   *
   * @param openMMDualTopologyEnergy OpenMMDualTopologyEnergy instance.
   */
  public OpenMMDualTopologySystem(OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    this.openMMDualTopologyEnergy = openMMDualTopologyEnergy;

    openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(0);
    openMMEnergy2 = openMMDualTopologyEnergy.getOpenMMEnergy(1);

    MolecularAssembly molecularAssembly = openMMDualTopologyEnergy.getMolecularAssembly(0);
    forceField = molecularAssembly.getForceField();

    // This array contains the shared and alchemical atoms for the dual topology system.
    atoms = openMMDualTopologyEnergy.getDualTopologyAtoms(0);

    // Load atoms.
    try {
      addAtoms();
    } catch (Exception e) {
      logger.severe(" Atom without mass encountered.");
    }

    logger.info(format("\n OpenMM dual-topology system created with %d atoms.", atoms.length));
  }

  /**
   * Add forces to the system.
   */
  @Override
  public void addForces() {
    setDefaultPeriodicBoxVectors();

    // Set up rigid constraints.
    // These flags need to be set before bonds and angles are created below.
    boolean rigidHydrogen = forceField.getBoolean("RIGID_HYDROGEN", false);
    boolean rigidBonds = forceField.getBoolean("RIGID_BONDS", false);
    boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    if (rigidHydrogen || rigidBonds || rigidHydrogenAngles) {
      logger.severe(" Dual Topology does not support rigid hydrogen atoms.");
    }

    logger.info("\n Bonded Terms\n");

    // Add Bond Force.
    bondForce = (BondForce) BondForce.constructForce(0, openMMDualTopologyEnergy);
    bondForce2 = (BondForce) BondForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(bondForce);
    addForce(bondForce2);

    VanDerWaals vdW1 = openMMEnergy.getVdwNode();
    VanDerWaals vdW2 = openMMEnergy2.getVdwNode();
    if (vdW1 != null || vdW2 != null) {
      logger.info("\n Non-Bonded Terms");

      if (vdW1 != null) {
        // Add the vdW Force for Topology 1.
        amoebaVDWForce = (AmoebaVdwForce) AmoebaVdwForce.constructForce(0, openMMDualTopologyEnergy);
        addForce(amoebaVDWForce);
      }
      if (vdW2 != null) {
        // Add the vdW Force for Topology 2.
        amoebaVDWForce2 = (AmoebaVdwForce) AmoebaVdwForce.constructForce(1, openMMDualTopologyEnergy);
        amoebaVDWForce2.setUseLambdaComplement(OpenMM_True);
        addForce(amoebaVDWForce2);
      }

      ParticleMeshEwald pme = openMMEnergy.getPmeNode();
      ParticleMeshEwald pme2 = openMMEnergy2.getPmeNode();
      if (pme != null) {
        amoebaMultipoleForce = (AmoebaMultipoleForce) AmoebaMultipoleForce.constructForce(0, openMMDualTopologyEnergy);
        addForce(amoebaMultipoleForce);
      }
      if (pme2 != null) {
        amoebaMultipoleForce2 = (AmoebaMultipoleForce) AmoebaMultipoleForce.constructForce(1, openMMDualTopologyEnergy);
        addForce(amoebaMultipoleForce2);
      }

    }
  }

  /**
   * Get the number of variables.
   *
   * @return the number of variables.
   */
  @Override
  public int getNumberOfVariables() {
    int nActive = 0;
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      if (atoms[i].isActive()) {
        nActive++;
      }
    }
    int ret = nActive * 3;
    return ret;
  }

  /**
   * Calculate the number of degrees of freedom.
   *
   * @return Number of degrees of freedom.
   */
  @Override
  public int calculateDegreesOfFreedom() {
    int dof = getNumberOfVariables();

    // Remove OpenMM constraints.
    dof = dof - getNumConstraints();

    // Remove center of mass motion.
    if (cmMotionRemover != null) {
      dof -= 3;
    }

    return dof;
  }

  /**
   * Destroy the dual-topology system.
   */
  @Override
  public void free() {
    if (getPointer() != null) {
      logger.fine(" Free OpenMM dual-topology system.");
      destroy();
      logger.fine(" Free OpenMM dual-topology system completed.");
    }
  }

  /**
   * Update parameters if the Use flags and/or Lambda value has changed.
   *
   * @param atoms Atoms in this list are considered.
   */
  @Override
  public void updateParameters(@Nullable Atom[] atoms) {
    VanDerWaals vanDerWaals = openMMEnergy.getVdwNode();
    if (vanDerWaals != null && vanDerWaals.getLambdaTerm()) {
      double lambdaVDW = vanDerWaals.getLambda();
      // Update the lambda value.
      openMMDualTopologyEnergy.getContext().setParameter("AmoebaVdwLambda", lambdaVDW);
    }

    if (amoebaVDWForce != null) {
      atoms = openMMEnergy.getAtomArray();
      amoebaVDWForce.updateForce(atoms, 0, openMMDualTopologyEnergy);
    }

    if (amoebaVDWForce2 != null) {
      atoms = openMMEnergy2.getAtomArray();
      amoebaVDWForce2.updateForce(atoms, 1, openMMDualTopologyEnergy);
    }

    if (amoebaMultipoleForce != null) {
      atoms = openMMEnergy.getAtomArray();
      amoebaMultipoleForce.updateForce(atoms, 0, openMMDualTopologyEnergy);
    }

    if (amoebaMultipoleForce2 != null) {
      atoms = openMMEnergy2.getAtomArray();
      amoebaMultipoleForce2.updateForce(atoms, 1, openMMDualTopologyEnergy);
    }
  }
}