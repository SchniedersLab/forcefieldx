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

import ffx.openmm.amoeba.TorsionTorsionForce;
import ffx.potential.MolecularAssembly;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;

import javax.annotation.Nullable;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
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
   * The ForceFieldEnergy instance for the first topology.
   */
  protected ForceFieldEnergy forceFieldEnergy;
  /**
   * The ForceFieldEnergy instance for the second topology.
   */
  protected ForceFieldEnergy forceFieldEnergy2;
  /**
   * OpenMM Custom Bond Force for topology 2.
   */
  protected BondForce bondForce2 = null;
  /**
   * OpenMM Custom Angle Force for topology 2.
   */
  protected AngleForce angleForce2 = null;
  /**
   * OpenMM Custom In-Plane Angle Force for topology 2.
   */
  protected InPlaneAngleForce inPlaneAngleForce2 = null;
  /**
   * OpenMM Custom Stetch-Bend Force for topology 2.
   */
  protected StretchBendForce stretchBendForce2 = null;
  /**
   * OpenMM Custom Urey-Bradley Force for topology 2.
   */
  protected UreyBradleyForce ureyBradleyForce2 = null;
  /**
   * OpenMM Custom Out-of-Plane Bend Force for topology 2.
   */
  protected OutOfPlaneBendForce outOfPlaneBendForce2 = null;
  /**
   * OpenMM Custom Pi-Orbital Torsion Force for topology 2.
   */
  protected PiOrbitalTorsionForce piOrbitalTorsionForce2 = null;
  /**
   * OpenMM Custom Torsion Force for topology 2.
   */
  protected TorsionForce torsionForce2 = null;
  /**
   * OpenMM Custom Improper Torsion Force for topology 2.
   */
  protected ImproperTorsionForce improperTorsionForce2 = null;
  /**
   * OpenMM Custom Stretch-Torsion Force for topology 2.
   */
  protected StretchTorsionForce stretchTorsionForce2 = null;
  /**
   * OpenMM Custom Angle-Torsion Force for topology 2.
   */
  protected AngleTorsionForce angleTorsionForce2 = null;
  /**
   * OpenMM Custom Torsion-Torsion Force for topology 2.
   * ToDo: There is no updateParametersInContext method for the AmoebaTorsionTorsionForce.
   * We assume that the AmoebaTorsionTorsionForce is constant along the alchemical path.
   */
  protected AmoebaTorsionTorsionForce amoebaTorsionTorsionForce2 = null;
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

    forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(0);
    forceFieldEnergy2 = openMMDualTopologyEnergy.getForceFieldEnergy(1);
    forceField = openMMDualTopologyEnergy.getMolecularAssembly(0).getForceField();

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
   * Get the Crystal instance.
   *
   * @return the Crystal instance.
   */
  @Override
  public Crystal getCrystal() {
    return forceFieldEnergy.getCrystal();
  }

  /**
   * Set the default values of the vectors defining the axes of the periodic box (measured in nm).
   *
   * <p>Any newly created Context will have its box vectors set to these. They will affect any
   * Force added to the System that uses periodic boundary conditions.
   *
   * <p>Triclinic boxes are supported, but the vectors must satisfy certain requirements. In
   * particular, a must point in the x direction, b must point "mostly" in the y direction, and c
   * must point "mostly" in the z direction. See the documentation for details.
   */
  @Override
  protected void setDefaultPeriodicBoxVectors() {
    Crystal crystal = forceFieldEnergy.getCrystal();
    if (!crystal.aperiodic()) {
      OpenMM_Vec3 a = new OpenMM_Vec3();
      OpenMM_Vec3 b = new OpenMM_Vec3();
      OpenMM_Vec3 c = new OpenMM_Vec3();
      double[][] Ai = crystal.Ai;
      a.x = Ai[0][0] * OpenMM_NmPerAngstrom;
      a.y = Ai[0][1] * OpenMM_NmPerAngstrom;
      a.z = Ai[0][2] * OpenMM_NmPerAngstrom;
      b.x = Ai[1][0] * OpenMM_NmPerAngstrom;
      b.y = Ai[1][1] * OpenMM_NmPerAngstrom;
      b.z = Ai[1][2] * OpenMM_NmPerAngstrom;
      c.x = Ai[2][0] * OpenMM_NmPerAngstrom;
      c.y = Ai[2][1] * OpenMM_NmPerAngstrom;
      c.z = Ai[2][2] * OpenMM_NmPerAngstrom;
      setDefaultPeriodicBoxVectors(a, b, c);
    }
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

    // Add Angle Force.
    angleForce = (AngleForce) AngleForce.constructForce(0, openMMDualTopologyEnergy);
    angleForce2 = (AngleForce) AngleForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(angleForce);
    addForce(angleForce2);

    // Add In-Plane Angle Force.
    inPlaneAngleForce = (InPlaneAngleForce) InPlaneAngleForce.constructForce(0, openMMDualTopologyEnergy);
    inPlaneAngleForce2 = (InPlaneAngleForce) InPlaneAngleForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(inPlaneAngleForce);
    addForce(inPlaneAngleForce2);

    // Add Stretch-Bend Force.
    stretchBendForce = (StretchBendForce) StretchBendForce.constructForce(0, openMMDualTopologyEnergy);
    stretchBendForce2 = (StretchBendForce) StretchBendForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(stretchBendForce);
    addForce(stretchBendForce2);

    // Add Urey-Bradley Force.
    ureyBradleyForce = (UreyBradleyForce) UreyBradleyForce.constructForce(0, openMMDualTopologyEnergy);
    ureyBradleyForce2 = (UreyBradleyForce) UreyBradleyForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(ureyBradleyForce);
    addForce(ureyBradleyForce2);

    // Add Out-of-Plane Bend Force.
    outOfPlaneBendForce = (OutOfPlaneBendForce) OutOfPlaneBendForce.constructForce(0, openMMDualTopologyEnergy);
    outOfPlaneBendForce2 = (OutOfPlaneBendForce) OutOfPlaneBendForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(outOfPlaneBendForce);
    addForce(outOfPlaneBendForce2);

    // Add Pi-Torsion Force.
    piOrbitalTorsionForce = (PiOrbitalTorsionForce) PiOrbitalTorsionForce.constructForce(0, openMMDualTopologyEnergy);
    piOrbitalTorsionForce2 = (PiOrbitalTorsionForce) PiOrbitalTorsionForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(piOrbitalTorsionForce);
    addForce(piOrbitalTorsionForce2);

    // Add Torsion Force.
    torsionForce = (TorsionForce) TorsionForce.constructForce(0, openMMDualTopologyEnergy);
    torsionForce2 = (TorsionForce) TorsionForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(torsionForce);
    addForce(torsionForce2);

    // Add Improper Torsion Force.
    improperTorsionForce = (ImproperTorsionForce) ImproperTorsionForce.constructForce(0, openMMDualTopologyEnergy);
    improperTorsionForce2 = (ImproperTorsionForce) ImproperTorsionForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(improperTorsionForce);
    addForce(improperTorsionForce2);

    // Add Stretch-Torsion coupling terms.
    stretchTorsionForce = (StretchTorsionForce) StretchTorsionForce.constructForce(0, openMMDualTopologyEnergy);
    stretchTorsionForce2 = (StretchTorsionForce) StretchTorsionForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(stretchTorsionForce);
    addForce(stretchTorsionForce2);

    // Add Angle-Torsion coupling terms.
    angleTorsionForce = (AngleTorsionForce) AngleTorsionForce.constructForce(0, openMMDualTopologyEnergy);
    angleTorsionForce2 = (AngleTorsionForce) AngleTorsionForce.constructForce(1, openMMDualTopologyEnergy);
    addForce(angleTorsionForce);
    addForce(angleTorsionForce2);

    // Add Torsion-Torsion Force.
    // ToDo: There is no updateParametersInContext method for the AmoebaTorsionTorsionForce.
    amoebaTorsionTorsionForce = (AmoebaTorsionTorsionForce) AmoebaTorsionTorsionForce.constructForce(0, openMMDualTopologyEnergy);
    addForce(amoebaTorsionTorsionForce);
    // amoebaTorsionTorsionForce2 = (AmoebaTorsionTorsionForce) AmoebaTorsionTorsionForce.constructForce(1, openMMDualTopologyEnergy);
    // addForce(amoebaTorsionTorsionForce2);

    // Check that the number of Torsion-Torsion forces in the two topologies is equal.
    TorsionTorsion[] torsionTorsion = openMMDualTopologyEnergy.getForceFieldEnergy(0).getTorsionTorsions();
    TorsionTorsion[] torsionTorsion2 = openMMDualTopologyEnergy.getForceFieldEnergy(1).getTorsionTorsions();
    int numTorsionTorsions = 0;
    int numTorsionTorsions2 = 0;
    if (torsionTorsion != null) {
      numTorsionTorsions = torsionTorsion.length;
    }
    if (torsionTorsion2 != null) {
      numTorsionTorsions2 = torsionTorsion2.length;
    }
    if (numTorsionTorsions != numTorsionTorsions2) {
      logger.severe(" The number of Torsion-Torsion forces in the two topologies do not match: "
          + numTorsionTorsions + " vs. " + numTorsionTorsions2);
    } else {
      // Check that the Torsion-Torsion instances match.
      for (int i = 0; i < numTorsionTorsions; i++) {
        TorsionTorsion tt1 = torsionTorsion[i];
        TorsionTorsion tt2 = torsionTorsion2[i];
        if (!tt1.equals(tt2)) {
          logger.severe(" The Torsion-Torsion terms in the two topologies do not match: " + tt1 + " vs. " + tt2);
        }
      }
    }

    VanDerWaals vdW1 = forceFieldEnergy.getVdwNode();
    VanDerWaals vdW2 = forceFieldEnergy2.getVdwNode();
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
        amoebaVDWForce2.setLambdaName("AmoebaVdwLambda2");
        addForce(amoebaVDWForce2);
      }

      ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
      ParticleMeshEwald pme2 = forceFieldEnergy2.getPmeNode();
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

//    if (updateBondedTerms) {
    if (bondForce != null) {
      bondForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (bondForce2 != null) {
      bondForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (angleForce != null) {
      angleForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (angleForce2 != null) {
      angleForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (inPlaneAngleForce != null) {
      inPlaneAngleForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (inPlaneAngleForce2 != null) {
      inPlaneAngleForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (stretchBendForce != null) {
      stretchBendForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (stretchBendForce2 != null) {
      stretchBendForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (ureyBradleyForce != null) {
      ureyBradleyForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (ureyBradleyForce2 != null) {
      ureyBradleyForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (outOfPlaneBendForce != null) {
      outOfPlaneBendForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (outOfPlaneBendForce2 != null) {
      outOfPlaneBendForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (piOrbitalTorsionForce != null) {
      piOrbitalTorsionForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (piOrbitalTorsionForce2 != null) {
      piOrbitalTorsionForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (torsionForce != null) {
      torsionForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (torsionForce2 != null) {
      torsionForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (improperTorsionForce != null) {
      improperTorsionForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (improperTorsionForce2 != null) {
      improperTorsionForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (stretchTorsionForce != null) {
      stretchTorsionForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (stretchTorsionForce2 != null) {
      stretchTorsionForce2.updateForce(1, openMMDualTopologyEnergy);
    }
    if (angleTorsionForce != null) {
      angleTorsionForce.updateForce(0, openMMDualTopologyEnergy);
    }
    if (angleTorsionForce2 != null) {
      angleTorsionForce2.updateForce(1, openMMDualTopologyEnergy);
    }

    // ToDo: there is no support in the OpenMM AmoebaTorsionTorsionForce to updateParametersInContext.

    if (amoebaVDWForce != null) {
      VanDerWaals vanDerWaals = forceFieldEnergy.getVdwNode();
      if (vanDerWaals.getLambdaTerm()) {
        double lambdaVDW = vanDerWaals.getLambda();
        openMMDualTopologyEnergy.getContext().setParameter("AmoebaVdwLambda", lambdaVDW);
      }
      atoms = forceFieldEnergy.getAtomArray();
      amoebaVDWForce.updateForce(atoms, 0, openMMDualTopologyEnergy);
    }
    if (amoebaVDWForce2 != null) {
      VanDerWaals vanDerWaals = forceFieldEnergy2.getVdwNode();
      if (vanDerWaals.getLambdaTerm()) {
        double lambdaVDW = vanDerWaals.getLambda();
        openMMDualTopologyEnergy.getContext().setParameter("AmoebaVdwLambda2", lambdaVDW);
      }
      atoms = forceFieldEnergy2.getAtomArray();
      amoebaVDWForce2.updateForce(atoms, 1, openMMDualTopologyEnergy);
    }
    if (amoebaMultipoleForce != null) {
      atoms = forceFieldEnergy.getAtomArray();
      amoebaMultipoleForce.updateForce(atoms, 0, openMMDualTopologyEnergy);
    }
    if (amoebaMultipoleForce2 != null) {
      atoms = forceFieldEnergy2.getAtomArray();
      amoebaMultipoleForce2.updateForce(atoms, 1, openMMDualTopologyEnergy);
    }
  }
}
