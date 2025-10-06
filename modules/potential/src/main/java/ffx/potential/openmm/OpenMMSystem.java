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

import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.openmm.AndersenThermostat;
import ffx.openmm.CMMotionRemover;
import ffx.openmm.Force;
import ffx.openmm.MonteCarloBarostat;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.nonbonded.implicit.ChandlerCavitation;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.terms.AnglePotentialEnergy;
import ffx.potential.terms.AngleTorsionPotentialEnergy;
import ffx.potential.terms.BondPotentialEnergy;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static ffx.potential.parameters.VDWType.VDW_TYPE.LENNARD_JONES;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * Create and manage an OpenMM System.
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
public class OpenMMSystem extends ffx.openmm.System {

  private static final Logger logger = Logger.getLogger(OpenMMSystem.class.getName());

  /**
   * The OpenMMEnergy instance.
   */
  private final OpenMMEnergy openMMEnergy;
  /**
   * The Force Field in use.
   */
  protected ForceField forceField;
  /**
   * Array of atoms in the system.
   */
  protected Atom[] atoms;

  /**
   * This flag indicates bonded force constants and equilibria are updated (e.g. during ManyBody
   * titration).
   */
  protected boolean updateBondedTerms = false;
  /**
   * Fixed charge softcore vdW force boolean.
   */
  protected boolean softcoreCreated = false;

  /**
   * OpenMM Custom Bond Force
   */
  protected BondForce bondForce = null;
  /**
   * OpenMM Custom Angle Force
   */
  protected AngleForce angleForce = null;
  /**
   * OpenMM Custom In-Plane Angle Force
   */
  protected InPlaneAngleForce inPlaneAngleForce = null;
  /**
   * OpenMM Custom Stretch-Bend Force
   */
  protected StretchBendForce stretchBendForce = null;
  /**
   * OpenMM Custom Urey-Bradley Force
   */
  protected UreyBradleyForce ureyBradleyForce = null;
  /**
   * OpenMM Custom Out-of-Plane Bend Force
   */
  protected OutOfPlaneBendForce outOfPlaneBendForce = null;
  /**
   * OpenMM Custom Pi-Torsion Force
   */
  protected PiOrbitalTorsionForce piOrbitalTorsionForce = null;
  /**
   * OpenMM AMOEBA Torsion Force.
   */
  protected TorsionForce torsionForce = null;
  /**
   * OpenMM Improper Torsion Force.
   */
  protected ImproperTorsionForce improperTorsionForce = null;
  /**
   * OpenMM Torsion-Torsion Force.
   */
  protected AmoebaTorsionTorsionForce amoebaTorsionTorsionForce = null;
  /**
   * OpenMM Stretch-Torsion Force.
   */
  protected StretchTorsionForce stretchTorsionForce = null;
  /**
   * OpenMM Angle-Torsion Force.
   */
  protected AngleTorsionForce angleTorsionForce = null;
  /**
   * OpenMM Restraint-Torsion Force.
   */
  protected RestrainTorsionsForce restrainTorsionsForce = null;
  /**
   * OpenMM Restrain-Position Force.
   */
  protected RestrainPositionsForce restrainPositionsForce = null;
  /**
   * OpenMM Restrain-Groups Force.
   */
  protected RestrainGroupsForce restrainGroupsForce = null;
  /**
   * OpenMM AMOEBA van der Waals Force.
   */
  protected AmoebaVdwForce amoebaVDWForce = null;
  /**
   * OpenMM AMOEBA Multipole Force.
   */
  protected AmoebaMultipoleForce amoebaMultipoleForce = null;
  /**
   * OpenMM Generalized Kirkwood Force.
   */
  protected AmoebaGeneralizedKirkwoodForce amoebaGeneralizedKirkwoodForce = null;
  /**
   * OpenMM AMOEBA WCA Dispersion Force.
   */
  protected AmoebaWcaDispersionForce amoebaWcaDispersionForce = null;
  /**
   * OpenMM AMOEBA WCA Cavitation Force.
   */
  protected AmoebaGKCavitationForce amoebaGKCavitationForce = null;
  /**
   * OpenMM Custom GB Force.
   */
  protected FixedChargeGBForce fixedChargeGBForce = null;
  /**
   * OpenMM Fixed Charge Non-Bonded Force.
   */
  protected FixedChargeNonbondedForce fixedChargeNonBondedForce = null;
  /**
   * Custom forces to handle alchemical transformations for fixed charge systems.
   */
  protected FixedChargeAlchemicalForces fixedChargeAlchemicalForces = null;

  /**
   * OpenMM thermostat. Currently, an Andersen thermostat is supported.
   */
  protected AndersenThermostat andersenThermostat = null;
  /**
   * Barostat to be added if NPT (isothermal-isobaric) dynamics is requested.
   */
  protected MonteCarloBarostat monteCarloBarostat = null;
  /**
   * OpenMM center-of-mass motion remover.
   */
  protected CMMotionRemover cmMotionRemover = null;

  /**
   * OpenMMSystem constructor.
   */
  public OpenMMSystem() {
    openMMEnergy = null;
    forceField = null;
    atoms = null;
  }

  /**
   * OpenMMSystem constructor.
   *
   * @param openMMEnergy ForceFieldEnergyOpenMM instance.
   */
  public OpenMMSystem(OpenMMEnergy openMMEnergy) {
    this.openMMEnergy = openMMEnergy;

    MolecularAssembly molecularAssembly = openMMEnergy.getMolecularAssembly();
    forceField = molecularAssembly.getForceField();
    atoms = molecularAssembly.getAtomArray();

    // Load atoms.
    try {
      addAtoms();
    } catch (Exception e) {
      logger.severe(" Atom without mass encountered.");
    }

    logger.info(format("\n OpenMM system created with %d atoms.", atoms.length));
  }

  /**
   * Get the Potential in use.
   *
   * @return The Potential.
   */
  public Potential getPotential() {
    return openMMEnergy;
  }

  /**
   * Get the atoms in the system.
   * @return Array of atoms in the system.
   */
  public Atom[] getAtoms() {
    return atoms;
  }

  /**
   * Add forces to the system.
   */
  public void addForces() {

    // Set up rigid constraints.
    // These flags need to be set before bonds and angles are created below.
    boolean rigidHydrogen = forceField.getBoolean("RIGID_HYDROGEN", false);
    boolean rigidBonds = forceField.getBoolean("RIGID_BONDS", false);
    boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    if (rigidHydrogen) {
      addHydrogenConstraints();
    }
    if (rigidBonds) {
      addUpBondConstraints();
    }
    if (rigidHydrogenAngles) {
      setUpHydrogenAngleConstraints();
    }

    logger.info("\n Bonded Terms\n");

    // Add Bond Force.
    if (rigidBonds) {
      logger.info(" Not creating AmoebaBondForce because bonds are constrained.");
    } else {
      bondForce = (BondForce) BondForce.constructForce(openMMEnergy);
      addForce(bondForce);
    }

    // Add Angle Force.
    angleForce = (AngleForce) AngleForce.constructForce(openMMEnergy);
    addForce(angleForce);

    // Add In-Plane Angle Force.
    inPlaneAngleForce = (InPlaneAngleForce) InPlaneAngleForce.constructForce(openMMEnergy);
    addForce(inPlaneAngleForce);

    // Add Stretch-Bend Force.
    stretchBendForce = (StretchBendForce) StretchBendForce.constructForce(openMMEnergy);
    addForce(stretchBendForce);

    // Add Urey-Bradley Force.
    ureyBradleyForce = (UreyBradleyForce) UreyBradleyForce.constructForce(openMMEnergy);
    addForce(ureyBradleyForce);

    // Out-of Plane Bend Force.
    outOfPlaneBendForce = (OutOfPlaneBendForce) OutOfPlaneBendForce.constructForce(openMMEnergy);
    addForce(outOfPlaneBendForce);

    // Add Pi-Torsion Force.
    piOrbitalTorsionForce = (PiOrbitalTorsionForce) PiOrbitalTorsionForce.constructForce(openMMEnergy);
    addForce(piOrbitalTorsionForce);

    // Add Torsion Force.
    torsionForce = (TorsionForce) TorsionForce.constructForce(openMMEnergy);
    addForce(torsionForce);

    // Add Improper Torsion Force.
    improperTorsionForce = (ImproperTorsionForce) ImproperTorsionForce.constructForce(openMMEnergy);
    addForce(improperTorsionForce);

    // Add Stretch-Torsion coupling terms.
    stretchTorsionForce = (StretchTorsionForce) StretchTorsionForce.constructForce(openMMEnergy);
    addForce(stretchTorsionForce);

    // Add Angle-Torsion coupling terms.
    angleTorsionForce = (AngleTorsionForce) AngleTorsionForce.constructForce(openMMEnergy);
    addForce(angleTorsionForce);

    // Add Torsion-Torsion Force.
    amoebaTorsionTorsionForce = (AmoebaTorsionTorsionForce) AmoebaTorsionTorsionForce.constructForce(openMMEnergy);
    addForce(amoebaTorsionTorsionForce);

    if (openMMEnergy.getRestrainMode() == ForceFieldEnergy.RestrainMode.ENERGY) {
      // Add Restrain Positions force.
      restrainPositionsForce = (RestrainPositionsForce) RestrainPositionsForce.constructForce(openMMEnergy);
      addForce(restrainPositionsForce);

      // Add Restrain Bonds force for each functional form.
      for (BondType.BondFunction function : BondType.BondFunction.values()) {
        RestrainDistanceForce restrainDistanceForce = (RestrainDistanceForce) RestrainDistanceForce.constructForce(function, openMMEnergy);
        addForce(restrainDistanceForce);
      }

      // Add Restrain Torsions force.
      restrainTorsionsForce = (RestrainTorsionsForce) RestrainTorsionsForce.constructForce(openMMEnergy);
      addForce(restrainTorsionsForce);
    }

    // Add Restrain Groups force.
    restrainGroupsForce = (RestrainGroupsForce) RestrainGroupsForce.constructForce(openMMEnergy);
    addForce(restrainGroupsForce);

    setDefaultPeriodicBoxVectors();

    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW != null) {
      logger.info("\n Non-Bonded Terms");
      VanDerWaalsForm vdwForm = vdW.getVDWForm();
      if (vdwForm.vdwType == LENNARD_JONES) {
        fixedChargeNonBondedForce = (FixedChargeNonbondedForce) FixedChargeNonbondedForce.constructForce(openMMEnergy);
        addForce(fixedChargeNonBondedForce);
        if (fixedChargeNonBondedForce != null) {
          GeneralizedKirkwood gk = openMMEnergy.getGK();
          if (gk != null) {
            fixedChargeGBForce = (FixedChargeGBForce) FixedChargeGBForce.constructForce(openMMEnergy);
            addForce(fixedChargeGBForce);
          }
        }
      } else {
        // Add vdW Force.
        amoebaVDWForce = (AmoebaVdwForce) AmoebaVdwForce.constructForce(openMMEnergy);
        addForce(amoebaVDWForce);

        // Add Multipole Force.
        amoebaMultipoleForce = (AmoebaMultipoleForce) AmoebaMultipoleForce.constructForce(openMMEnergy);
        addForce(amoebaMultipoleForce);

        // Add Generalized Kirkwood Force.
        GeneralizedKirkwood gk = openMMEnergy.getGK();
        if (gk != null) {
          amoebaGeneralizedKirkwoodForce = (AmoebaGeneralizedKirkwoodForce) AmoebaGeneralizedKirkwoodForce.constructForce(openMMEnergy);
          addForce(amoebaGeneralizedKirkwoodForce);

          // Add WCA Dispersion Force.
          DispersionRegion dispersionRegion = gk.getDispersionRegion();
          if (dispersionRegion != null) {
            amoebaWcaDispersionForce = (AmoebaWcaDispersionForce) AmoebaWcaDispersionForce.constructForce(openMMEnergy);
            addForce(amoebaWcaDispersionForce);
          }

          // Add a GaussVol Cavitation Force.
          ChandlerCavitation chandlerCavitation = gk.getChandlerCavitation();
          if (chandlerCavitation != null && chandlerCavitation.getGaussVol() != null) {
            amoebaGKCavitationForce = (AmoebaGKCavitationForce) AmoebaGKCavitationForce.constructForce(openMMEnergy);
            addForce(amoebaGKCavitationForce);
          }
        }
      }
    }
  }

  /**
   * Remove a force from the OpenMM System.
   * The OpenMM memory associated with the removed Force object is deleted.
   */
  public void removeForce(Force force) {
    if (force != null) {
      removeForce(force.getForceIndex());
    }
  }

  /**
   * Add an Andersen thermostat to the system.
   *
   * @param targetTemp Target temperature in Kelvins.
   */
  public void addAndersenThermostatForce(double targetTemp) {
    double collisionFreq = forceField.getDouble("COLLISION_FREQ", 1.0);
    addAndersenThermostatForce(targetTemp, collisionFreq);
  }

  /**
   * Add an Andersen thermostat to the system.
   *
   * @param targetTemp    Target temperature in Kelvins.
   * @param collisionFreq Collision frequency in 1/psec.
   */
  public void addAndersenThermostatForce(double targetTemp, double collisionFreq) {
    if (andersenThermostat != null) {
      removeForce(andersenThermostat);
      andersenThermostat = null;
    }
    andersenThermostat = new AndersenThermostat(targetTemp, collisionFreq);
    addForce(andersenThermostat);
    logger.info("\n Adding an Andersen thermostat");
    logger.info(format("  Target Temperature:   %6.2f (K)", targetTemp));
    logger.info(format("  Collision Frequency:  %6.2f (1/psec)", collisionFreq));
  }

  /**
   * Adds a force that removes center-of-mass motion.
   */
  public void addCOMMRemoverForce() {
    if (cmMotionRemover != null) {
      removeForce(cmMotionRemover);
      cmMotionRemover = null;
    }
    int frequency = forceField.getInteger("REMOVE-COM-FREQUENCY", 100);
    cmMotionRemover = new CMMotionRemover(frequency);
    int forceGroup = forceField.getInteger("COMM_FORCE_GROUP", 1);
    cmMotionRemover.setForceGroup(forceGroup);
    addForce(cmMotionRemover);
    logger.info("\n Added a Center of Mass Motion Remover");
    logger.info(format("  Frequency:            %6d", frequency));
    logger.info(format("  Force Group:          %6d", forceGroup));
  }

  /**
   * Add a Monte Carlo Barostat to the system.
   *
   * @param targetPressure The target pressure (in atm).
   * @param targetTemp     The target temperature.
   * @param frequency      The frequency to apply the barostat.
   */
  public void addMonteCarloBarostatForce(double targetPressure, double targetTemp, int frequency) {
    if (monteCarloBarostat != null) {
      removeForce(monteCarloBarostat);
      monteCarloBarostat = null;
    }
    double pressureInBar = targetPressure * Constants.ATM_TO_BAR;
    monteCarloBarostat = new MonteCarloBarostat(pressureInBar, targetTemp, frequency);
    CompositeConfiguration properties = openMMEnergy.getMolecularAssembly().getProperties();
    if (properties.containsKey("barostat-seed")) {
      int randomSeed = properties.getInt("barostat-seed", 0);
      logger.info(format(" Setting random seed %d for Monte Carlo Barostat", randomSeed));
      monteCarloBarostat.setRandomNumberSeed(randomSeed);
    }
    addForce(monteCarloBarostat);
    logger.info("\n Added a Monte Carlo Barostat");
    logger.info(format("  Target Pressure:      %6.2f (atm)", targetPressure));
    logger.info(format("  Target Temperature:   %6.2f (K)", targetTemp));
    logger.info(format("  MC Move Frequency:    %6d", frequency));
  }

  /**
   * Calculate the number of degrees of freedom.
   *
   * @return Number of degrees of freedom.
   */
  public int calculateDegreesOfFreedom() {
    // Begin from the 3 times the number of active atoms.
    int dof = openMMEnergy.getNumberOfVariables();
    // Remove OpenMM constraints.
    dof = dof - getNumConstraints();
    // Remove center of mass motion.
    if (cmMotionRemover != null) {
      dof -= 3;
    }
    return dof;
  }

  public double getTemperature(double kineticEnergy) {
    double dof = calculateDegreesOfFreedom();
    return 2.0 * kineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);
  }

  /**
   * Get the ForceField in use.
   *
   * @return ForceField instance.
   */
  public ForceField getForceField() {
    return forceField;
  }

  /**
   * Get the Crystal instance.
   *
   * @return the Crystal instance.
   */
  public Crystal getCrystal() {
    return openMMEnergy.getCrystal();
  }

  /**
   * Get the number of variables.
   *
   * @return the number of variables.
   */
  public int getNumberOfVariables() {
    return openMMEnergy.getNumberOfVariables();
  }

  /**
   * Destroy the system.
   */
  public void free() {
    if (getPointer() != null) {
      logger.fine(" Free OpenMM system.");
      destroy();
      logger.fine(" Free OpenMM system completed.");
    }
  }

  public void setUpdateBondedTerms(boolean updateBondedTerms) {
    this.updateBondedTerms = updateBondedTerms;
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
  protected void setDefaultPeriodicBoxVectors() {
    Crystal crystal = openMMEnergy.getCrystal();
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
   * Update parameters if the Use flags and/or Lambda value has changed.
   *
   * @param atoms Atoms in this list are considered.
   */
  public void updateParameters(@Nullable Atom[] atoms) {
    VanDerWaals vanDerWaals = openMMEnergy.getVdwNode();
    if (vanDerWaals != null) {
      boolean vdwLambdaTerm = vanDerWaals.getLambdaTerm();
      if (vdwLambdaTerm) {
        double lambdaVDW = vanDerWaals.getLambda();
        if (fixedChargeNonBondedForce != null) {
          if (!softcoreCreated) {
            fixedChargeAlchemicalForces = new FixedChargeAlchemicalForces(openMMEnergy, fixedChargeNonBondedForce);
            addForce(fixedChargeAlchemicalForces.getFixedChargeSoftcoreForce());
            addForce(fixedChargeAlchemicalForces.getAlchemicalAlchemicalStericsForce());
            addForce(fixedChargeAlchemicalForces.getNonAlchemicalAlchemicalStericsForce());
            // Re-initialize the context.
            openMMEnergy.getContext().reinitialize(OpenMM_True);
            softcoreCreated = true;
          }
          // Update the lambda value.
          openMMEnergy.getContext().setParameter("vdw_lambda", lambdaVDW);
        } else if (amoebaVDWForce != null) {
          // Update the lambda value.
          openMMEnergy.getContext().setParameter("AmoebaVdwLambda", lambdaVDW);
          if (softcoreCreated) {
            ParticleMeshEwald pme = openMMEnergy.getPmeNode();
            // Avoid any updateParametersInContext calls if vdwLambdaTerm is true, but not other alchemical terms.
            if (pme == null || !pme.getLambdaTerm()) {
              return;
            }
          } else {
            softcoreCreated = true;
          }
        }
      }
    }

    // Note Stretch-Torsion and Angle-Torsion terms (for nucleic acids)
    // and Torsion-Torsion terms (for protein backbones) are not updated yet.
    if (updateBondedTerms) {
      if (bondForce != null) {
        bondForce.updateForce(openMMEnergy);
      }
      if (angleForce != null) {
        angleForce.updateForce(openMMEnergy);
      }
      if (inPlaneAngleForce != null) {
        inPlaneAngleForce.updateForce(openMMEnergy);
      }
      if (stretchBendForce != null) {
        stretchBendForce.updateForce(openMMEnergy);
      }
      if (ureyBradleyForce != null) {
        ureyBradleyForce.updateForce(openMMEnergy);
      }
      if (outOfPlaneBendForce != null) {
        outOfPlaneBendForce.updateForce(openMMEnergy);
      }
      if (piOrbitalTorsionForce != null) {
        piOrbitalTorsionForce.updateForce(openMMEnergy);
      }
      if (torsionForce != null) {
        torsionForce.updateForce(openMMEnergy);
      }
      if (improperTorsionForce != null) {
        improperTorsionForce.updateForce(openMMEnergy);
      }
    }

    if (restrainTorsionsForce != null) {
      restrainTorsionsForce.updateForce(openMMEnergy);
    }

    if (atoms == null || atoms.length == 0) {
      return;
    }

    // Update fixed charge non-bonded parameters.
    if (fixedChargeNonBondedForce != null) {
      // Need to pass all atoms due to non-bonded exceptions.
      fixedChargeNonBondedForce.updateForce(this.atoms, openMMEnergy);
    }

    // Update fixed charge GB parameters.
    if (fixedChargeGBForce != null) {
      fixedChargeGBForce.updateForce(atoms, openMMEnergy);
    }

    // Update AMOEBA vdW parameters.
    if (amoebaVDWForce != null) {
      amoebaVDWForce.updateForce(atoms, openMMEnergy);
    }

    // Update AMOEBA polarizable multipole parameters.
    if (amoebaMultipoleForce != null) {
      amoebaMultipoleForce.updateForce(atoms, openMMEnergy);
    }

    // Update GK force.
    if (amoebaGeneralizedKirkwoodForce != null) {
      amoebaGeneralizedKirkwoodForce.updateForce(atoms, openMMEnergy);
    }

    // Update WCA Force.
    if (amoebaWcaDispersionForce != null) {
      amoebaWcaDispersionForce.updateForce(atoms, openMMEnergy);
    }

    // Update WCA Force.
    if (amoebaGKCavitationForce != null) {
      amoebaGKCavitationForce.updateForce(atoms, openMMEnergy);
    }
  }

  /**
   * Adds atoms from the molecular assembly to the OpenMM System and reports to the user the number
   * of particles added.
   */
  protected void addAtoms() throws Exception {
    for (Atom atom : atoms) {
      double mass = atom.getMass();
      if (mass < 0.0) {
        throw new Exception(" Atom with mass less than 0.");
      }
      if (mass == 0.0) {
        logger.info(format(" Atom %s has zero mass.", atom));
      }
      addParticle(mass);
    }
  }

  /**
   * This method sets the mass of inactive atoms to zero.
   *
   * @return Returns true if any inactive atoms were found and set to zero mass.
   */
  public boolean updateAtomMass() {
    int index = 0;
    int inactiveCount = 0;
    for (Atom atom : atoms) {
      double mass = 0.0;
      if (atom.isActive()) {
        mass = atom.getMass();
      } else {
        inactiveCount++;
      }
      setParticleMass(index++, mass);
    }
    if (inactiveCount > 0) {
      logger.fine(format(" Inactive atoms (%d) set to zero mass.", inactiveCount));
      return true;
    }
    return false;
  }

  public boolean hasAmoebaCavitationForce() {
    return amoebaGKCavitationForce != null;
  }

  /**
   * Add a constraint to every bond.
   */
  protected void addUpBondConstraints() {
    BondPotentialEnergy bondPotentialEnergy = openMMEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return;
    }
    Bond[] bonds = bondPotentialEnergy.getBondArray();
    logger.info(" Adding constraints for all bonds.");
    for (Bond bond : bonds) {
      Atom atom1 = bond.getAtom(0);
      Atom atom2 = bond.getAtom(1);
      int iAtom1 = atom1.getXyzIndex() - 1;
      int iAtom2 = atom2.getXyzIndex() - 1;
      addConstraint(iAtom1, iAtom2, bond.bondType.distance * OpenMM_NmPerAngstrom);
    }
  }

  /**
   * Add a constraint to every bond that includes a hydrogen atom.
   */
  protected void addHydrogenConstraints() {
    BondPotentialEnergy bondPotentialEnergy = openMMEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return;
    }
    Bond[] bonds = bondPotentialEnergy.getBondArray();
    logger.info(" Adding constraints for hydrogen bonds.");
    for (Bond bond : bonds) {
      Atom atom1 = bond.getAtom(0);
      Atom atom2 = bond.getAtom(1);
      if (atom1.isHydrogen() || atom2.isHydrogen()) {
        BondType bondType = bond.bondType;
        int iAtom1 = atom1.getXyzIndex() - 1;
        int iAtom2 = atom2.getXyzIndex() - 1;
        addConstraint(iAtom1, iAtom2, bondType.distance * OpenMM_NmPerAngstrom);
      }
    }
  }

  /**
   * Add a constraint to every angle that includes two hydrogen atoms.
   */
  protected void setUpHydrogenAngleConstraints() {
    AnglePotentialEnergy anglePotentialEnergy = openMMEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return;
    }
    Angle[] angles = anglePotentialEnergy.getAngleArray();
    logger.info(" Adding hydrogen angle constraints.");
    for (Angle angle : angles) {
      if (isHydrogenAngle(angle)) {
        Atom atom1 = angle.getAtom(0);
        Atom atom3 = angle.getAtom(2);

        // Calculate a "false bond" length between atoms 1 and 3 to constrain the angle using the
        // law of cosines.
        Bond bond1 = angle.getBond(0);
        double distance1 = bond1.bondType.distance;

        Bond bond2 = angle.getBond(1);
        double distance2 = bond2.bondType.distance;

        // Equilibrium angle value in degrees.
        double angleVal = angle.angleType.angle[angle.nh];

        // Law of cosines.
        double falseBondLength = sqrt(
            distance1 * distance1 + distance2 * distance2 - 2.0 * distance1 * distance2 * cos(
                toRadians(angleVal)));

        int iAtom1 = atom1.getXyzIndex() - 1;
        int iAtom3 = atom3.getXyzIndex() - 1;
        addConstraint(iAtom1, iAtom3, falseBondLength * OpenMM_NmPerAngstrom);
      }
    }
  }

  /**
   * Check to see if an angle is a hydrogen angle. This method only returns true for hydrogen
   * angles that are less than 160 degrees.
   *
   * @param angle Angle to check.
   * @return boolean indicating whether an angle is a hydrogen angle that is less than 160 degrees.
   */
  protected boolean isHydrogenAngle(Angle angle) {
    if (angle.containsHydrogen()) {
      // Equilibrium angle value in degrees.
      double angleVal = angle.angleType.angle[angle.nh];
      // Make sure angle is less than 160 degrees because greater than 160 degrees will not be
      // constrained
      // well using the law of cosines.
      if (angleVal < 160.0) {
        Atom atom1 = angle.getAtom(0);
        Atom atom2 = angle.getAtom(1);
        Atom atom3 = angle.getAtom(2);
        // Setting constraints only on angles where atom1 or atom3 is a hydrogen while atom2 is
        // not a hydrogen.
        return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
      }
    }
    return false;
  }
}
