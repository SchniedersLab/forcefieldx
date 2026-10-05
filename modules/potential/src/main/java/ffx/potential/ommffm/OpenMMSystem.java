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
package ffx.potential.ommffm;

import ffx.crystal.Crystal;
import ffx.numerics.Potential;
import ffx.openmm.ffm.AndersenThermostat;
import ffx.openmm.ffm.CMMotionRemover;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.MonteCarloBarostat;
import ffx.openmm.ffm.System;
import ffx.openmm.ffm.Vec3;
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
import ffx.potential.terms.BondPotentialEnergy;
import ffx.utilities.Constants;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static ffx.potential.parameters.VDWType.VDW_TYPE.LENNARD_JONES;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;
import static java.lang.Math.toRadians;
import static java.lang.String.format;

/**
 * Create and manage an OpenMM System using the FFM backend.
 */
public class OpenMMSystem extends System {

  private static final Logger logger = Logger.getLogger(OpenMMSystem.class.getName());

  private final OpenMMEnergy openMMEnergy;
  protected ForceField forceField;
  protected Atom[] atoms;
  protected boolean updateBondedTerms = false;
  protected BondForce bondForce = null;
  protected AngleForce angleForce = null;
  protected InPlaneAngleForce inPlaneAngleForce = null;
  protected StretchBendForce stretchBendForce = null;
  protected UreyBradleyForce ureyBradleyForce = null;
  protected OutOfPlaneBendForce outOfPlaneBendForce = null;
  protected PiOrbitalTorsionForce piOrbitalTorsionForce = null;
  protected TorsionForce torsionForce = null;
  protected ImproperTorsionForce improperTorsionForce = null;
  protected StretchTorsionForce stretchTorsionForce = null;
  protected AngleTorsionForce angleTorsionForce = null;
  protected AmoebaTorsionTorsionForce amoebaTorsionTorsionForce = null;
  protected RestrainPositionsForce restrainPositionsForce = null;
  protected RestrainTorsionsForce restrainTorsionsForce = null;
  protected RestrainGroupsForce restrainGroupsForce = null;
  protected FixedChargeGBForce fixedChargeGBForce = null;
  protected FixedChargeNonbondedForce fixedChargeNonBondedForce = null;
  protected FixedChargeAlchemicalForces fixedChargeAlchemicalForces = null;
  protected AmoebaVdwForce amoebaVDWForce = null;
  protected AmoebaMultipoleForce amoebaMultipoleForce = null;
  protected AmoebaGeneralizedKirkwoodForce amoebaGeneralizedKirkwoodForce = null;
  protected AmoebaWcaDispersionForce amoebaWcaDispersionForce = null;
  protected AmoebaGKCavitationForce amoebaGKCavitationForce = null;
  protected AndersenThermostat andersenThermostat = null;
  protected MonteCarloBarostat monteCarloBarostat = null;
  protected CMMotionRemover cmMotionRemover = null;
  private boolean softcoreCreated = false;

  public OpenMMSystem() {
    openMMEnergy = null;
    forceField = null;
    atoms = null;
  }

  public OpenMMSystem(OpenMMEnergy openMMEnergy) {
    this.openMMEnergy = openMMEnergy;

    MolecularAssembly molecularAssembly = openMMEnergy.getMolecularAssembly();
    forceField = molecularAssembly.getForceField();
    atoms = molecularAssembly.getAtomArray();

    try {
      addAtoms();
    } catch (Exception e) {
      logger.severe(" Atom without mass encountered.");
    }

    logger.info(format("\n OpenMM system created with %d atoms.", atoms.length));
  }

  public Potential getPotential() {
    return openMMEnergy;
  }

  public Atom[] getAtoms() {
    return atoms;
  }

  public void addForces() {
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

    if (rigidBonds) {
      logger.info(" Not creating AmoebaBondForce because bonds are constrained.");
    } else {
      bondForce = (BondForce) BondForce.constructForce(openMMEnergy);
      if (bondForce != null) {
        addForce(bondForce);
      }
    }

    // Add Angle Force.
    angleForce = (AngleForce) AngleForce.constructForce(openMMEnergy);
    if (angleForce != null) {
      addForce(angleForce);
    }

    // Add In-Plane Angle Force.
    inPlaneAngleForce = (InPlaneAngleForce) InPlaneAngleForce.constructForce(openMMEnergy);
    if (inPlaneAngleForce != null) {
      addForce(inPlaneAngleForce);
    }

    // Add Stretch-Bend Force.
    stretchBendForce = (StretchBendForce) StretchBendForce.constructForce(openMMEnergy);
    if (stretchBendForce != null) {
      addForce(stretchBendForce);
    }

    // Add Urey-Bradley Force.
    ureyBradleyForce = (UreyBradleyForce) UreyBradleyForce.constructForce(openMMEnergy);
    if (ureyBradleyForce != null) {
      addForce(ureyBradleyForce);
    }

    // Out-of-Plane Bend Force.
    outOfPlaneBendForce = (OutOfPlaneBendForce) OutOfPlaneBendForce.constructForce(openMMEnergy);
    if (outOfPlaneBendForce != null) {
      addForce(outOfPlaneBendForce);
    }

    // Add Pi-Torsion Force.
    piOrbitalTorsionForce = (PiOrbitalTorsionForce) PiOrbitalTorsionForce.constructForce(openMMEnergy);
    if (piOrbitalTorsionForce != null) {
      addForce(piOrbitalTorsionForce);
    }

    // Add Torsion Force.
    torsionForce = (TorsionForce) TorsionForce.constructForce(openMMEnergy);
    if (torsionForce != null) {
      addForce(torsionForce);
    }

    // Add Improper Torsion Force.
    improperTorsionForce = (ImproperTorsionForce) ImproperTorsionForce.constructForce(openMMEnergy);
    if (improperTorsionForce != null) {
      addForce(improperTorsionForce);
    }

    // Add Stretch-Torsion coupling terms.
    stretchTorsionForce = (StretchTorsionForce) StretchTorsionForce.constructForce(openMMEnergy);
    if (stretchTorsionForce != null) {
      addForce(stretchTorsionForce);
    }

    // Add Angle-Torsion coupling terms.
    angleTorsionForce = (AngleTorsionForce) AngleTorsionForce.constructForce(openMMEnergy);
    if (angleTorsionForce != null) {
      addForce(angleTorsionForce);
    }

    // Add Torsion-Torsion Force.
    amoebaTorsionTorsionForce = (AmoebaTorsionTorsionForce) AmoebaTorsionTorsionForce.constructForce(openMMEnergy);
    if (amoebaTorsionTorsionForce != null) {
      addForce(amoebaTorsionTorsionForce);
    }

    if (openMMEnergy.getRestrainMode() == ForceFieldEnergy.RestrainMode.ENERGY) {
      // Add Restrain Positions force.
      restrainPositionsForce = (RestrainPositionsForce) RestrainPositionsForce.constructForce(openMMEnergy);
      if (restrainPositionsForce != null) {
        addForce(restrainPositionsForce);
      }

      // Add Restrain Bonds force for each functional form.
      for (BondType.BondFunction function : BondType.BondFunction.values()) {
        RestrainDistanceForce restrainDistanceForce = (RestrainDistanceForce) RestrainDistanceForce.constructForce(function, openMMEnergy);
        if (restrainDistanceForce != null) {
          addForce(restrainDistanceForce);
        }
      }

      // Add Restrain Torsions force.
      restrainTorsionsForce = (RestrainTorsionsForce) RestrainTorsionsForce.constructForce(openMMEnergy);
      if (restrainTorsionsForce != null) {
        addForce(restrainTorsionsForce);
      }
    }

    // Add Restrain Groups force.
    restrainGroupsForce = (RestrainGroupsForce) RestrainGroupsForce.constructForce(openMMEnergy);
    if (restrainGroupsForce != null) {
      addForce(restrainGroupsForce);
    }

    setDefaultPeriodicBoxVectors();

    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW != null) {
      logger.info("\n Non-Bonded Terms");
      VanDerWaalsForm vdwForm = vdW.getVDWForm();
      if (vdwForm.vdwType == LENNARD_JONES) {
        fixedChargeNonBondedForce = (FixedChargeNonbondedForce) FixedChargeNonbondedForce.constructForce(openMMEnergy);
        if (fixedChargeNonBondedForce != null) {
          addForce(fixedChargeNonBondedForce);
          GeneralizedKirkwood gk = openMMEnergy.getGK();
          if (gk != null) {
            fixedChargeGBForce = (FixedChargeGBForce) FixedChargeGBForce.constructForce(openMMEnergy);
            if (fixedChargeGBForce != null) {
              addForce(fixedChargeGBForce);
            }
          }
        }
      } else {
        // Add vdW Force.
        amoebaVDWForce = (AmoebaVdwForce) AmoebaVdwForce.constructForce(openMMEnergy);
        if (amoebaVDWForce != null) {
          addForce(amoebaVDWForce);
        }

        // Add Multipole Force.
        amoebaMultipoleForce = (AmoebaMultipoleForce) AmoebaMultipoleForce.constructForce(openMMEnergy);
        if (amoebaMultipoleForce != null) {
          addForce(amoebaMultipoleForce);
        }

        // Add Generalized Kirkwood Force.
        GeneralizedKirkwood gk = openMMEnergy.getGK();
        if (gk != null) {
          amoebaGeneralizedKirkwoodForce = (AmoebaGeneralizedKirkwoodForce) AmoebaGeneralizedKirkwoodForce.constructForce(openMMEnergy);
          if (amoebaGeneralizedKirkwoodForce != null) {
            addForce(amoebaGeneralizedKirkwoodForce);
          }

          // Add WCA Dispersion Force.
          DispersionRegion dispersionRegion = gk.getDispersionRegion();
          if (dispersionRegion != null) {
            amoebaWcaDispersionForce = (AmoebaWcaDispersionForce) AmoebaWcaDispersionForce.constructForce(openMMEnergy);
            if (amoebaWcaDispersionForce != null) {
              addForce(amoebaWcaDispersionForce);
            }
          }

          // Add a GaussVol Cavitation Force.
          ChandlerCavitation chandlerCavitation = gk.getChandlerCavitation();
          if (chandlerCavitation != null && chandlerCavitation.getGaussVol() != null) {
            amoebaGKCavitationForce = (AmoebaGKCavitationForce) AmoebaGKCavitationForce.constructForce(openMMEnergy);
            if (amoebaGKCavitationForce != null) {
              addForce(amoebaGKCavitationForce);
            }
          }
        }
      }
    }
  }

  public void removeForce(Force force) {
    if (force != null) {
      removeForce(force.getForceIndex());
    }
  }

  public void addAndersenThermostatForce(double targetTemp) {
    double collisionFreq = forceField.getDouble("COLLISION_FREQ", 1.0);
    addAndersenThermostatForce(targetTemp, collisionFreq);
  }

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
    logger.info("\n Added a center-of-mass motion remover");
    logger.info(format("  Frequency:            %6d", frequency));
  }

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

  public int calculateDegreesOfFreedom() {
    int dof = openMMEnergy.getNumberOfVariables();
    dof = dof - getNumConstraints();
    if (cmMotionRemover != null) {
      dof -= 3;
    }
    return dof;
  }

  public double getTemperature(double kineticEnergy) {
    double dof = calculateDegreesOfFreedom();
    return 2.0 * kineticEnergy * Constants.KCAL_TO_GRAM_ANG2_PER_PS2 / (Constants.kB * dof);
  }

  public ForceField getForceField() {
    return forceField;
  }

  public Crystal getCrystal() {
    return openMMEnergy.getCrystal();
  }

  public int getNumberOfVariables() {
    return openMMEnergy.getNumberOfVariables();
  }

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

  protected void setDefaultPeriodicBoxVectors() {
    Crystal crystal = openMMEnergy.getCrystal();
    if (!crystal.aperiodic()) {
      double[][] Ai = crystal.Ai;
      Vec3 a = new Vec3(Ai[0][0] * NM_PER_ANGSTROM, Ai[0][1] * NM_PER_ANGSTROM, Ai[0][2] * NM_PER_ANGSTROM);
      Vec3 b = new Vec3(Ai[1][0] * NM_PER_ANGSTROM, Ai[1][1] * NM_PER_ANGSTROM, Ai[1][2] * NM_PER_ANGSTROM);
      Vec3 c = new Vec3(Ai[2][0] * NM_PER_ANGSTROM, Ai[2][1] * NM_PER_ANGSTROM, Ai[2][2] * NM_PER_ANGSTROM);
      setDefaultPeriodicBoxVectors(a, b, c);
    }
  }

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
            openMMEnergy.getContext().reinitialize(true);
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
      if (stretchTorsionForce != null) {
        stretchTorsionForce.updateForce(openMMEnergy);
      }
      if (angleTorsionForce != null) {
        angleTorsionForce.updateForce(openMMEnergy);
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

    // Update Cavitation Force.
    if (amoebaGKCavitationForce != null) {
      amoebaGKCavitationForce.updateForce(atoms, openMMEnergy);
    }
  }

  public FixedChargeNonbondedForce getFixedChargeNonbondedForce() {
    return fixedChargeNonBondedForce;
  }

  public FixedChargeGBForce getFixedChargeGBForce() {
    return fixedChargeGBForce;
  }

  public AmoebaVdwForce getAmoebaVdwForce() {
    return amoebaVDWForce;
  }

  public AmoebaMultipoleForce getAmoebaMultipoleForce() {
    return amoebaMultipoleForce;
  }

  public AmoebaGeneralizedKirkwoodForce getAmoebaGeneralizedKirkwoodForce() {
    return amoebaGeneralizedKirkwoodForce;
  }

  public AmoebaWcaDispersionForce getAmoebaWcaDispersionForce() {
    return amoebaWcaDispersionForce;
  }

  public FixedChargeAlchemicalForces getFixedChargeAlchemicalForces() {
    return fixedChargeAlchemicalForces;
  }

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
      addConstraint(iAtom1, iAtom2, bond.bondType.distance * NM_PER_ANGSTROM);
    }
  }

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
        addConstraint(iAtom1, iAtom2, bondType.distance * NM_PER_ANGSTROM);
      }
    }
  }

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
        addConstraint(iAtom1, iAtom3, falseBondLength * NM_PER_ANGSTROM);
      }
    }
  }

  public boolean hasAmoebaCavitationForce() {
    return amoebaGKCavitationForce != null;
  }

  protected boolean isHydrogenAngle(Angle angle) {
    if (angle.containsHydrogen()) {
      double angleVal = angle.angleType.angle[angle.nh];
      if (angleVal < 160.0) {
        Atom atom1 = angle.getAtom(0);
        Atom atom2 = angle.getAtom(1);
        Atom atom3 = angle.getAtom(2);
        return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
      }
    }
    return false;
  }
}
