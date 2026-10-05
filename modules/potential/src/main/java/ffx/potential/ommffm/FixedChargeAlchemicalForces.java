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
import ffx.openmm.ffm.CustomBondForce;
import ffx.openmm.ffm.CustomNonbondedForce;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.IntSet;
import ffx.openmm.ffm.NonbondedForce;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.ForceField;

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static ffx.potential.parameters.VDWType.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.parameters.VDWType.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.parameters.VDWType.VDW_TYPE.LENNARD_JONES;
import static java.lang.String.format;

/**
 * Fixed Charge Alchemical Forces backed by FFM {@link CustomNonbondedForce} and {@link CustomBondForce}.
 * <p>
 * 1. Handle interactions between non-alchemical atoms with our default OpenMM NonBondedForce.
 * Note that alchemical atoms must have eps=0 to turn them off in this force.
 * <p>
 * 2. Handle interactions between alchemical atoms and mixed non-alchemical with alchemical
 * interactions with an OpenMM CustomNonBondedForce.
 */
public class FixedChargeAlchemicalForces {

  private static final Logger logger = Logger.getLogger(FixedChargeAlchemicalForces.class.getName());

  private CustomNonbondedForce fixedChargeSoftcoreForce;
  private CustomBondForce alchemicalAlchemicalStericsForce;
  private CustomBondForce nonAlchemicalAlchemicalStericsForce;

  public FixedChargeAlchemicalForces(OpenMMEnergy openMMEnergy,
                                     FixedChargeNonbondedForce fixedChargeNonBondedForce) {

    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      return;
    }

    /*
     Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
     for epsilon is supported.
    */
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    if (vdwForm.vdwType != LENNARD_JONES || vdwForm.radiusRule != ARITHMETIC
        || vdwForm.epsilonRule != GEOMETRIC) {
      logger.info(format(" VDW Type:         %s", vdwForm.vdwType));
      logger.info(format(" VDW Radius Rule:  %s", vdwForm.radiusRule));
      logger.info(format(" VDW Epsilon Rule: %s", vdwForm.epsilonRule));
      logger.log(Level.SEVERE, " Unsupported van der Waals functional form.");
      return;
    }

    // Sterics mixing rules.
    String stericsMixingRules = " epsilon = sqrt(epsilon1*epsilon2);";
    stericsMixingRules += " rmin = 0.5 * (sigma1 + sigma2) * 1.122462048309372981;";

    // Softcore Lennard-Jones, with a form equivalent to that used in FFX VanDerWaals class.
    String stericsEnergyExpression = "(vdw_lambda^beta)*epsilon*x*(x-2.0);";
    // Effective softcore distance for sterics.
    stericsEnergyExpression += " x = 1.0 / (alpha*(1.0-vdw_lambda)^2.0 + (r/rmin)^6.0);";
    // Define energy expression for sterics.
    String energyExpression = stericsEnergyExpression + stericsMixingRules;

    fixedChargeSoftcoreForce = new CustomNonbondedForce(energyExpression);

    double alpha = vdW.getAlpha();
    double beta = vdW.getBeta();

    fixedChargeSoftcoreForce.addGlobalParameter("vdw_lambda", 1.0);
    fixedChargeSoftcoreForce.addGlobalParameter("alpha", alpha);
    fixedChargeSoftcoreForce.addGlobalParameter("beta", beta);
    fixedChargeSoftcoreForce.addPerParticleParameter("sigma");
    fixedChargeSoftcoreForce.addPerParticleParameter("epsilon");

    // Add particles.
    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    try (IntSet alchemicalGroup = new IntSet();
         IntSet nonAlchemicalGroup = new IntSet();
         DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Atom atom : atoms) {
        if (atom.applyLambda()) {
          alchemicalGroup.insert(index);
        } else {
          nonAlchemicalGroup.insert(index);
        }

        NonbondedForce.ParticleParameters particleParameters = fixedChargeNonBondedForce.getParticleParameters(index);
        double sigmaValue = particleParameters.sigma();
        double epsValue = particleParameters.epsilon();

        // Handle cases where sigma is 0.0; for example Amber99 tyrosine hydrogen atoms.
        if (sigmaValue == 0.0) {
          sigmaValue = 1.0;
          epsValue = 0.0;
        }

        parameters.append(sigmaValue);
        parameters.append(epsValue);
        fixedChargeSoftcoreForce.addParticle(parameters);
        parameters.resize(0);
        index++;
      }
      fixedChargeSoftcoreForce.addInteractionGroup(alchemicalGroup, alchemicalGroup);
      fixedChargeSoftcoreForce.addInteractionGroup(alchemicalGroup, nonAlchemicalGroup);
    }

    Crystal crystal = openMMEnergy.getCrystal();
    if (crystal.aperiodic()) {
      fixedChargeSoftcoreForce.setNonbondedMethod(0); // NoCutoff
    } else {
      fixedChargeSoftcoreForce.setNonbondedMethod(2); // CutoffPeriodic
      NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
      double off = nonbondedCutoff.off;
      double cut = nonbondedCutoff.cut;
      fixedChargeSoftcoreForce.setCutoffDistance(NM_PER_ANGSTROM * off);
      fixedChargeSoftcoreForce.setUseSwitchingFunction(true);
      if (cut == off) {
        logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
        if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
          logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
          cut = 1E40;
        } else {
          logger.info(format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut, off));
          cut *= 0.99;
        }
      }
      fixedChargeSoftcoreForce.setSwitchingDistance(NM_PER_ANGSTROM * cut);
    }

    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    int forceGroup = forceField.getInteger("VDW_LAMBDA_FORCE_GROUP", 0);
    fixedChargeSoftcoreForce.setForceGroup(forceGroup);

    // Create custom bond forces to handle scaled 1-4 interactions.
    alchemicalAlchemicalStericsForce = new CustomBondForce(stericsEnergyExpression + stericsMixingRules);
    alchemicalAlchemicalStericsForce.addPerBondParameter("rmin");
    alchemicalAlchemicalStericsForce.addPerBondParameter("epsilon");
    alchemicalAlchemicalStericsForce.addGlobalParameter("vdw_lambda", 1.0);
    alchemicalAlchemicalStericsForce.addGlobalParameter("alpha", alpha);
    alchemicalAlchemicalStericsForce.addGlobalParameter("beta", beta);

    nonAlchemicalAlchemicalStericsForce = new CustomBondForce(stericsEnergyExpression + stericsMixingRules);
    nonAlchemicalAlchemicalStericsForce.addPerBondParameter("rmin");
    nonAlchemicalAlchemicalStericsForce.addPerBondParameter("epsilon");
    nonAlchemicalAlchemicalStericsForce.addGlobalParameter("vdw_lambda", 1.0);
    nonAlchemicalAlchemicalStericsForce.addGlobalParameter("alpha", alpha);
    nonAlchemicalAlchemicalStericsForce.addGlobalParameter("beta", beta);

    int range = fixedChargeNonBondedForce.getNumExceptions();
    int[][] torsionMask = vdW.getMask14();

    try (DoubleArray bondParameters = new DoubleArray(0)) {
      for (int i = 0; i < range; i++) {
        NonbondedForce.ExceptionParameters exceptionParameters = fixedChargeNonBondedForce.getExceptionParameters(i);
        int atomi = exceptionParameters.particle1();
        int atomj = exceptionParameters.particle2();
        double sigmaVal = exceptionParameters.sigma();
        double epsVal = exceptionParameters.epsilon();

        // Omit both Exclusions (1-2, 1-3) and Exceptions (scaled 1-4) from the
        // CustomNonbondedForce.
        fixedChargeSoftcoreForce.addExclusion(atomi, atomj);

        // Deal with scaled 1-4 torsions using the CustomBondForce
        int[] maskI = torsionMask[atomi];
        boolean epsException = false;

        for (int mask : maskI) {
          if (mask == atomj) {
            epsException = true;
            break;
          }
        }

        if (epsException) {
          Atom atom1 = atoms[atomi];
          Atom atom2 = atoms[atomj];
          boolean bothAlchemical = false;
          boolean oneAlchemical = false;
          if (atom1.applyLambda() && atom2.applyLambda()) {
            bothAlchemical = true;
          } else if ((atom1.applyLambda() && !atom2.applyLambda()) || (!atom1.applyLambda() && atom2.applyLambda())) {
            oneAlchemical = true;
          }
          if (bothAlchemical || oneAlchemical) {
            bondParameters.append(sigmaVal * 1.122462048309372981);
            bondParameters.append(epsVal);
            if (bothAlchemical) {
              alchemicalAlchemicalStericsForce.addBond(atomi, atomj, bondParameters);
            } else {
              nonAlchemicalAlchemicalStericsForce.addBond(atomi, atomj, bondParameters);
            }
            bondParameters.resize(0);
          }
        }
      }
    }
    alchemicalAlchemicalStericsForce.setForceGroup(forceGroup);
    nonAlchemicalAlchemicalStericsForce.setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Added fixed charge softcore force \t%d", forceGroup));
    logger.log(Level.INFO, format("   Alpha = %8.6f and beta = %8.6f", alpha, beta));
  }

  public CustomNonbondedForce getFixedChargeSoftcoreForce() {
    return fixedChargeSoftcoreForce;
  }

  public CustomBondForce getAlchemicalAlchemicalStericsForce() {
    return alchemicalAlchemicalStericsForce;
  }

  public CustomBondForce getNonAlchemicalAlchemicalStericsForce() {
    return nonAlchemicalAlchemicalStericsForce;
  }

}
