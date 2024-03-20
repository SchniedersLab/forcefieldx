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
package ffx.potential.openmm;

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import edu.uiowa.jopenmm.OpenMMLibrary;
import ffx.crystal.Crystal;
import ffx.openmm.BondArray;
import ffx.openmm.Force;
import ffx.openmm.NonbondedForce;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.VDWType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_SIZE.RADIUS;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_TYPE.R_MIN;
import static ffx.potential.nonbonded.VanDerWaalsForm.VDW_TYPE.LENNARD_JONES;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.abs;

/**
 * Define a fixed charge non-bonded force.
 * <p>
 * Uses arithmetic mean to define sigma and geometric mean for epsilon.
 */
public class FixedChargeNonbondedForce extends NonbondedForce {

  private static final Logger logger = Logger.getLogger(FixedChargeNonbondedForce.class.getName());

  /**
   * Boolean array, holds charge exclusion list.
   */
  private boolean[] chargeExclusion;
  /**
   * Boolean array, holds van Der Waals exclusion list.
   */
  private boolean[] vdWExclusion;
  /**
   * Double array, holds charge quantity value for exceptions.
   */
  private double[] exceptionChargeProd;
  /**
   * Double array, holds epsilon quantity value for exceptions.
   */
  private double[] exceptionEps;

  public FixedChargeNonbondedForce(OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      // Free the memory created by the OpenMMNonbondedForce constructor.
      destroy();
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

    // OpenMM vdW force requires a diameter (i.e. not radius).
    double radScale = 1.0;
    if (vdwForm.radiusSize == RADIUS) {
      radScale = 2.0;
    }

    // OpenMM vdw force requires atomic sigma values (i.e. not r-min).
    if (vdwForm.radiusType == R_MIN) {
      radScale /= 1.122462048309372981;
    }

    // Add particles.
    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    for (Atom atom : atoms) {
      VDWType vdwType = atom.getVDWType();
      double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
      double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
      double charge = 0.0;
      MultipoleType multipoleType = atom.getMultipoleType();
      if (multipoleType != null && atom.getElectrostatics()) {
        charge = multipoleType.charge;
      }
      addParticle(charge, sigma, eps);
    }

    // Define 1-4 scale factors.
    double lj14Scale = vdwForm.getScale14();
    double coulomb14Scale = 1.0 / 1.2;
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme != null) {
      coulomb14Scale = pme.getScale14();
    }
    Bond[] bonds = openMMEnergy.getBonds();
    if (bonds != null && bonds.length > 0) {
      BondArray bondArray = new BondArray(0);
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getXyzIndex() - 1;
        int i2 = bond.getAtom(1).getXyzIndex() - 1;
        bondArray.append(i1, i2);
      }
      createExceptionsFromBonds(bondArray, coulomb14Scale, lj14Scale);
      bondArray.destroy();

      int num = getNumExceptions();
      chargeExclusion = new boolean[num];
      vdWExclusion = new boolean[num];
      exceptionChargeProd = new double[num];
      exceptionEps = new double[num];
      IntByReference particle1 = new IntByReference();
      IntByReference particle2 = new IntByReference();
      DoubleByReference chargeProd = new DoubleByReference();
      DoubleByReference sigma = new DoubleByReference();
      DoubleByReference eps = new DoubleByReference();
      for (int i = 0; i < num; i++) {
        getExceptionParameters(i, particle1, particle2, chargeProd, sigma, eps);
        if (abs(chargeProd.getValue()) > 0.0) {
          chargeExclusion[i] = false;
          exceptionChargeProd[i] = chargeProd.getValue();
        } else {
          exceptionChargeProd[i] = 0.0;
          chargeExclusion[i] = true;
        }
        if (abs(eps.getValue()) > 0.0) {
          vdWExclusion[i] = false;
          exceptionEps[i] = eps.getValue();
        } else {
          vdWExclusion[i] = true;
          exceptionEps[i] = 0.0;
        }
      }
    }

    Crystal crystal = openMMEnergy.getCrystal();
    if (crystal.aperiodic()) {
      setNonbondedMethod(OpenMMLibrary.OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_NoCutoff);
    } else {
      setNonbondedMethod(OpenMMLibrary.OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_PME);
      if (pme != null) {
        // Units of the Ewald coefficient are A^-1; Multiply by AngstromsPerNM to convert to
        // (Nm^-1).
        double aEwald = OpenMM_AngstromsPerNm * pme.getEwaldCoefficient();
        int nx = pme.getReciprocalSpace().getXDim();
        int ny = pme.getReciprocalSpace().getYDim();
        int nz = pme.getReciprocalSpace().getZDim();
        setPMEParameters(aEwald, nx, ny, nz);
      }

      NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
      double off = nonbondedCutoff.off;
      double cut = nonbondedCutoff.cut;
      setCutoffDistance(OpenMM_NmPerAngstrom * off);
      setUseSwitchingFunction(OpenMM_True);
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
      setSwitchingDistance(OpenMM_NmPerAngstrom * cut);
    }

    setUseDispersionCorrection(OpenMM_False);

    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
    int pmeGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
    if (forceGroup != pmeGroup) {
      logger.severe(format(" ERROR: VDW-FORCE-GROUP is %d while PME-FORCE-GROUP is %d. "
              + "This is invalid for fixed-charge force fields with combined non-bonded forces.",
          forceGroup, pmeGroup));
    }

    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Fixed charge non-bonded force \t%1d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Non-Bonded Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the Non-Bonded Force information.
   * @return An OpenMM Non-Bonded Force, or null if there is no vdW information.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      // Free the memory created by the OpenMMNonbondedForce constructor.
      return null;
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
      return null;
    }

    return new FixedChargeNonbondedForce(openMMEnergy);
  }

  /**
   * Update an existing non-bonded force for the OpenMM System.
   *
   * @param atoms        The Atom array.
   * @param openMMEnergy The OpenMM Energy instance that contains the non-bonded force information.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    // Only 6-12 LJ with arithmetic mean to define sigma and geometric mean for epsilon is
    // supported.
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    if (vdwForm.vdwType != LENNARD_JONES || vdwForm.radiusRule != ARITHMETIC
        || vdwForm.epsilonRule != GEOMETRIC) {
      logger.log(Level.SEVERE, " Unsupported van der Waals functional form.");
      return;
    }

    // OpenMM vdW force requires a diameter (i.e. not radius).
    double radScale = 1.0;
    if (vdwForm.radiusSize == RADIUS) {
      radScale = 2.0;
    }

    // OpenMM vdw force requires atomic sigma values (i.e. not r-min).
    if (vdwForm.radiusType == R_MIN) {
      radScale /= 1.122462048309372981;
    }

    // Update parameters.
    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      boolean applyLambda = atom.applyLambda();

      double charge = Double.MIN_VALUE;
      MultipoleType multipoleType = atom.getMultipoleType();
      if (multipoleType != null && atom.getElectrostatics()) {
        charge = multipoleType.charge;
      }

      VDWType vdwType = atom.getVDWType();
      double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
      double eps = OpenMM_KJPerKcal * vdwType.wellDepth;

      if (applyLambda) {
        OpenMMSystem system = openMMEnergy.getSystem();
        // If we're using vdwLambdaTerm, this atom's vdW interactions are handled by the Custom
        // Non-Bonded force.
        if (system.getVdwLambdaTerm()) {
          eps = 0.0;
        }
        // Always scale the charge by lambdaElec
        charge *= system.getLambdaElec();
      }

      if (!atom.getUse()) {
        eps = 0.0;
        charge = 0.0;
      }

      setParticleParameters(index, charge, sigma, eps);
    }

    // Update Exceptions.
    IntByReference particle1 = new IntByReference();
    IntByReference particle2 = new IntByReference();
    DoubleByReference chargeProd = new DoubleByReference();
    DoubleByReference sigma = new DoubleByReference();
    DoubleByReference eps = new DoubleByReference();

    int numExceptions = getNumExceptions();
    for (int i = 0; i < numExceptions; i++) {

      // Only update exceptions.
      if (chargeExclusion[i] && vdWExclusion[i]) {
        continue;
      }

      getExceptionParameters(i, particle1, particle2, chargeProd, sigma, eps);

      int i1 = particle1.getValue();
      int i2 = particle2.getValue();

      double qq = exceptionChargeProd[i];
      double epsilon = exceptionEps[i];

      Atom atom1 = atoms[i1];
      Atom atom2 = atoms[i2];

      /*
      Note that the minimum epsilon value cannot be zero, or OpenMM may
      report an error that the number of Exceptions has changed.
      */
      double minEpsilon = 1.0e-12;

      OpenMMSystem system = openMMEnergy.getSystem();
      double lambdaValue = system.getLambdaElec();
      boolean vdwLambdaTerm = system.getVdwLambdaTerm();

      if (lambdaValue < minEpsilon) {
        lambdaValue = minEpsilon;
      }

      if (atom1.applyLambda()) {
        qq *= lambdaValue;
        if (vdwLambdaTerm) {
          epsilon = minEpsilon;
        }
      }
      if (atom2.applyLambda()) {
        qq *= lambdaValue;
        if (vdwLambdaTerm) {
          epsilon = minEpsilon;
        }
      }
      if (!atom1.getUse() || !atom2.getUse()) {
        qq = minEpsilon;
        epsilon = minEpsilon;
      }

      setExceptionParameters(i, i1, i2, qq, sigma.getValue(), epsilon);
    }

    // Update the parameters in the OpenMM Context.
    updateParametersInContext(openMMEnergy.getContext());
  }

}
