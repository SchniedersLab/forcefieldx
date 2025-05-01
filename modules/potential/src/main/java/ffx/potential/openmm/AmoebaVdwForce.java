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

import ffx.crystal.Crystal;
import ffx.openmm.Force;
import ffx.openmm.IntArray;
import ffx.openmm.amoeba.VdwForce;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWPairType;
import ffx.potential.parameters.VDWType;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_AlchemicalMethod.OpenMM_AmoebaVdwForce_Annihilate;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_AlchemicalMethod.OpenMM_AmoebaVdwForce_Decouple;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_CutoffPeriodic;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static java.lang.String.format;

/**
 * The Amoeba vdW Force.
 */
public class AmoebaVdwForce extends VdwForce {

  private static final Logger logger = Logger.getLogger(AmoebaVdwForce.class.getName());

  /**
   * The vdW class used to specify no vdW interactions for an atom will be Zero
   * if all atom classes are greater than zero.
   * <p>
   * Otherwise:
   * vdWClassForNoInteraction = min(atomClass) - 1
   */
  private int vdWClassForNoInteraction = 0;

  /**
   * A map from vdW class values to OpenMM vdW types.
   */
  private final Map<Integer, Integer> vdwClassToOpenMMType = new HashMap<>();

  /**
   * The Amoeba vdW Force constructor.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the vdW parameters.
   */
  public AmoebaVdwForce(OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      destroy();
      return;
    }

    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
    Crystal crystal = openMMEnergy.getCrystal();

    double radScale = 1.0;
    if (vdwForm.radiusSize == VDWType.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();

    Map<String, VDWType> vdwTypes = forceField.getVDWTypes();

    for (VDWType vdwType : vdwTypes.values()) {
      int atomClass = vdwType.atomClass;
      if (!vdwClassToOpenMMType.containsKey(atomClass)) {
        double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
        double rad = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
        // OpenMM AMOEBA vdW class does not allow a radius of 0.
        if (rad == 0) {
          rad = OpenMM_NmPerAngstrom * radScale;
        }
        int type = addParticleType(rad, eps);
        vdwClassToOpenMMType.put(atomClass, type);
        if (atomClass <= vdWClassForNoInteraction) {
          vdWClassForNoInteraction = atomClass - 1;
        }
      }
    }

    // Add a special vdW type for zero vdW energy and forces (e.g. to support the FFX "use" flag).
    int type = addParticleType(OpenMM_NmPerAngstrom, 0.0);
    vdwClassToOpenMMType.put(vdWClassForNoInteraction, type);

    Map<String, VDWPairType> vdwPairTypeMap = forceField.getVDWPairTypes();
    for (VDWPairType vdwPairType : vdwPairTypeMap.values()) {
      int c1 = vdwPairType.atomClasses[0];
      int c2 = vdwPairType.atomClasses[1];
      int type1 = vdwClassToOpenMMType.get(c1);
      int type2 = vdwClassToOpenMMType.get(c2);
      double rMin = vdwPairType.radius * OpenMM_NmPerAngstrom;
      double eps = vdwPairType.wellDepth * OpenMM_KJPerKcal;
      addTypePair(type1, type2, rMin, eps);
      addTypePair(type2, type1, rMin, eps);
    }

    ExtendedSystem extendedSystem = vdW.getExtendedSystem();
    double[] vdwPrefactorAndDerivs = new double[3];

    int[] ired = vdW.getReductionIndex();
    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      VDWType vdwType = atom.getVDWType();
      int atomClass = vdwType.atomClass;
      type = vdwClassToOpenMMType.get(atomClass);
      int isAlchemical = atom.applyLambda() ? 1 : 0;
      double scaleFactor = 1.0;
      if (extendedSystem != null) {
        extendedSystem.getVdwPrefactor(i, vdwPrefactorAndDerivs);
        scaleFactor = vdwPrefactorAndDerivs[0];
      }
      addParticle_1(ired[i], type, vdwType.reductionFactor, isAlchemical, scaleFactor);
    }

    setCutoffDistance(nonbondedCutoff.off * OpenMM_NmPerAngstrom);
    if (vdW.getDoLongRangeCorrection()) {
      setUseDispersionCorrection(OpenMM_True);
    } else {
      setUseDispersionCorrection(OpenMM_False);
    }

    if (crystal.aperiodic()) {
      setNonbondedMethod(OpenMM_AmoebaVdwForce_NoCutoff);
    } else {
      setNonbondedMethod(OpenMM_AmoebaVdwForce_CutoffPeriodic);
    }

    if (vdW.getLambdaTerm()) {
      boolean annihilate = vdW.getIntramolecularSoftcore();
      if (annihilate) {
        setAlchemicalMethod(OpenMM_AmoebaVdwForce_Annihilate);
      } else {
        setAlchemicalMethod(OpenMM_AmoebaVdwForce_Decouple);
      }
      setSoftcoreAlpha(vdW.getAlpha());
      setSoftcorePower((int) vdW.getBeta());
    }

    int[][] bondMask = vdW.getMask12();
    int[][] angleMask = vdW.getMask13();

    // Create exclusion lists.
    IntArray exclusions = new IntArray(0);
    for (int i = 0; i < nAtoms; i++) {
      exclusions.append(i);
      final int[] bondMaski = bondMask[i];
      for (int value : bondMaski) {
        exclusions.append(value);
      }
      final int[] angleMaski = angleMask[i];
      for (int value : angleMaski) {
        exclusions.append(value);
      }
      setParticleExclusions(i, exclusions);
      exclusions.resize(0);
    }
    exclusions.destroy();

    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
    setForceGroup(forceGroup);

    logger.log(Level.INFO, vdW.toString());
    logger.log(Level.FINE, format("   Force group:\t\t%d\n", forceGroup));
  }

  /**
   * Convenience method to construct an AMOEBA vdW force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the vdW information.
   * @return An AMOEBA vdW Force, or null if there are no vdW interactions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    if (vdW == null) {
      return null;
    }
    return new AmoebaVdwForce(openMMEnergy);
  }

  /**
   * Update the vdW force.
   *
   * @param atoms        The atoms to update.
   * @param openMMEnergy The OpenMM Energy instance that contains the vdW parameters.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    VanDerWaals vdW = openMMEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    double radScale = 1.0;
    if (vdwForm.radiusSize == VDWType.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    ExtendedSystem extendedSystem = vdW.getExtendedSystem();
    double[] vdwPrefactorAndDerivs = new double[3];

    int[] ired = vdW.getReductionIndex();
    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      VDWType vdwType = atom.getVDWType();

      // Get the OpenMM index for this vdW type.
      int type = vdwClassToOpenMMType.get(vdwType.atomClass);
      if (!atom.getUse()) {
        // Get the OpenMM index for a special vdW type that has no interactions.
        type = vdwClassToOpenMMType.get(vdWClassForNoInteraction);
      }
      int isAlchemical = atom.applyLambda() ? 1 : 0;
      double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
      double rad = OpenMM_NmPerAngstrom * vdwType.radius * radScale;

      double scaleFactor = 1.0;
      if (extendedSystem != null) {
        extendedSystem.getVdwPrefactor(index, vdwPrefactorAndDerivs);
        scaleFactor = vdwPrefactorAndDerivs[0];
      }

      setParticleParameters(index, ired[index], rad, eps, vdwType.reductionFactor, isAlchemical, type, scaleFactor);
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

}
