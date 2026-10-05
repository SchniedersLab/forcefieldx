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
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.IntArray;
import ffx.openmm.ffm.amoeba.VdwForce;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWPairType;
import ffx.potential.parameters.VDWType;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static java.lang.Math.sqrt;

/**
 * The AMOEBA vdW Force backed by FFM {@link VdwForce}.
 */
public class AmoebaVdwForce extends VdwForce {

  private static final Logger logger = Logger.getLogger(AmoebaVdwForce.class.getName());

  public static final int NONBONDED_METHOD_NO_CUTOFF = 0;
  public static final int NONBONDED_METHOD_CUTOFF_PERIODIC = 1;

  public static final int ALCHEMICAL_METHOD_NONE = 0;
  public static final int ALCHEMICAL_METHOD_DECOUPLE = 1;
  public static final int ALCHEMICAL_METHOD_ANNIHILATE = 2;

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

    // Configure the Amoeba vdW Force.
    configureForce(openMMEnergy);

    // Add the particles.
    ExtendedSystem extendedSystem = vdW.getExtendedSystem();
    double[] vdwPrefactorAndDerivs = new double[3];

    int[] ired = vdW.getReductionIndex();
    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      VDWType vdwType = atom.getVDWType();
      int atomClass = vdwType.atomClass;
      int type = vdwClassToOpenMMType.get(atomClass);
      boolean isAlchemical = atom.applyLambda();
      double scaleFactor = 1.0;
      if (extendedSystem != null) {
        extendedSystem.getVdwPrefactor(i, vdwPrefactorAndDerivs);
        scaleFactor = vdwPrefactorAndDerivs[0];
      }
      addParticle(ired[i], type, vdwType.reductionFactor, isAlchemical, scaleFactor);
    }

    // Create exclusion lists.
    int[][] bondMask = vdW.getMask12();
    int[][] angleMask = vdW.getMask13();
    try (IntArray exclusions = new IntArray(0)) {
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
    }
  }

  /**
   * The Dual-Topology Amoeba vdW Force constructor.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public AmoebaVdwForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    VanDerWaals vdW = forceFieldEnergy.getVdwNode();
    if (vdW == null) {
      destroy();
      return;
    }

    double scale = sqrt(openMMDualTopologyEnergy.getTopologyScale(topology));

    // Configure the Amoeba vdW Force.
    configureForce(forceFieldEnergy);

    // Add the particles.
    ExtendedSystem extendedSystem = vdW.getExtendedSystem();
    if (extendedSystem != null) {
      logger.severe(" Extended system is not supported for dual-topology simulations.");
    }

    int nAtoms = openMMDualTopologyEnergy.getNumberOfAtoms();
    int[] ired = vdW.getReductionIndex();

    // Add a particle for each atom in the dual topology.
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = openMMDualTopologyEnergy.getDualTopologyAtom(topology, i);
      int top = atom.getTopologyIndex();
      if (top == topology) {
        int index = atom.getArrayIndex();
        int ir = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, ired[index]);
        VDWType vdwType = atom.getVDWType();
        int atomClass = vdwType.atomClass;
        int type = vdwClassToOpenMMType.get(atomClass);
        boolean isAlchemical = atom.applyLambda();
        addParticle(ir, type, vdwType.reductionFactor, isAlchemical, scale);
      } else {
        // Add a fake particle for an atom not in this topology.
        int index = atom.getTopologyAtomIndex();
        int type = vdwClassToOpenMMType.get(vdWClassForNoInteraction);
        boolean isAlchemical = true;
        double scaleFactor = 0.0;
        addParticle(index, type, 1.0, isAlchemical, scaleFactor);
      }
    }

    // Create exclusion lists only for the atoms in this topology.
    int[][] bondMask = vdW.getMask12();
    int[][] angleMask = vdW.getMask13();
    try (IntArray exclusions = new IntArray(0)) {
      for (int index = 0; index < nAtoms; index++) {
        Atom atom = openMMDualTopologyEnergy.getDualTopologyAtom(topology, index);
        if (atom.getTopologyIndex() != topology) {
          continue; // Skip atoms not in this topology.
        }
        exclusions.append(index);

        // Exclude 1-2 interactions.
        final int[] bondMaski = bondMask[atom.getArrayIndex()];
        for (int value : bondMaski) {
          value = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, value);
          exclusions.append(value);
        }

        // Exclude 1-3 interactions.
        final int[] angleMaski = angleMask[atom.getArrayIndex()];
        for (int value : angleMaski) {
          value = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, value);
          exclusions.append(value);
        }
        setParticleExclusions(index, exclusions);
        exclusions.resize(0);
      }
    }
  }

  /**
   * Configuration of the AMOEBA vdW force that is used for single-topology simulations.
   *
   * @param forceFieldEnergy The ForceFieldEnergy instance that contains the vdW parameters.
   */
  private void configureForce(ForceFieldEnergy forceFieldEnergy) {
    VanDerWaals vdW = forceFieldEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();

    double radScale = 1.0;
    if (vdwForm.radiusSize == VDWType.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    ForceField forceField = forceFieldEnergy.getMolecularAssembly().getForceField();
    Map<String, VDWType> vdwTypes = forceField.getVDWTypes();
    for (VDWType vdwType : vdwTypes.values()) {
      int atomClass = vdwType.atomClass;
      if (!vdwClassToOpenMMType.containsKey(atomClass)) {
        double eps = KJ_PER_KCAL * vdwType.wellDepth;
        double rad = NM_PER_ANGSTROM * vdwType.radius * radScale;
        // OpenMM AMOEBA vdW class does not allow a radius of 0.
        if (rad == 0) {
          rad = NM_PER_ANGSTROM * radScale;
        }
        int type = addParticleType(rad, eps);
        vdwClassToOpenMMType.put(atomClass, type);
        if (atomClass <= vdWClassForNoInteraction) {
          vdWClassForNoInteraction = atomClass - 1;
        }
      }
    }

    // Add a special vdW type for zero vdW energy and forces (e.g. to support the FFX "use" flag).
    int type = addParticleType(NM_PER_ANGSTROM, 0.0);
    vdwClassToOpenMMType.put(vdWClassForNoInteraction, type);

    Map<String, VDWPairType> vdwPairTypeMap = forceField.getVDWPairTypes();
    for (VDWPairType vdwPairType : vdwPairTypeMap.values()) {
      int c1 = vdwPairType.atomClasses[0];
      int c2 = vdwPairType.atomClasses[1];
      int type1 = vdwClassToOpenMMType.get(c1);
      int type2 = vdwClassToOpenMMType.get(c2);
      double rMin = vdwPairType.radius * NM_PER_ANGSTROM;
      double eps = vdwPairType.wellDepth * KJ_PER_KCAL;
      addTypePair(type1, type2, rMin, eps);
      addTypePair(type2, type1, rMin, eps);
    }

    // Set the nonbonded cutoff and dispersion correction.
    NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
    setCutoffDistance(nonbondedCutoff.off * NM_PER_ANGSTROM);
    setUseDispersionCorrection(vdW.getDoLongRangeCorrection());

    // Set the nonbonded method based on the crystal periodicity.
    Crystal crystal = forceFieldEnergy.getCrystal();
    if (crystal.aperiodic()) {
      setNonbondedMethod(NONBONDED_METHOD_NO_CUTOFF);
    } else {
      setNonbondedMethod(NONBONDED_METHOD_CUTOFF_PERIODIC);
    }

    // Set the alchemical method if the vdW force has a lambda term.
    if (vdW.getLambdaTerm()) {
      boolean annihilate = vdW.getIntramolecularSoftcore();
      if (annihilate) {
        setAlchemicalMethod(ALCHEMICAL_METHOD_ANNIHILATE);
      } else {
        setAlchemicalMethod(ALCHEMICAL_METHOD_DECOUPLE);
      }
      setSoftcoreAlpha(vdW.getAlpha());
      setSoftcorePower((int) vdW.getBeta());
    }

    int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
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
   * Convenience method to construct a Dual Topology AMOEBA vdW force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An AMOEBA vdW Force, or null if there are no vdW interactions.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    VanDerWaals vdW = forceFieldEnergy.getVdwNode();
    if (vdW == null) {
      return null;
    }
    return new AmoebaVdwForce(topology, openMMDualTopologyEnergy);
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
      int index = atom.getArrayIndex();
      VDWType vdwType = atom.getVDWType();

      // Get the OpenMM index for this vdW type.
      int type = vdwClassToOpenMMType.get(vdwType.atomClass);
      if (!atom.getUse()) {
        // Get the OpenMM index for a special vdW type that has no interactions.
        type = vdwClassToOpenMMType.get(vdWClassForNoInteraction);
      }
      boolean isAlchemical = atom.applyLambda();
      double eps = KJ_PER_KCAL * vdwType.wellDepth;
      double rad = NM_PER_ANGSTROM * vdwType.radius * radScale;

      double scaleFactor = 1.0;
      if (extendedSystem != null) {
        extendedSystem.getVdwPrefactor(index, vdwPrefactorAndDerivs);
        scaleFactor = vdwPrefactorAndDerivs[0];
      }

      setParticleParameters(index, ired[index], rad, eps, vdwType.reductionFactor, isAlchemical, type, scaleFactor);
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update existing AMOEBA vdW force for the Dual-Topology OpenMM System.
   *
   * @param atoms                    The atoms to update.
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(Atom[] atoms, int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    double scale = sqrt(openMMDualTopologyEnergy.getTopologyScale(topology));

    VanDerWaals vdW = forceFieldEnergy.getVdwNode();
    VanDerWaalsForm vdwForm = vdW.getVDWForm();
    double radScale = 1.0;
    if (vdwForm.radiusSize == VDWType.RADIUS_SIZE.DIAMETER) {
      radScale = 0.5;
    }

    // Remap the reduction index to the dual topology index.
    int[] ired = vdW.getReductionIndex();

    for (Atom atom : atoms) {
      if (atom.getTopologyIndex() != topology) {
        // Skip atoms not in this topology.
        continue;
      }

      // Get the dual topology index for this atom.
      int indexDT = atom.getTopologyAtomIndex();
      // Map the reduction index for this atom from single to dual topology.
      int ir = ired[atom.getArrayIndex()];
      ir = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, ir);

      VDWType vdwType = atom.getVDWType();
      // Get the OpenMM index for this vdW type.
      int type = vdwClassToOpenMMType.get(vdwType.atomClass);
      if (!atom.getUse()) {
        // Get the OpenMM index for a special vdW type that has no interactions.
        type = vdwClassToOpenMMType.get(vdWClassForNoInteraction);
      }
      boolean isAlchemical = atom.applyLambda();
      double eps = KJ_PER_KCAL * vdwType.wellDepth;
      double rad = NM_PER_ANGSTROM * vdwType.radius * radScale;
      setParticleParameters(indexDT, ir, rad, eps, vdwType.reductionFactor, isAlchemical, type, scale);
    }
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
