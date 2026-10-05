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

import ffx.openmm.ffm.CustomCompoundBondForce;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.IntArray;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.terms.StretchTorsionPotentialEnergy;

import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static java.lang.String.format;

/**
 * OpenMM Stretch-Torsion Force backed by FFM {@link CustomCompoundBondForce}.
 */
public class StretchTorsionForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(StretchTorsionForce.class.getName());

  /**
   * Create an OpenMM Stretch-Torsion Force.
   *
   * @param stretchTorsionPotentialEnergy The StretchTorsionPotentialEnergy instance that contains the stretch-torsions.
   */
  public StretchTorsionForce(StretchTorsionPotentialEnergy stretchTorsionPotentialEnergy) {
    super(4, StretchTorsion.stretchTorsionForm());
    StretchTorsion[] stretchTorsions = stretchTorsionPotentialEnergy.getStretchTorsionArray();
    addGlobalParameter("phi1", 0);
    addGlobalParameter("phi2", Math.PI);
    addGlobalParameter("phi3", 0);
    for (int m = 1; m < 4; m++) {
      for (int n = 1; n < 4; n++) {
        addPerBondParameter(format("k%d%d", m, n));
      }
    }
    for (int m = 1; m < 4; m++) {
      addPerBondParameter(format("b%d", m));
    }

    final double unitConv = KJ_PER_KCAL / NM_PER_ANGSTROM;

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      for (StretchTorsion stretchTorsion : stretchTorsions) {
        double[] constants = stretchTorsion.getConstants();
        for (int m = 0; m < 3; m++) {
          for (int n = 0; n < 3; n++) {
            int index = (3 * m) + n;
            parameters.append(constants[index] * unitConv);
          }
        }
        parameters.append(stretchTorsion.bondType1.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType2.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType3.distance * NM_PER_ANGSTROM);

        Atom[] atoms = stretchTorsion.getAtomArray(true);
        for (int i = 0; i < 4; i++) {
          particles.append(atoms[i].getArrayIndex());
        }

        addBond(particles, parameters);
        parameters.resize(0);
        particles.resize(0);
      }
    }

    int forceGroup = stretchTorsionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Stretch-Torsions:                  %10d", stretchTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create a Dual Topology OpenMM Stretch-Torsion Force.
   *
   * @param stretchPotentialEnergy   The StretchTorsionPotentialEnergy instance that contains the stretch-torsions.
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public StretchTorsionForce(StretchTorsionPotentialEnergy stretchPotentialEnergy,
                             int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(4, StretchTorsion.stretchTorsionForm());
    StretchTorsion[] stretchTorsions = stretchPotentialEnergy.getStretchTorsionArray();
    addGlobalParameter("phi1", 0);
    addGlobalParameter("phi2", Math.PI);
    addGlobalParameter("phi3", 0);
    for (int m = 1; m < 4; m++) {
      for (int n = 1; n < 4; n++) {
        addPerBondParameter(format("k%d%d", m, n));
      }
    }
    for (int m = 1; m < 4; m++) {
      addPerBondParameter(format("b%d", m));
    }

    final double unitConv = KJ_PER_KCAL / NM_PER_ANGSTROM;
    double scaleDT = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      for (StretchTorsion stretchTorsion : stretchTorsions) {
        double scale = 1.0;
        // Don't apply lambda scale to alchemical stretch-torsion
        if (!stretchTorsion.applyLambda()) {
          scale = scaleDT;
        }
        double[] constants = stretchTorsion.getConstants();
        for (int m = 0; m < 3; m++) {
          for (int n = 0; n < 3; n++) {
            int index = (3 * m) + n;
            parameters.append(constants[index] * unitConv * scale);
          }
        }
        parameters.append(stretchTorsion.bondType1.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType2.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType3.distance * NM_PER_ANGSTROM);

        Atom[] atoms = stretchTorsion.getAtomArray(true);
        for (int i = 0; i < 4; i++) {
          int atomIndex = atoms[i].getArrayIndex();
          atomIndex = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, atomIndex);
          particles.append(atomIndex);
        }

        addBond(particles, parameters);
        parameters.resize(0);
        particles.resize(0);
      }
    }

    int forceGroup = stretchPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Stretch-Torsions:                  %10d", stretchTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Stretch-Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the stretch-torsions.
   * @return A Stretch-Torsion Force, or null if there are no stretch-torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    StretchTorsionPotentialEnergy stretchTorsionPotentialEnergy =
        openMMEnergy.getStretchTorsionPotentialEnergy();
    if (stretchTorsionPotentialEnergy == null) {
      return null;
    }
    return new StretchTorsionForce(stretchTorsionPotentialEnergy);
  }

  /**
   * Convenience method to construct a Dual Topology OpenMM Stretch-Torsion Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Stretch-Torsion Force, or null if there are no stretch-torsions.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    StretchTorsionPotentialEnergy stretchTorsionPotentialEnergy =
        forceFieldEnergy.getStretchTorsionPotentialEnergy();
    if (stretchTorsionPotentialEnergy == null) {
      return null;
    }
    return new StretchTorsionForce(stretchTorsionPotentialEnergy, topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the Stretch-Torsion parameters in the OpenMM Context.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the stretch-torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    StretchTorsionPotentialEnergy stretchTorsionPotentialEnergy =
        openMMEnergy.getStretchTorsionPotentialEnergy();
    if (stretchTorsionPotentialEnergy == null) {
      return;
    }
    StretchTorsion[] stretchTorsions = stretchTorsionPotentialEnergy.getStretchTorsionArray();
    final double unitConv = KJ_PER_KCAL / NM_PER_ANGSTROM;

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      int bondIndex = 0;
      for (StretchTorsion stretchTorsion : stretchTorsions) {
        double[] constants = stretchTorsion.getConstants();
        for (int m = 0; m < 3; m++) {
          for (int n = 0; n < 3; n++) {
            int index = (3 * m) + n;
            parameters.append(constants[index] * unitConv);
          }
        }
        parameters.append(stretchTorsion.bondType1.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType2.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType3.distance * NM_PER_ANGSTROM);

        Atom[] atoms = stretchTorsion.getAtomArray(true);
        for (int i = 0; i < 4; i++) {
          particles.append(atoms[i].getArrayIndex());
        }

        setBondParameters(bondIndex++, particles, parameters);
        parameters.resize(0);
        particles.resize(0);
      }
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the Dual Topology Stretch-Torsion Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    StretchTorsionPotentialEnergy stretchTorsionPotentialEnergy =
        forceFieldEnergy.getStretchTorsionPotentialEnergy();
    if (stretchTorsionPotentialEnergy == null) {
      return;
    }
    StretchTorsion[] stretchTorsions = stretchTorsionPotentialEnergy.getStretchTorsionArray();
    final double unitConv = KJ_PER_KCAL / NM_PER_ANGSTROM;
    double scaleDT = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      int stIndex = 0;
      for (StretchTorsion stretchTorsion : stretchTorsions) {
        double scale = 1.0;
        // Don't apply lambda scale to alchemical stretch-torsion
        if (!stretchTorsion.applyLambda()) {
          scale = scaleDT;
        }
        double[] constants = stretchTorsion.getConstants();
        for (int m = 0; m < 3; m++) {
          for (int n = 0; n < 3; n++) {
            int index = (3 * m) + n;
            parameters.append(constants[index] * unitConv * scale);
          }
        }
        parameters.append(stretchTorsion.bondType1.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType2.distance * NM_PER_ANGSTROM);
        parameters.append(stretchTorsion.bondType3.distance * NM_PER_ANGSTROM);

        Atom[] atoms = stretchTorsion.getAtomArray(true);
        for (int i = 0; i < 4; i++) {
          int atomIndex = atoms[i].getArrayIndex();
          atomIndex = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, atomIndex);
          particles.append(atomIndex);
        }

        setBondParameters(stIndex++, particles, parameters);
        parameters.resize(0);
        particles.resize(0);
      }
    }

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }
}
