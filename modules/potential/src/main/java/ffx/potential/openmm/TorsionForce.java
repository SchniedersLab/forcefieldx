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

import ffx.openmm.Force;
import ffx.openmm.PeriodicTorsionForce;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Torsion;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static java.lang.String.format;

public class TorsionForce extends PeriodicTorsionForce {

  private static final Logger logger = Logger.getLogger(TorsionForce.class.getName());

  private final boolean manyBodyTitration;

  public TorsionForce(OpenMMEnergy openMMEnergy) {
    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);
    Torsion[] torsions = openMMEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return;
    }

    for (Torsion torsion : torsions) {
      int a1 = torsion.getAtom(0).getArrayIndex();
      int a2 = torsion.getAtom(1).getArrayIndex();
      int a3 = torsion.getAtom(2).getArrayIndex();
      int a4 = torsion.getAtom(3).getArrayIndex();
      TorsionType torsionType = torsion.torsionType;
      int nTerms = torsionType.phase.length;
      for (int j = 0; j < nTerms; j++) {
        double k = torsion.getTorsionScale() * torsionType.torsionUnit * torsionType.amplitude[j];
        addTorsion(a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree, OpenMM_KJPerKcal * k);
      }
      // Enforce 6-fold torsions since TorsionType instances can have different lengths
      // when side-chain protonation changes.
      if (manyBodyTitration) {
        for (int j = nTerms; j < 6; j++) {
          addTorsion(a1, a2, a3, a4, j + 1, 0.0, 0.0);
        }
      }
    }

    int forceGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Torsions:                          %10d", torsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Torsion Force constructor for Dual Topology.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public TorsionForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ForceField forceField = forceFieldEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);
    Torsion[] torsions = forceFieldEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return;
    }

    if (manyBodyTitration) {
      logger.severe("OpenMM Dual Topology does not suppport many body titration.");
    }

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    for (Torsion torsion : torsions) {
      int a1 = torsion.getAtom(0).getArrayIndex();
      int a2 = torsion.getAtom(1).getArrayIndex();
      int a3 = torsion.getAtom(2).getArrayIndex();
      int a4 = torsion.getAtom(3).getArrayIndex();
      a1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a1);
      a2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a2);
      a3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a3);
      a4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a4);
      TorsionType torsionType = torsion.torsionType;
      int nTerms = torsionType.phase.length;
      for (int j = 0; j < nTerms; j++) {
        double k = torsion.getTorsionScale() * torsionType.torsionUnit * torsionType.amplitude[j];
        // Don't apply lambda scale to alchemical torsion
        if (!torsion.applyLambda()) {
          k = k * scale;
        }
        addTorsion(a1, a2, a3, a4, j + 1,
                torsionType.phase[j] * OpenMM_RadiansPerDegree, OpenMM_KJPerKcal * k);
      }
    }

    int forceGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Torsions:                          %10d", torsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the torsions.
   * @return A Torsion Force, or null if there are no torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    Torsion[] torsions = openMMEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return null;
    }
    return new TorsionForce(openMMEnergy);
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Torsion Force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return A Torsion Force, or null if there are no torsions.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    Torsion[] torsions = forceFieldEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return null;
    }
    return new TorsionForce(topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the Torsion force.
   *
   * @param openMMEnergy The OpenMM Energy that contains the torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    // Check if this system has torsions.
    Torsion[] torsions = openMMEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return;
    }

    int index = 0;
    for (Torsion torsion : torsions) {
      TorsionType torsionType = torsion.torsionType;
      int nTerms = torsionType.phase.length;
      int a1 = torsion.getAtom(0).getArrayIndex();
      int a2 = torsion.getAtom(1).getArrayIndex();
      int a3 = torsion.getAtom(2).getArrayIndex();
      int a4 = torsion.getAtom(3).getArrayIndex();
      for (int j = 0; j < nTerms; j++) {
        double k = torsion.getTorsionScale() * torsionType.torsionUnit * torsionType.amplitude[j];
        setTorsionParameters(index++, a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree, OpenMM_KJPerKcal * k);
      }
      // Enforce 6-fold torsions since TorsionType instances can have different lengths
      // when side-chain protonation changes.
      if (manyBodyTitration) {
        for (int j = nTerms; j < 6; j++) {
          setTorsionParameters(index++, a1, a2, a3, a4, j + 1, 0.0, 0.0);
        }
      }
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the Torsion force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    // Check if this system has torsions.
    Torsion[] torsions = forceFieldEnergy.getTorsions();
    if (torsions == null || torsions.length < 1) {
      return;
    }

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    int index = 0;
    for (Torsion torsion : torsions) {
      TorsionType torsionType = torsion.torsionType;
      int nTerms = torsionType.phase.length;
      int a1 = torsion.getAtom(0).getArrayIndex();
      int a2 = torsion.getAtom(1).getArrayIndex();
      int a3 = torsion.getAtom(2).getArrayIndex();
      int a4 = torsion.getAtom(3).getArrayIndex();
      a1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a1);
      a2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a2);
      a3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a3);
      a4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a4);
      for (int j = 0; j < nTerms; j++) {
        double k = torsion.getTorsionScale() * torsionType.torsionUnit * torsionType.amplitude[j];
        // Don't apply lambda scale to alchemical torsion
        if (!torsion.applyLambda()) {
          k = k * scale;
        }
        setTorsionParameters(index++, a1, a2, a3, a4, j + 1,
                torsionType.phase[j] * OpenMM_RadiansPerDegree, OpenMM_KJPerKcal * k);
      }
    }
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
