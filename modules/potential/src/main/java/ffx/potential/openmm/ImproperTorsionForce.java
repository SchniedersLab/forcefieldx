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
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.parameters.ImproperTorsionType;
import ffx.potential.terms.ImproperTorsionPotentialEnergy;

import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static java.lang.String.format;

/**
 * OpenMM Improper Torsion Force.
 */
public class ImproperTorsionForce extends PeriodicTorsionForce {

  private static final Logger logger = Logger.getLogger(ImproperTorsionForce.class.getName());

  /**
   * Create an OpenMM Improper Torsion Force.
   *
   * @param improperTorsionPotentialEnergy The ImproperTorsionPotentialEnergy instance that contains the improper torsions.
   */
  public ImproperTorsionForce(ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy) {
    ImproperTorsion[] improperTorsions = improperTorsionPotentialEnergy.getImproperTorsionArray();
    for (ImproperTorsion improperTorsion : improperTorsions) {
      int a1 = improperTorsion.getAtom(0).getArrayIndex();
      int a2 = improperTorsion.getAtom(1).getArrayIndex();
      int a3 = improperTorsion.getAtom(2).getArrayIndex();
      int a4 = improperTorsion.getAtom(3).getArrayIndex();
      ImproperTorsionType type = improperTorsion.improperType;
      double forceConstant = OpenMM_KJPerKcal * type.impTorUnit * improperTorsion.scaleFactor * type.k;
      addTorsion(a1, a2, a3, a4, type.periodicity, type.phase * OpenMM_RadiansPerDegree, forceConstant);
    }
    int forceGroup = improperTorsionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Improper Torsions:                 %10d", improperTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create a Dual Topology OpenMM Improper Torsion Force.
   *
   * @param improperTorsionPotentialEnergy The ImproperTorsionPotentialEnergy instance that contains the improper torsions.
   * @param topology                       The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy       The OpenMMDualTopologyEnergy instance.
   */
  public ImproperTorsionForce(ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy,
                              int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ImproperTorsion[] improperTorsions = improperTorsionPotentialEnergy.getImproperTorsionArray();
    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    for (ImproperTorsion improperTorsion : improperTorsions) {
      int a1 = improperTorsion.getAtom(0).getArrayIndex();
      int a2 = improperTorsion.getAtom(1).getArrayIndex();
      int a3 = improperTorsion.getAtom(2).getArrayIndex();
      int a4 = improperTorsion.getAtom(3).getArrayIndex();
      a1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a1);
      a2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a2);
      a3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a3);
      a4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a4);
      ImproperTorsionType type = improperTorsion.improperType;
      double forceConstant = OpenMM_KJPerKcal * type.impTorUnit * improperTorsion.scaleFactor * type.k;
      // Don't apply lambda scale to alchemical improper torsion
      if (!improperTorsion.applyLambda()) {
        forceConstant *= scale;
      }
      addTorsion(a1, a2, a3, a4, type.periodicity, type.phase * OpenMM_RadiansPerDegree, forceConstant);
    }

    int forceGroup = improperTorsionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Improper Torsions:                 %10d", improperTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Improper Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the improper torsions.
   * @return An Improper Torsion Force, or null if there are no improper torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy = openMMEnergy.getImproperTorsionPotentialEnergy();
    if (improperTorsionPotentialEnergy == null) {
      return null;
    }
    return new ImproperTorsionForce(improperTorsionPotentialEnergy);
  }

  /**
   * Convenience method to construct a Dual Topology OpenMM Improper Torsion Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return A Torsion Force, or null if there are no torsions.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy = forceFieldEnergy.getImproperTorsionPotentialEnergy();
    if (improperTorsionPotentialEnergy == null) {
      return null;
    }
    return new ImproperTorsionForce(improperTorsionPotentialEnergy, topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the Improper Torsion force.
   *
   * @param openMMEnergy The OpenMM Energy that contains the improper torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy = openMMEnergy.getImproperTorsionPotentialEnergy();
    if (improperTorsionPotentialEnergy == null) {
      return;
    }
    ImproperTorsion[] improperTorsions = improperTorsionPotentialEnergy.getImproperTorsionArray();
    int nImproperTorsions = improperTorsions.length;
    for (int i = 0; i < nImproperTorsions; i++) {
      ImproperTorsion improperTorsion = improperTorsions[i];
      int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
      ImproperTorsionType type = improperTorsion.improperType;
      double forceConstant = OpenMM_KJPerKcal * type.impTorUnit * improperTorsion.scaleFactor * type.k;
      setTorsionParameters(i, a1, a2, a3, a4, type.periodicity, type.phase * OpenMM_RadiansPerDegree, forceConstant);
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the Dual Topology Improper Torsion force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy = forceFieldEnergy.getImproperTorsionPotentialEnergy();
    if (improperTorsionPotentialEnergy == null) {
      return;
    }
    ImproperTorsion[] improperTorsions = improperTorsionPotentialEnergy.getImproperTorsionArray();
    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    int nImproperTorsions = improperTorsions.length;
    for (int i = 0; i < nImproperTorsions; i++) {
      ImproperTorsion improperTorsion = improperTorsions[i];
      int a1 = improperTorsion.getAtom(0).getArrayIndex();
      int a2 = improperTorsion.getAtom(1).getArrayIndex();
      int a3 = improperTorsion.getAtom(2).getArrayIndex();
      int a4 = improperTorsion.getAtom(3).getArrayIndex();
      a1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a1);
      a2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a2);
      a3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a3);
      a4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a4);
      ImproperTorsionType type = improperTorsion.improperType;
      double forceConstant = OpenMM_KJPerKcal * type.impTorUnit * improperTorsion.scaleFactor * type.k;
      // Don't apply lambda scale to alchemical improper torsion
      if (!improperTorsion.applyLambda()) {
        forceConstant *= scale;
      }
      setTorsionParameters(i, a1, a2, a3, a4, type.periodicity, type.phase * OpenMM_RadiansPerDegree, forceConstant);
    }

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
