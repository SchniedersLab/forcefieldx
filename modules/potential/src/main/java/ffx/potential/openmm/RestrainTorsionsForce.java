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
import ffx.potential.parameters.TorsionType;
import ffx.potential.terms.RestrainTorsionPotentialEnergy;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static java.lang.String.format;

/**
 * Restrain Torsions Force.
 */
public class RestrainTorsionsForce extends PeriodicTorsionForce {

  private static final Logger logger = Logger.getLogger(RestrainTorsionsForce.class.getName());

  /**
   * Restrain Torsion Force constructor.
   *
   * @param restrainTorsionPotentialEnergy The RestrainTorsionPotentialEnergy instance that contains the torsions.
   */
  public RestrainTorsionsForce(RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy) {
    Torsion[] restrainTorsions = restrainTorsionPotentialEnergy.getRestrainTorsionArray();
    for (Torsion restrainTorsion : restrainTorsions) {
      int a1 = restrainTorsion.getAtom(0).getArrayIndex();
      int a2 = restrainTorsion.getAtom(1).getArrayIndex();
      int a3 = restrainTorsion.getAtom(2).getArrayIndex();
      int a4 = restrainTorsion.getAtom(3).getArrayIndex();
      TorsionType torsionType = restrainTorsion.torsionType;
      int nTerms = torsionType.phase.length;
      for (int j = 0; j < nTerms; j++) {
        addTorsion(a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree,
            OpenMM_KJPerKcal * torsionType.torsionUnit * torsionType.amplitude[j]);
      }
    }
    int forceGroup = restrainTorsionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Restrain-Torsions \t%6d\t\t%1d", restrainTorsions.length, forceGroup));
  }

  /**
   * Restrain Torsion Force constructor for Dual Topology.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public RestrainTorsionsForce(RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy, int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    Torsion[] restrainTorsions = restrainTorsionPotentialEnergy.getRestrainTorsionArray();

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    for (Torsion restrainTorsion : restrainTorsions) {
      int a1 = restrainTorsion.getAtom(0).getArrayIndex();
      int a2 = restrainTorsion.getAtom(1).getArrayIndex();
      int a3 = restrainTorsion.getAtom(2).getArrayIndex();
      int a4 = restrainTorsion.getAtom(3).getArrayIndex();
      a1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a1);
      a2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a2);
      a3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a3);
      a4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a4);
      TorsionType torsionType = restrainTorsion.torsionType;
      int nTerms = torsionType.phase.length;
      for (int j = 0; j < nTerms; j++) {
        double k = restrainTorsion.getTorsionScale() * torsionType.torsionUnit * torsionType.amplitude[j];
        k = k * (1 - scale); // multiply force constant by 1 minus lambda
        addTorsion(a1, a2, a3, a4, j + 1,
                torsionType.phase[j] * OpenMM_RadiansPerDegree,
                OpenMM_KJPerKcal * k);
      }
    }

    int forceGroup = forceFieldEnergy.getMolecularAssembly().getForceField().getInteger("RESTRAIN_TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Restrain-Torsions:                 %10d", restrainTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the torsions.
   * @return A Torsion Force, or null if there are no torsions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy =
        openMMEnergy.getRestrainTorsionPotentialEnergy();
    if (restrainTorsionPotentialEnergy == null) {
      return null;
    }
    return new RestrainTorsionsForce(restrainTorsionPotentialEnergy);
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
    RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy = forceFieldEnergy.getRestrainTorsionPotentialEnergy();
    if (restrainTorsionPotentialEnergy == null) {
      return null;
    }
    return new RestrainTorsionsForce(restrainTorsionPotentialEnergy, topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the Restraint-Torsion force.
   *
   * @param openMMEnergy The OpenMM Energy that contains the restraint-torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    // Check if this system has restraintTorsions.
    RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy =
        openMMEnergy.getRestrainTorsionPotentialEnergy();
    if (restrainTorsionPotentialEnergy == null) {
      return;
    }
    Torsion[] restrainTorsions = restrainTorsionPotentialEnergy.getRestrainTorsionArray();
    int index = 0;
    for (Torsion restrainTorsion : restrainTorsions) {
      TorsionType torsionType = restrainTorsion.torsionType;
      int nTerms = torsionType.phase.length;
      int a1 = restrainTorsion.getAtom(0).getArrayIndex();
      int a2 = restrainTorsion.getAtom(1).getArrayIndex();
      int a3 = restrainTorsion.getAtom(2).getArrayIndex();
      int a4 = restrainTorsion.getAtom(3).getArrayIndex();
      for (int j = 0; j < nTerms; j++) {
        double forceConstant = OpenMM_KJPerKcal * torsionType.torsionUnit * torsionType.amplitude[j];
        setTorsionParameters(index++, a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree, forceConstant);
      }
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the Restraint-Torsion force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy = forceFieldEnergy.getRestrainTorsionPotentialEnergy();
    if (restrainTorsionPotentialEnergy == null) {
      return;
    }
    // Check if this system has restraintTorsions.
    Torsion[] restrainTorsions = restrainTorsionPotentialEnergy.getRestrainTorsionArray();

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    int index = 0;
    for (Torsion restrainTorsion : restrainTorsions) {
      TorsionType torsionType = restrainTorsion.torsionType;
      int nTerms = torsionType.phase.length;
      int a1 = restrainTorsion.getAtom(0).getArrayIndex();
      int a2 = restrainTorsion.getAtom(1).getArrayIndex();
      int a3 = restrainTorsion.getAtom(2).getArrayIndex();
      int a4 = restrainTorsion.getAtom(3).getArrayIndex();
      a1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a1);
      a2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a2);
      a3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a3);
      a4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, a4);
      for (int j = 0; j < nTerms; j++) {
        double forceConstant = restrainTorsion.getTorsionScale() * torsionType.torsionUnit * torsionType.amplitude[j];
        forceConstant = forceConstant * (1 - scale); // multiply force constant by 1 minus lambda
        setTorsionParameters(index++, a1, a2, a3, a4, j + 1,
                torsionType.phase[j] * OpenMM_RadiansPerDegree, OpenMM_KJPerKcal * forceConstant);
      }
    }
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
