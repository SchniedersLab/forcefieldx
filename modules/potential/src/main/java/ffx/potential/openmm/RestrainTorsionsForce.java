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

import ffx.openmm.Force;
import ffx.openmm.PeriodicTorsionForce;
import ffx.potential.bonded.RestraintTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.parameters.TorsionType;

import java.util.List;
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
   * @param openMMEnergy The OpenMM Energy that contains the restraint-torsions.
   */
  public RestrainTorsionsForce(OpenMMEnergy openMMEnergy) {
    List<RestraintTorsion> restraintTorsions = openMMEnergy.getRestraintTorsions();
    if (restraintTorsions == null || restraintTorsions.isEmpty()) {
      return;
    }

    for (RestraintTorsion restraintTorsion : restraintTorsions) {
      int a1 = restraintTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = restraintTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = restraintTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = restraintTorsion.getAtom(3).getXyzIndex() - 1;
      TorsionType torsionType = restraintTorsion.torsionType;
      int nTerms = torsionType.phase.length;
      for (int j = 0; j < nTerms; j++) {
        addTorsion(a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree,
            OpenMM_KJPerKcal * torsionType.torsionUnit * torsionType.amplitude[j]);
      }
    }

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("RESTRAINT_TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Restraint-Torsions \t%6d\t\t%1d", restraintTorsions.size(), forceGroup));
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
    return new RestrainTorsionsForce(openMMEnergy);
  }

  /**
   * Update the Restraint-Torsion force.
   *
   * @param openMMEnergy The OpenMM Energy that contains the restraint-torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    // Check if this system has restraintTorsions.
    List<RestraintTorsion> restraintTorsions = openMMEnergy.getRestraintTorsions();
    if (restraintTorsions == null || restraintTorsions.isEmpty()) {
      return;
    }

    int index = 0;
    for (RestraintTorsion restraintTorsion : restraintTorsions) {
      TorsionType torsionType = restraintTorsion.torsionType;
      int nTerms = torsionType.phase.length;
      int a1 = restraintTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = restraintTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = restraintTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = restraintTorsion.getAtom(3).getXyzIndex() - 1;
      for (int j = 0; j < nTerms; j++) {
        double forceConstant = OpenMM_KJPerKcal * torsionType.torsionUnit * torsionType.amplitude[j];
        setTorsionParameters(index++, a1, a2, a3, a4, j + 1,
            torsionType.phase[j] * OpenMM_RadiansPerDegree, forceConstant);
      }
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

}
