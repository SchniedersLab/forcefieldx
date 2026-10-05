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

import ffx.openmm.ffm.CustomExternalForce;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.RestrainPosition;
import ffx.potential.terms.RestrainPositionPotentialEnergy;

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static java.lang.String.format;

/**
 * Restrain Positions Force backed by FFM {@link CustomExternalForce}.
 */
public class RestrainPositionsForce extends CustomExternalForce {

  private static final Logger logger = Logger.getLogger(RestrainPositionsForce.class.getName());

  /**
   * Restrain Positions Force constructor.
   *
   * @param restrainPositionPotentialEnergy RestrainPositionPotentialEnergy instance.
   */
  public RestrainPositionsForce(RestrainPositionPotentialEnergy restrainPositionPotentialEnergy) {
    super(RestrainPositionPotentialEnergy.getRestrainPositionEnergyString());
    RestrainPosition[] restrainPositions = restrainPositionPotentialEnergy.getRestrainPositionArray();
    // Define per-particle parameters.
    addPerParticleParameter("k0");
    addPerParticleParameter("x0");
    addPerParticleParameter("y0");
    addPerParticleParameter("z0");

    int nRestraints = restrainPositions.length;
    double convert = KJ_PER_KCAL / (NM_PER_ANGSTROM * NM_PER_ANGSTROM);

    try (DoubleArray parameters = new DoubleArray(4)) {
      for (RestrainPosition restrainPosition : restrainPositions) {
        double forceConstant = restrainPosition.getForceConstant() * convert;
        Atom[] restrainPositionAtoms = restrainPosition.getAtoms();
        int numAtoms = restrainPosition.getNumAtoms();
        double[][] equilibriumCoordinates = restrainPosition.getEquilibriumCoordinates();
        for (int i = 0; i < numAtoms; i++) {
          equilibriumCoordinates[i][0] *= NM_PER_ANGSTROM;
          equilibriumCoordinates[i][1] *= NM_PER_ANGSTROM;
          equilibriumCoordinates[i][2] *= NM_PER_ANGSTROM;
        }

        for (int i = 0; i < numAtoms; i++) {
          int index = restrainPositionAtoms[i].getXyzIndex() - 1;
          parameters.set(0, forceConstant);
          for (int j = 0; j < 3; j++) {
            parameters.set(j + 1, equilibriumCoordinates[i][j]);
          }
          addParticle(index, parameters);
        }
      }
    }

    int forceGroup = restrainPositionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Restrain Positions\t%6d\t\t%d", nRestraints, forceGroup));
  }

  /**
   * Add a Restrain-Position force to the OpenMM System.
   *
   * @param openMMEnergy The OpenMM Energy instance.
   * @return Force instance or null if no restrain positions potential energy exists.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    RestrainPositionPotentialEnergy restrainPositionPotentialEnergy =
        openMMEnergy.getRestrainPositionPotentialEnergy();
    if (restrainPositionPotentialEnergy == null) {
      return null;
    }
    return new RestrainPositionsForce(restrainPositionPotentialEnergy);
  }
}
