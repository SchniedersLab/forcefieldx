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

import ffx.openmm.CustomExternalForce;
import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.RestrainPosition;
import ffx.potential.terms.RestrainPositionPotentialEnergy;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;

/**
 * Restrain Positions Force.
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
    // Define per particle parameters.
    addPerParticleParameter("k0");
    addPerParticleParameter("x0");
    addPerParticleParameter("y0");
    addPerParticleParameter("z0");

    int nRestraints = restrainPositions.length;
    double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    for (RestrainPosition restrainPosition : restrainPositions) {
      double forceConstant = restrainPosition.getForceConstant() * convert;
      Atom[] restrainPositionAtoms = restrainPosition.getAtoms();
      int numAtoms = restrainPosition.getNumAtoms();
      double[][] equilibriumCoordinates = restrainPosition.getEquilibriumCoordinates();
      for (int i = 0; i < numAtoms; i++) {
        equilibriumCoordinates[i][0] *= OpenMM_NmPerAngstrom;
        equilibriumCoordinates[i][1] *= OpenMM_NmPerAngstrom;
        equilibriumCoordinates[i][2] *= OpenMM_NmPerAngstrom;
      }

      DoubleArray parameters = new DoubleArray(4);
      for (int i = 0; i < numAtoms; i++) {
        int index = restrainPositionAtoms[i].getXyzIndex() - 1;
        parameters.set(0, forceConstant);
        for (int j = 0; j < 3; j++) {
          parameters.set(j + 1, equilibriumCoordinates[i][j]);
        }
        addParticle(index, parameters);
      }
      parameters.destroy();
    }

    int forceGroup = restrainPositionPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Restrain Positions\t%6d\t\t%d", nRestraints, forceGroup));
  }

  /**
   * Add a Restrain-Position force to the OpenMM System.
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
