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
package ffx.potential.openmm;

import ffx.openmm.CustomCentroidBondForce;
import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.RestrainGroups;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static java.lang.String.format;

/**
 * Restrain Groups Force.
 */
public class RestrainGroupsForce extends CustomCentroidBondForce {

  private static final Logger logger = Logger.getLogger(RestrainGroupsForce.class.getName());

  private static final String energy = "step(distance(g1,g2)-u)*k*(distance(g1,g2)-u)^2+step(l-distance(g1,g2))*k*(distance(g1,g2)-l)^2";

  /**
   * Restrain Groups Force constructor.
   *
   * @param openMMEnergy The OpenMM Energy.
   */
  public RestrainGroupsForce(OpenMMEnergy openMMEnergy) {
    super(2, energy);
    RestrainGroups restrainGroups = openMMEnergy.getRestrainGroups();
    if (restrainGroups == null) {
      destroy();
      return;
    }

    addPerBondParameter("k");
    addPerBondParameter("l");
    addPerBondParameter("u");

    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();

    // Create the Restrain Groups.
    int nGroups = restrainGroups.getNumberOfGroups();
    IntArray group = new IntArray(0);
    DoubleArray weight = new DoubleArray(0);
    for (int j = 0; j < nGroups; j++) {
      int[] groupMembers = restrainGroups.getGroupMembers(j);
      for (int i : groupMembers) {
        group.append(i);
        weight.append(atoms[i].getMass());
      }
      addGroup(group, weight);
      group.resize(0);
      weight.resize(0);
    }
    group.destroy();
    weight.destroy();

    // Add the restraints between groups.
    double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    int nRestraints = restrainGroups.getNumberOfRestraints();
    int[] group1 = restrainGroups.getGroup1();
    int[] group2 = restrainGroups.getGroup2();
    double[] forceConstants = restrainGroups.getForceConstants();
    double[] smallerDistance = restrainGroups.getSmallerDistance();
    double[] largerDistance = restrainGroups.getLargerDistance();
    group = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (int i = 0; i < nRestraints; i++) {
      group.append(group1[i]);
      group.append(group2[i]);
      parameters.append(forceConstants[i] * convert);
      parameters.append(smallerDistance[i] * OpenMM_NmPerAngstrom);
      parameters.append(largerDistance[i] * OpenMM_NmPerAngstrom);
      addBond(group, parameters);
      group.resize(0);
      parameters.resize(0);
    }
    group.destroy();
    parameters.destroy();

    // Add the constraint force.
    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("RESTRAIN_GROUPS_FORCE_GROUP", 0);
    setForceGroup(forceGroup);

    if (openMMEnergy.getCrystal().aperiodic()) {
      setUsesPeriodicBoundaryConditions(OpenMM_False);
    } else {
      setUsesPeriodicBoundaryConditions(OpenMM_True);
    }
    logger.log(Level.INFO, format("  Restrain Groups \t%6d\t\t%1d", nRestraints, forceGroup));
  }

  /**
   * Add a Restrain-Groups force to the OpenMM System.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    RestrainGroups restrainGroups = openMMEnergy.getRestrainGroups();
    if (restrainGroups == null) {
      return null;
    }
    return new RestrainGroupsForce(openMMEnergy);
  }
}
