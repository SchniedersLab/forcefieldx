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

import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;
import ffx.openmm.CustomCompoundBondForce;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static java.lang.String.format;

/**
 * OpenMM Angle-Torsion Force.
 */
public class AngleTorsionForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(AngleTorsionForce.class.getName());

  /**
   * Create an OpenMM Angle-Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the Angle-Torsions.
   */
  public AngleTorsionForce(OpenMMEnergy openMMEnergy) {
    super(4, AngleTorsion.angleTorsionForm());
    AngleTorsion[] angleTorsions = openMMEnergy.getAngleTorsions();
    if (angleTorsions == null || angleTorsions.length < 1) {
      // Free the memory created by the call to super.
      destroy();
      return;
    }
    addGlobalParameter("phi1", 0);
    addGlobalParameter("phi2", Math.PI);
    addGlobalParameter("phi3", 0);
    for (int m = 1; m < 3; m++) {
      for (int n = 1; n < 4; n++) {
        addPerBondParameter(format("k%d%d", m, n));
      }
    }
    for (int m = 1; m < 3; m++) {
      addPerBondParameter(format("a%d", m));
    }
    for (AngleTorsion angleTorsion : angleTorsions) {
      double[] constants = angleTorsion.getConstants();
      DoubleArray parameters = new DoubleArray(0);
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 3; n++) {
          int index = (3 * m) + n;
          parameters.append(constants[index] * OpenMM_KJPerKcal);
        }
      }
      Atom[] atoms = angleTorsion.getAtomArray(true);
      parameters.append(angleTorsion.angleType1.angle[0] * OpenMM_RadiansPerDegree);
      parameters.append(angleTorsion.angleType2.angle[0] * OpenMM_RadiansPerDegree);

      IntArray particles = new IntArray(0);
      for (int i = 0; i < 4; i++) {
        particles.append(atoms[i].getArrayIndex());
      }

      addBond(particles, parameters);
      parameters.destroy();
      particles.destroy();
    }

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("ANGLE_TORSION_FORCE_GROUP", 0);

    setForceGroup(forceGroup);
    logger.info(format("  Angle-Torsions:                    %10d", angleTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create a Dual Topology OpenMM Angle-Torsion Force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public AngleTorsionForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(4, AngleTorsion.angleTorsionForm());

    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    AngleTorsion[] angleTorsions = forceFieldEnergy.getAngleTorsions();
    if (angleTorsions == null || angleTorsions.length < 1) {
      // Free the memory created by the call to super.
      destroy();
      return;
    }
    addGlobalParameter("phi1", 0);
    addGlobalParameter("phi2", Math.PI);
    addGlobalParameter("phi3", 0);
    for (int m = 1; m < 3; m++) {
      for (int n = 1; n < 4; n++) {
        addPerBondParameter(format("k%d%d", m, n));
      }
    }
    for (int m = 1; m < 3; m++) {
      addPerBondParameter(format("a%d", m));
    }

    double scaleDT = openMMDualTopologyEnergy.getTopologyScale(topology);

    for (AngleTorsion angleTorsion : angleTorsions) {
      double scale = 1.0;
      // Don't apply lambda scale to alchemical stretch-torsion
      if (!angleTorsion.applyLambda()) {
        scale = scaleDT;
      }
      double[] constants = angleTorsion.getConstants();
      DoubleArray parameters = new DoubleArray(0);
      for (int m = 0; m < 2; m++) {
        for (int n = 0; n < 3; n++) {
          int index = (3 * m) + n;
          parameters.append(constants[index] * OpenMM_KJPerKcal * scale);
        }
      }
      Atom[] atoms = angleTorsion.getAtomArray(true);
      parameters.append(angleTorsion.angleType1.angle[0] * OpenMM_RadiansPerDegree);
      parameters.append(angleTorsion.angleType2.angle[0] * OpenMM_RadiansPerDegree);

      IntArray particles = new IntArray(0);
      for (int i = 0; i < 4; i++) {
        int atomIndex = atoms[i].getArrayIndex();
        atomIndex = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, atomIndex);
        particles.append(atomIndex);
      }

      addBond(particles, parameters);
      parameters.destroy();
      particles.destroy();
    }

    int forceGroup = forceFieldEnergy.getMolecularAssembly().getForceField().getInteger("ANGLE_TORSION_FORCE_GROUP", 0);

    setForceGroup(forceGroup);
    logger.info(format("  Angle-Torsions:                    %10d", angleTorsions.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Angle-Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angle-torsions.
   * @return An OpenMM Stretch-Bend Force, or null if there are no angle-torsion.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    AngleTorsion[] angleTorsions = openMMEnergy.getAngleTorsions();
    if (angleTorsions == null || angleTorsions.length < 1) {
      return null;
    }
    return new AngleTorsionForce(openMMEnergy);
  }

  /**
   * Convenience method to construct a Dual Topology OpenMM Angle-Torsion Force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Stretch-Bend Force, or null if there are no angle-torsion.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    AngleTorsion[] angleTorsions = forceFieldEnergy.getAngleTorsions();
    if (angleTorsions == null || angleTorsions.length < 1) {
      return null;
    }
    return new AngleTorsionForce(topology, openMMDualTopologyEnergy);
  }
}
