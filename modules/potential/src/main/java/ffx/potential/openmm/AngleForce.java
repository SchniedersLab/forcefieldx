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
import ffx.openmm.CustomAngleForce;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static java.lang.String.format;

/**
 * OpenMM Angle Force.
 */
public class AngleForce extends CustomAngleForce {

  private static final Logger logger = Logger.getLogger(AngleForce.class.getName());

  private int nAngles = 0;
  private final boolean manyBodyTitration;
  private final boolean rigidHydrogenAngles;

  /**
   * Create an OpenMM Angle Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   */
  public AngleForce(OpenMMEnergy openMMEnergy) {
    super(openMMEnergy.getAngleEnergyString());
    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);
    rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    addPerAngleParameter("theta0");
    addPerAngleParameter("k");
    setName("Angle");

    DoubleArray parameters = new DoubleArray(0);
    Angle[] angles = openMMEnergy.getAngles();
    for (Angle angle : angles) {
      AngleType angleType = angle.getAngleType();
      AngleType.AngleMode angleMode = angleType.angleMode;
      if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
        // Skip In-Plane angles unless this is ManyBody Titration.
      } else if (isHydrogenAngle(angle) && rigidHydrogenAngles) {
        logger.log(Level.INFO, " Constrained angle %s was not added the AngleForce.", angle);
      } else {
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;

        double theta0 = angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angleType.angleUnit * angleType.forceConstant;
        if (angleMode == AngleType.AngleMode.IN_PLANE) {
          // This is a place-holder Angle, in case the In-Plane Angle is swtiched to a
          // Normal Angle during in the udpateAngleForce.
          k = 0.0;
        }
        parameters.append(theta0);
        parameters.append(k);
        addAngle(i1, i2, i3, parameters);
        nAngles++;
        parameters.resize(0);
      }
    }
    parameters.destroy();

    if (nAngles > 0) {
      int forceGroup = forceField.getInteger("ANGLE_FORCE_GROUP", 0);
      setForceGroup(forceGroup);
      logger.info(format("  Angles:                            %10d", nAngles));
      logger.fine(format("   Force Group:                      %10d", forceGroup));
    }
  }

  /**
   * Convenience method to construct an OpenMM Angle Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   * @return An Angle Force, or null if there are no angles.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return null;
    }
    AngleForce angleForce = new AngleForce(openMMEnergy);
    if (angleForce.nAngles > 0) {
      return angleForce;
    }
    return null;
  }

  /**
   * Update an existing angle force for the OpenMM System.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return;
    }

    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;
      if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
        // Skip In-Plane angles unless this is ManyBody Titration.
      }
      // Update angles that do not involve rigid hydrogen atoms.
      else if (!rigidHydrogenAngles || !isHydrogenAngle(angle)) {
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        if (angleMode == AngleType.AngleMode.IN_PLANE) {
          // Zero the force constant for In-Plane Angles.
          k = 0.0;
        }
        parameters.append(theta0);
        parameters.append(k);
        setAngleParameters(index++, i1, i2, i3, parameters);
        parameters.resize(0);
      }
    }
    parameters.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Check to see if an angle is a hydrogen angle. This method only returns true for hydrogen
   * angles that are less than 160 degrees.
   *
   * @param angle Angle to check.
   * @return boolean indicating whether an angle is a hydrogen angle that is less than 160 degrees.
   */
  private boolean isHydrogenAngle(Angle angle) {
    if (angle.containsHydrogen()) {
      // Equilibrium angle value in degrees.
      double angleVal = angle.angleType.angle[angle.nh];
      // Make sure angle is less than 160 degrees because greater than 160 degrees will not be
      // constrained
      // well using the law of cosines.
      if (angleVal < 160.0) {
        Atom atom1 = angle.getAtom(0);
        Atom atom2 = angle.getAtom(1);
        Atom atom3 = angle.getAtom(2);
        // Setting constraints only on angles where atom1 or atom3 is a hydrogen while atom2 is
        // not a hydrogen.
        return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
      }
    }
    return false;
  }

}
