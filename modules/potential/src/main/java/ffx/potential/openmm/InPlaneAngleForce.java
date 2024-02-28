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

import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;
import ffx.openmm.CustomCompoundBondForce;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static java.lang.String.format;

/**
 * OpenMM In-Plane Angle Force.
 */
public class InPlaneAngleForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(InPlaneAngleForce.class.getName());

  private int nAngles = 0;
  private final boolean manyBodyTitration;

  /**
   * Create an OpenMM Angle Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   */
  public InPlaneAngleForce(OpenMMEnergy openMMEnergy) {
    super(4, openMMEnergy.getInPlaneAngleEnergyString());
    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      // Clean up the Memory allocated by the OpenMMCustomCompoundBondForce constructor.
      destroy();
      return;
    }
    addPerBondParameter("theta0");
    addPerBondParameter("k");
    setName("InPlaneAngle");

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;

      if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
        // Skip Normal angles unless this is ManyBody Titration.
      } else {
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        int i4 = 0;
        if (angleMode == AngleType.AngleMode.NORMAL) {
          // This is a place-holder Angle, in case the Normal Angle is switched to an
          // In-Plane Angle during in the udpateInPlaneAngleForce.
          k = 0.0;
          Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
          if (fourthAtom != null) {
            i4 = fourthAtom.getXyzIndex() - 1;
          } else {
            while (i1 == i4 || i2 == i4 || i3 == i4) {
              i4++;
            }
          }
        } else {
          i4 = angle.getAtom4().getXyzIndex() - 1;
        }
        particles.append(i1);
        particles.append(i2);
        particles.append(i3);
        particles.append(i4);
        parameters.append(theta0);
        parameters.append(k);
        addBond(particles, parameters);
        nAngles++;
        particles.resize(0);
        parameters.resize(0);
      }
    }
    particles.destroy();
    parameters.destroy();

    if (nAngles > 0) {
      int forceGroup = forceField.getInteger("IN_PLANE_ANGLE_FORCE_GROUP", 0);
      setForceGroup(forceGroup);
      logger.log(Level.INFO, format("  In-Plane Angles \t%6d\t\t%1d", nAngles, forceGroup));
    }
  }

  /**
   * Convenience method to construct an OpenMM In-Plane Angle Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   * @return An OpenMM Angle Force, or null if there are no angles.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    Angle[] angles = openMMEnergy.getAngles();
    if (angles == null || angles.length < 1) {
      return null;
    }
    InPlaneAngleForce angleForce = new InPlaneAngleForce(openMMEnergy);
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
    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (Angle angle : angles) {
      AngleType.AngleMode angleMode = angle.angleType.angleMode;
      if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
        // Skip Normal angles unless this is ManyBody Titration.
      } else {
        double theta0 = angle.angleType.angle[angle.nh];
        double k = OpenMM_KJPerKcal * angle.angleType.angleUnit * angle.angleType.forceConstant;
        int i1 = angle.getAtom(0).getXyzIndex() - 1;
        int i2 = angle.getAtom(1).getXyzIndex() - 1;
        int i3 = angle.getAtom(2).getXyzIndex() - 1;
        // There is no 4th atom for normal angles, so set the index to first atom.
        int i4 = 0;
        if (angleMode == AngleType.AngleMode.NORMAL) {
          // Zero the force constant for Normal Angles.
          k = 0.0;
          Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
          if (fourthAtom != null) {
            i4 = fourthAtom.getXyzIndex() - 1;
          } else {
            while (i1 == i4 || i2 == i4 || i3 == i4) {
              i4++;
            }
          }
        } else {
          i4 = angle.getAtom4().getXyzIndex() - 1;
        }
        particles.append(i1);
        particles.append(i2);
        particles.append(i3);
        particles.append(i4);
        parameters.append(theta0);
        parameters.append(k);
        setBondParameters(index++, particles, parameters);
        particles.resize(0);
        parameters.resize(0);
      }
    }
    particles.destroy();
    parameters.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }
}
