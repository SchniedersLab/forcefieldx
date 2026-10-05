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

import ffx.openmm.ffm.CustomAngleForce;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;
import ffx.potential.terms.AnglePotentialEnergy;

import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static java.lang.String.format;

/**
 * OpenMM Angle Force backed by FFM {@link CustomAngleForce}.
 */
public class AngleForce extends CustomAngleForce {

  private static final Logger logger = Logger.getLogger(AngleForce.class.getName());

  private int nAngles = 0;
  private final boolean manyBodyTitration;
  private final boolean rigidHydrogenAngles;

  /**
   * Create an OpenMM Angle Force.
   *
   * @param anglePotentialEnergy The AnglePotentialEnergy instance that contains the angles.
   * @param openMMEnergy         The OpenMM Energy instance.
   */
  public AngleForce(AnglePotentialEnergy anglePotentialEnergy, OpenMMEnergy openMMEnergy) {
    super(anglePotentialEnergy.getAngleEnergyString());
    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);
    rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    addPerAngleParameter("theta0");
    addPerAngleParameter("k");
    setName("Angle");

    Angle[] angles = anglePotentialEnergy.getAngleArray();
    try (DoubleArray parameters = new DoubleArray(0)) {
      for (Angle angle : angles) {
        AngleType angleType = angle.getAngleType();
        AngleType.AngleMode angleMode = angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
          // Skip In-Plane angles unless this is ManyBody Titration.
        } else if (isHydrogenAngle(angle) && rigidHydrogenAngles) {
          logger.log(Level.INFO, " Constrained angle %s was not added to the AngleForce.", angle);
        } else {
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();

          double theta0 = angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angleType.angleUnit * angleType.forceConstant;
          if (angleMode == AngleType.AngleMode.IN_PLANE) {
            // This is a placeholder Angle, in case the In-Plane Angle is switched to a
            // Normal Angle during updateAngleForce.
            k = 0.0;
          }
          parameters.append(theta0);
          parameters.append(k);
          addAngle(i1, i2, i3, parameters);
          nAngles++;
          parameters.resize(0);
        }
      }
    }

    if (nAngles > 0) {
      int forceGroup = anglePotentialEnergy.getForceGroup();
      setForceGroup(forceGroup);
      logger.info(format("  Angles:                            %10d", nAngles));
      logger.fine(format("   Force Group:                      %10d", forceGroup));
    }
  }

  /**
   * Create an OpenMM Angle Force for Dual Topology.
   *
   * @param anglePotentialEnergy      The AnglePotentialEnergy instance that contains the angles.
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public AngleForce(AnglePotentialEnergy anglePotentialEnergy, int topology,
                    OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(anglePotentialEnergy.getAngleEnergyString());
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ForceField forceField = forceFieldEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);
    rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
    addPerAngleParameter("theta0");
    addPerAngleParameter("k");
    setName("Angle");

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    Angle[] angles = anglePotentialEnergy.getAngleArray();
    try (DoubleArray parameters = new DoubleArray(0)) {
      for (Angle angle : angles) {
        AngleType angleType = angle.getAngleType();
        AngleType.AngleMode angleMode = angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
          // Skip In-Plane angles unless this is ManyBody Titration.
        } else if (isHydrogenAngle(angle) && rigidHydrogenAngles) {
          logger.log(Level.INFO, " Constrained angle %s was not added to the AngleForce.", angle);
        } else {
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();

          double theta0 = angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angleType.angleUnit * angleType.forceConstant;
          if (angleMode == AngleType.AngleMode.IN_PLANE) {
            k = 0.0;
          }
          // Don't apply lambda scale to alchemical angle
          if (!angle.applyLambda()) {
            k = k * scale;
          }
          i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
          i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
          i3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i3);
          parameters.append(theta0);
          parameters.append(k);
          addAngle(i1, i2, i3, parameters);
          nAngles++;
          parameters.resize(0);
        }
      }
    }

    if (nAngles > 0) {
      int forceGroup = anglePotentialEnergy.getForceGroup();
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
    AnglePotentialEnergy anglePotentialEnergy = openMMEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return null;
    }
    AngleForce angleForce = new AngleForce(anglePotentialEnergy, openMMEnergy);
    if (angleForce.nAngles > 0) {
      return angleForce;
    }
    angleForce.destroy();
    return null;
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Angle Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Angle Force, or null if there are no angles.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    AnglePotentialEnergy anglePotentialEnergy = forceFieldEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return null;
    }
    AngleForce angleForce = new AngleForce(anglePotentialEnergy, topology, openMMDualTopologyEnergy);
    if (angleForce.nAngles > 0) {
      return angleForce;
    }
    angleForce.destroy();
    return null;
  }

  /**
   * Update an existing angle force for the OpenMM System.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    AnglePotentialEnergy anglePotentialEnergy = openMMEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return;
    }
    Angle[] angles = anglePotentialEnergy.getAngleArray();

    try (DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Angle angle : angles) {
        AngleType.AngleMode angleMode = angle.angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
          // Skip In-Plane angles unless this is ManyBody Titration.
        } else if (!rigidHydrogenAngles || !isHydrogenAngle(angle)) {
          // Update angles that do not involve rigid hydrogen atoms.
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();
          double theta0 = angle.angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angle.angleType.angleUnit * angle.angleType.forceConstant;
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
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update existing angle force for the Dual-Topology OpenMM System.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    AnglePotentialEnergy anglePotentialEnergy = forceFieldEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return;
    }
    Angle[] angles = anglePotentialEnergy.getAngleArray();
    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Angle angle : angles) {
        AngleType.AngleMode angleMode = angle.angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleType.AngleMode.IN_PLANE) {
          // Skip In-Plane angles unless this is ManyBody Titration.
        } else if (!rigidHydrogenAngles || !isHydrogenAngle(angle)) {
          // Update angles that do not involve rigid hydrogen atoms.
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();
          double theta0 = angle.angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angle.angleType.angleUnit * angle.angleType.forceConstant;
          if (angleMode == AngleType.AngleMode.IN_PLANE) {
            // Zero the force constant for In-Plane Angles.
            k = 0.0;
          }
          // Don't apply lambda scale to alchemical angle
          if (!angle.applyLambda()) {
            k = k * scale;
          }
          i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
          i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
          i3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i3);
          parameters.append(theta0);
          parameters.append(k);
          setAngleParameters(index++, i1, i2, i3, parameters);
          parameters.resize(0);
        }
      }
    }
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
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
      double angleVal = angle.angleType.angle[angle.nh];
      if (angleVal < 160.0) {
        Atom atom1 = angle.getAtom(0);
        Atom atom2 = angle.getAtom(1);
        Atom atom3 = angle.getAtom(2);
        return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
      }
    }
    return false;
  }
}
