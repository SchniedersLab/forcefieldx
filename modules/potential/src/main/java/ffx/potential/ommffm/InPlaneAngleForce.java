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

import ffx.openmm.ffm.CustomCompoundBondForce;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.IntArray;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.ForceField;
import ffx.potential.terms.AnglePotentialEnergy;

import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static java.lang.String.format;

/**
 * OpenMM In-Plane Angle Force backed by FFM {@link CustomCompoundBondForce}.
 */
public class InPlaneAngleForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(InPlaneAngleForce.class.getName());

  private int nAngles = 0;
  private final boolean manyBodyTitration;

  /**
   * Create an OpenMM In-Plane Angle Force.
   *
   * @param anglePotentialEnergy The AnglePotentialEnergy that contains the angles.
   * @param openMMEnergy         The OpenMM Energy instance that contains the angles.
   */
  public InPlaneAngleForce(AnglePotentialEnergy anglePotentialEnergy, OpenMMEnergy openMMEnergy) {
    super(4, anglePotentialEnergy.getInPlaneAngleEnergyString());
    Angle[] angles = anglePotentialEnergy.getAngleArray();
    addPerBondParameter("theta0");
    addPerBondParameter("k");
    setName("InPlaneAngle");

    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      for (Angle angle : angles) {
        AngleType.AngleMode angleMode = angle.angleType.angleMode;

        if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
          // Skip Normal angles unless this is ManyBody Titration.
        } else {
          double theta0 = angle.angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angle.angleType.angleUnit * angle.angleType.forceConstant;
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();
          int i4 = 0;
          if (angleMode == AngleType.AngleMode.NORMAL) {
            // This is a placeholder Angle, in case the Normal Angle is switched to an
            // In-Plane Angle during updateInPlaneAngleForce.
            k = 0.0;
            Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
            if (fourthAtom != null) {
              i4 = fourthAtom.getArrayIndex();
            } else {
              while (i1 == i4 || i2 == i4 || i3 == i4) {
                i4++;
              }
            }
          } else {
            i4 = angle.getAtom4().getArrayIndex();
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
    }

    if (nAngles > 0) {
      int forceGroup = anglePotentialEnergy.getForceGroup();
      setForceGroup(forceGroup);
      logger.info(format("  In-Plane Angles:                   %10d", nAngles));
      logger.fine(format("   Force Group:                      %10d", forceGroup));
    }
  }

  /**
   * Create an OpenMM In-Plane Angle Force for Dual Topology.
   *
   * @param anglePotentialEnergy      The AnglePotentialEnergy that contains the angles.
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public InPlaneAngleForce(AnglePotentialEnergy anglePotentialEnergy, int topology,
                           OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(4, anglePotentialEnergy.getInPlaneAngleEnergyString());
    Angle[] angles = anglePotentialEnergy.getAngleArray();
    addPerBondParameter("theta0");
    addPerBondParameter("k");
    setName("InPlaneAngle");

    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ForceField forceField = forceFieldEnergy.getMolecularAssembly().getForceField();
    manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      for (Angle angle : angles) {
        AngleType.AngleMode angleMode = angle.angleType.angleMode;

        if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
          // Skip Normal angles unless this is ManyBody Titration.
        } else {
          double theta0 = angle.angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angle.angleType.angleUnit * angle.angleType.forceConstant;
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();
          int i4 = 0;
          if (angleMode == AngleType.AngleMode.NORMAL) {
            k = 0.0;
            Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
            if (fourthAtom != null) {
              i4 = fourthAtom.getArrayIndex();
            } else {
              while (i1 == i4 || i2 == i4 || i3 == i4) {
                i4++;
              }
            }
          } else {
            i4 = angle.getAtom4().getArrayIndex();
          }
          // Don't apply lambda scale to alchemical in-plane angle
          if (!angle.applyLambda()) {
            k = k * scale;
          }
          i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
          i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
          i3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i3);
          i4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i4);
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
    }

    if (nAngles > 0) {
      int forceGroup = anglePotentialEnergy.getForceGroup();
      setForceGroup(forceGroup);
      logger.info(format("  In-Plane Angles:                   %10d", nAngles));
      logger.fine(format("   Force Group:                      %10d", forceGroup));
    }
  }

  /**
   * Convenience method to construct an OpenMM In-Plane Angle Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   * @return An In-Plane Angle Force, or null if there are no In-Plane angles.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    AnglePotentialEnergy anglePotentialEnergy = openMMEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return null;
    }
    InPlaneAngleForce inPlaneAngleForce = new InPlaneAngleForce(anglePotentialEnergy, openMMEnergy);
    if (inPlaneAngleForce.nAngles > 0) {
      return inPlaneAngleForce;
    }
    inPlaneAngleForce.destroy();
    return null;
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM In-Plane Angle Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An In-Plane Angle Force, or null if there are no In-Plane angles.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    AnglePotentialEnergy anglePotentialEnergy = forceFieldEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return null;
    }
    InPlaneAngleForce inPlaneAngleForce = new InPlaneAngleForce(anglePotentialEnergy, topology, openMMDualTopologyEnergy);
    if (inPlaneAngleForce.nAngles > 0) {
      return inPlaneAngleForce;
    }
    inPlaneAngleForce.destroy();
    return null;
  }

  /**
   * Update the In-Plane Angle force.
   *
   * @param openMMEnergy The OpenMM Energy instance.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    AnglePotentialEnergy anglePotentialEnergy = openMMEnergy.getAnglePotentialEnergy();
    if (anglePotentialEnergy == null) {
      return;
    }
    Angle[] angles = anglePotentialEnergy.getAngleArray();

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Angle angle : angles) {
        AngleType.AngleMode angleMode = angle.angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
          // Skip Normal angles unless this is ManyBody Titration.
        } else {
          double theta0 = angle.angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angle.angleType.angleUnit * angle.angleType.forceConstant;
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();
          int i4 = 0;
          if (angleMode == AngleType.AngleMode.NORMAL) {
            k = 0.0;
            Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
            if (fourthAtom != null) {
              i4 = fourthAtom.getArrayIndex();
            } else {
              while (i1 == i4 || i2 == i4 || i3 == i4) {
                i4++;
              }
            }
          } else {
            i4 = angle.getAtom4().getArrayIndex();
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
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the Dual-Topology In-Plane Angle force.
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

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Angle angle : angles) {
        AngleType.AngleMode angleMode = angle.angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleType.AngleMode.NORMAL) {
          // Skip Normal angles unless this is ManyBody Titration.
        } else {
          double theta0 = angle.angleType.angle[angle.nh];
          double k = KJ_PER_KCAL * angle.angleType.angleUnit * angle.angleType.forceConstant;
          int i1 = angle.getAtom(0).getArrayIndex();
          int i2 = angle.getAtom(1).getArrayIndex();
          int i3 = angle.getAtom(2).getArrayIndex();
          int i4 = 0;
          if (angleMode == AngleType.AngleMode.NORMAL) {
            k = 0.0;
            Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
            if (fourthAtom != null) {
              i4 = fourthAtom.getArrayIndex();
            } else {
              while (i1 == i4 || i2 == i4 || i3 == i4) {
                i4++;
              }
            }
          } else {
            i4 = angle.getAtom4().getArrayIndex();
          }
          // Don't apply lambda scale to alchemical in-plane angle
          if (!angle.applyLambda()) {
            k = k * scale;
          }
          i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
          i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
          i3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i3);
          i4 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i4);
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
    }

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }
}
