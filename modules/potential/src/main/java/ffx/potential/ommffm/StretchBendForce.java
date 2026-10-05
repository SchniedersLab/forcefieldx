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
import ffx.potential.bonded.StretchBend;
import ffx.potential.terms.StretchBendPotentialEnergy;

import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static ffx.openmm.ffm.OpenMMUnits.RADIANS_PER_DEGREE;
import static java.lang.String.format;

/**
 * OpenMM Stretch-Bend Force backed by FFM {@link CustomCompoundBondForce}.
 */
public class StretchBendForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(StretchBendForce.class.getName());

  /**
   * Create an OpenMM Stretch-Bend Force.
   *
   * @param stretchBendPotentialEnergy The StretchBendPotentialEnergy instance that contains the stretch-bends.
   */
  public StretchBendForce(StretchBendPotentialEnergy stretchBendPotentialEnergy) {
    super(3, StretchBendPotentialEnergy.getStretchBendEnergyString());
    StretchBend[] stretchBends = stretchBendPotentialEnergy.getStretchBendArray();
    addPerBondParameter("r12");
    addPerBondParameter("r23");
    addPerBondParameter("theta0");
    addPerBondParameter("k1");
    addPerBondParameter("k2");
    setName("AmoebaStretchBend");

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      for (StretchBend stretchBend : stretchBends) {
        int i1 = stretchBend.getAtom(0).getArrayIndex();
        int i2 = stretchBend.getAtom(1).getArrayIndex();
        int i3 = stretchBend.getAtom(2).getArrayIndex();
        double r12 = stretchBend.bond0Eq * NM_PER_ANGSTROM;
        double r23 = stretchBend.bond1Eq * NM_PER_ANGSTROM;
        double theta0 = stretchBend.angleEq * RADIANS_PER_DEGREE;
        double k1 = stretchBend.force0 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        double k2 = stretchBend.force1 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        particles.append(i1);
        particles.append(i2);
        particles.append(i3);
        parameters.append(r12);
        parameters.append(r23);
        parameters.append(theta0);
        parameters.append(k1);
        parameters.append(k2);
        addBond(particles, parameters);
        particles.resize(0);
        parameters.resize(0);
      }
    }

    int forceGroup = stretchBendPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Stretch-Bends:                     %10d", stretchBends.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Stretch-Bend Force.
   *
   * @param stretchBendPotentialEnergy The StretchBendPotentialEnergy instance that contains the stretch-bends.
   * @param topology                   The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy   The OpenMMDualTopologyEnergy instance.
   */
  public StretchBendForce(StretchBendPotentialEnergy stretchBendPotentialEnergy,
                          int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(3, StretchBendPotentialEnergy.getStretchBendEnergyString());
    StretchBend[] stretchBends = stretchBendPotentialEnergy.getStretchBendArray();
    addPerBondParameter("r12");
    addPerBondParameter("r23");
    addPerBondParameter("theta0");
    addPerBondParameter("k1");
    addPerBondParameter("k2");
    setName("AmoebaStretchBend");

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      for (StretchBend stretchBend : stretchBends) {
        int i1 = stretchBend.getAtom(0).getArrayIndex();
        int i2 = stretchBend.getAtom(1).getArrayIndex();
        int i3 = stretchBend.getAtom(2).getArrayIndex();
        double r12 = stretchBend.bond0Eq * NM_PER_ANGSTROM;
        double r23 = stretchBend.bond1Eq * NM_PER_ANGSTROM;
        double theta0 = stretchBend.angleEq * RADIANS_PER_DEGREE;
        double k1 = stretchBend.force0 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        double k2 = stretchBend.force1 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        // Don't apply lambda scale to alchemical stretch bend
        if (!stretchBend.applyLambda()) {
          k1 = k1 * scale;
          k2 = k2 * scale;
        }
        i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
        i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
        i3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i3);
        particles.append(i1);
        particles.append(i2);
        particles.append(i3);
        parameters.append(r12);
        parameters.append(r23);
        parameters.append(theta0);
        parameters.append(k1);
        parameters.append(k2);
        addBond(particles, parameters);
        particles.resize(0);
        parameters.resize(0);
      }
    }

    int forceGroup = stretchBendPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Stretch-Bends:                     %10d", stretchBends.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Stretch-Bend Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the stretch-bends.
   * @return A Stretch-Bend Force, or null if there are no stretch-bends.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    StretchBendPotentialEnergy stretchBendPotentialEnergy =
        openMMEnergy.getStretchBendPotentialEnergy();
    if (stretchBendPotentialEnergy == null) {
      return null;
    }
    return new StretchBendForce(stretchBendPotentialEnergy);
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Stretch-Bend Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Stretch-Bend Force, or null if there are no stretch-bends.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    StretchBendPotentialEnergy stretchBendPotentialEnergy = forceFieldEnergy.getStretchBendPotentialEnergy();
    if (stretchBendPotentialEnergy == null) {
      return null;
    }
    return new StretchBendForce(stretchBendPotentialEnergy, topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the Stretch-Bend parameters in the OpenMM Context.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the stretch-bends.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    StretchBendPotentialEnergy stretchBendPotentialEnergy =
        openMMEnergy.getStretchBendPotentialEnergy();
    if (stretchBendPotentialEnergy == null) {
      return;
    }
    StretchBend[] stretchBends = stretchBendPotentialEnergy.getStretchBendArray();

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (StretchBend stretchBend : stretchBends) {
        int i1 = stretchBend.getAtom(0).getArrayIndex();
        int i2 = stretchBend.getAtom(1).getArrayIndex();
        int i3 = stretchBend.getAtom(2).getArrayIndex();
        double r12 = stretchBend.bond0Eq * NM_PER_ANGSTROM;
        double r23 = stretchBend.bond1Eq * NM_PER_ANGSTROM;
        double theta0 = stretchBend.angleEq * RADIANS_PER_DEGREE;
        double k1 = stretchBend.force0 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        double k2 = stretchBend.force1 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        particles.append(i1);
        particles.append(i2);
        particles.append(i3);
        parameters.append(r12);
        parameters.append(r23);
        parameters.append(theta0);
        parameters.append(k1);
        parameters.append(k2);
        setBondParameters(index++, particles, parameters);
        particles.resize(0);
        parameters.resize(0);
      }
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update existing Stretch-Bend Force for the Dual-Topology OpenMM System.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    StretchBendPotentialEnergy stretchBendPotentialEnergy = forceFieldEnergy.getStretchBendPotentialEnergy();
    if (stretchBendPotentialEnergy == null) {
      return;
    }
    StretchBend[] stretchBends = stretchBendPotentialEnergy.getStretchBendArray();

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (IntArray particles = new IntArray(0);
         DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (StretchBend stretchBend : stretchBends) {
        int i1 = stretchBend.getAtom(0).getArrayIndex();
        int i2 = stretchBend.getAtom(1).getArrayIndex();
        int i3 = stretchBend.getAtom(2).getArrayIndex();
        double r12 = stretchBend.bond0Eq * NM_PER_ANGSTROM;
        double r23 = stretchBend.bond1Eq * NM_PER_ANGSTROM;
        double theta0 = stretchBend.angleEq * RADIANS_PER_DEGREE;
        double k1 = stretchBend.force0 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        double k2 = stretchBend.force1 * KJ_PER_KCAL / NM_PER_ANGSTROM;
        // Don't apply lambda scale to alchemical stretch bend
        if (!stretchBend.applyLambda()) {
          k1 = k1 * scale;
          k2 = k2 * scale;
        }
        i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
        i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
        i3 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i3);
        particles.append(i1);
        particles.append(i2);
        particles.append(i3);
        parameters.append(r12);
        parameters.append(r23);
        parameters.append(theta0);
        parameters.append(k1);
        parameters.append(k2);
        setBondParameters(index++, particles, parameters);
        particles.resize(0);
        parameters.resize(0);
      }
    }

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }
}
