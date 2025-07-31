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
import ffx.potential.bonded.StretchBend;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static java.lang.String.format;

/**
 * OpenMM Stretch-Bend Force.
 */
public class StretchBendForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(StretchBendForce.class.getName());

  /**
   * Create an OpenMM Stretch-Bend Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the Stretch-Bends.
   */
  public StretchBendForce(OpenMMEnergy openMMEnergy) {
    super(3, openMMEnergy.getStretchBendEnergyString());
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return;
    }
    addPerBondParameter("r12");
    addPerBondParameter("r23");
    addPerBondParameter("theta0");
    addPerBondParameter("k1");
    addPerBondParameter("k2");
    setName("AmoebaStretchBend");

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (StretchBend stretchBend : stretchBends) {
      int i1 = stretchBend.getAtom(0).getArrayIndex();
      int i2 = stretchBend.getAtom(1).getArrayIndex();
      int i3 = stretchBend.getAtom(2).getArrayIndex();
      double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
      double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
      double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
      double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
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
    particles.destroy();
    parameters.destroy();

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("STRETCH_BEND_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Stretch-Bends:                     %10d", stretchBends.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create an OpenMM Stretch-Bend Force for Dual Topology.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public StretchBendForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(3, openMMDualTopologyEnergy.getOpenMMEnergy(topology).getStretchBendEnergyString());

    OpenMMEnergy openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(topology);
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return;
    }
    addPerBondParameter("r12");
    addPerBondParameter("r23");
    addPerBondParameter("theta0");
    addPerBondParameter("k1");
    addPerBondParameter("k2");
    setName("AmoebaStretchBend");

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (StretchBend stretchBend : stretchBends) {
      int i1 = stretchBend.getAtom(0).getArrayIndex();
      int i2 = stretchBend.getAtom(1).getArrayIndex();
      int i3 = stretchBend.getAtom(2).getArrayIndex();
      double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
      double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
      double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
      double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      // Don't apply lambda scale to alchemical stretch bend
      if (!stretchBend.applyLambda()) {
        k1 = k1 * scale; // todo how to do scale - scale both by lambda or square root of lambda - look at equation
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
    particles.destroy();
    parameters.destroy();

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("STRETCH_BEND_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Stretch-Bends:                     %10d", stretchBends.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Stretch-Bend Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the stretch-bends.
   * @return An OpenMM Stretch-Bend Force, or null if there are no stretch-bends.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return null;
    }
    return new StretchBendForce(openMMEnergy);
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Stretch-Bend Force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Stretch-Bend Force, or null if there are no stretch-bends.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    OpenMMEnergy openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(topology);
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return null;
    }
    return new StretchBendForce(topology, openMMDualTopologyEnergy);
  }

  /**
   * Update this Stretch-Bend Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the stretch-bends.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return;
    }

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (StretchBend stretchBend : stretchBends) {
      int i1 = stretchBend.getAtom(0).getArrayIndex();
      int i2 = stretchBend.getAtom(1).getArrayIndex();
      int i3 = stretchBend.getAtom(2).getArrayIndex();
      double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
      double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
      double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
      double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
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
    particles.destroy();
    parameters.destroy();

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update existing Stretch-Bend Force for the Dual-Topology OpenMM System.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    OpenMMEnergy openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(topology);
    StretchBend[] stretchBends = openMMEnergy.getStretchBends();
    if (stretchBends == null || stretchBends.length < 1) {
      return;
    }

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (StretchBend stretchBend : stretchBends) {
      int i1 = stretchBend.getAtom(0).getArrayIndex();
      int i2 = stretchBend.getAtom(1).getArrayIndex();
      int i3 = stretchBend.getAtom(2).getArrayIndex();
      double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
      double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
      double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
      double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
      // Don't apply lambda scale to alchemical stretch bend
      if (!stretchBend.applyLambda()) {
        k1 = k1 * scale; // todo - how to apply?
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
    particles.destroy();
    parameters.destroy();

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }
}
