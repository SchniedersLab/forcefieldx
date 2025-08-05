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
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.parameters.OutOfPlaneBendType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static java.lang.String.format;

/**
 * OpenMM Out-of-Plane Bend Force.
 */
public class OutOfPlaneBendForce extends CustomCompoundBondForce {

  private static final Logger logger = Logger.getLogger(OutOfPlaneBendForce.class.getName());

  /**
   * Create an Out-of-Plane Bend Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the out-of-plane bends.
   */
  public OutOfPlaneBendForce(OpenMMEnergy openMMEnergy) {
    super(4, openMMEnergy.getOutOfPlaneEnergyString());
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return;
    }

    addPerBondParameter("k");
    setName("OutOfPlaneBend");

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
      OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
      int i1 = outOfPlaneBend.getAtom(0).getArrayIndex();
      int i2 = outOfPlaneBend.getAtom(1).getArrayIndex();
      int i3 = outOfPlaneBend.getAtom(2).getArrayIndex();
      int i4 = outOfPlaneBend.getAtom(3).getArrayIndex();
      double k = OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * outOfPlaneBendType.opBendUnit;
      particles.append(i1);
      particles.append(i2);
      particles.append(i3);
      particles.append(i4);
      parameters.append(k);
      addBond(particles, parameters);
      particles.resize(0);
      parameters.resize(0);
    }
    particles.destroy();
    parameters.destroy();
    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("OUT_OF_PLANE_BEND_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Out-of-Plane Bends:                %10d", outOfPlaneBends.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create an Out-of-Plane Bend Force for Dual Topology.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public OutOfPlaneBendForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(4, openMMDualTopologyEnergy.getOpenMMEnergy(topology).getOutOfPlaneEnergyString());

    OpenMMEnergy openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(topology);
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return;
    }

    addPerBondParameter("k");
    setName("OutOfPlaneBend");

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
      OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
      int i1 = outOfPlaneBend.getAtom(0).getArrayIndex();
      int i2 = outOfPlaneBend.getAtom(1).getArrayIndex();
      int i3 = outOfPlaneBend.getAtom(2).getArrayIndex();
      int i4 = outOfPlaneBend.getAtom(3).getArrayIndex();
      double k = OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * outOfPlaneBendType.opBendUnit;
      // Don't apply lambda scale to alchemcial out-of-plane bend - todo not sure if this is required for this bonded force
      if (!outOfPlaneBend.applyLambda()) {
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
      parameters.append(k);
      addBond(particles, parameters);
      particles.resize(0);
      parameters.resize(0);
    }
    particles.destroy();
    parameters.destroy();
    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("OUT_OF_PLANE_BEND_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Out-of-Plane Bends:                %10d", outOfPlaneBends.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Out-of-Plane Bend Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the out-of-plane bends.
   * @return An OpenMM Out-of-Plane Bend Force, or null if there are no out-of-plane bends.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return null;
    }
    return new OutOfPlaneBendForce(openMMEnergy);
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Out-of-Plane Bend Force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Out-of-Plane Bend Force, or null if there are no out-of-plane bends.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    OpenMMEnergy openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(topology);
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return null;
    }
    return new OutOfPlaneBendForce(topology, openMMDualTopologyEnergy);
  }

  /**
   * Update an existing angle force for the OpenMM System.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return;
    }

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
      OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
      int i1 = outOfPlaneBend.getAtom(0).getArrayIndex();
      int i2 = outOfPlaneBend.getAtom(1).getArrayIndex();
      int i3 = outOfPlaneBend.getAtom(2).getArrayIndex();
      int i4 = outOfPlaneBend.getAtom(3).getArrayIndex();
      double k = OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * outOfPlaneBendType.opBendUnit;
      particles.append(i1);
      particles.append(i2);
      particles.append(i3);
      particles.append(i4);
      parameters.append(k);
      setBondParameters(index++, particles, parameters);
      particles.resize(0);
      parameters.resize(0);
    }
    particles.destroy();
    parameters.destroy();

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update an existing angle force for the Dual-Topology OpenMM System.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    OpenMMEnergy openMMEnergy = openMMDualTopologyEnergy.getOpenMMEnergy(topology);
    OutOfPlaneBend[] outOfPlaneBends = openMMEnergy.getOutOfPlaneBends();
    if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
      return;
    }

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    IntArray particles = new IntArray(0);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
      OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
      int i1 = outOfPlaneBend.getAtom(0).getArrayIndex();
      int i2 = outOfPlaneBend.getAtom(1).getArrayIndex();
      int i3 = outOfPlaneBend.getAtom(2).getArrayIndex();
      int i4 = outOfPlaneBend.getAtom(3).getArrayIndex();
      double k = OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * outOfPlaneBendType.opBendUnit;
      // Don't apply lambda scale to alchemcial out-of-plane bend - todo not sure if this is required for this bonded force
      if (!outOfPlaneBend.applyLambda()) {
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
      parameters.append(k);
      setBondParameters(index++, particles, parameters);
      particles.resize(0);
      parameters.resize(0);
    }
    particles.destroy();
    parameters.destroy();

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }
}
