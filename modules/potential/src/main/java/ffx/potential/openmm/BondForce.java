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

import ffx.openmm.CustomBondForce;
import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.BondType;
import ffx.potential.terms.BondPotentialEnergy;

import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;

/**
 * Bond Force.
 */
public class BondForce extends CustomBondForce {

  private static final Logger logger = Logger.getLogger(BondForce.class.getName());

  /**
   * Bond Force constructor.
   *
   * @param bondPotentialEnergy BondPotentialEnergy that contains the Bond instances
   */
  public BondForce(BondPotentialEnergy bondPotentialEnergy) {
    super(bondPotentialEnergy.getBondEnergyString());
    Bond[] bonds = bondPotentialEnergy.getBondArray();
    addPerBondParameter("r0");
    addPerBondParameter("k");
    setName("AmoebaBond");

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    DoubleArray parameters = new DoubleArray(0);
    for (Bond bond : bonds) {
      int i1 = bond.getAtom(0).getArrayIndex();
      int i2 = bond.getAtom(1).getArrayIndex();
      BondType bondType = bond.bondType;
      double r0 = bondType.distance * OpenMM_NmPerAngstrom;
      double k = kParameterConversion * bondType.forceConstant * bond.bondType.bondUnit;
      parameters.append(r0);
      parameters.append(k);
      addBond(i1, i2, parameters);
      parameters.resize(0);
    }
    parameters.destroy();

    int forceGroup = bondPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Bonds:                             %10d", bonds.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Bond Force constructor.
   *
   * @param bondPotentialEnergy      BondPotentialEnergy that contains the Bond instances
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public BondForce(BondPotentialEnergy bondPotentialEnergy,
                   int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(bondPotentialEnergy.getBondEnergyString());
    Bond[] bonds = bondPotentialEnergy.getBondArray();
    addPerBondParameter("r0");
    addPerBondParameter("k");
    setName("AmoebaBond");

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    DoubleArray parameters = new DoubleArray(0);
    for (Bond bond : bonds) {
      int i1 = bond.getAtom(0).getArrayIndex();
      int i2 = bond.getAtom(1).getArrayIndex();
      BondType bondType = bond.bondType;
      double r0 = bondType.distance * OpenMM_NmPerAngstrom;
      double k = kParameterConversion * bondType.forceConstant * bond.bondType.bondUnit;
      // Don't apply lambda scale to alchemical bond
      if (!bond.applyLambda()) {
        k = k * scale;
      }
      parameters.append(r0);
      parameters.append(k);
      i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
      i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
      addBond(i1, i2, parameters);
      parameters.resize(0);
    }
    parameters.destroy();

    int forceGroup = bondPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Bonds:                             %10d", bonds.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Creat a bond force for the OpenMM System.
   *
   * @param openMMEnergy OpenMM Energy that contains the Bond instances.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    BondPotentialEnergy bondPotentialEnergy = openMMEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return null;
    }
    return new BondForce(bondPotentialEnergy);
  }

  /**
   * Add a bond force to the OpenMM System.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    BondPotentialEnergy bondPotentialEnergy = forceFieldEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return null;
    }
    return new BondForce(bondPotentialEnergy, topology, openMMDualTopologyEnergy);
  }

  /**
   * Update an existing bond force for the OpenMM System.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    BondPotentialEnergy bondPotentialEnergy = openMMEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return;
    }
    Bond[] bonds = bondPotentialEnergy.getBondArray();

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (Bond bond : bonds) {
      int i1 = bond.getAtom(0).getArrayIndex();
      int i2 = bond.getAtom(1).getArrayIndex();
      BondType bondType = bond.bondType;
      double r0 = bondType.distance * OpenMM_NmPerAngstrom;
      double k = kParameterConversion * bondType.forceConstant * bondType.bondUnit;
      parameters.append(r0);
      parameters.append(k);
      setBondParameters(index++, i1, i2, parameters);
      parameters.resize(0);
    }
    parameters.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update an existing bond force for the OpenMM System.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    BondPotentialEnergy bondPotentialEnergy = forceFieldEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return;
    }
    Bond[] bonds = bondPotentialEnergy.getBondArray();

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    DoubleArray parameters = new DoubleArray(0);
    int index = 0;
    for (Bond bond : bonds) {
      int i1 = bond.getAtom(0).getArrayIndex();
      int i2 = bond.getAtom(1).getArrayIndex();
      BondType bondType = bond.bondType;
      double r0 = bondType.distance * OpenMM_NmPerAngstrom;
      double k = kParameterConversion * bondType.forceConstant * bondType.bondUnit;
      // Don't apply lambda scale to alchemical bond
      if (!bond.applyLambda()) {
        k = k * scale;
      }
      parameters.append(r0);
      parameters.append(k);
      i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
      i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
      setBondParameters(index++, i1, i2, parameters);
      parameters.resize(0);
    }
    parameters.destroy();
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
