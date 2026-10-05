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

import ffx.openmm.ffm.CustomBondForce;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.BondType;
import ffx.potential.terms.BondPotentialEnergy;

import java.util.logging.Logger;

import static ffx.openmm.ffm.OpenMMUnits.KJ_PER_KCAL;
import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static java.lang.String.format;

/**
 * OpenMM Bond Force backed by FFM {@link CustomBondForce}.
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

    double kParameterConversion = KJ_PER_KCAL / (NM_PER_ANGSTROM * NM_PER_ANGSTROM);
    try (DoubleArray parameters = new DoubleArray(0)) {
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getArrayIndex();
        int i2 = bond.getAtom(1).getArrayIndex();
        BondType bondType = bond.bondType;
        double r0 = bondType.distance * NM_PER_ANGSTROM;
        double k = kParameterConversion * bondType.forceConstant * bond.bondType.bondUnit;
        parameters.append(r0);
        parameters.append(k);
        addBond(i1, i2, parameters);
        parameters.resize(0);
      }
    }

    int forceGroup = bondPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Bonds:                             %10d", bonds.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create a Bond Force for Dual Topology.
   *
   * @param bondPotentialEnergy      BondPotentialEnergy that contains the Bond instances.
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public BondForce(BondPotentialEnergy bondPotentialEnergy, int topology,
                   OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    super(bondPotentialEnergy.getBondEnergyString());
    Bond[] bonds = bondPotentialEnergy.getBondArray();
    addPerBondParameter("r0");
    addPerBondParameter("k");
    setName("AmoebaBond");

    double kParameterConversion = KJ_PER_KCAL / (NM_PER_ANGSTROM * NM_PER_ANGSTROM);
    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (DoubleArray parameters = new DoubleArray(0)) {
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getArrayIndex();
        int i2 = bond.getAtom(1).getArrayIndex();
        BondType bondType = bond.bondType;
        double r0 = bondType.distance * NM_PER_ANGSTROM;
        double k = kParameterConversion * bondType.forceConstant * bond.bondType.bondUnit;
        // Don't apply lambda scale to alchemical bond
        if (!bond.applyLambda()) {
          k = k * scale;
        }
        i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
        i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
        parameters.append(r0);
        parameters.append(k);
        addBond(i1, i2, parameters);
        parameters.resize(0);
      }
    }

    int forceGroup = bondPotentialEnergy.getForceGroup();
    setForceGroup(forceGroup);
    logger.info(format("  Bonds:                             %10d", bonds.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Create a bond force for the OpenMM System.
   *
   * @param openMMEnergy OpenMM Energy that contains the Bond instances.
   * @return Force instance or null if no bond potential energy exists.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    BondPotentialEnergy bondPotentialEnergy = openMMEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return null;
    }
    return new BondForce(bondPotentialEnergy);
  }

  /**
   * Convenience method to construct a Dual-Topology OpenMM Bond Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   * @return An OpenMM Bond Force, or null if there are no bonds.
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
   *
   * @param openMMEnergy OpenMM Energy that contains the Bond instances.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    BondPotentialEnergy bondPotentialEnergy = openMMEnergy.getBondPotentialEnergy();
    if (bondPotentialEnergy == null) {
      return;
    }
    Bond[] bonds = bondPotentialEnergy.getBondArray();

    double kParameterConversion = KJ_PER_KCAL / (NM_PER_ANGSTROM * NM_PER_ANGSTROM);
    try (DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getArrayIndex();
        int i2 = bond.getAtom(1).getArrayIndex();
        BondType bondType = bond.bondType;
        double r0 = bondType.distance * NM_PER_ANGSTROM;
        double k = kParameterConversion * bondType.forceConstant * bondType.bondUnit;
        parameters.append(r0);
        parameters.append(k);
        setBondParameters(index++, i1, i2, parameters);
        parameters.resize(0);
      }
    }
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update existing bond force for the Dual-Topology OpenMM System.
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

    double kParameterConversion = KJ_PER_KCAL / (NM_PER_ANGSTROM * NM_PER_ANGSTROM);
    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    try (DoubleArray parameters = new DoubleArray(0)) {
      int index = 0;
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getArrayIndex();
        int i2 = bond.getAtom(1).getArrayIndex();
        BondType bondType = bond.bondType;
        double r0 = bondType.distance * NM_PER_ANGSTROM;
        double k = kParameterConversion * bondType.forceConstant * bondType.bondUnit;
        // Don't apply lambda scale to alchemical bond
        if (!bond.applyLambda()) {
          k = k * scale;
        }
        i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
        i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
        parameters.append(r0);
        parameters.append(k);
        setBondParameters(index++, i1, i2, parameters);
        parameters.resize(0);
      }
    }
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }
}
