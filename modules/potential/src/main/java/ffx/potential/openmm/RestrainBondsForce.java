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
import ffx.openmm.CustomBondForce;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.RestrainDistance;
import ffx.potential.parameters.BondType;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;

/**
 * Restrain Bonds Force.
 */
public class RestrainBondsForce extends CustomBondForce {

  private static final Logger logger = Logger.getLogger(RestrainBondsForce.class.getName());

  /**
   * Restrain Bond Force constructor.
   *
   * @param bondFunction The bond function.
   * @param openMMEnergy The OpenMM Energy.
   */
  public RestrainBondsForce(BondType.BondFunction bondFunction, OpenMMEnergy openMMEnergy) {
    super(bondFunction.toMathematicalForm());

    List<RestrainDistance> restrainDistances = openMMEnergy.getRestrainDistances(bondFunction);
    if (restrainDistances == null || restrainDistances.isEmpty()) {
      destroy();
      return;
    }

    addPerBondParameter("k");
    addPerBondParameter("r0");
    if (bondFunction.hasFlatBottom()) {
      addPerBondParameter("fb");
    }

    BondType bondType = restrainDistances.getFirst().bondType;
    switch (bondFunction) {
      case QUARTIC, FLAT_BOTTOM_QUARTIC -> {
        addGlobalParameter("cubic", bondType.cubic / OpenMM_NmPerAngstrom);
        addGlobalParameter("quartic", bondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
      }
    }

    // OpenMM's HarmonicBondForce class uses k, not 1/2*k as does FFX.
    double forceConvert = 2.0 * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
    DoubleArray parameters = new DoubleArray(0);
    for (RestrainDistance restrainDistance : restrainDistances) {
      bondType = restrainDistance.bondType;
      double forceConstant = bondType.forceConstant * bondType.bondUnit * forceConvert;
      double distance = bondType.distance * OpenMM_NmPerAngstrom;
      Atom[] atoms = restrainDistance.getAtomArray();
      int i1 = atoms[0].getXyzIndex() - 1;
      int i2 = atoms[1].getXyzIndex() - 1;
      parameters.append(forceConstant);
      parameters.append(distance);
      if (bondFunction.hasFlatBottom()) {
        parameters.append(bondType.flatBottomRadius * OpenMM_NmPerAngstrom);
      }
      addBond(i1, i2, parameters);
      parameters.destroy();
    }
    parameters.destroy();

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("BOND_RESTRAINT_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Restraint bonds force \t%6d\t%d", restrainDistances.size(), forceGroup));
  }

  /**
   * Add a Restrain-Bond force to the OpenMM System.
   *
   * @param bondFunction The bond function.
   * @param openMMEnergy The OpenMM Energy.
   */
  public static Force constructForce(BondType.BondFunction bondFunction, OpenMMEnergy openMMEnergy) {
    List<RestrainDistance> restrainDistances = openMMEnergy.getRestrainDistances(bondFunction);
    if (restrainDistances == null || restrainDistances.isEmpty()) {
      return null;
    }
    return new RestrainBondsForce(bondFunction, openMMEnergy);
  }
}
