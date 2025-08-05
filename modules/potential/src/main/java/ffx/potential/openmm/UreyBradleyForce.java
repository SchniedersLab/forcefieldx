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

import ffx.openmm.Force;
import ffx.openmm.HarmonicBondForce;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.parameters.UreyBradleyType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;

/**
 * Urey-Bradley Force.
 */
public class UreyBradleyForce extends HarmonicBondForce {

  private static final Logger logger = Logger.getLogger(UreyBradleyForce.class.getName());

  /**
   * Urey-Bradly Force constructor.
   *
   * @param openMMEnergy The OpenMMEnergy instance that contains the Urey-Bradley terms.
   */
  public UreyBradleyForce(OpenMMEnergy openMMEnergy) {
    UreyBradley[] ureyBradleys = openMMEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return;
    }

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    for (UreyBradley ureyBradley : ureyBradleys) {
      int i1 = ureyBradley.getAtom(0).getArrayIndex();
      int i2 = ureyBradley.getAtom(2).getArrayIndex();
      UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
      double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
      // The implementation of UreyBradley in FFX & Tinker: k x^2
      // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
      double k = 2.0 * ureyBradleyType.forceConstant * ureyBradleyType.ureyUnit * kParameterConversion;
      addBond(i1, i2, length, k);
    }

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("UREY_BRADLEY_FORCE", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Urey-Bradleys:                     %10d", ureyBradleys.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Urey-Bradly Force constructor.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public UreyBradleyForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    UreyBradley[] ureyBradleys = forceFieldEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return;
    }

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    for (UreyBradley ureyBradley : ureyBradleys) {
      int i1 = ureyBradley.getAtom(0).getArrayIndex();
      int i2 = ureyBradley.getAtom(2).getArrayIndex();
      UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
      double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
      // The implementation of UreyBradley in FFX & Tinker: k x^2
      // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
      double k = 2.0 * ureyBradleyType.forceConstant * ureyBradleyType.ureyUnit * kParameterConversion;
      // Don't apply lambda scale to alchemical Urey-Bradley
      if (!ureyBradley.applyLambda()) {
        k = k * scale;
      }
      i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
      i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
      addBond(i1, i2, length, k);
    }

    int forceGroup = forceFieldEnergy.getMolecularAssembly().getForceField().getInteger("UREY_BRADLEY_FORCE", 0);
    setForceGroup(forceGroup);
    logger.info(format("  Urey-Bradleys:                     %10d", ureyBradleys.length));
    logger.fine(format("   Force Group:                      %10d", forceGroup));
  }

  /**
   * Convenience method to construct an OpenMM Urey-Bradley Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the angles.
   * @return An OpenMM Urey-Bradley Force, or null if there are no Urey-Bradley.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    UreyBradley[] ureyBradleys = openMMEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return null;
    }
    return new UreyBradleyForce(openMMEnergy);
  }

  /**
   * Convenience method to construct an OpenMM Urey-Bradley Force.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);

    UreyBradley[] ureyBradleys = forceFieldEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return null;
    }
    return new UreyBradleyForce(topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the Urey-Bradley parameters in the OpenMM Context.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the Urey-Bradley terms.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    UreyBradley[] ureyBradleys = openMMEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return;
    }

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    int index = 0;
    for (UreyBradley ureyBradley : ureyBradleys) {
      int i1 = ureyBradley.getAtom(0).getArrayIndex();
      int i2 = ureyBradley.getAtom(2).getArrayIndex();
      UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
      double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
      // The implementation of UreyBradley in FFX & Tinker: k x^2
      // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
      double k = 2.0 * ureyBradleyType.forceConstant * ureyBradleyType.ureyUnit * kParameterConversion;
      setBondParameters(index++, i1, i2, length, k);
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the Urey-Bradley parameters in the OpenMM Context.
   *
   * @param topology The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    UreyBradley[] ureyBradleys = forceFieldEnergy.getUreyBradleys();
    if (ureyBradleys == null || ureyBradleys.length < 1) {
      return;
    }

    double scale = openMMDualTopologyEnergy.getTopologyScale(topology);

    double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

    int index = 0;
    for (UreyBradley ureyBradley : ureyBradleys) {
      int i1 = ureyBradley.getAtom(0).getArrayIndex();
      int i2 = ureyBradley.getAtom(2).getArrayIndex();
      UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
      double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
      // The implementation of UreyBradley in FFX & Tinker: k x^2
      // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
      double k = 2.0 * ureyBradleyType.forceConstant * ureyBradleyType.ureyUnit * kParameterConversion;
      // Don't apply lambda scale to alchemical Urey-Bradley
      if (!ureyBradley.applyLambda()) {
        k = k * scale;
      }
      i1 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i1);
      i2 = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, i2);
      setBondParameters(index++, i1, i2, length, k);
    }

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
