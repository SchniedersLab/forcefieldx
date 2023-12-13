// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.parameters.ImproperTorsionType;

import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static java.lang.String.format;

/**
 * OpenMM Improper Torsion Force.
 */
public class ImproperTorsionForce extends OpenMMPeriodicTorsionForce {

  private static final Logger logger = Logger.getLogger(ImproperTorsionForce.class.getName());

  private double lambdaTorsion = 1.0;

  /**
   * Create an OpenMM Improper Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the improper torsions.
   */
  public ImproperTorsionForce(OpenMMEnergy openMMEnergy) {
    ImproperTorsion[] improperTorsions = openMMEnergy.getImproperTorsions();
    if (improperTorsions == null || improperTorsions.length < 1) {
      // Clean up the memory allocated by the OpenMMPeriodicTorsionForce constructor.
      destroy();
      return;
    }

    for (ImproperTorsion improperTorsion : improperTorsions) {
      int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
      ImproperTorsionType type = improperTorsion.improperType;
      double forceConstant = OpenMM_KJPerKcal * type.impTorUnit * improperTorsion.scaleFactor * type.k;
      addTorsion(a1, a2, a3, a4, type.periodicity, type.phase * OpenMM_RadiansPerDegree, forceConstant);
    }

    int forceGroup = openMMEnergy.getMolecularAssembly().getForceField().getInteger("IMPROPER_TORSION_FORCE_GROUP", 0);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  Improper Torsions \t%6d\t\t%1d", improperTorsions.length, forceGroup));
  }

  /**
   * Set the lambda torsion scale factor.
   *
   * @param lambdaTorsion The lambda torsion scale factor.
   */
  public void setLambdaTorsion(double lambdaTorsion) {
    this.lambdaTorsion = lambdaTorsion;
  }

  /**
   * Convenience method to construct an OpenMM Improper Torsion Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the improper torsions.
   * @return An Improper Torsion Force, or null if there are no improper torsions.
   */
  public static OpenMMForce constructForce(OpenMMEnergy openMMEnergy) {
    ImproperTorsion[] improperTorsions = openMMEnergy.getImproperTorsions();
    if (improperTorsions == null || improperTorsions.length < 1) {
      return null;
    }
    return new ImproperTorsionForce(openMMEnergy);
  }

  /**
   * Update the Improper Torsion force.
   *
   * @param openMMEnergy The OpenMM Energy that contains the improper torsions.
   */
  public void updateForce(OpenMMEnergy openMMEnergy) {
    ImproperTorsion[] improperTorsions = openMMEnergy.getImproperTorsions();
    if (improperTorsions == null || improperTorsions.length < 1) {
      return;
    }

    int nImproperTorsions = improperTorsions.length;
    for (int i = 0; i < nImproperTorsions; i++) {
      ImproperTorsion improperTorsion = improperTorsions[i];
      int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
      int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
      int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
      int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
      ImproperTorsionType type = improperTorsion.improperType;
      double forceConstant = OpenMM_KJPerKcal * type.impTorUnit * improperTorsion.scaleFactor * type.k * lambdaTorsion;
      setTorsionParameters(i, a1, a2, a3, a4, type.periodicity, type.phase * OpenMM_RadiansPerDegree, forceConstant);
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

}
