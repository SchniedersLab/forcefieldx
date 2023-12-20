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

import ffx.openmm.CustomIntegrator;

import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.copyOfRange;

/**
 * OpenMM Custom MTS Integrator.
 */
public class CustomMTSIntegrator extends CustomIntegrator {

  private static final Logger logger = Logger.getLogger(CustomMTSIntegrator.class.getName());

  public CustomMTSIntegrator(double dt, double constraintTolerance, boolean hasAmoebaCavitationForce) {
    super(dt);
    setConstraintTolerance(constraintTolerance);

    int n = 4;
    // Force group 1 contains slowly varying forces.
    // Force group 0 contains the fast varying forces.
    int[] forceGroups = {1, 0};
    // There will be 1 force evaluation per outer step, and 4 per inner step.
    int[] subSteps = {1, 4};
    if (hasAmoebaCavitationForce) {
      n = 8;
      // Force group 2 contains the cavitation force.
      // Force group 1 contains slowly varying forces.
      // Force group 0 contains the fast varying forces.
      forceGroups = new int[]{2, 1, 0};
      // There will be 1 force evaluation per outer step.
      // There will be 2 force evaluations per middle step.
      // There will be 8 force evaluations per inner step.
      subSteps = new int[]{1, 2, 8};
    }

    addPerDofVariable("x1", 0.0);
    addUpdateContextState();
    createMTSSubStep(1, forceGroups, subSteps);
    addConstrainVelocities();
    logger.info("  Custom MTS Integrator");
    logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
    logger.info(format("  Inner Time step:      %6.2f (fsec)", dt / n * 1000));
  }

  /**
   * Create substeps for the MTS CustomIntegrator.
   *
   * @param parentSubsteps The number of substeps for the previous force group.
   * @param forceGroups    The force groups to be evaluated.
   * @param subSteps       The number of substeps for each force group.
   */
  public void createMTSSubStep(int parentSubsteps, int[] forceGroups, int[] subSteps) {
    int forceGroup = forceGroups[0];
    int steps = subSteps[0];
    int stepsPerParentStep = steps / parentSubsteps;
    if (stepsPerParentStep < 1 || steps % parentSubsteps != 0) {
      throw new IllegalArgumentException(
          "The number for substeps for each group must be a multiple of the number for the previous group");
    }
    if (forceGroup < 0 || forceGroup > 31) {
      throw new IllegalArgumentException("Force group must be between 0 and 31");
    }
    for (int i = 0; i < stepsPerParentStep; i++) {
      addComputePerDof("v", "v+0.5*(dt/" + steps + ")*f" + forceGroup + "/m");
      if (forceGroups.length == 1) {
        addComputePerDof("x", "x+(dt/" + steps + ")*v");
        addComputePerDof("x1", "x");
        addConstrainPositions();
        addComputePerDof("v", "v+(x-x1)/(dt/" + steps + ")");
        addConstrainVelocities();
      } else {
        createMTSSubStep(steps, copyOfRange(forceGroups, 1, forceGroups.length),
            copyOfRange(subSteps, 1, subSteps.length));
      }
      addComputePerDof("v", "v+0.5*(dt/" + steps + ")*f" + forceGroup + "/m");
    }
  }
}
