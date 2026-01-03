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
package ffx.potential;

import ffx.numerics.Potential;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.parameters.ForceField;

/**
 * FiniteDifference
 *
 * @version 1.0
 */
public class FiniteDifferenceUtils {

  /**
   * Compute dE/dL using finite differences.
   *
   * @param potential       The potential to compute the energy of.
   * @param lambdaInterface The lambda interface to use.
   * @param forceField      The force field to use.
   * @return The computed dE/dL.
   */
  public static double computedEdL(Potential potential, LambdaInterface lambdaInterface, ForceField forceField) {
    // Small optimization to only create the x array once.
    int n = potential.getNumberOfVariables();
    double[] x = new double[n];
    potential.getCoordinates(x);

    double finiteDifferenceStepSize = forceField.getDouble("FD_DLAMBDA", 0.001);
    boolean twoSidedFiniteDifference = forceField.getBoolean("FD_TWO_SIDED", true);

    double currentLambda = lambdaInterface.getLambda();
    double width = finiteDifferenceStepSize;
    double ePlus;
    double eMinus;

    if (twoSidedFiniteDifference) {
      if (currentLambda + finiteDifferenceStepSize > 1.0) {
        lambdaInterface.setLambda(currentLambda - finiteDifferenceStepSize);
        eMinus = potential.energy(x);
        lambdaInterface.setLambda(currentLambda);
        ePlus = potential.energy(x);
      } else if (currentLambda - finiteDifferenceStepSize < 0.0) {
        lambdaInterface.setLambda(currentLambda + finiteDifferenceStepSize);
        ePlus = potential.energy(x);
        lambdaInterface.setLambda(currentLambda);
        eMinus = potential.energy(x);
      } else {
        // Two sided finite difference estimate of dE/dL.
        lambdaInterface.setLambda(currentLambda + finiteDifferenceStepSize);
        ePlus = potential.energy(x);
        lambdaInterface.setLambda(currentLambda - finiteDifferenceStepSize);
        eMinus = potential.energy(x);
        width *= 2.0;
        lambdaInterface.setLambda(currentLambda);
      }
    } else {
      // One-sided finite difference estimates of dE/dL
      if (currentLambda + finiteDifferenceStepSize > 1.0) {
        lambdaInterface.setLambda(currentLambda - finiteDifferenceStepSize);
        eMinus = potential.energy(x);
        lambdaInterface.setLambda(currentLambda);
        ePlus = potential.energy(x);
      } else {
        lambdaInterface.setLambda(currentLambda + finiteDifferenceStepSize);
        ePlus = potential.energy(x);
        lambdaInterface.setLambda(currentLambda);
        eMinus = potential.energy(x);
      }
    }

    // Compute the finite difference derivative.
    double dEdL = (ePlus - eMinus) / width;

    // logger.info(format(" getdEdL currentLambda: CL=%8.6f L=%8.6f dEdL=%12.6f", currentLambda,
    // lambda, dEdL));
    return dEdL;

  }

}
