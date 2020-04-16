// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.optimize;

import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.algorithms.AlgorithmListener;
import ffx.numerics.Potential;
import ffx.numerics.optimization.LineSearch;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.ForceFieldEnergyOpenMM;
import ffx.potential.ForceFieldEnergyOpenMM.Context;
import ffx.potential.ForceFieldEnergyOpenMM.State;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import java.util.logging.Logger;

/**
 * OpenMM accelerated L-BFGS minimization.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MinimizeOpenMM extends Minimize {

  private static final Logger logger = Logger.getLogger(MinimizeOpenMM.class.getName());

  public MinimizeOpenMM(MolecularAssembly molecularAssembly) {
    super(molecularAssembly, molecularAssembly.getPotentialEnergy(), null);
  }

  public MinimizeOpenMM(
      MolecularAssembly molecularAssembly, ForceFieldEnergyOpenMM forceFieldEnergyOpenMM) {
    super(molecularAssembly, forceFieldEnergyOpenMM, null);
  }

  public MinimizeOpenMM(
      MolecularAssembly molecularAssembly,
      ForceFieldEnergyOpenMM forceFieldEnergyOpenMM,
      AlgorithmListener algorithmListener) {
    super(molecularAssembly, forceFieldEnergyOpenMM, algorithmListener);
  }

  /**
   * Note the OpenMM L-BFGS minimizer does not accept the parameter "m" for the number of previous
   * steps used to estimate the Hessian.
   *
   * @param m The number of previous steps used to estimate the Hessian (ignored).
   * @param eps The convergence criteria.
   * @param maxIterations The maximum number of iterations.
   * @return The potential.
   */
  @Override
  public Potential minimize(int m, double eps, int maxIterations) {
    return minimize(eps, maxIterations);
  }

  /**
   * minimize
   *
   * @param eps The convergence criteria.
   * @param maxIterations The maximum number of iterations.
   * @return a {@link ffx.numerics.Potential} object.
   */
  @Override
  public Potential minimize(double eps, int maxIterations) {

    ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();

    if (forceFieldEnergy instanceof ForceFieldEnergyOpenMM) {
      time = -System.nanoTime();
      ForceFieldEnergyOpenMM openMM = (ForceFieldEnergyOpenMM) forceFieldEnergy;

      // Respect the use flag, and lambda state.
      Atom[] atoms = molecularAssembly.getAtomArray();
      openMM.updateParameters(atoms);

      // Respect (in)active atoms.
      openMM.setActiveAtoms();

      // Get the coordinates to start from.
      openMM.getCoordinates(x);

      // Calculate the starting energy before optimization.
      double e = openMM.energy(x);
      logger.info(format("\n Initial energy:                 %12.6f (kcal/mol)", e));

      // Run the minimization in the current OpenMM Context.
      Context context = openMM.getContext();
      context.optimize(eps, maxIterations);

      // Get the minimized coordinates, forces and potential energy back from OpenMM.
      State state = openMM.createState(true, true, true, false);
      energy = state.potentialEnergy;
      state.getPositions(x);
      state.getGradient(grad);
      state.free();

      // Compute the RMS gradient.
      int index = 0;
      double grad2 = 0;
      for (Atom atom : atoms) {
        if (atom.isActive()) {
          double fx = grad[index++];
          double fy = grad[index++];
          double fz = grad[index++];
          grad2 += fx * fx + fy * fy + fz * fz;
        }
      }
      rmsGradient = sqrt(grad2 / n);

      double[] ffxGrad = new double[n];
      openMM.getCoordinates(x);
      double ffxEnergy = openMM.energyAndGradientFFX(x, ffxGrad);
      double grmsFFX = 0.0;
      for (int i = 0; i < n; i++) {
        double gi = ffxGrad[i];
        if (isNaN(gi) || isInfinite(gi)) {
          String message = format(" The gradient of variable %d is %8.3f.", i, gi);
          logger.warning(message);
        }
        grmsFFX += gi * gi;
      }
      grmsFFX = sqrt(grmsFFX / n);

      time += System.nanoTime();
      logger.info(
          format(
              " Final energy for OpenMM         %12.6f vs. FFX %12.6f in %8.3f (sec).",
              energy, ffxEnergy, time * 1.0e-9));
      logger.info(
          format(
              " Convergence criteria for OpenMM %12.6f vs. FFX %12.6f (kcal/mol/A).",
              rmsGradient, grmsFFX));
    }

    if (algorithmListener != null) {
      algorithmListener.algorithmUpdate(molecularAssembly);
    }

    return forceFieldEnergy;
  }

  /**
   * MinimizeOpenMM does not support the OptimizationListener interface.
   *
   * @since 1.0
   */
  @Override
  public boolean optimizationUpdate(
      int iteration,
      int functionEvaluations,
      double rmsGradient,
      double rmsCoordinateChange,
      double energy,
      double energyChange,
      double angle,
      LineSearch.LineSearchResult lineSearchResult) {
    logger.warning(" MinimizeOpenMM does not support updates at each optimization step.");
    return false;
  }
}
