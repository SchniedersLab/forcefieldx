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
package ffx.potential.nonbonded.pme;

import static ffx.numerics.math.ScalarMath.binomial;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import ffx.potential.parameters.ForceField;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;

@SuppressWarnings("deprecation")
public class SCFPredictorParameters {

  private static final Logger logger = Logger.getLogger(SCFPredictorParameters.class.getName());

  /** Induced dipole predictor order. */
  public int predictorOrder;
  /** Induced dipole predictor index. */
  public int predictorStartIndex;
  /** Induced dipole predictor count. */
  public int predictorCount;
  /** Dimensions of [mode][predictorOrder][nAtoms][3] */
  public double[][][][] predictorInducedDipole;
  /** Dimensions of [mode][predictorOrder][nAtoms][3] */
  public double[][][][] predictorInducedDipoleCR;

  private LeastSquaresPredictor leastSquaresPredictor;
  public LevenbergMarquardtOptimizer leastSquaresOptimizer;

  public final SCFPredictor scfPredictor;
  private final int nAtoms;

  public SCFPredictorParameters(SCFPredictor scfPredictor, int nAtoms) {
    this.nAtoms = nAtoms;
    this.scfPredictor = scfPredictor;
  }

  /** Always-stable predictor-corrector for the mutual induced dipoles. */
  public void aspcPredictor(LambdaMode lambdaMode,
      double[][][] inducedDipole, double[][][] inducedDipoleCR) {

    if (predictorCount < 6) {
      return;
    }

    int mode;
    switch (lambdaMode) {
      case CONDENSED_NO_LIGAND:
        mode = 1;
        break;
      case VAPOR:
        mode = 2;
        break;
      default:
        mode = 0;
    }

    final double[] aspc = {
        22.0 / 7.0, -55.0 / 14.0, 55.0 / 21.0, -22.0 / 21.0, 5.0 / 21.0, -1.0 / 42.0
    };

    // Initialize a pointer into predictor induced dipole array.
    int index = predictorStartIndex;

    // Expansion loop.
    for (int k = 0; k < 6; k++) {

      // Set the current predictor coefficient.
      double c = aspc[k];
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
          inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
        }
      }
      index++;
      if (index >= predictorOrder) {
        index = 0;
      }
    }
  }

  public void init(ForceField forceField) {
    predictorCount = 0;
    int defaultOrder = 6;
    predictorOrder = forceField.getInteger("SCF_PREDICTOR_ORDER", defaultOrder);
    if (scfPredictor == SCFPredictor.LS) {
      leastSquaresPredictor = new LeastSquaresPredictor();
      double eps = 1.0e-4;
      leastSquaresOptimizer =
          new LevenbergMarquardtOptimizer(
              new org.apache.commons.math3.optimization.SimpleVectorValueChecker(eps, eps));
    } else if (scfPredictor == SCFPredictor.ASPC) {
      predictorOrder = 6;
    }
    predictorStartIndex = 0;
  }

  /**
   * The least-squares predictor with induced dipole information from 8-10 previous steps reduces the
   * number SCF iterations by ~50%.
   */
  public void leastSquaresPredictor(LambdaMode lambdaMode,
      double[][][] inducedDipole, double[][][] inducedDipoleCR) {
    if (predictorCount < 2) {
      return;
    }
    try {
      /*
       The Jacobian and target do not change during the LS optimization,
       so it's most efficient to update them once before the
       Least-Squares optimizer starts.
      */
      leastSquaresPredictor.lambdaMode = lambdaMode;
      leastSquaresPredictor.updateJacobianAndTarget();
      int maxEvals = 100;
      fill(leastSquaresPredictor.initialSolution, 0.0);
      leastSquaresPredictor.initialSolution[0] = 1.0;
      org.apache.commons.math3.optimization.PointVectorValuePair optimum =
          leastSquaresOptimizer.optimize(
              maxEvals,
              leastSquaresPredictor,
              leastSquaresPredictor.calculateTarget(),
              leastSquaresPredictor.weights,
              leastSquaresPredictor.initialSolution);
      double[] optimalValues = optimum.getPoint();
      if (logger.isLoggable(Level.FINEST)) {
        logger.finest(format("\n LS RMS:            %10.6f", leastSquaresOptimizer.getRMS()));
        logger.finest(format(" LS Iterations:     %10d", leastSquaresOptimizer.getEvaluations()));
        logger.finest(
            format(" Jacobian Evals:    %10d", leastSquaresOptimizer.getJacobianEvaluations()));
        logger.finest(format(" Chi Square:        %10.6f", leastSquaresOptimizer.getChiSquare()));
        logger.finest(" LS Coefficients");
        for (int i = 0; i < predictorOrder - 1; i++) {
          logger.finest(format(" %2d  %10.6f", i + 1, optimalValues[i]));
        }
      }

      int mode;
      switch (lambdaMode) {
        case CONDENSED_NO_LIGAND:
          mode = 1;
          break;
        case VAPOR:
          mode = 2;
          break;
        default:
          mode = 0;
      }

      // Initialize a pointer into predictor induced dipole array.
      int index = predictorStartIndex;

      // Apply the LS coefficients in order to provide an initial guess at the converged induced
      // dipoles.
      for (int k = 0; k < predictorOrder - 1; k++) {

        // Set the current coefficient.
        double c = optimalValues[k];
        for (int i = 0; i < nAtoms; i++) {
          for (int j = 0; j < 3; j++) {
            inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
            inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
          }
        }
        index++;
        if (index >= predictorOrder) {
          index = 0;
        }
      }
    } catch (Exception e) {
      logger.log(Level.WARNING, " Exception computing predictor coefficients", e);
    }
  }

  /** Polynomial predictor for the mutual induced dipoles. */
  public void polynomialPredictor(LambdaMode lambdaMode,
      double[][][] inducedDipole, double[][][] inducedDipoleCR) {

    if (predictorCount == 0) {
      return;
    }

    int mode;
    switch (lambdaMode) {
      case CONDENSED_NO_LIGAND:
        mode = 1;
        break;
      case VAPOR:
        mode = 2;
        break;
      default:
        mode = 0;
    }

    // Check the number of previous induced dipole vectors available.
    int n = predictorOrder;
    if (predictorCount < predictorOrder) {
      n = predictorCount;
    }

    // Initialize a pointer into predictor induced dipole array.
    int index = predictorStartIndex;

    // Initialize the sign of the polynomial expansion.
    double sign = -1.0;

    // Expansion loop.
    for (int k = 0; k < n; k++) {

      // Set the current predictor sign and coefficient.
      sign *= -1.0;
      double c = sign * binomial(n, k);
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < 3; j++) {
          inducedDipole[0][i][j] += c * predictorInducedDipole[mode][index][i][j];
          inducedDipoleCR[0][i][j] += c * predictorInducedDipoleCR[mode][index][i][j];
        }
      }
      index++;
      if (index >= predictorOrder) {
        index = 0;
      }
    }
  }

  /** Save the current converged mutual induced dipoles. */
  public void saveMutualInducedDipoles(LambdaMode lambdaMode,
      double[][][] inducedDipole, double[][][] inducedDipoleCR,
      double[][] directDipole, double[][] directDipoleCR) {

    int mode;
    switch (lambdaMode) {
      case CONDENSED_NO_LIGAND:
        mode = 1;
        break;
      case VAPOR:
        mode = 2;
        break;
      default:
        mode = 0;
    }

    // Current induced dipoles are saved before those from the previous step.
    predictorStartIndex--;
    if (predictorStartIndex < 0) {
      predictorStartIndex = predictorOrder - 1;
    }

    if (predictorCount < predictorOrder) {
      predictorCount++;
    }

    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        predictorInducedDipole[mode][predictorStartIndex][i][j] =
            inducedDipole[0][i][j] - directDipole[i][j];
        predictorInducedDipoleCR[mode][predictorStartIndex][i][j] =
            inducedDipoleCR[0][i][j] - directDipoleCR[i][j];
      }
    }
  }

  private class LeastSquaresPredictor implements DifferentiableMultivariateVectorFunction {

    double[] weights;
    double[] target;
    double[] values;
    double[][] jacobian;
    double[] initialSolution;
    private final MultivariateMatrixFunction multivariateMatrixFunction = this::jacobian;

    public LambdaMode lambdaMode = null;

    LeastSquaresPredictor() {
      weights = new double[2 * nAtoms * 3];
      target = new double[2 * nAtoms * 3];
      values = new double[2 * nAtoms * 3];
      jacobian = new double[2 * nAtoms * 3][predictorOrder - 1];
      initialSolution = new double[predictorOrder - 1];
      fill(weights, 1.0);
      initialSolution[0] = 1.0;
    }

    @Override
    public MultivariateMatrixFunction jacobian() {
      return multivariateMatrixFunction;
    }

    @Override
    public double[] value(double[] variables) {
      int mode;
      switch (lambdaMode) {
        case CONDENSED_NO_LIGAND:
          mode = 1;
          break;
        case VAPOR:
          mode = 2;
          break;
        default:
          mode = 0;
      }

      for (int i = 0; i < nAtoms; i++) {
        int index = 6 * i;
        values[index] = 0;
        values[index + 1] = 0;
        values[index + 2] = 0;
        values[index + 3] = 0;
        values[index + 4] = 0;
        values[index + 5] = 0;
        int pi = predictorStartIndex + 1;
        if (pi >= predictorOrder) {
          pi = 0;
        }
        for (int j = 0; j < predictorOrder - 1; j++) {
          values[index] += variables[j] * predictorInducedDipole[mode][pi][i][0];
          values[index + 1] += variables[j] * predictorInducedDipole[mode][pi][i][1];
          values[index + 2] += variables[j] * predictorInducedDipole[mode][pi][i][2];
          values[index + 3] += variables[j] * predictorInducedDipoleCR[mode][pi][i][0];
          values[index + 4] += variables[j] * predictorInducedDipoleCR[mode][pi][i][1];
          values[index + 5] += variables[j] * predictorInducedDipoleCR[mode][pi][i][2];
          pi++;
          if (pi >= predictorOrder) {
            pi = 0;
          }
        }
      }
      return values;
    }

    double[] calculateTarget() {
      return target;
    }

    void updateJacobianAndTarget() {
      int mode;
      switch (lambdaMode) {
        case CONDENSED_NO_LIGAND:
          mode = 1;
          break;
        case VAPOR:
          mode = 2;
          break;
        default:
          mode = 0;
      }

      // Update the target.
      int index = 0;
      for (int i = 0; i < nAtoms; i++) {
        target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][0];
        target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][1];
        target[index++] = predictorInducedDipole[mode][predictorStartIndex][i][2];
        target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][0];
        target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][1];
        target[index++] = predictorInducedDipoleCR[mode][predictorStartIndex][i][2];
      }

      // Update the Jacobian.
      index = predictorStartIndex + 1;
      if (index >= predictorOrder) {
        index = 0;
      }
      for (int j = 0; j < predictorOrder - 1; j++) {
        int ji = 0;
        for (int i = 0; i < nAtoms; i++) {
          jacobian[ji++][j] = predictorInducedDipole[mode][index][i][0];
          jacobian[ji++][j] = predictorInducedDipole[mode][index][i][1];
          jacobian[ji++][j] = predictorInducedDipole[mode][index][i][2];
          jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][0];
          jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][1];
          jacobian[ji++][j] = predictorInducedDipoleCR[mode][index][i][2];
        }
        index++;
        if (index >= predictorOrder) {
          index = 0;
        }
      }
    }

    private double[][] jacobian(double[] variables) {
      return jacobian;
    }
  }
}
