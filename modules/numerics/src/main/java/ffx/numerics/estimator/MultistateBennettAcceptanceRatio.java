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
package ffx.numerics.estimator;

import ffx.numerics.OptimizationInterface;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch;
import ffx.numerics.optimization.OptimizationListener;
import ffx.utilities.Constants;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.MathArrays;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;

import static ffx.numerics.estimator.EstimateBootstrapper.getBootstrapIndices;
import static ffx.numerics.estimator.Zwanzig.Directionality.BACKWARDS;
import static ffx.numerics.estimator.Zwanzig.Directionality.FORWARDS;
import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.stream;
import static org.apache.commons.lang3.ArrayFill.fill;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The MultistateBennettAcceptanceRatio class defines a statistical estimator based on a generalization
 * to the Bennett Acceptance Ratio (BAR) method for multiple lambda windows. It requires an input of
 * K X N array of energies (every window at every snap at every lambda value). No support for different
 * number of snapshots at each window. This will be caught by the filter, but not by the Harmonic Oscillators
 * testcase.
 * <p>
 * This class implements the method discussed in:
 * Shirts, M. R. and Chodera, J. D. (2008) Statistically optimal analysis of snaps from multiple equilibrium
 * states. J. Chem. Phys. 129, 124105. doi:10.1063/1.2978177
 * <p>
 * This class is based heavily on the pymbar code, which is available at:
 * https://github.com/choderalab/pymbar/tree/master
 *
 * @author Matthew J. Speranza
 * @since 1.0
 */
public class MultistateBennettAcceptanceRatio extends SequentialEstimator implements BootstrappableEstimator, OptimizationInterface {
  private static final Logger logger = Logger.getLogger(MultistateBennettAcceptanceRatio.class.getName());

  /**
   * Default MBAR convergence tolerance.
   */
  private static final double DEFAULT_TOLERANCE = 1.0E-7;
  /**
   * Number of free of differences between simulation windows. Calculated at the very end.
   */
  private final int nFreeEnergyDiffs;
  /**
   * MBAR free-energy difference estimates (nFreeEnergyDiffs values).
   */
  private final double[] mbarFEDifferenceEstimates;
  /**
   * Number of lamda states (basically nFreeEnergyDiffs + 1).
   */
  private final int nLambdaStates;
  /**
   * MBAR free-energy estimates at each lambda value (nStates values). The first value is defined as 0 throughout to
   * promote stability. This estimate is novel to the MBAR method and is not seen in BAR or Zwanzig. Only the differences
   * between these values have physical significance.
   */
  double[] mbarFEEstimates;
  /**
   * MBAR free-energy difference uncertainties.
   */
  private double[] mbarUncertainties;
  /**
   * Matrix of free-energy difference uncertainties between all i & j
   */
  private double[][] diffMatrix;
  /**
   * MBAR convergence tolerance.
   */
  private final double tolerance;
  /**
   * Random number generator used for bootstrapping.
   */
  private final Random random;
  /**
   * Total MBAR free-energy difference estimate.
   */
  private double totalMBAREstimate;
  /**
   * Total MBAR free-energy difference uncertainty.
   */
  private double totalMBARUncertainty;
  /**
   * MBAR Enthalpy estimates
   */
  private final double[] mbarEnthalpy;

  /**
   * "Reduced" potential energies. -ln(exp(beta * -U)) or more practically U * (1 / RT).
   * Has shape (nLambdaStates, numSnaps * nLambdaStates)
   */
  private double[][] reducedPotentials;
  /**
   * Seed MBAR calculation with another free energy estimation (BAR,ZWANZIG) or zeros
   */
  private SeedType seedType;
  /**
   * Enum of MBAR seed types.
   */
  public enum SeedType {BAR, ZWANZIG, ZEROS}

  /**
   * Constructor for MBAR estimator.
   *
   * @param lambdaValues array of lambda values
   * @param energiesAll  array of energies at each lambda value
   * @param temperature  array of temperatures
   */
  public MultistateBennettAcceptanceRatio(double[] lambdaValues, double[][][] energiesAll, double[] temperature) {
    this(lambdaValues, energiesAll, temperature, DEFAULT_TOLERANCE, SeedType.ZWANZIG);
  }

  /**
   * Constructor for MBAR estimator.
   *
   * @param lambdaValues array of lambda values
   * @param energiesAll  array of energies at each lambda value
   * @param temperature  array of temperatures
   * @param tolerance    convergence tolerance
   * @param seedType     seed type for MBAR
   */
  public MultistateBennettAcceptanceRatio(double[] lambdaValues, double[][][] energiesAll, double[] temperature,
                                          double tolerance, SeedType seedType) {
    super(lambdaValues, energiesAll, temperature);
    this.tolerance = tolerance;
    this.seedType = seedType;

    // MBAR calculates free energy at each lambda value (only the differences between them have physical significance)
    nLambdaStates = lambdaValues.length;
    mbarFEEstimates = new double[nLambdaStates];

    nFreeEnergyDiffs = lambdaValues.length - 1;
    mbarFEDifferenceEstimates = new double[nFreeEnergyDiffs];
    mbarUncertainties = new double[nFreeEnergyDiffs];
    mbarEnthalpy = new double[nFreeEnergyDiffs];
    random = new Random();
    estimateDG();
  }

  /**
   * Set the MBAR seed energies using BAR, Zwanzig or zeros.
   */
  private void seedEnergies() {
    switch (seedType) {
      case BAR:
        try {
          SequentialEstimator barEstimator = new BennettAcceptanceRatio(lamValues, eLow, eAt, eHigh, temperatures);
          mbarFEEstimates[0] = 0.0;
          double[] barEstimates = barEstimator.getBinEnergies();
          for (int i = 0; i < nFreeEnergyDiffs; i++) {
            mbarFEEstimates[i + 1] = mbarFEEstimates[i] + barEstimates[i];
          }
          break;
        } catch (IllegalArgumentException e) {
          logger.warning(" BAR failed to converge. Zwanzig will be used for seed energies.");
          seedType = SeedType.ZWANZIG;
          seedEnergies();
          return;
        }
      case ZWANZIG:
        try{
          Zwanzig forwardsFEP = new Zwanzig(lamValues, eLow, eAt, eHigh, temperatures, FORWARDS);
          Zwanzig backwardsFEP = new Zwanzig(lamValues, eLow, eAt, eHigh, temperatures, BACKWARDS);
          double[] forwardZwanzig = forwardsFEP.getBinEnergies();
          double[] backwardZwanzig = backwardsFEP.getBinEnergies();
          mbarFEEstimates[0] = 0.0;
          for (int i = 0; i < nFreeEnergyDiffs; i++) {
            mbarFEEstimates[i + 1] = mbarFEEstimates[i] + .5 * (forwardZwanzig[i] + backwardZwanzig[i]);
          }
          if (stream(mbarFEEstimates).anyMatch(Double::isInfinite) || stream(mbarFEEstimates).anyMatch(Double::isNaN)) {
            throw new IllegalArgumentException("MBAR contains NaNs or Infs after seeding.");
          }
        break;
        } catch (IllegalArgumentException e) {
          logger.warning(" Zwanzig failed to converge. Zeros will be used for seed energies.");
          seedType = SeedType.ZEROS;
          seedEnergies();
          return;
        }
      case SeedType.ZEROS:
        fill(mbarFEEstimates, 0.0);
        break;
      default:
        throw new IllegalArgumentException("Seed type not supported");
    }
  }

  /**
   * Get the MBAR free-energy estimates at each lambda value.
   */
  @Override
  public void estimateDG() {
    estimateDG(false);
  }

  /**
   * Implementation of MBAR solved with self-consistent iteration and L-BFGS optimization.
   */
  @Override
  public void estimateDG(boolean randomSamples) {
    // Bootstrap needs resetting to zeros
    fill(mbarFEEstimates, 0.0);
    seedEnergies();
    if (stream(mbarFEEstimates).anyMatch(Double::isInfinite) || stream(mbarFEEstimates).anyMatch(Double::isNaN)) {
      seedType = SeedType.ZEROS;
      seedEnergies();
    }

    // Precompute beta for each state.
    double[] rtValues = new double[nLambdaStates];
    double[] invRTValues = new double[nLambdaStates];
    for (int i = 0; i < nLambdaStates; i++) {
      rtValues[i] = Constants.R * temperatures[i];
      invRTValues[i] = 1.0 / rtValues[i];
    }
    int numEvaluations = eAllFlat[0].length;

    // Sample random snapshots from each window.
    int[][] indices = new int[nLambdaStates][numEvaluations];
    if (randomSamples) {
      int[] randomIndices = getBootstrapIndices(numEvaluations, random);
      for (int i = 0; i < nLambdaStates; i++) {
        // Use the same random indices across all lambda values
        indices[i] = randomIndices;
      }
    } else {
      for (int i = 0; i < numEvaluations; i++) {
        for (int j = 0; j < nLambdaStates; j++) {
          indices[j][i] = i;
        }
      }
    }

    // Precompute reducedPotentials since it doesn't change
    reducedPotentials = new double[nLambdaStates][numEvaluations];
    double minPotential = Double.POSITIVE_INFINITY;
    for (int state = 0; state < nLambdaStates; state++) { // For each lambda value
      for (int n = 0; n < numEvaluations; n++) {
        reducedPotentials[state][n] = eAllFlat[state][indices[state][n]] * invRTValues[state];
        if (reducedPotentials[state][n] < minPotential) {
          minPotential = reducedPotentials[state][n];
        }
      }
    }

    // Subtract the minimum potential from all potentials (we are calculating relative free energies anyway)
    for (int state = 0; state < nLambdaStates; state++) {
      for (int n = 0; n < numEvaluations; n++) {
        reducedPotentials[state][n] -= minPotential;
      }
    }

    // Remove reduced potential arrays where snaps are zero since they cause issues for N.R. and SCI
    // i.e. no trajectories for that lambda were generated/sampled, but other trajectories had potentials evaluated at that lambda
    ArrayList<Integer> zeroSnapLambdas = new ArrayList<>();
    ArrayList<Integer> sampledLambdas = new ArrayList<>();
    for(int i = 0; i < nLambdaStates; i++) {
      if (snaps[i] == 0) {
        zeroSnapLambdas.add(i);
      } else {
        sampledLambdas.add(i);
      }
    }
    int nLambdaStatesTemp = nLambdaStates - zeroSnapLambdas.size();
    double[][] reducedPotentialsTemp = new double[nLambdaStates - zeroSnapLambdas.size()][numEvaluations];
    double[] mbarFEEstimatesTemp = new double[nLambdaStates - zeroSnapLambdas.size()];
    int[] snapsTemp = new int[nLambdaStates - zeroSnapLambdas.size()];
    if (!zeroSnapLambdas.isEmpty()) {
      int index = 0;
      for (int i = 0; i < nLambdaStates; i++) {
        if (!zeroSnapLambdas.contains(i)) {
          reducedPotentialsTemp[index] = reducedPotentials[i];
          mbarFEEstimatesTemp[index] = mbarFEEstimates[i];
          snapsTemp[index] = snaps[i];
          index++;
        }
      }
    } else { // If there aren't any zero snap lambdas, just use the original arrays
      reducedPotentialsTemp = reducedPotentials;
      mbarFEEstimatesTemp = mbarFEEstimates;
      snapsTemp = snaps;
    }

    // SCI iterations used to start optimization of MBAR objective function.
    // Optimizers can struggle when starting too far from the minimum, but SCI doesn't.
    double[] prevMBAR;
    double omega = 1.5; // Parameter chosen empirically to work with most systems (> 2 works but not always).
    for (int i = 0; i < 10; i++) {
      prevMBAR = copyOf(mbarFEEstimatesTemp, nLambdaStatesTemp);
      mbarFEEstimatesTemp = selfConsistentUpdate(reducedPotentialsTemp, snapsTemp, mbarFEEstimatesTemp);
      for (int j = 0; j < nLambdaStatesTemp; j++) { // SOR
        mbarFEEstimatesTemp[j] = omega * mbarFEEstimatesTemp[j] + (1 - omega) * prevMBAR[j];
      }
      if (stream(mbarFEEstimatesTemp).anyMatch(Double::isInfinite) || stream(mbarFEEstimatesTemp).anyMatch(Double::isNaN)) {
        throw new IllegalArgumentException("MBAR contains NaNs or Infs during startup SCI ");
      }
    }

    try {
      if (nLambdaStatesTemp > 100) { // L-BFGS optimization for high granularity windows where hessian is expensive
        int mCorrections = 5;
        double[] x = new double[nLambdaStatesTemp];
        arraycopy(mbarFEEstimatesTemp, 0, x, 0, nLambdaStatesTemp);
        double[] grad = mbarGradient(reducedPotentialsTemp, snapsTemp, mbarFEEstimatesTemp);
        double eps = 1.0E-4; // Gradient tolarance -> chosen since L-BFGS seems unstable with tight tolerances
        OptimizationListener listener = getOptimizationListener();
        LBFGS.minimize(nLambdaStatesTemp, mCorrections, x, mbarObjectiveFunction(reducedPotentialsTemp, snapsTemp, mbarFEEstimatesTemp),
            grad, eps, 1000, this, listener);
        arraycopy(x, 0, mbarFEEstimatesTemp, 0, nLambdaStatesTemp);
      } else { // Newton optimization if hessian inversion isn't too expensive
        mbarFEEstimatesTemp = newton(mbarFEEstimatesTemp, reducedPotentialsTemp, snapsTemp, 1.0, 100, tolerance);
      }

      // Update the FE estimates with the optimized values from derivative-based optimization
      int count = 0;
      for(Integer i : sampledLambdas){
        mbarFEEstimates[i] = mbarFEEstimatesTemp[count];
        count++;
      }
    } catch (Exception e) {
      logger.warning(" L-BFGS/Newton failed to converge. Finishing w/ self-consistent iteration.");
      logger.warning(e.getMessage());
    }

    // Self-consistent iteration is used to finish off optimization of MBAR objective function
    int sciIter = 0;
    do {
      prevMBAR = copyOf(mbarFEEstimates, nLambdaStates);
      mbarFEEstimates = selfConsistentUpdate(reducedPotentials, snaps, mbarFEEstimates);
      for (int i = 0; i < nLambdaStates; i++) { // SOR
        mbarFEEstimates[i] = omega * mbarFEEstimates[i] + (1 - omega) * prevMBAR[i];
      }
      if (stream(mbarFEEstimates).anyMatch(Double::isInfinite) || stream(mbarFEEstimates).anyMatch(Double::isNaN)) {
        throw new IllegalArgumentException("MBAR estimate contains NaNs or Infs after iteration " + sciIter);
      }
      sciIter++;
    } while (!converged(prevMBAR));

    logger.fine(" MBAR converged after " + sciIter + " iterations with omega " + omega + ".");

    // Zero out the first term
    double f0 = mbarFEEstimates[0];
    for (int i = 0; i < nLambdaStates; i++) {
      mbarFEEstimates[i] -= f0;
    }

    // Calculate uncertainties
    mbarUncertainties = mbarUncertaintyCalc(reducedPotentials, snaps, mbarFEEstimates);
    totalMBARUncertainty = mbarTotalUncertaintyCalc(reducedPotentials, snaps, mbarFEEstimates);
    diffMatrix = diffMatrixCalculation(reducedPotentials, snaps, mbarFEEstimates);

    // Convert to kcal/mol & calculate differences/sums
    for (int i = 0; i < nLambdaStates; i++) {
      mbarFEEstimates[i] = mbarFEEstimates[i] * rtValues[i];
    }
    for (int i = 0; i < nFreeEnergyDiffs; i++) {
      mbarFEDifferenceEstimates[i] = mbarFEEstimates[i + 1] - mbarFEEstimates[i];
    }

    totalMBAREstimate = stream(mbarFEDifferenceEstimates).sum();
  }

  /**
   * Checks if the MBAR free energy estimates have converged by comparing the difference
   * between the previous and current free energies. The tolerance is set by the user.
   *
   * @param prevMBAR previous MBAR free energy estimates.
   * @return true if converged, false otherwise
   */
  private boolean converged(double[] prevMBAR) {
    double[] differences = new double[prevMBAR.length];
    for (int i = 0; i < prevMBAR.length; i++) {
      differences[i] = abs(prevMBAR[i] - mbarFEEstimates[i]);
    }
    return stream(differences).allMatch(d -> d < tolerance);
  }

  //////// Methods for calculating MBAR variables, vectors, and matrices. ////////

  /**
   * MBAR objective function. This is used for L-BFGS optimization.
   *
   * @param reducedPotentials -ln(boltzmann weights)
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return The objective function value.
   */
  private static double mbarObjectiveFunction(double[][] reducedPotentials, int[] snapsPerLambda, double[] freeEnergyEstimates) {
    if (stream(freeEnergyEstimates).anyMatch(Double::isInfinite) || stream(freeEnergyEstimates).anyMatch(Double::isNaN)) {
      throw new IllegalArgumentException("MBAR contains NaNs or Infs.");
    }
    int nStates = freeEnergyEstimates.length;
    double[] log_denom_n = new double[reducedPotentials[0].length];
    for (int i = 0; i < reducedPotentials[0].length; i++) {
      double[] temp = new double[nStates];
      double maxTemp = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < nStates; j++) {
        temp[j] = freeEnergyEstimates[j] - reducedPotentials[j][i];
        if (temp[j] > maxTemp) {
          maxTemp = temp[j];
        }
      }
      log_denom_n[i] = logSumExp(temp, snapsPerLambda, maxTemp);
    }
    double[] dotNkFk = new double[snapsPerLambda.length];
    for (int i = 0; i < snapsPerLambda.length; i++) {
      dotNkFk[i] = snapsPerLambda[i] * freeEnergyEstimates[i];
    }
    return stream(log_denom_n).sum() - stream(dotNkFk).sum();
  }

  /**
   * Gradient of the MBAR objective function. This is used for L-BFGS optimization.
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return Gradient for the mbar objective function.
   */
  private static double[] mbarGradient(double[][] reducedPotentials, int[] snapsPerLambda, double[] freeEnergyEstimates) {
    int nStates = freeEnergyEstimates.length;
    double[] log_num_k = new double[nStates];
    double[] log_denom_n = new double[reducedPotentials[0].length];
    double[][] logDiff = new double[reducedPotentials.length][reducedPotentials[0].length];
    double maxLogDiff = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < reducedPotentials[0].length; i++) {
      double[] temp = new double[nStates];
      double maxTemp = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < nStates; j++) {
        temp[j] = freeEnergyEstimates[j] - reducedPotentials[j][i];
        if (temp[j] > maxTemp) {
          maxTemp = temp[j];
        }
      }
      log_denom_n[i] = logSumExp(temp, snapsPerLambda, maxTemp);
      for (int j = 0; j < nStates; j++) {
        logDiff[j][i] = -log_denom_n[i] - reducedPotentials[j][i];
        if (logDiff[j][i] > maxLogDiff) {
          maxLogDiff = logDiff[j][i];
        }
      }
    }
    for (int i = 0; i < nStates; i++) {
      log_num_k[i] = logSumExp(logDiff[i], maxLogDiff);
    }
    double[] grad = new double[nStates];
    for (int i = 0; i < nStates; i++) {
      grad[i] = -1.0 * snapsPerLambda[i] * (1.0 - exp(freeEnergyEstimates[i] + log_num_k[i]));
    }
    return grad;
  }

  /**
   * Hessian of the MBAR objective function. This is used for Newton optimization.
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return Hessian for the mbar objective function.
   */
  private static double[][] mbarHessian(double[][] reducedPotentials, int[] snapsPerLambda, double[] freeEnergyEstimates) {
    int nStates = freeEnergyEstimates.length;
    double[][] W = mbarW(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    // h = dot(W.T, W) * snapsPerLambda * snapsPerLambda[:, newaxis] - diag(W.sum(0) * snapsPerLambda)
    double[][] hessian = new double[nStates][nStates];
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < nStates; j++) {
        double sum = 0.0;
        for (int k = 0; k < reducedPotentials[0].length; k++) {
          sum += W[i][k] * W[j][k];
        }
        hessian[i][j] = sum * snapsPerLambda[i] * snapsPerLambda[j];
      }
      double wSum = 0.0;
      for (int k = 0; k < W[i].length; k++) {
        wSum += W[i][k];
      }
      hessian[i][i] -= wSum * snapsPerLambda[i];
    }
    // h = -h
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < nStates; j++) {
        hessian[i][j] = -hessian[i][j];
      }
    }
    return hessian;
  }

  /**
   * W = exp(freeEnergyEstimates - reducedPotentials.T - log_denominator_n[:, newaxis])
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return W matrix.
   */
  private static double[][] mbarW(double[][] reducedPotentials, int[] snapsPerLambda, double[] freeEnergyEstimates) {
    int nStates = freeEnergyEstimates.length;
    double[] log_denom_n = new double[reducedPotentials[0].length];
    double[][] logDiff = new double[reducedPotentials.length][reducedPotentials[0].length];
    double maxLogDiff = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < reducedPotentials[0].length; i++) {
      double[] temp = new double[nStates];
      double maxTemp = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < nStates; j++) {
        temp[j] = freeEnergyEstimates[j] - reducedPotentials[j][i];
        if (temp[j] > maxTemp) {
          maxTemp = temp[j];
        }
      }
      log_denom_n[i] = logSumExp(temp, snapsPerLambda, maxTemp);
      for (int j = 0; j < nStates; j++) {
        logDiff[j][i] = -log_denom_n[i] - reducedPotentials[j][i];
        if (logDiff[j][i] > maxLogDiff) {
          maxLogDiff = logDiff[j][i];
        }
      }
    }
    // logW = freeEnergyEstimates - reducedPotentials.T - log_denominator_n[:, newaxis]
    double[][] W = new double[nStates][reducedPotentials[0].length];
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < reducedPotentials[0].length; j++) {
        W[i][j] = exp(freeEnergyEstimates[i] - reducedPotentials[i][j] - log_denom_n[j]);
      }
    }
    return W;
  }

  /**
   * logW = freeEnergyEstimates - reducedPotentials.T - log_denominator_n[:, newaxis]
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return logW matrix.
   */
  private static double[][] mbarLogW(double[][] reducedPotentials, int[] snapsPerLambda, double[] freeEnergyEstimates) {
    int nStates = freeEnergyEstimates.length;
    // double[] log_num_k = new double[nStates];
    double[] log_denom_n = new double[reducedPotentials[0].length];
    double[][] logDiff = new double[reducedPotentials.length][reducedPotentials[0].length];
    double maxLogDiff = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < reducedPotentials[0].length; i++) {
      double[] temp = new double[nStates];
      double maxTemp = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < nStates; j++) {
        temp[j] = freeEnergyEstimates[j] - reducedPotentials[j][i];
        if (temp[j] > maxTemp) {
          maxTemp = temp[j];
        }
      }
      log_denom_n[i] = logSumExp(temp, snapsPerLambda, maxTemp);
      for (int j = 0; j < nStates; j++) {
        logDiff[j][i] = -log_denom_n[i] - reducedPotentials[j][i];
        if (logDiff[j][i] > maxLogDiff) {
          maxLogDiff = logDiff[j][i];
        }
      }
    }
    // logW = freeEnergyEstimates - reducedPotentials.T - log_denominator_n[:, newaxis]
    double[][] logW = new double[nStates][reducedPotentials[0].length];
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < reducedPotentials[0].length; j++) {
        logW[i][j] = freeEnergyEstimates[i] - reducedPotentials[i][j] - log_denom_n[j];
      }
    }
    return logW;
  }

  /**
   * Theta = W.T @ (I - W @ diag(snapsPerState) @ W.T)^-1 @ W.
   * <p>
   * Requires calculation and inversion of W matrix.
   * D4 from supp info of MBAR paper used instead.
   *
   * @param reducedPotentials energies
   * @param snapsPerState  number of snaps per state
   * @param snapsPerLambda  free energies
   * @return Theta matrix.
   */
  private static double[][] mbarTheta(double[][] reducedPotentials, int[] snapsPerState, double[] snapsPerLambda) {
    // SVD of W
    double[][] W = mbarW(reducedPotentials, snapsPerState, snapsPerLambda);
    RealMatrix WMatrix = MatrixUtils.createRealMatrix(W).transpose();
    RealMatrix I = MatrixUtils.createRealIdentityMatrix(snapsPerLambda.length);
    RealMatrix NkMatrix = MatrixUtils.createRealDiagonalMatrix(stream(snapsPerState).mapToDouble(i -> i).toArray());
    SingularValueDecomposition svd = new SingularValueDecomposition(WMatrix);
    RealMatrix V = svd.getV();
    RealMatrix S = MatrixUtils.createRealDiagonalMatrix(svd.getSingularValues());

    // W.T @ (I - W @ diag(snapsPerState) @ W.T)^-1 @ W
    // = V @ S @ (I - S @ V.T @ diag(snapsPerState) @ V @ S)^-1 @ S @ V.T
    RealMatrix theta = S.multiply(V.transpose());
    theta = theta.multiply(NkMatrix).multiply(V).multiply(S);
    theta = I.subtract(theta);
    theta = new SingularValueDecomposition(theta).getSolver().getInverse(); // pinv equivalent
    theta = V.multiply(S).multiply(theta).multiply(S).multiply(V.transpose());

    return theta.getData();
  }

  /**
   * MBAR uncertainty calculation.
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return Uncertainties for the MBAR free energy estimates.
   */
  private static double[] mbarUncertaintyCalc(double[][] reducedPotentials, int[] snapsPerLambda,
                                              double[] freeEnergyEstimates) {
    double[][] theta = mbarTheta(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[] uncertainties = new double[freeEnergyEstimates.length - 1];
    // del(dFij) = Theta[i,i] - 2 * Theta[i,j] + Theta[j,j]
    for (int i = 0; i < freeEnergyEstimates.length - 1; i++) {
      // TODO: Figure out why this is happening (likely due to theta calculation differing from pymbar's)
      double variance = theta[i][i] - 2 * theta[i][i + 1] + theta[i + 1][i + 1];
      if (variance < 0) {
        logger.warning(" Negative variance detected in MBAR uncertainty calculation. " +
                "Multiplying by -1 to get real value. Check diff matrix to see which variances were negative. " +
                "They should be NaN.");
        variance *= -1;
      }
      uncertainties[i] = sqrt(variance);
    }
    return uncertainties;
  }

  /**
   * MBAR total uncertainty calculation.
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return Total uncertainty for the MBAR free energy estimates.
   */
  private static double mbarTotalUncertaintyCalc(double[][] reducedPotentials, int[] snapsPerLambda,
                                                 double[] freeEnergyEstimates) {
    double[][] theta = mbarTheta(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    int nStates = freeEnergyEstimates.length;
    return sqrt(theta[0][0] - 2 * theta[0][nStates - 1] + theta[nStates - 1][nStates - 1]);
  }

  /**
   * MBAR diff Matrix calculation.
   *
   * @param reducedPotential energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return Diff matrix for the MBAR free energy estimates.
   */
  private static double[][] diffMatrixCalculation(double[][] reducedPotential, int[] snapsPerLambda,
                                                  double[] freeEnergyEstimates) {
    double[][] theta = mbarTheta(reducedPotential, snapsPerLambda, freeEnergyEstimates);
    double[][] diffMatrix = new double[freeEnergyEstimates.length][freeEnergyEstimates.length];
    for (int i = 0; i < freeEnergyEstimates.length; i++) {
      for (int j = 0; j < freeEnergyEstimates.length; j++) {
        diffMatrix[i][j] = sqrt(theta[i][i] - 2 * theta[i][j] + theta[j][j]);
      }
    }
    return diffMatrix;
  }

  //////// Methods for solving MBAR with self-consistent iteration, L-BFGS optimization, and Newton-Raphson. ////////

  /**
   * Self-consistent iteration to update free energies.
   *
   * @param reducedPotential energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return updated free energies
   */
  private static double[] selfConsistentUpdate(double[][] reducedPotential, int[] snapsPerLambda,
                                               double[] freeEnergyEstimates) {
    int nStates = freeEnergyEstimates.length;
    double[] updatedF_k = new double[nStates];
    double[] log_denom_n = new double[reducedPotential[0].length];
    double[][] logDiff = new double[reducedPotential.length][reducedPotential[0].length];
    double[] maxLogDiff = new double[nStates];
    fill(maxLogDiff, Double.NEGATIVE_INFINITY);
    for (int i = 0; i < reducedPotential[0].length; i++) {
      double[] temp = new double[nStates];
      double maxTemp = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < nStates; j++) {
        temp[j] = freeEnergyEstimates[j] - reducedPotential[j][i];
        if (temp[j] > maxTemp) {
          maxTemp = temp[j];
        }
      }
      log_denom_n[i] = logSumExp(temp, snapsPerLambda, maxTemp);
      for (int j = 0; j < nStates; j++) {
        logDiff[j][i] = -log_denom_n[i] - reducedPotential[j][i];
        if (logDiff[j][i] > maxLogDiff[j]) {
          maxLogDiff[j] = logDiff[j][i];
        }
      }
    }

    for (int i = 0; i < nStates; i++) {
      updatedF_k[i] = -1.0 * logSumExp(logDiff[i], maxLogDiff[i]);
    }

    // Constrain f1=0 over the course of iterations to prevent uncontrolled growth in magnitude
    double norm = updatedF_k[0];
    updatedF_k[0] = 0.0;
    for (int i = 1; i < nStates; i++) {
      updatedF_k[i] = updatedF_k[i] - norm;
    }

    return updatedF_k;
  }

  /**
   * Newton-Raphson step for MBAR optimization.
   *
   * @param n        current free energies.
   * @param grad     gradient of the objective function.
   * @param hessian  hessian of the objective function.
   * @param stepSize step size for the Newton-Raphson step.
   * @return updated free energies.
   */
  private static double[] newtonStep(double[] n, double[] grad, double[][] hessian, double stepSize) {
    double[] nPlusOne = new double[n.length];
    RealMatrix hessianInverse = MatrixUtils.inverse(MatrixUtils.createRealMatrix(hessian));
    double[] step = hessianInverse.preMultiply(grad);
    // Zero out the first term of the step
    double temp = step[0];
    step[0] = 0.0;
    for (int i = 1; i < step.length; i++) {
      step[i] -= temp;
    }
    for (int i = 0; i < n.length; i++) {
      nPlusOne[i] = n[i] - step[i] * stepSize;
    }
    return nPlusOne;
  }

  /**
   * Newton-Raphson optimization for MBAR.
   *
   * @param freeEnergyEstimates       free energies.
   * @param reducedPotentials      energies.
   * @param snapsPerLambda       number of snaps per state.
   * @param stepSize  step size for the Newton-Raphson step.
   * @param maxIter   maximum number of iterations.
   * @param tolerance convergence tolerance.
   * @return updated free energies.
   */
  private static double[] newton(double[] freeEnergyEstimates, double[][] reducedPotentials,
                                 int[] snapsPerLambda, double stepSize, int maxIter, double tolerance) {
    double[] grad = mbarGradient(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[][] hessian = mbarHessian(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[] f_kPlusOne = newtonStep(freeEnergyEstimates, grad, hessian, stepSize);
    int iter = 1;
    while (iter < maxIter && MathArrays.distance1(freeEnergyEstimates, f_kPlusOne) > tolerance) {
      freeEnergyEstimates = f_kPlusOne;
      grad = mbarGradient(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
      hessian = mbarHessian(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
      f_kPlusOne = newtonStep(freeEnergyEstimates, grad, hessian, stepSize);
      iter++;
    }

    logger.fine(" Newton converged after " + iter + " iterations.");

    return f_kPlusOne;
  }

  /**
   * Calculates the log of the sum of the exponentials of the given values.
   * <p>
   * The max value is subtracted from each value in the array before exponentiation to prevent overflow.
   *
   * @param values The values to exponentiate and sum.
   * @param max    The max value is subtracted from each value in the array prior to exponentiation.
   * @return the sum
   */
  private static double logSumExp(double[] values, double max) {
    int[] b = fill(new int[values.length], 1);
    return logSumExp(values, b, max);
  }

  /**
   * Calculates the log of the sum of the exponentials of the given values.
   * <p>
   * The max value is subtracted from each value in the array before exponentiation to prevent overflow.
   *
   * @param values The values to exponentiate and sum.
   * @param max    The max value is subtracted from each value in the array prior to exponentiation.
   * @param b      Weights for each value in the array.
   * @return the sum
   */
  private static double logSumExp(double[] values, int[] b, double max) {
    // ChatGPT mostly wrote this and I tweaked it to match more closely with scipy's logsumexp implementation
    // Find the maximum value in the array.
    assert values.length == b.length : "values and b must be the same length";

    // Subtract the maximum value from each value in the array, exponentiate the result, and add up these values.
    double sum = 0.0;
    for (int i = 0; i < values.length; i++) {
      sum += b[i] * exp(values[i] - max);
    }

    // Take the natural logarithm of the sum and add the maximum value back in.
    return max + log(sum);
  }

  /**
   * TODO: Log out the MBAR optimization progress.
   *
   * @return
   */
  private OptimizationListener getOptimizationListener() {
    return new OptimizationListener() {
      @Override
      public boolean optimizationUpdate(int iter, int nBFGS, int nFunctionEvals, double gradientRMS,
                                        double coordinateRMS, double f, double df, double angle,
                                        LineSearch.LineSearchResult info) {
        return true;
      }
    };
  }

  /**
   * MBAR objective function evaluation at a given free energy estimate for L-BFGS optimization.
   *
   * @param x Input parameters.
   * @return The objective function value at the given parameters.
   */
  @Override
  public double energy(double[] x) {
    // Zero out the first term
    double tempO = x[0];
    x[0] = 0.0;
    for (int i = 1; i < x.length; i++) {
      x[i] -= tempO;
    }
    return mbarObjectiveFunction(reducedPotentials, snaps, x);
  }

  /**
   * MBAR objective function evaluation and gradient at a given free energy estimate for L-BFGS optimization.
   *
   * @param x Input parameters.
   * @param g The gradient with respect to each parameter.
   * @return The objective function value at the given parameters.
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    double tempO = x[0];
    x[0] = 0.0;
    for (int i = 1; i < x.length; i++) {
      x[i] -= tempO;
    }
    double[] tempG = mbarGradient(reducedPotentials, snaps, x);
    arraycopy(tempG, 0, g, 0, g.length);
    return mbarObjectiveFunction(reducedPotentials, snaps, x);
  }

  @Override
  public double[] getCoordinates(double[] parameters) {
    return new double[0];
  }

  @Override
  public int getNumberOfVariables() {
    return 0;
  }

  @Override
  public double[] getScaling() {
    return null;
  }

  @Override
  public void setScaling(double[] scaling) {
  }

  @Override
  public double getTotalEnergy() {
    return 0;
  }

  //////// Getters and setters ////////
  public BennettAcceptanceRatio getBAR() {
    return new BennettAcceptanceRatio(lamValues, eLow, eAt, eHigh, temperatures);
  }

  @Override
  public MultistateBennettAcceptanceRatio copyEstimator() {
    return new MultistateBennettAcceptanceRatio(lamValues, eAll, temperatures, tolerance, seedType);
  }

  @Override
  public double[] getBinEnergies() {
    return mbarFEDifferenceEstimates;
  }

  public double[] getMBARFreeEnergies() {
    return mbarFEEstimates;
  }

  @Override
  public double[] getBinUncertainties() {
    return mbarUncertainties;
  }

  public double[][] getDiffMatrix() {
    return diffMatrix;
  }

  @Override
  public double getFreeEnergy() {
    return totalMBAREstimate;
  }

  @Override
  public double getUncertainty() {
    return totalMBARUncertainty;
  }

  @Override
  public int numberOfBins() {
    return nFreeEnergyDiffs;
  }

  @Override
  public double[] getBinEnthalpies() {
    return mbarEnthalpy;
  }

  /**
   * Harmonic oscillators test case generates data for testing the MBAR implementation
   */
  public static class HarmonicOscillatorsTestCase {

    /**
     * Inverse temperature.
     */
    private final double beta;
    /**
     * Equilibrium positions.
     */
    private final double[] O_k;
    /**
     * Number of states.
     */
    private final int n_states;
    /**
     * Spring constants.
     */
    private final double[] K_k;

    /**
     * Constructor for HarmonicOscillatorsTestCase
     *
     * @param O_k  array of equilibrium positions
     * @param K_k  array of spring constants
     * @param beta inverse temperature
     */
    public HarmonicOscillatorsTestCase(double[] O_k, double[] K_k, double beta) {
      this.beta = beta;
      this.O_k = O_k;
      this.n_states = O_k.length;
      this.K_k = K_k;

      if (this.K_k.length != this.n_states) {
        throw new IllegalArgumentException("Lengths of K_k and O_k should be equal");
      }
    }

    public double[] analyticalMeans() {
      return O_k;
    }

    public double[] analyticalVariances() {
      double[] variances = new double[n_states];
      for (int i = 0; i < n_states; i++) {
        variances[i] = 1.0 / (beta * K_k[i]);
      }
      return variances;
    }

    public double[] analyticalStandardDeviations() {
      double[] deviations = new double[n_states];
      for (int i = 0; i < n_states; i++) {
        deviations[i] = Math.sqrt(1.0 / (beta * K_k[i]));
      }
      return deviations;
    }

    public double[] analyticalObservable(String observable) {
      double[] result = new double[n_states];

      switch (observable) {
        case "position" -> {
          return analyticalMeans();
        }
        case "potential energy" -> {
          for (int i = 0; i < n_states; i++) {
            result[i] = 0.5 / beta;
          }
        }
        case "position^2" -> {
          for (int i = 0; i < n_states; i++) {
            result[i] = 1.0 / (beta * K_k[i]) + Math.pow(O_k[i], 2);
          }
        }
        case "RMS displacement" -> {
          return analyticalStandardDeviations();
        }
      }

      return result;
    }

    public double[] analyticalFreeEnergies() {
      int subtractComponentIndex = 0;
      double[] fe = new double[n_states];
      double subtract = 0.0;
      for (int i = 0; i < n_states; i++) {
        fe[i] = -0.5 * Math.log(2 * Math.PI / (beta * K_k[i]));
        if (i == 0) {
          subtract = fe[subtractComponentIndex];
        }
        fe[i] -= subtract;
      }
      return fe;
    }

    public double[] analyticalEntropies(int subtractComponent) {
      double[] entropies = new double[n_states];
      double[] potentialEnergy = analyticalObservable("analytical entropy");
      double[] freeEnergies = analyticalFreeEnergies();

      for (int i = 0; i < n_states; i++) {
        entropies[i] = potentialEnergy[i] - freeEnergies[i];
      }

      return entropies;
    }

    /**
     * Sample from harmonic oscillator w/ gaussian & std
     *
     * @param N_k  number of snaps per state
     * @param mode only u_kn -> return K x N_tot matrix where u_kn[k,n] is reduced potential of sample n evaluated at state k
     * @return u_kn[k, n] is reduced potential of sample n evaluated at state k
     */
    public Object[] sample(int[] N_k, String mode, Long seed) {
      Random random = new Random(seed);

      int N_max = 0;
      for (int N : N_k) {
        if (N > N_max) {
          N_max = N;
        }
      }

      int N_tot = 0;
      for (int N : N_k) {
        N_tot += N;
      }

      double[][] x_kn = new double[n_states][N_max];
      double[][] u_kn = new double[n_states][N_tot];
      double[][][] u_kln = new double[n_states][n_states][N_max];
      double[] x_n = new double[N_tot];
      int[] s_n = new int[N_tot];

      // Sample harmonic oscillators
      int index = 0;
      for (int k = 0; k < n_states; k++) {
        double x0 = O_k[k];
        double sigma = Math.sqrt(1.0 / (beta * K_k[k]));

        // Number of snaps
        for (int n = 0; n < N_k[k]; n++) {
          double x = x0 + random.nextGaussian() * sigma;
          x_kn[k][n] = x;
          x_n[index] = x;
          s_n[index] = k;
          // Potential energy evaluations
          for (int l = 0; l < n_states; l++) {
            double u = beta * 0.5 * K_k[l] * Math.pow(x - O_k[l], 2.0);
            u_kln[k][l][n] = u;
            u_kn[l][index] = u;
          }
          index++;
        }
        // Set the rest of the array to NaN
        for(int n = N_k[k]; n < N_max; n++){
          for ( int l =0; l < n_states; l++) {
            u_kln[k][l][n] = Double.NaN;
            u_kn[l][index] = Double.NaN;
          }
        }
      }

      // Setting corrections
      if ("u_kn".equals(mode)) {
        return new Object[]{x_n, u_kn, N_k, s_n};
      } else if ("u_kln".equals(mode)) {
        return new Object[]{x_n, u_kln, N_k, s_n};
      } else {
        throw new IllegalArgumentException("Unknown mode: " + mode);
      }
    }

    public static Object[] evenlySpacedOscillators(
        int n_states, int n_samplesPerState, double lower_O_k, double upper_O_k,
        double lower_K_k, double upper_K_k, Long seed) {
      // Random random = new Random(seed);

      double[] O_k = new double[n_states];
      double[] K_k = new double[n_states];
      int[] N_k = new int[n_states];

      double stepO_k = (upper_O_k - lower_O_k) / (n_states - 1);
      double stepK_k = (upper_K_k - lower_K_k) / (n_states - 1);

      for (int i = 0; i < n_states; i++) {
        O_k[i] = lower_O_k + i * stepO_k;
        K_k[i] = lower_K_k + i * stepK_k;
        N_k[i] = n_samplesPerState;
      }

      HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(O_k, K_k, 1.0);
      Object[] result = testCase.sample(N_k, "u_kn", System.currentTimeMillis());

      return new Object[]{testCase, result[0], result[1], result[2], result[3]};
    }

    public static void main(String[] args) {
      // Example parameters
      double[] O_k = {0, 1, 2, 3, 4};
      double[] K_k = {1, 2, 4, 8, 16};
      double beta = 1.0;
      System.out.println("Beta: " + beta);

      // Create an instance of HarmonicOscillatorsTestCase
      HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(O_k, K_k, beta);

      // Print results of various functions
      System.out.println("Analytical Means: " + Arrays.toString(testCase.analyticalMeans()));
      System.out.println("Analytical Variances: " + Arrays.toString(testCase.analyticalVariances()));
      System.out.println("Analytical Standard Deviations: " + Arrays.toString(testCase.analyticalStandardDeviations()));
      System.out.println("Analytical Free Energies: " + Arrays.toString(testCase.analyticalFreeEnergies()));

      // Example usage of sample function with u_kn mode
      int[] N_k = {10, 20, 30, 40, 50};
      String setting = "u_kln";
      Object[] sampleResult = testCase.sample(N_k, setting, System.currentTimeMillis());

      System.out.println("Sample x_n: " + Arrays.toString((double[]) sampleResult[0]));
      if ("u_kn".equals(setting)) {
        System.out.println("Sample u_kn: " + Arrays.deepToString((double[][]) sampleResult[1]));
      } else {
        System.out.println("Sample u_kln: " + Arrays.deepToString((double[][][]) sampleResult[1]));
      }
      System.out.println("Sample N_k: " + Arrays.toString((int[]) sampleResult[2]));
      System.out.println("Sample s_n: " + Arrays.toString((int[]) sampleResult[3]));
    }
  }

  public static void writeFile(double[][] energies, File file, double temperature) {
    try (FileWriter fw = new FileWriter(file);
         BufferedWriter bw = new BufferedWriter(fw)) {
      // Write the number of snapshots and the temperature on the first line
      bw.write(energies[0].length + " " + temperature);
      bw.newLine();

      // Write the energies
      StringBuilder sb = new StringBuilder();
      for (int i = 0; i < energies[0].length; i++) {
        sb.append("     ").append(i).append(" "); // Write the index of the snapshot
        for (int j = 0; j < energies.length; j++) {
          sb.append("    ").append(energies[j][i]).append(" ");
        }
        sb.append("\n");
        bw.write(sb.toString());
        sb = new StringBuilder(); // Very important
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    double[] O_k = {0, 1, 2, 3, 4}; // Equilibrium positions
    double[] K_k = {1, 3, 7, 10, 15}; // Spring constants
    int[] N_k = {0, 500, 300, 0, 10000}; // No support for different number of snapshots
    double beta = 1.0;

    // Create an instance of HarmonicOscillatorsTestCase
    HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(O_k, K_k, beta);

    // Generate sample data
    String setting = "u_kln";
    System.out.print("Generating sample data... ");
    Object[] sampleResult = testCase.sample(N_k, setting, (long) 2); // Set seed to fixed value for reproducibility
    System.out.println("done. \n");
    double[][][] u_kln = (double[][][]) sampleResult[1];
    double[] temps = {1 / Constants.R};

    // Write file for comparison with pymbar
    // Output to forcefieldx/testing/mbar/data/harmonic_oscillators/mbarFiles/energies_{i}.mbar
    // Get absolute path to root of project

    /*
    String rootPath = new File("").getAbsolutePath();
    File outputPath = new File(rootPath + "/testing/mbar/data/harmonic_oscillators/mbarFiles");
    if (!outputPath.exists() && !outputPath.mkdirs()) {
      throw new RuntimeException("Failed to create directory: " + outputPath);
    }

    double[] temperatures = new double[O_k.length];
    Arrays.fill(temperatures, temps[0]);
    for (int i = 0; i < u_kln.length; i++) {
      File file = new File(outputPath, "energies_" + i + ".mbar");
      writeFile(u_kln[i], file, temperatures[i]);
    } */

    // Create an instance of MultistateBennettAcceptanceRatio
    System.out.print("Creating MBAR instance and estimateDG() with standard tol & Zeros seeding.");
    //MBARFilter mbarFilter = new MBARFilter(new File("/Users/matthewsperanza/Programs/forcefieldx/testing/mbar/LYS_Umod/zeroSampleTest"));
    MultistateBennettAcceptanceRatio mbar = new MultistateBennettAcceptanceRatio(O_k, u_kln, temps, 1e-7, SeedType.ZEROS);
    double[] mbarFEEstimates = Arrays.copyOf(mbar.mbarFEEstimates, mbar.mbarFEEstimates.length);
    double[] mbarUncertainties = Arrays.copyOf(mbar.mbarUncertainties, mbar.mbarUncertainties.length);
    double[][] mbarDiffMatrix = Arrays.copyOf(mbar.diffMatrix, mbar.diffMatrix.length);

    EstimateBootstrapper bootstrapper = new EstimateBootstrapper(mbar);
    bootstrapper.bootstrap(50);
    System.out.println("done. \n");

    // Get the analytical free energy differences
    double[] analyticalFreeEnergies = testCase.analyticalFreeEnergies();
    // Calculate the error
    double[] error = new double[analyticalFreeEnergies.length];
    for (int i = 0; i < error.length; i++) {
      error[i] = -mbarFEEstimates[i] + analyticalFreeEnergies[i];
    }

    // Compare the calculated free energy differences with the analytical ones
    System.out.println("Analytical Free Energies: " + Arrays.toString(analyticalFreeEnergies));
    System.out.println("MBAR Free Energies:       " + Arrays.toString(mbarFEEstimates));
    System.out.println("MBAR Uncertainties:       " + Arrays.toString(mbarUncertainties));
    System.out.println("Free Energy Error:        " + Arrays.toString(error));
    System.out.println();
    System.out.println("Diff Matrix: ");
    for (double[] matrix : mbarDiffMatrix) {
      System.out.println(Arrays.toString(matrix));
    }
    System.out.println("\n\n");

    // Get the calculated free energy differences
    double[] mbarBootstrappedEstimates = bootstrapper.getFE();
    double[] mbarBootstrappedFE = new double[mbarBootstrappedEstimates.length + 1];
    for (int i = 0; i < mbarBootstrappedEstimates.length; i++) {
      mbarBootstrappedFE[i + 1] = mbarBootstrappedEstimates[i] + mbarBootstrappedFE[i];
    }
    mbarUncertainties = bootstrapper.getUncertainty();
    // Calculate the error
    double[] errors = new double[mbarBootstrappedFE.length];
    for (int i = 0; i < errors.length; i++) {
      errors[i] = -mbarBootstrappedFE[i] + analyticalFreeEnergies[i];
    }

    System.out.println("Analytical Estimates:         " + Arrays.toString(analyticalFreeEnergies));
    System.out.println("MBAR Bootstrapped Estimates:  " + Arrays.toString(mbarBootstrappedFE));
    System.out.println("MBAR Bootstrap Uncertainties: " + Arrays.toString(mbarUncertainties));
    System.out.println("Bootstrap Free Energy Error:  " + Arrays.toString(errors));
  }
}
