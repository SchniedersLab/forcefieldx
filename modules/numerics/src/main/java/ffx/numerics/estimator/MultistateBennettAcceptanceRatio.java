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
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

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
import static org.apache.commons.math3.util.FastMath.*;

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
  private double[] mbarFEEstimates;
  /**
   * MBAR observable ensemble estimates.
   */
  private double[] mbarObservableEnsembleAverages;
  private double[] mbarObservableEnsembleAverageUncertainties;
  /**
   * MBAR free-energy difference uncertainties.
   */
  private double[] mbarUncertainties;
  /**
   * Matrix of free-energy difference uncertainties between all i & j
   */
  private double[][] uncertaintyMatrix;
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
  private double[] mbarEnthalpy;

  /**
   * MBAR Entropy estimates
   */
  private double[] mbarEntropy;

  public double[] rtValues;

  /**
   * "Reduced" potential energies. -ln(exp(beta * -U)) or more practically U * (1 / RT).
   * Has shape (nLambdaStates, numSnaps * nLambdaStates)
   */
  private double[][] reducedPotentials;

  private double[][] oAllFlat;
  private double[][] biasFlat;
  /**
   * Seed MBAR calculation with another free energy estimation (BAR,ZWANZIG) or zeros
   */
  private SeedType seedType;

  /**
   * Enum of MBAR seed types.
   */
  public enum SeedType {BAR, ZWANZIG, ZEROS}
  public static boolean FORCE_ZEROS_SEED = false;
  public static boolean VERBOSE = false;

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
    mbarEntropy = new double[nFreeEnergyDiffs];
    random = new Random();
    estimateDG();
  }

  /**
   * Set the MBAR seed energies using BAR, Zwanzig, or zeros.
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
   * MBAR solved with self-consistent iteration and Newton/L-BFGS optimization.
   */
  @Override
  public void estimateDG(boolean randomSamples) {
    if (MultistateBennettAcceptanceRatio.VERBOSE){
      logger.setLevel(java.util.logging.Level.FINE);
    }

    // Bootstrap needs resetting to zeros
    fill(mbarFEEstimates, 0.0);
    if (FORCE_ZEROS_SEED) {
      seedType = SeedType.ZEROS;
    }
    seedEnergies();
    if (stream(mbarFEEstimates).anyMatch(Double::isInfinite) || stream(mbarFEEstimates).anyMatch(Double::isNaN)) {
      seedType = SeedType.ZEROS;
      seedEnergies();
    }
    if(MultistateBennettAcceptanceRatio.VERBOSE) {
      logger.info(" Seed Type: " + seedType);
      logger.info(" MBAR FE Estimates after seeding: " + Arrays.toString(mbarFEEstimates));
    }

    // Precompute beta for each state.
    rtValues = new double[nLambdaStates];
    double[] invRTValues = new double[nLambdaStates];
    for (int i = 0; i < nLambdaStates; i++) {
      rtValues[i] = Constants.R * temperatures[i];
      invRTValues[i] = 1.0 / rtValues[i];
    }

    // Find repeated snapshots if from continuous lambda
    int numEvaluations = eAllFlat[0].length;

    // Sample random snapshots from each window.
    int[][] indices = new int[nLambdaStates][numEvaluations];
    if (randomSamples) {
      // Build random indices vector maintaining snapshot nums!
      int[] randomIndices = new int[numEvaluations];
      int sum = 0;
      for (int snap : snaps) {
        System.arraycopy(getBootstrapIndices(snap, random), 0, randomIndices, sum, snap);
        sum += snap;
      }
      for (int i = 0; i < nLambdaStates; i++) {
        // Use the same random indices across lambda values
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
    for (int state = 0; state < eAllFlat.length; state++) { // For each lambda value
      for (int n = 0; n < eAllFlat[0].length; n++) {
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
      logger.info(" Sampled Lambdas: " + sampledLambdas);
      logger.info(" Zero Snap Lambdas: " + zeroSnapLambdas);
    } else { // If there aren't any zero snap lambdas, just use the original arrays
      reducedPotentialsTemp = reducedPotentials;
      mbarFEEstimatesTemp = mbarFEEstimates;
      snapsTemp = snaps;
    }

    // SCI iterations used to start optimization of MBAR objective function.
    // Optimizers can struggle when starting too far from the minimum, but SCI doesn't.
    double[] prevMBAR = copyOf(mbarFEEstimatesTemp, nLambdaStatesTemp);;
    double omega = 1.5; // Parameter chosen empirically to work with most systems (> 2 works but not always).
    for (int i = 0; i < 10; i++) {
      prevMBAR = copyOf(mbarFEEstimatesTemp, nLambdaStatesTemp);
      mbarFEEstimatesTemp = mbarSelfConsistentUpdate(reducedPotentialsTemp, snapsTemp, mbarFEEstimatesTemp);
      for (int j = 0; j < nLambdaStatesTemp; j++) { // SOR
        mbarFEEstimatesTemp[j] = omega * mbarFEEstimatesTemp[j] + (1 - omega) * prevMBAR[j];
      }
      if (stream(mbarFEEstimatesTemp).anyMatch(Double::isInfinite) || stream(mbarFEEstimatesTemp).anyMatch(Double::isNaN)) {
        throw new IllegalArgumentException("MBAR contains NaNs or Infs during startup SCI ");
      }
      if(converged(prevMBAR)) {
        break;
      }
    }
    if (MultistateBennettAcceptanceRatio.VERBOSE) {
      logger.info(" Omega for SCI w/ relaxation: " + omega);
      logger.info(" MBAR FE Estimates after 10 SCI iterations: " + Arrays.toString(mbarFEEstimatesTemp));
    }

    try {
      if (nLambdaStatesTemp > 100 && !converged(prevMBAR)) { // L-BFGS optimization for high granularity windows where hessian^-1 is expensive
        if (MultistateBennettAcceptanceRatio.VERBOSE) {
          logger.info(" L-BFGS optimization started.");
        }
        int mCorrections = 5;
        double[] x = new double[nLambdaStatesTemp];
        arraycopy(mbarFEEstimatesTemp, 0, x, 0, nLambdaStatesTemp);
        double[] grad = mbarGradient(reducedPotentialsTemp, snapsTemp, mbarFEEstimatesTemp);
        double eps = 1.0E-4; // Gradient tolarance -> chosen since L-BFGS seems unstable with tight tolerances
        OptimizationListener listener = getOptimizationListener();
        LBFGS.minimize(nLambdaStatesTemp, mCorrections, x, mbarObjectiveFunction(reducedPotentialsTemp, snapsTemp, mbarFEEstimatesTemp),
            grad, eps, 1000, this, listener);
        arraycopy(x, 0, mbarFEEstimatesTemp, 0, nLambdaStatesTemp);
      } else if (!converged(prevMBAR)){ // Newton optimization if hessian inversion isn't too expensive
        if (MultistateBennettAcceptanceRatio.VERBOSE) {
          logger.info(" Newton optimization started.");
        }
        mbarFEEstimatesTemp = newton(mbarFEEstimatesTemp, reducedPotentialsTemp, snapsTemp, tolerance);
      }
    } catch (Exception e) {
      logger.warning(" L-BFGS/Newton failed to converge. Finishing w/ self-consistent iteration. Message: " +
              e.getMessage());
    }
    if(MultistateBennettAcceptanceRatio.VERBOSE) {
      logger.info(" MBAR FE Estimates after gradient optimization: " + Arrays.toString(mbarFEEstimatesTemp));
    }

    // Update the FE estimates with the optimized values from derivative-based optimization
    int count = 0;
    for(Integer i : sampledLambdas){
      if (!Double.isNaN(mbarFEEstimatesTemp[count])) { // Should be !NaN
        mbarFEEstimates[i] = mbarFEEstimatesTemp[count];
      }
      count++;
    }

    // Self-consistent iteration is used to finish off optimization of MBAR objective function
    int sciIter = 0;
    while (!converged(prevMBAR) && sciIter < 1000) {
      prevMBAR = copyOf(mbarFEEstimates, nLambdaStates);
      mbarFEEstimates = mbarSelfConsistentUpdate(reducedPotentials, snaps, mbarFEEstimates);
      for (int i = 0; i < nLambdaStates; i++) { // SOR for acceleration
        mbarFEEstimates[i] = omega * mbarFEEstimates[i] + (1 - omega) * prevMBAR[i];
      }
      if (stream(mbarFEEstimates).anyMatch(Double::isInfinite) || stream(mbarFEEstimates).anyMatch(Double::isNaN)) {
        throw new IllegalArgumentException("MBAR estimate contains NaNs or Infs after iteration " + sciIter);
      }
      sciIter++;
    }
    if (MultistateBennettAcceptanceRatio.VERBOSE) {
      logger.info(" SCI iterations (max 1000): " + sciIter);
    }

    // Calculate uncertainties
    double[][] theta = mbarTheta(reducedPotentials, snaps, mbarFEEstimates); // Quite expensive
    mbarUncertainties = mbarUncertaintyCalc(theta);
    totalMBARUncertainty = mbarTotalUncertaintyCalc(theta);
    uncertaintyMatrix = diffMatrixCalculation(theta);
    if (!randomSamples && MultistateBennettAcceptanceRatio.VERBOSE) { // Never log for bootstrapping
      logWeights();
    }

    // Convert to kcal/mol & calculate differences/sums
    for (int i = 0; i < nLambdaStates; i++) {
      mbarFEEstimates[i] = mbarFEEstimates[i] * rtValues[i];
    }
    for (int i = 0; i < nFreeEnergyDiffs; i++) {
      mbarFEDifferenceEstimates[i] = mbarFEEstimates[i + 1] - mbarFEEstimates[i];
    }

    mbarEnthalpy = mbarEnthalpyCalc(eAllFlat, mbarFEEstimates);
    mbarEntropy = mbarEntropyCalc(mbarEnthalpy, mbarFEEstimates);

    totalMBAREstimate = stream(mbarFEDifferenceEstimates).sum();
  }

  //////// Misc. Methods ////////////

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

  /**
   * Print out, for each FE expectation, the sum of the weights for each trajectory. This
   * gives an array of length nLambdaStates, where each element is the sum of the weights
   * coming from the trajectory sampled at the lambda value corresponding to that index.
   *
   * <p>i.e. collapsedW[0][0] is the sum of the weights in W[0] from the trajectory sampled at
   * lambda 0. The diagonal of this matrix should be larger than all other values if that
   * window had proper sampling.
   *
   */
  private void logWeights() {
    logger.info(" MBAR Weight Matrix Information Collapsed:");
    double[][] W = mbarW(reducedPotentials, snaps, mbarFEEstimates);
    double[][] collapsedW = new double[W.length][W.length]; // Collapse W trajectory-wise (into K x K)
    for (int i = 0; i < snaps.length; i++) {
      for (int j = 0; j < W.length; j++) {
        int start = 0;
        for(int k = 0; k < i; k++) {
          start += snaps[k];
        }
        for(int k = 0; k < snaps[i]; k++) {
          collapsedW[j][i] += W[j][start + k];
        }
      }
    }
    for(int i = 0; i < W.length; i++) {
      logger.info( "\n Estimation " + i + ": " + Arrays.toString(collapsedW[i]));
    }
    double[] rowSum = new double[W.length];
    for(int i = 0; i < collapsedW[0].length; i++) {
      for (double[] trajectory : collapsedW) {
        rowSum[i] += trajectory[i];
      }
    }
    softMax(rowSum);
    logger.info("\n Softmax of trajectory weight: " + Arrays.toString(rowSum));
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
   * Gradient of the MBAR objective function. C6 in Shirts and Chodera 2008.
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
    double[] maxLogDiff = new double[nStates];
    Arrays.fill(maxLogDiff, Double.NEGATIVE_INFINITY);
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
        if (logDiff[j][i] > maxLogDiff[j]) {
          maxLogDiff[j] = logDiff[j][i];
        }
      }
    }
    for (int i = 0; i < nStates; i++) {
      log_num_k[i] = logSumExp(logDiff[i], maxLogDiff[i]);
    }
    double[] grad = new double[nStates];
    for (int i = 0; i < nStates; i++) {
      grad[i] = -1.0 * snapsPerLambda[i] * (1.0 - exp(freeEnergyEstimates[i] + log_num_k[i]));
    }
    return grad;
  }

  /**
   * Hessian of the MBAR objective function. C9 in Shirts and Chodera 2008.
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
   * Eq. 9 in Shirts and Chodera 2008.
   *
   * @param reducedPotentials energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return W matrix.
   */
  private static double[][] mbarW(double[][] reducedPotentials, int[] snapsPerLambda, double[] freeEnergyEstimates) {
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
      // log_denom_n = calculates log(sumOverStates(N_k * exp(FE[j] - reducedPotentials[j][i])))
      log_denom_n[i] = logSumExp(temp, snapsPerLambda, maxTemp);
    }
    // logW = freeEnergyEstimates - reducedPotentials.T - log_denominator_n[:, newaxis]
    // freeEnergyEstimates[i] = log(ck / ci) --> ratio of normalization constants
    double[][] W = new double[nStates][reducedPotentials[0].length];
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < reducedPotentials[0].length; j++) {
        W[i][j] = exp(freeEnergyEstimates[i] - reducedPotentials[i][j] - log_denom_n[j]);
      }
    }
    return W;
  }

  private double[] mbarEnthalpyCalc(double[][] reducedPotentials, double[] mbarFEEstimates) {
    double[] enthalpy = new double[mbarFEEstimates.length - 1];
    double[] averagePotential = new double[mbarFEEstimates.length];
    for(int i = 0; i < reducedPotentials.length; i++) {
      averagePotential[i] = computeExpectations(eAllFlat[i])[i]; // average potential of ith lambda
    }
    for(int i = 0; i < enthalpy.length; i++) {
      enthalpy[i] = averagePotential[i + 1] - averagePotential[i];
    }
    return enthalpy;
  }

  private double[] mbarEntropyCalc(double[] mbarEnthalpy, double[] mbarFEEstimates) {
    double[] entropy = new double[mbarFEEstimates.length - 1];
    for(int i = 0; i < entropy.length; i++) {
      entropy[i] = mbarEnthalpy[i] - mbarFEDifferenceEstimates[i]; // dG = dH - TdS || TdS = dH - dG
    }
    return entropy;
  }

  /**
   * Weight observable by exp(bias/RT) prior to computing expectation when set.
   *
   * @param biasAll
   * @param multiDataObservable
   */
  public void setBiasData(double[][][] biasAll, boolean multiDataObservable) {
    biasFlat = new double[biasAll.length][biasAll.length * biasAll[0][0].length];
    if(multiDataObservable){ // Flatten data
      int[] snapsT = new int[biasAll.length];
      int[] nanCount = new int[biasAll.length];
      for (int i = 0; i < biasAll.length; i++) {
        ArrayList<Double> temp = new ArrayList<>();
        double maxBias = Double.NEGATIVE_INFINITY;
        for(int j = 0; j < biasAll.length; j++) {
          int count = 0;
          int countNaN = 0;
          for(int k = 0; k < biasAll[j][i].length; k++) {
            // Don't include NaN values
            if (!Double.isNaN(biasAll[j][i][k])) {
              temp.add(biasAll[j][i][k]);
              if(biasAll[j][i][k] > maxBias){
                maxBias = biasAll[j][i][k];
              }
              count++;
            } else {
              countNaN++;
            }
          }
          snapsT[j] = count;
          nanCount[j] = countNaN;
        }
        biasFlat[i] = temp.stream().mapToDouble(Double::doubleValue).toArray();
        // Regularize bias for this lambda
        for(int j = 0; j < biasFlat[i].length; j++){
          biasFlat[i][j] -= maxBias;
        }
      }
    } else { // Put relevant data into the 0th index
      int count = 0;
      double maxBias = Double.NEGATIVE_INFINITY;
      for (int i = 0; i < biasAll.length; i++){
        for(int j = 0; j < biasAll[0][0].length; j++){
          if(!Double.isNaN(biasAll[i][i][j])){
            biasFlat[0][count] = biasAll[i][i][j];
            if(biasAll[i][i][j] > maxBias){
              maxBias = biasAll[i][i][j];
            }
            count++;
          }
        }
      }
      // Regularize bias for this lambda
      for(int i = 0; i < biasFlat[0].length; i++){
        biasFlat[0][i] -= maxBias;
      }
    }
  }

  public void setObservableData(double[][][] oAll, boolean multiDataObservable, boolean uncertainties) {
    oAllFlat = new double[oAll.length][oAll.length * oAll[0][0].length];
    if(multiDataObservable){ // Flatten data
      int[] snapsT = new int[oAll.length];
      int[] nanCount = new int[oAll.length];
      for (int i = 0; i < oAll.length; i++) {
        ArrayList<Double> temp = new ArrayList<>();
        for(int j = 0; j < oAll.length; j++) {
          int count = 0;
          int countNaN = 0;
          for(int k = 0; k < oAll[j][i].length; k++) {
            // Don't include NaN values
            if (!Double.isNaN(oAll[j][i][k])) {
              temp.add(oAll[j][i][k]);
              count++;
            } else {
              countNaN++;
            }
          }
          snapsT[j] = count;
          nanCount[j] = countNaN;
        }
        oAllFlat[i] = temp.stream().mapToDouble(Double::doubleValue).toArray();
      }
    } else { // Put relevant data into the 0th index
      int count = 0;
      for (int i = 0; i < oAll.length; i++){
        for(int j = 0; j < oAll[0][0].length; j++){
          if(!Double.isNaN(oAll[i][i][j])){
            oAllFlat[0][count] = oAll[i][i][j]; // Note [i][i] indexing
            count++;
          }
        }
      }
    }
    // OST Data
    if (biasFlat != null) {
      for(int i = 0; i< oAllFlat.length; i++) {
        for (int j = 0; j < oAllFlat[i].length; j++) {
          oAllFlat[i][j] *= exp(biasFlat[i][j]/rtValues[i]);
        }
      }
    }
    this.fillObservationExpectations(multiDataObservable, uncertainties);
  }

  /**
   * Calculate expectation of samples from W matrix. Optionally calculate the uncertainty with
   * augmented W matrix (incurs a significant computational cost ~10-20x MBAR calculation).
   *
   * @return Uncertainty of the observable.
   */
  private void fillObservationExpectations(boolean multiData, boolean uncertainties){
    if(multiData){
      mbarObservableEnsembleAverages = new double[oAllFlat.length];
      mbarObservableEnsembleAverageUncertainties = new double[oAllFlat.length];
      for(int i = 0; i < oAllFlat.length; i++){
        mbarObservableEnsembleAverages[i] = computeExpectations(oAllFlat[i])[i];
        if (uncertainties) {
          mbarObservableEnsembleAverageUncertainties[i] = computeExpectationStd(oAllFlat[i])[i];
        }
      }
    } else {
      mbarObservableEnsembleAverages = computeExpectations(oAllFlat[0]);
      if (uncertainties) {
        mbarObservableEnsembleAverageUncertainties = computeExpectationStd(oAllFlat[0]);
      }
    }
  }

  /**
   * Compute the MBAR expectation of a given observable (1xN) for each K. This observable
   * could be something like x, x^2 (where x is equilibrium for a harmonic oscillator),
   * or some other function of the configuration X like RMSD from a target conformation.
   * Additionally, it could be evaluations of some potential at a specific lambda value.
   * Each trajectory snap should have a corresponding observable value (or evaluation).
   *
   * @param samples
   * @return
   */
  private double[] computeExpectations(double[] samples){
    double[][] W = mbarW(reducedPotentials, snaps, mbarFEEstimates);
    double[] expectation = new double[W.length];
    for(int i = 0; i < W.length; i++){
      for(int j = 0; j < W[i].length; j++){
        expectation[i] += W[i][j] * samples[j];
      }
    }
    return expectation;
  }

  /**
   * Eq. 13-15 in Shirts and Chodera (2008) for the MBAR observable uncertainty calculation.
   * Originally implemented as seen in paper, but switched to logsumexp version because of
   * Inf/NaN issues for large values captured in samples (i.e. potential energies).
   *
   * @return WnA matrix.
   */
  private double[][] mbarAugmentedW(double[] samples) {
    int nStates = mbarFEEstimates.length;
    // Enforce positivity of samples --> from pymbar
    double minSample = stream(samples).min().getAsDouble() - 3*java.lang.Math.ulp(1.0); // ulp to avoid zeros
    if (minSample < 0) {
      for (int i = 0; i < samples.length; i++) {
        samples[i] -= minSample;
      }
    }
    // Eq. 14 in Shirts and Chodera (2008)
    double[][] logCATerms = new double[nStates][reducedPotentials[0].length];
    double[] maxLogCATerm = new double[reducedPotentials[0].length];
    Arrays.fill(maxLogCATerm, Double.NEGATIVE_INFINITY);
    double[] logCA = new double[nStates];
    double[] log_denom_n = new double[reducedPotentials[0].length];
    for (int i = 0; i < reducedPotentials[0].length; i++) {
      double[] temp = new double[nStates];
      double maxTemp = Double.NEGATIVE_INFINITY;
      for (int j = 0; j < nStates; j++) {
        temp[j] = mbarFEEstimates[j] - reducedPotentials[j][i];
        if (temp[j] > maxTemp) {
          maxTemp = temp[j];
        }
      }
      log_denom_n[i] = logSumExp(temp, snaps, maxTemp);
      for(int j = 0; j < nStates; j++){
        logCATerms[j][i] = log(samples[i]) - reducedPotentials[j][i] - log_denom_n[i];
        if(logCATerms[j][i] > maxLogCATerm[i]){
          maxLogCATerm[j] = logCATerms[j][i];
        }
      }
    }
    for(int i = 0; i < nStates; i++){
      logCA[i] = logSumExp(logCATerms[i], maxLogCATerm[i]);
    }
    // Eq. 13 in Shirts and Chodera (2008)
    double[][] WnA = new double[nStates][reducedPotentials[0].length];
    double[][] Wna = new double[nStates][reducedPotentials[0].length]; // normal W matrix
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < reducedPotentials[0].length; j++) {
        WnA[i][j] = samples[j] * exp(-logCA[i] - reducedPotentials[i][j] - log_denom_n[j]);
        Wna[i][j] = exp(-mbarFEEstimates[i] - reducedPotentials[i][j] - log_denom_n[j]);
      }
    }
    if (minSample < 0) { // reset samples
      for (int i = 0; i < samples.length; i++) {
        samples[i] += minSample;
      }
    }
    double[][] augmentedW = new double[nStates * 2][reducedPotentials[0].length];
    for (int i = 0; i < augmentedW.length; i++) {
      augmentedW[i] = i < nStates ? Wna[i] : WnA[(i-nStates)];
    }
    return augmentedW;
  }

  /**
   * Compute the MBAR uncertainty of an observable. The equations for this are not clear,
   * but we append an augmented weight matrix (calculated by multiplying the observed values
   * into the W matrix calculation) to the original W matrix. This is then used to calculate
   * theta.
   *
   * @param samples
   * @return
   */
  private double[] computeExpectationStd(double[] samples){
    int[] extendedSnaps = new int[snaps.length * 2];
    System.arraycopy(snaps, 0, extendedSnaps, 0, snaps.length);
    RealMatrix theta = MatrixUtils.createRealMatrix(mbarTheta(extendedSnaps, mbarAugmentedW(samples)));
    double[] expectations = computeExpectations(samples);
    double[] diag = new double[expectations.length*2];
    for(int i = 0; i < expectations.length; i++){
      diag[i] = expectations[i];
      diag[i+expectations.length] = expectations[i];
    }
    RealMatrix diagMatrix = MatrixUtils.createRealDiagonalMatrix(diag);
    theta = diagMatrix.multiply(theta).multiply(diagMatrix);
    RealMatrix ul = theta.getSubMatrix(0, expectations.length-1, 0, expectations.length-1);
    RealMatrix ur = theta.getSubMatrix(0, expectations.length-1, expectations.length, expectations.length*2-1);
    RealMatrix ll = theta.getSubMatrix(expectations.length, expectations.length*2-1, 0, expectations.length-1);
    RealMatrix lr = theta.getSubMatrix(expectations.length, expectations.length*2-1, expectations.length, expectations.length*2-1);
    double[][] covA = ul.add(lr).subtract(ur).subtract(ll).getData(); // Loose precision here
    double[] sigma = new double[covA.length];
    for(int i = 0; i < covA.length; i++){
      sigma[i] = sqrt(abs(covA[i][i]));
    }
    return sigma;
  }

  /**
   * MBAR uncertainty calculation.
   *
   * @return Uncertainties for the MBAR free energy estimates.
   */
  private static double[] mbarUncertaintyCalc(double[][] theta) {
    double[] uncertainties = new double[theta.length - 1];
    // del(dFij) = Theta[i,i] - 2 * Theta[i,j] + Theta[j,j]
    for (int i = 0; i < theta.length - 1; i++) {
      // TODO: Figure out why negative var is happening (likely due to theta calculation differing from pymbar's)
      double variance = theta[i][i] - 2 * theta[i][i + 1] + theta[i + 1][i + 1];
      if (variance < 0) {
        if (MultistateBennettAcceptanceRatio.VERBOSE) {
          logger.warning(" Negative variance detected in MBAR uncertainty calculation. " +
                  "Multiplying by -1 to get real value. Check diff matrix to see which variances were negative. " +
                  "They should be NaN.");
        }
        variance *= -1;
      }
      uncertainties[i] = sqrt(variance);
    }
    return uncertainties;
  }

  /**
   * MBAR total uncertainty calculation. Eq 12 in Shirts and Chodera (2008).
   *
   * @param theta matrix of covariances
   * @return Total uncertainty for the MBAR free energy estimates.
   */
  private static double mbarTotalUncertaintyCalc(double[][] theta) {
    int nStates = theta.length;
    return sqrt(abs(theta[0][0] - 2 * theta[0][nStates - 1] + theta[nStates - 1][nStates - 1]));
  }

  /**
   * Theta = W.T @ (I - W @ diag(snapsPerState) @ W.T)^-1 @ W.
   * <p>
   * Requires calculation and inversion of W matrix.
   * D4 from supp info of MBAR paper used instead to reduce storage and comp. complexity.
   *
   * @param reducedPotentials energies
   * @param snapsPerState  number of snaps per state
   * @param freeEnergies  free energies
   * @return Theta matrix.
   */
  private static double[][] mbarTheta(double[][] reducedPotentials, int[] snapsPerState, double[] freeEnergies) {
    return mbarTheta(snapsPerState, mbarW(reducedPotentials, snapsPerState, freeEnergies));
  }

  /**
   * Compute theta with a given W matrix.
   *
   * @param snapsPerState
   * @param W
   * @return
   */
  private static double[][] mbarTheta(int[] snapsPerState, double[][] W) {
    RealMatrix WMatrix = MatrixUtils.createRealMatrix(W).transpose();
    RealMatrix I = MatrixUtils.createRealIdentityMatrix(snapsPerState.length);
    RealMatrix NkMatrix = MatrixUtils.createRealDiagonalMatrix(stream(snapsPerState).mapToDouble(i -> i).toArray());
    SingularValueDecomposition svd = new SingularValueDecomposition(WMatrix);
    RealMatrix V = svd.getV();
    RealMatrix S = MatrixUtils.createRealDiagonalMatrix(svd.getSingularValues());

    // W.T @ (I - W @ diag(snapsPerState) @ W.T)^-1 @ W
    // = V @ S @ (I - S @ V.T @ diag(snapsPerState) @ V @ S)^-1 @ S @ V.T
    RealMatrix theta = S.multiply(V.transpose());
    theta = theta.multiply(NkMatrix).multiply(V).multiply(S);
    theta = I.subtract(theta);
    theta = MatrixUtils.inverse(theta); // pinv equivalent
    theta = V.multiply(S).multiply(theta).multiply(S).multiply(V.transpose());

    return theta.getData();
  }

  /**
   * MBAR uncertainty matrix calculation. diff[i][j] gives FE uncertainty of moving between
   * lambda i-> j.
   *
   * @param theta matrix of covariances
   * @return Diff matrix for the MBAR free energy estimates.
   */
  private static double[][] diffMatrixCalculation(double[][] theta) {
    double[][] diffMatrix = new double[theta.length][theta.length];
    for (int i = 0; i < diffMatrix.length; i++) {
      for (int j = 0; j < diffMatrix.length; j++) {
        diffMatrix[i][j] = sqrt(theta[i][i] - 2 * theta[i][j] + theta[j][j]);
      }
    }
    return diffMatrix;
  }

  //////// Methods for solving MBAR with self-consistent iteration, L-BFGS optimization, and Newton-Raphson. ////////

  /**
   * Self-consistent iteration to update free energies. Eq. 11 from Shirts and Chodera (2008).
   *
   * @param reducedPotential energies
   * @param snapsPerLambda  number of snaps per state
   * @param freeEnergyEstimates  free energies
   * @return updated free energies
   */
  private static double[] mbarSelfConsistentUpdate(double[][] reducedPotential, int[] snapsPerLambda,
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
   * Newton-Raphson step for MBAR optimization. Falls back to the steepest descent if hessian is singular.
   *
   * The matrix can come back from being singular after several iterations, so it isn't worth moving to L-BFGS.
   *
   * @param n        current free energies.
   * @param grad     gradient of the objective function.
   * @param hessian  hessian of the objective function.
   * @param stepSize step size for the Newton-Raphson step.
   * @return updated free energies.
   */
  private static double[] newtonStep(double[] n, double[] grad, double[][] hessian, double stepSize) {
    double[] nPlusOne = new double[n.length];
    double[] step;
    try {
      RealMatrix hessianInverse = MatrixUtils.inverse(MatrixUtils.createRealMatrix(hessian));
      step = hessianInverse.preMultiply(grad);
    } catch (IllegalArgumentException e){
        if(MultistateBennettAcceptanceRatio.VERBOSE) {
          logger.info(" Singular matrix detected in MBAR Newton-Raphson step. Performing steepest descent step.");
        }
        step = grad;
        stepSize = 1e-5;
    }
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
   * @param freeEnergyEstimates free energies.
   * @param reducedPotentials   energies.
   * @param snapsPerLambda      number of snaps per state.
   * @param tolerance           convergence tolerance.
   * @return updated free energies.
   */
  private static double[] newton(double[] freeEnergyEstimates, double[][] reducedPotentials,
                                 int[] snapsPerLambda, double tolerance) {
    double[] grad = mbarGradient(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[][] hessian = mbarHessian(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[] f_kPlusOne = newtonStep(freeEnergyEstimates, grad, hessian, 1.0);
    int iter = 1;
    while (iter < 15) { // Quadratic convergence is expected, SCI will run anyway
      freeEnergyEstimates = f_kPlusOne;
      grad = mbarGradient(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
      hessian = mbarHessian(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
      // Catches singular matrices and performs steepest descent
      f_kPlusOne = newtonStep(freeEnergyEstimates, grad, hessian, 1.0);
      double eps = 0.0;
      for(int i = 0; i < freeEnergyEstimates.length; i++){
        eps += abs(grad[i]);
      }
      if (eps < tolerance) {
        break;
      }
      iter++;
    }
    if(MultistateBennettAcceptanceRatio.VERBOSE) {
      logger.info(" Newton iterations (max 15): " + iter);
    }

    return f_kPlusOne;
  }

  /**
   * Calculates the log of the sum of the exponential of the given values.
   * <p>
   * The max value is subtracted from each value in the array before exponentiation to prevent overflow.
   *
   * @param values The values to exponential and sum.
   * @param max    The max value is subtracted from each value in the array prior to exponentiation.
   * @return the sum
   */
  private static double logSumExp(double[] values, double max) {
    int[] b = fill(new int[values.length], 1);
    return logSumExp(values, b, max);
  }

  /**
   * Calculates the log of the sum of the exponential of the given values.
   * <p>
   * The max value is subtracted from each value in the array before exponentiation to prevent overflow.
   * MBAR calculation is easiest to do in log terms, only exponentiating when required. Prevents zeros
   * in the denominator.
   *
   * @param values The values to exponential and sum.
   * @param max    The max value is subtracted from each value in the array prior to exponentiation.
   * @param b      Weights for each value in the array.
   * @return the sum
   */
  private static double logSumExp(double[] values, int[] b, double max) {
    // ChatGPT mostly wrote this and I tweaked it to match more closely with scipy's log-sum-exp implementation
    // Find the maximum value in the array.
    assert values.length == b.length : "values and b must be the same length";

    // Subtract the maximum value from each value in the array, exponential the result, and add up these values.
    double sum = 0.0;
    for (int i = 0; i < values.length; i++) {
      sum += b[i] * exp(values[i] - max);
    }

    // Take the natural logarithm of the sum and add the maximum value back in.
    return max + log(sum);
  }

  /**
   * Turns vector into probability distribution.
   *
   * @param values
   */
  private static void softMax(double[] values){
    double max = stream(values).max().getAsDouble();
    double sum = 0.0;
    for(int i = 0; i < values.length; i++){
      values[i] = exp(values[i]-max);
      sum += values[i];
    }
    for(int i = 0; i < values.length; i++){
      values[i] /= sum;
    }
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

  public double[][] getReducedPotentials() {
    return reducedPotentials;
  }

  public int[] getSnaps() {
    return snaps;
  }

  @Override
  public double[] getBinUncertainties() {
    return mbarUncertainties;
  }

  public double[] getObservationEnsembleAverages() {
    return mbarObservableEnsembleAverages;
  }

  public double[] getObservationEnsembleUncertainties() {
    return mbarObservableEnsembleAverageUncertainties;
  }

  public double[][] getUncertaintyMatrix() {
    return uncertaintyMatrix;
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

  public double[] getBinEntropies() {
    return mbarEntropy;
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

  /**
   * Test all MBAR methods individually with a simple Harmonic Oscillator test case with an
   * excess of samples. "PASS" indicates that the test passed, while "FAIL" followed by the
   * method name indicates that the test failed.
   *
   * Last updated - 06/11/2024
   *
   * @return array of test results
   */
  public static String[] testMBARMethods(){
    // Set up highly converged test case
    double[] O_k = {1, 2, 3, 4};
    double[] K_k = {.5, 1.0, 1.5, 2};
    int[] N_k = {100000, 100000, 100000, 100000};
    double beta = 1.0;
    HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(O_k, K_k, beta);
    String setting = "u_kln";
    Object[] sampleResult = testCase.sample(N_k, setting, (long) 0);
    double[][][] u_kln = (double[][][]) sampleResult[1];
    double[] temps = {1 / Constants.R};
    MultistateBennettAcceptanceRatio mbar = new MultistateBennettAcceptanceRatio(O_k, u_kln, temps, 1.0E-7, MultistateBennettAcceptanceRatio.SeedType.ZEROS);
    MultistateBennettAcceptanceRatio mbarHigherTol = new MultistateBennettAcceptanceRatio(O_k, u_kln, temps, 1.0, MultistateBennettAcceptanceRatio.SeedType.ZEROS);
    String[] results = new String[7];
    // Get required information for all methods
    double[][] reducedPotentials = mbar.getReducedPotentials();
    double[] freeEnergyEstimates = mbar.getMBARFreeEnergies();
    double[] highTolFEEstimates = mbarHigherTol.getMBARFreeEnergies();
    double[] zeros = new double[freeEnergyEstimates.length];
    int[] snapsPerLambda = mbar.getSnaps();

    // getMBARFreeEnergies()
    double[] expectedFEEstimates = new double[]{0.0, 0.3474485596619945, 0.5460865684340613, 0.6866650788765148};
    boolean pass = normDiff(freeEnergyEstimates, expectedFEEstimates) < 1e-5;
    expectedFEEstimates = new double[]{0.0, 0.35798124225733474, 0.44721370511807645, 0.477203739646745};
    pass = normDiff(highTolFEEstimates, expectedFEEstimates) < 1e-5 && pass;
    results[0] = pass ? "PASS" : "FAIL getMBARFreeEnergies()";

    // mbarObjectiveFunction()
    double objectiveFunction = mbarObjectiveFunction(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    pass = !(abs(objectiveFunction - 4786294.2692739945) > 1e-5);
    objectiveFunction = mbarObjectiveFunction(reducedPotentials, snapsPerLambda, highTolFEEstimates);
    pass = !(abs(objectiveFunction - 4787001.700838844) > 1e-5) && pass;
    objectiveFunction = mbarObjectiveFunction(reducedPotentials, snapsPerLambda, zeros);
    pass = !(abs(objectiveFunction - 4792767.352152844) > 1e-5) && pass;
    results[1] = pass ? "PASS" : "FAIL mbarObjectiveFunction()";

    // mbarGradient()
    double[] gradient = mbarGradient(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[] expected = new double[]{6.067113034191607E-4, -8.777718552011038E-4, 8.210768953631487E-4, -5.500246369471995E-4};
    pass = !(normDiff(gradient, expected) > 4e-5);
    gradient = mbarGradient(reducedPotentials, snapsPerLambda, highTolFEEstimates);
    expected = new double[]{1969.705314577408, 5108.841258429764, -1072.9526887468976, -6005.593884267446};
    pass = !(normDiff(gradient, expected) > 4e-5) && pass;
    gradient = mbarGradient(reducedPotentials, snapsPerLambda, zeros);
    expected = new double[]{22797.82037585665, -3273.72282675803, -8859.999065013779, -10664.098484078011};
    pass = !(normDiff(gradient, expected) > 4e-5) && pass;
    results[2] = pass ? "PASS" : "FAIL mbarGradient()";

    pass = true;
    // mbarHessian()
    double[][] hessian = mbarHessian(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[][] expected2d = new double[][]{{47600.586808418964, -29977.008359691405, -12870.425573135915, -4753.1528755909385},
            {-29977.008359691405, 63767.745823769576, -24597.198354108747, -9193.539109971487},
            {-12870.425573135915, -24597.198354108747, 64584.87112481013, -27117.247197561417},
            {-4753.1528755909385, -9193.539109971487, -27117.247197561417, 41063.93918312612}};
    pass = !(normDiff(hessian, expected2d) > 16e-5);
    hessian = mbarHessian(reducedPotentials, snapsPerLambda, highTolFEEstimates);
    expected2d = new double[][]{{49168.30161780381, -31256.519016487477, -12983.708230229113, -4928.074371082683},
            {-31256.519016487477, 66075.94621325849, -25339.462656640117, -9479.964540130917},
            {-12983.708230229113, -25339.462656640117, 64308.30940252403, -25985.13851565483},
            {-4928.074371082683, -9479.964540130917, -25985.13851565483, 40393.1774268678}};
    pass = !(normDiff(hessian, expected2d) > 16e-5) && pass;
    hessian = mbarHessian(reducedPotentials, snapsPerLambda, zeros);
    expected2d = new double[][]{{56125.271437145464, -33495.87894376072, -15738.011263498352, -6891.381229885624},
            {-33495.87894376072, 64613.515110188295, -21970.091845920833, -9147.544320511564},
            {-15738.011263498352, -21970.091845920833, 61407.66256511316, -23699.55945569241},
            {-6891.381229885624, -9147.544320511564, -23699.55945569241, 39738.48500608951}};
    pass = !(normDiff(hessian, expected2d) > 16e-5) && pass;
    results[3] = pass ? "PASS" : "FAIL mbarHessian()";

    pass = true;
    // mbarTheta() --> Checked by diffMatrix
    double[][] theta = mbarTheta(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    double[][] diff = diffMatrixCalculation(theta);
    expected2d = new double[][]{{0.0, 0.001953125, 0.003400485419234404, 0.004858337095247168},
            {0.0020716018980074633, 0.0, 0.002042627017905458, 0.004055968683065466},
            {0.003435363105339426, 0.002042627017905458, 0.0, 0.002560568476977909},
            {0.0048828125, 0.004055968683065466, 0.0025135815773894045, 0.0}};
    pass = !(normDiff(diff, expected2d) > 16e-5);
    results[4] = pass ? "PASS" : "FAIL mbarTheta() or diffMatrixCalculation()" ;

    pass = true;
    // selfConsistentUpdate()
    double[] updatedF_k = mbarSelfConsistentUpdate(reducedPotentials, snapsPerLambda, freeEnergyEstimates);
    expected = new double[]{0.0, 0.3474485745068261, 0.5460865662904055, 0.6866650904438742};
    pass = !(normDiff(updatedF_k, expected) > 1e-5);
    updatedF_k = mbarSelfConsistentUpdate(reducedPotentials, snapsPerLambda, highTolFEEstimates);
    expected = new double[]{0.0, 0.327660608017009, 0.4775067849198251, 0.5586442310038073};
    pass = !(normDiff(updatedF_k, expected) > 1e-5) && pass;
    updatedF_k = mbarSelfConsistentUpdate(reducedPotentials, snapsPerLambda, zeros);
    expected = new double[]{0.0, 0.23865416150488983, 0.29814247007871764, 0.31813582643116334};
    pass = !(normDiff(updatedF_k, expected) > 1e-5) && pass;
    results[5] = pass ? "PASS" : "FAIL mbarSelfConsistentUpdate()";

    pass = true;
    // newton()
    updatedF_k = newton(highTolFEEstimates, reducedPotentials, snapsPerLambda, 1e-7);
    pass = !(normDiff(updatedF_k, freeEnergyEstimates) > 1e-5);
    updatedF_k = newton(zeros, reducedPotentials, snapsPerLambda, 1e-7);
    pass = !(normDiff(updatedF_k, freeEnergyEstimates) > 1e-5) && pass;
    results[6] = pass ? "PASS" : "FAIL newton()";

    return results;
  }

  private static double normDiff(double[] a, double[] b){
    double sum = 0.0;
    for(int i = 0; i < a.length; i++){
      sum += abs(a[i] - b[i]);
    }
    return sum;
  }

  private static double normDiff(double[][] a, double[][] b){
    double sum = 0.0;
    for(int i = 0; i < a.length; i++){
      for(int j = 0; j < a[i].length; j++){
        sum += abs(a[i][j] - b[i][j]);
      }
    }
    return sum;
  }

  /**
   * Example MBAR code usage and comparison with analytic answers for Harmonic Oscillators.
   *
   * @param args
   */
  public static void main(String[] args) {
    // Generate sample data
    double[] equilPositions = {1, 2, 3, 4}; // Equilibrium positions
    double[] springConstants = {.5, 1.0, 1.5, 2}; // Spring constants
    int[] samples = {100000, 100000, 100000, 100000}; // Samples per state
    double beta = 1.0; // 1 / (kB * T) equivalent
    HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(equilPositions, springConstants, beta);
    String setting = "u_kln";
    System.out.print("Generating sample data... ");
    Object[] sampleResult = testCase.sample(samples, setting, (long) 0); // Set seed to fixed value for reproducibility
    System.out.println("done. \n");
    double[] x_n = (double[]) sampleResult[0];
    double[][][] u_kln = (double[][][]) sampleResult[1];
    double[] temps = {1 / Constants.R}; // To be passed into MBAR to cancel out beta within calculation

    // Write file for comparison with pymbar
    // Output to forcefieldx/testing/mbar/data/harmonic_oscillators/mbarFiles/energies_{i}.mbar
    // Get absolute path to root of project
    String rootPath = new File("").getAbsolutePath();
    File outputPath = new File(rootPath + "/testing/mbar/data/harmonic_oscillators/mbarFiles");
    if (!outputPath.exists() && !outputPath.mkdirs()) {
      throw new RuntimeException("Failed to create directory: " + outputPath);
    }

    double[] temperatures = new double[equilPositions.length];
    Arrays.fill(temperatures, temps[0]);
    for (int i = 0; i < u_kln.length; i++) {
      File file = new File(outputPath, "energies_" + i + ".mbar");
      writeFile(u_kln[i], file, temperatures[i]);
    }

    // Create an instance of MultistateBennettAcceptanceRatio
    System.out.print("Creating MBAR instance and .estimateDG(false) with standard tolerance & zeros seeding...");
    //MultistateBennettAcceptanceRatio.VERBOSE = true; // Log Newton/SCI iters and other relevant information
    MultistateBennettAcceptanceRatio mbar = new MultistateBennettAcceptanceRatio(equilPositions, u_kln, temps, 1e-7, SeedType.ZEROS);
    System.out.println("done! \n\n");
    double[] mbarFEEstimates = Arrays.copyOf(mbar.mbarFEEstimates, mbar.mbarFEEstimates.length);
    double[] mbarEnthalpyDiff = Arrays.copyOf(mbar.mbarEnthalpy, mbar.mbarEnthalpy.length);
    double[] mbarEntropyDiff = Arrays.copyOf(mbar.mbarEntropy, mbar.mbarEntropy.length);
    double[] mbarUncertainties = Arrays.copyOf(mbar.mbarUncertainties, mbar.mbarUncertainties.length);
    double[][] mbarDiffMatrix = Arrays.copyOf(mbar.uncertaintyMatrix, mbar.uncertaintyMatrix.length);

    // Analytical free energies and entropies
    double[] analyticalFreeEnergies = testCase.analyticalFreeEnergies();
    double[] error = new double[analyticalFreeEnergies.length];
    for (int i = 0; i < error.length; i++) {
      error[i] = analyticalFreeEnergies[i]-mbarFEEstimates[i];
    }
    double[] temp = testCase.analyticalEntropies(0);
    double[] analyticEntropyDiff = new double[temp.length - 1];
    double[] errorEntropy = new double[temp.length - 1];
    for(int i = 0; i < analyticEntropyDiff.length; i++){
      analyticEntropyDiff[i] = temp[i+1] - temp[i];
      errorEntropy[i] = analyticEntropyDiff[i] -mbarEntropyDiff[i];
    }

    // Compare the calculated free energy differences with the analytical ones
    System.out.println("STANDARD THERMODYNAMIC CALCULATIONS: \n");
    System.out.println("Analytical Free Energies: " + Arrays.toString(analyticalFreeEnergies));
    System.out.println("MBAR Free Energies:       " + Arrays.toString(mbarFEEstimates));
    System.out.println("Free Energy Error:        " + Arrays.toString(error));
    System.out.println();
    System.out.println("MBAR dG:                  " + Arrays.toString(mbar.mbarFEDifferenceEstimates));
    System.out.println("MBAR Uncertainties:       " + Arrays.toString(mbarUncertainties));
    System.out.println("MBAR Enthalpy Changes:    " + Arrays.toString(mbarEnthalpyDiff));
    System.out.println();
    System.out.println("MBAR Entropy Changes:     " + Arrays.toString(mbarEntropyDiff));
    System.out.println("Analytic Entropy Changes: " + Arrays.toString(analyticEntropyDiff));
    System.out.println("Entropy Error:            " + Arrays.toString(errorEntropy));
    System.out.println();
    System.out.println("Uncertainty Diff Matrix: ");
    for (double[] matrix : mbarDiffMatrix) {
      System.out.println(Arrays.toString(matrix));
    }
    System.out.println("\n\n");

    // Observables
    System.out.println("MBAR DERIVED OBSERVABLES: \n");
    mbar.setObservableData(u_kln, true, true);
    double[] mbarObservableEnsembleAverages = Arrays.copyOf(mbar.mbarObservableEnsembleAverages,
            mbar.mbarObservableEnsembleAverages.length);
    double[] mbarObservableEnsembleAverageUncertainties = Arrays.copyOf(mbar.mbarObservableEnsembleAverageUncertainties,
            mbar.mbarObservableEnsembleAverageUncertainties.length);
    System.out.println("Multi-Data Observable Example u_kln:");
    System.out.println("MBAR Observable Ensemble Averages (Potential):              " + Arrays.toString(mbarObservableEnsembleAverages));
    System.out.println("Analytical Observable Ensemble Averages (Potential):        " + Arrays.toString(testCase.analyticalObservable("potential energy")));
    System.out.println("MBAR Observable Ensemble Average Uncertainties (Potential): " + Arrays.toString(mbarObservableEnsembleAverageUncertainties));
    System.out.println();

    // Reads data from xAll[0]
    double[][][] xAll = new double[equilPositions.length][equilPositions.length][x_n.length];
    for(int i = 0; i < xAll[0].length; i++){
      for(int j = 0; j < xAll[0][0].length; j++){
        // Copy data multiple times into same window
        xAll[0][i][j] = x_n[j];
      }
    }
    mbar.setObservableData(xAll, false, true);
    mbarObservableEnsembleAverages = Arrays.copyOf(mbar.mbarObservableEnsembleAverages,
            mbar.mbarObservableEnsembleAverages.length);
    mbarObservableEnsembleAverageUncertainties = Arrays.copyOf(mbar.mbarObservableEnsembleAverageUncertainties,
            mbar.mbarObservableEnsembleAverageUncertainties.length);
    System.out.println("Single-Data Observable Example x_n:");
    System.out.println("MBAR Observable Ensemble Averages (Position):              " + Arrays.toString(mbarObservableEnsembleAverages));
    System.out.println("Analytical Observable Ensemble Averages (Position):        " + Arrays.toString(testCase.analyticalMeans()));
    System.out.println("MBAR Observable Ensemble Average Uncertainties (Position): " + Arrays.toString(mbarObservableEnsembleAverageUncertainties));
    System.out.println();
  }

  /**
   * Harmonic oscillators test case generates data for testing the MBAR implementation
   */
  public static class HarmonicOscillatorsTestCase {
    private final double beta;
    private final double[] equilPositions;
    private final int n_states;
    private final double[] springConstants;
    public HarmonicOscillatorsTestCase(double[] O_k, double[] K_k, double beta) {
      this.beta = beta;
      this.equilPositions = O_k;
      this.n_states = O_k.length;
      this.springConstants = K_k;

      if (this.springConstants.length != this.n_states) {
        throw new IllegalArgumentException("Lengths of K_k and O_k should be equal");
      }
    }

    public double[] analyticalMeans() {
      return equilPositions;
    }

    public double[] analyticalStandardDeviations() {
      double[] deviations = new double[n_states];
      for (int i = 0; i < n_states; i++) {
        deviations[i] = Math.sqrt(1.0 / (beta * springConstants[i]));
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
            result[i] = 1.0 / (beta * springConstants[i]) + Math.pow(equilPositions[i], 2);
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
        fe[i] = -0.5 * Math.log(2 * Math.PI / (beta * springConstants[i]));
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
        double x0 = equilPositions[k];
        double sigma = Math.sqrt(1.0 / (beta * springConstants[k]));

        // Number of snaps
        for (int n = 0; n < N_k[k]; n++) {
          double x = x0 + random.nextGaussian() * sigma;
          x_kn[k][n] = x;
          x_n[index] = x;
          s_n[index] = k;
          // Potential energy evaluations
          for (int l = 0; l < n_states; l++) {
            double u = beta * 0.5 * springConstants[l] * Math.pow(x - equilPositions[l], 2.0);
            u_kln[k][l][n] = u;
            u_kn[l][index] = u;
          }
          index++;
        }
        // Set the rest of the array to NaN
        for(int n = N_k[k]; n < N_max; n++){
          for ( int l =0; l < n_states; l++) {
            u_kln[k][l][n] = Double.NaN;
          }
        }
      }

      // Setting corrections
      if ("u_kn".equals(mode)) {
        return new Object[]{x_n, u_kn, N_k, s_n};
      } else if ("u_kln".equals(mode)) {
        return new Object[]{x_n, u_kln, N_k, s_n, u_kn};
      } else {
        throw new IllegalArgumentException("Unknown mode: " + mode);
      }
    }
  }
}
