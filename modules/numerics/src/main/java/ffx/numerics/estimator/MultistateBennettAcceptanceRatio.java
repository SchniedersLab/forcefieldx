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
import ffx.numerics.integrate.DataSet;
import ffx.numerics.integrate.DoublesDataSet;
import ffx.numerics.integrate.Integrate1DNumeric;
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
    if (MultistateBennettAcceptanceRatio.VERBOSE){
      logger.setLevel(java.util.logging.Level.FINE);
    }

    // Find repeated snapshots if from continuous lambda
    int numEvaluations = eAllFlat[0].length;
    int firstTrajLength = snaps[0];
    if(stream(snaps).allMatch(s -> s == firstTrajLength)) {
      double[][] eAllFlatTemp = new double[nLambdaStates][snaps[0]];
      for (int i = 0; i < nLambdaStates; i++) {
        for (int j = 0; j < snaps[0]; j++) {
          eAllFlatTemp[i][j] = eAllFlat[i][j];
        }
      }
      // Check if there are repeated snapshots
      int repeatedCount = 0;
      for (int i = 0; i < nLambdaStates; i++) {
        for (int j = snaps[0]; j < numEvaluations; j++) {
          if (eAllFlatTemp[i][j % (snaps[0])] - eAllFlat[i][j] < 1.0E-6) {
            repeatedCount++;
          }
        }
      }
      int expectedRepeats = eAllFlat.length * eAllFlat[0].length - eAllFlatTemp.length * eAllFlatTemp[0].length;
      if (repeatedCount == expectedRepeats) {
        logger.warning(" Repeated snapshots detected. MBAR may not converge.");
        eAllFlat = eAllFlatTemp;
        int reduction = numEvaluations / eAllFlat[0].length;
        numEvaluations = eAllFlat[0].length;
        for (int i = 0; i < snaps.length; i++) {
          snaps[i] /= reduction;
        }
      }
    }

    applyBiasCorrection();

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

    // Precompute beta for each state.
    rtValues = new double[nLambdaStates];
    double[] invRTValues = new double[nLambdaStates];
    for (int i = 0; i < nLambdaStates; i++) {
      rtValues[i] = Constants.R * temperatures[i];
      invRTValues[i] = 1.0 / rtValues[i];
    }

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
        mbarFEEstimatesTemp = newton(mbarFEEstimatesTemp, reducedPotentialsTemp, snapsTemp, tolerance);
      }
    } catch (Exception e) {
      logger.warning(" L-BFGS/Newton failed to converge. Finishing w/ self-consistent iteration. Message: " +
              e.getMessage());
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
    do {
      prevMBAR = copyOf(mbarFEEstimates, nLambdaStates);
      mbarFEEstimates = selfConsistentUpdate(reducedPotentials, snaps, mbarFEEstimates);
      for (int i = 0; i < nLambdaStates; i++) { // SOR for acceleration
        mbarFEEstimates[i] = omega * mbarFEEstimates[i] + (1 - omega) * prevMBAR[i];
      }
      if (stream(mbarFEEstimates).anyMatch(Double::isInfinite) || stream(mbarFEEstimates).anyMatch(Double::isNaN)) {
        throw new IllegalArgumentException("MBAR estimate contains NaNs or Infs after iteration " + sciIter);
      }
      sciIter++;
    } while (!converged(prevMBAR) && sciIter < 1000);
    if (MultistateBennettAcceptanceRatio.VERBOSE) {
      logger.info(" SCI iterations: " + sciIter);
    }

    // Calculate uncertainties
    double[][] theta = mbarTheta(reducedPotentials, snaps, mbarFEEstimates); // Quite expensive
    mbarUncertainties = mbarUncertaintyCalc(theta);
    totalMBARUncertainty = mbarTotalUncertaintyCalc(theta);
    diffMatrix = diffMatrixCalculation(theta);
    if (!randomSamples && MultistateBennettAcceptanceRatio.VERBOSE) { // Don't log for bootstrapping
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

  private void applyBiasCorrection() {
    if(biasFlat != null) {
      for (int i = 0; i < eAllFlat.length; i++) {
        for (int j = 0; j < eAllFlat[0].length; j++) {
          eAllFlat[i][j] += biasFlat[i][j];
        }
      }
    }
  }

  //////// Misc. Methods ////////////

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
   * Gradient of the MBAR objective function. This is used for L-BFGS & Newton optimization.
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

  private double[] mbarEntropyCalc(double[] mbarEnthalpy, double[] mbarFEEstimates) {
    double[] entropy = new double[mbarFEEstimates.length - 1];
    for(int i = 0; i < entropy.length; i++) {
      entropy[i] = mbarEnthalpy[i] - mbarFEDifferenceEstimates[i]; // dG = dH - TdS | TdS = dH - dG
    }
    return entropy;
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

  private void fillObservationExpectations(boolean multiData){
    if(multiData){
      mbarObservableEnsembleAverages = new double[oAllFlat.length];
      mbarObservableEnsembleAverageUncertainties = new double[oAllFlat.length];
      for(int i = 0; i < oAllFlat.length; i++){
        mbarObservableEnsembleAverages[i] = computeExpectations(oAllFlat[i])[i];
        mbarObservableEnsembleAverageUncertainties[i] = computeExpectationStd(oAllFlat[i])[i];
      }
    } else {
      mbarObservableEnsembleAverages = computeExpectations(oAllFlat[0]);
      mbarObservableEnsembleAverageUncertainties = computeExpectationStd(oAllFlat[0]);
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
    // Subtract min sample value --> pymbar does this and says there's not a diff (but I think it helps)
    double minSample = stream(samples).min().getAsDouble();
    samples = stream(samples).map(d -> d - minSample).toArray();
    double[] expectations = computeExpectations(samples);
    samples = stream(samples).map(d -> d + minSample).toArray(); // Don't alter values in samples
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
   * Theta = W.T @ (I - W @ diag(snapsPerState) @ W.T)^-1 @ W.
   * <p>
   * Requires calculation and inversion of W matrix.
   * D4 from supp info of MBAR paper used instead to reduce complexity to K^3.
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
   * @param theta matrix of covariances
   * @return Total uncertainty for the MBAR free energy estimates.
   */
  private static double mbarTotalUncertaintyCalc(double[][] theta) {
    int nStates = theta.length;
    return sqrt(abs(theta[0][0] - 2 * theta[0][nStates - 1] + theta[nStates - 1][nStates - 1]));
  }

  /**
   * MBAR uncertainty diff Matrix calculation.
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
   * Newton-Raphson step for MBAR optimization. Falls back to the steepest descent if hessian is singular.
   *
   * The matrix can come back from being singular after iterations, so it isn't worth moving to L-BFGS.
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
   * The logSumExp operation itself prevents causing 0 values from appearing due to large denominators.
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

  double[] getBinEntropies() {
    return mbarEntropy;
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

  public void setObservableData(double[][][] oAll, boolean multiDataObservable) {
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
            oAllFlat[0][count] = oAll[i][i][j];
            count++;
          }
        }
      }
    }
    if (biasFlat != null) {
      for(int i = 0; i< oAllFlat.length; i++) {
        for (int j = 0; j < oAllFlat[i].length; j++) {
          oAllFlat[i][j] *= exp(-2*biasFlat[i][j]);
        }
      }
    }
    this.fillObservationExpectations(multiDataObservable);
  }

  public void setBiasData(double[][][] biasAll, boolean multiDataObservable) {
    biasFlat = new double[biasAll.length][biasAll.length * biasAll[0][0].length];
    if(multiDataObservable){ // Flatten data
      int[] snapsT = new int[biasAll.length];
      int[] nanCount = new int[biasAll.length];
      for (int i = 0; i < biasAll.length; i++) {
        ArrayList<Double> temp = new ArrayList<>();
        for(int j = 0; j < biasAll.length; j++) {
          int count = 0;
          int countNaN = 0;
          for(int k = 0; k < biasAll[j][i].length; k++) {
            // Don't include NaN values
            if (!Double.isNaN(biasAll[j][i][k])) {
              temp.add(biasAll[j][i][k]);
              count++;
            } else {
              countNaN++;
            }
          }
          snapsT[j] = count;
          nanCount[j] = countNaN;
        }
        biasFlat[i] = temp.stream().mapToDouble(Double::doubleValue).toArray();
      }
    } else { // Put relevant data into the 0th index
      int count = 0;
      for (int i = 0; i < biasAll.length; i++){
        for(int j = 0; j < biasAll[0][0].length; j++){
          if(!Double.isNaN(biasAll[i][i][j])){
            biasFlat[0][count] = biasAll[i][i][j];
            count++;
          }
        }
      }
    }
  }

  public static void main(String[] args) {
    double[] O_k = {0, .1, .7, 3, 4}; // Equilibrium positions
    double[] K_k = {1, 2, 3, 5, 6}; // Spring constants
    int[] N_k = {10000, 10000, 10000, 10000, 10000};
    double beta = 1.0;

    // Create an instance of HarmonicOscillatorsTestCase
    HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(O_k, K_k, beta);

    // Generate sample data
    String setting = "u_kln";
    System.out.print("Generating sample data... ");
    Object[] sampleResult = testCase.sample(N_k, setting, System.currentTimeMillis()); // Set seed to fixed value for reproducibility
    System.out.println("done. \n");
    double[][] u_n = ((double[][]) sampleResult[4]);
    double[] x_n = (double[]) sampleResult[0];
    double[][][] u_kln = (double[][][]) sampleResult[1];
    double[] temps = {1 / Constants.R};

    // Write file for comparison with pymbar
    // Output to forcefieldx/testing/mbar/data/harmonic_oscillators/mbarFiles/energies_{i}.mbar
    // Get absolute path to root of project

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
    }

    // Create an instance of MultistateBennettAcceptanceRatio
    System.out.print("Creating MBAR instance and estimateDG() with standard tol & Zeros seeding...");
    File mbarParentFile = new File("/Users/matthewsperanza/Programs/forcefieldx/testing/mbar/hxacan/mbarBiasOST");
    MBARFilter mbarFilter = new MBARFilter(mbarParentFile);
    MultistateBennettAcceptanceRatio.VERBOSE = true;
    MultistateBennettAcceptanceRatio mbar = mbarFilter.getMBAR(SeedType.ZEROS, 1e-7);
    mbarFilter.readObservableData(true, true, false);
    mbar.estimateDG(); // Second run
    mbarFilter.readObservableData(true, false, true);
    double[] mbarObservableEnsembleAverages = Arrays.copyOf(mbar.mbarObservableEnsembleAverages,
            mbar.mbarObservableEnsembleAverages.length);
    double[] mbarObservableEnsembleAverageUncertainties = Arrays.copyOf(mbar.mbarObservableEnsembleAverageUncertainties,
            mbar.mbarObservableEnsembleAverageUncertainties.length);
    mbarFilter.readObservableData(false, false, true);
    double[] mbarObservableEnsembleAveragesSingle = Arrays.copyOf(mbar.mbarObservableEnsembleAverages,
            mbar.mbarObservableEnsembleAverages.length);
    double[] mbarObservableEnsembleAverageUncertaintiesSingle = Arrays.copyOf(mbar.mbarObservableEnsembleAverageUncertainties,
            mbar.mbarObservableEnsembleAverageUncertainties.length);

    //MultistateBennettAcceptanceRatio mbar = new MultistateBennettAcceptanceRatio(O_k, u_kln, temps, 1e-7, SeedType.ZEROS);
    //mbar.setObservableData(u_kln, true);

    double[] lambdas = new double[19];
    for(int i = 0; i < lambdas.length; i++){
      lambdas[i] = i/18.0;
    }
    double[] mbarFEEstimates = Arrays.copyOf(mbar.mbarFEEstimates, mbar.mbarFEEstimates.length);
    double[] mbarEnthalpy = Arrays.copyOf(mbar.mbarEnthalpy, mbar.mbarEnthalpy.length);
    double[] mbarEntropy = Arrays.copyOf(mbar.mbarEntropy, mbar.mbarEntropy.length);
    double[] mbarUncertainties = Arrays.copyOf(mbar.mbarUncertainties, mbar.mbarUncertainties.length);
    double[][] mbarDiffMatrix = Arrays.copyOf(mbar.diffMatrix, mbar.diffMatrix.length);

    EstimateBootstrapper bootstrapper = new EstimateBootstrapper(mbar);
    bootstrapper.bootstrap(0);
    System.out.println("done! \n");

    System.out.println("Lambdas:                           " + Arrays.toString(lambdas));
    System.out.println("MBAR Observable Ensemble Averages: " + Arrays.toString(mbarObservableEnsembleAverages));
    DataSet dSet = new DoublesDataSet(Integrate1DNumeric.generateXPoints(0,1, mbarObservableEnsembleAverages.length, false),
            mbarObservableEnsembleAverages, false);
    double integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.TRAPEZOIDAL);
    System.out.println("Integral Simp: " + integral);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.BOOLE);
    System.out.println("Integral Boole: " + integral);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.SIMPSONS);
    System.out.println("Integral Trap: " + integral);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.RECTANGULAR);
    System.out.println("Integral Rect: " + integral);
    System.out.println("Total FE: " + mbar.totalMBAREstimate);
    System.out.println("MBAR Ob  Ensemble Averages Single: " + Arrays.toString(mbarObservableEnsembleAveragesSingle));
    dSet = new DoublesDataSet(Integrate1DNumeric.generateXPoints(0,1, mbarObservableEnsembleAveragesSingle.length, false),
            mbarObservableEnsembleAveragesSingle, false);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.SIMPSONS);
    System.out.println("Integral Simp: " + integral);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.BOOLE);
    System.out.println("Integral Boole: " + integral);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.SIMPSONS);
    System.out.println("Integral Trap: " + integral);
    integral = Integrate1DNumeric.integrateData(dSet, Integrate1DNumeric.IntegrationSide.LEFT, Integrate1DNumeric.IntegrationType.RECTANGULAR);
    System.out.println("Integral Rect: " + integral);
    System.out.println("Total FE: " + mbar.totalMBAREstimate);
    System.out.println("MBAR Observable Ensemble Average Uncertainties: " + Arrays.toString(mbarObservableEnsembleAverageUncertainties));
    System.out.println("MBAR Ob  Ensemble Average Uncertainties Single: " + Arrays.toString(mbarObservableEnsembleAverageUncertaintiesSingle));
    System.out.println();

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
    System.out.println("MBAR dG:                  " + Arrays.toString(mbar.mbarFEDifferenceEstimates));
    System.out.println("MBAR Enthalpy Changes:    " + Arrays.toString(mbarEnthalpy));
    System.out.println("MBAR Entropy Changes:     " + Arrays.toString(mbarEntropy));
    double[] temp = testCase.analyticalEntropies(0);
    double[] temp2 = new double[temp.length - 1];
    for(int i = 0; i < temp2.length; i++){
      temp2[i] = temp[i+1] - temp[i];
    }
    System.out.println("Analytic Entropy Changes: " + Arrays.toString(temp2));
    System.out.println();
    System.out.println("Uncertainty Diff Matrix: ");
    for (double[] matrix : mbarDiffMatrix) {
      System.out.println(Arrays.toString(matrix));
    }
    System.out.println("\n");

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
