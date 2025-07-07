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
package ffx.numerics.estimator;

import ffx.numerics.math.RunningStatistics;
import ffx.numerics.math.SummaryStatistics;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.stream;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Bootstrap Free Energy Estimate.
 */
public class EstimateBootstrapper {

  private static final Logger logger = Logger.getLogger(EstimateBootstrapper.class.getName());
  private static final int DEFAULT_LOG_INTERVAL = 100;

  private final BootstrappableEstimator estimate;
  private final int nWindows;
  private final SummaryStatistics[] freeEnergyDifferenceResults;
  private final SummaryStatistics[] enthalpyResults;

  /**
   * Constructor.
   *
   * @param estimator Estimator to bootstrap.
   */
  public EstimateBootstrapper(BootstrappableEstimator estimator) {
    this.estimate = estimator;
    nWindows = estimate.getNumberOfBins();
    freeEnergyDifferenceResults = new SummaryStatistics[nWindows];
    enthalpyResults = new SummaryStatistics[nWindows];
  }

  /**
   * Gets randomized bootstrap indices; ensures there are at least two distinct indices.
   *
   * @param length Number of random indices to generate in range [0,length)
   * @return Randomized indices.
   */
  public static int[] getBootstrapIndices(int length) {
    return getBootstrapIndices(length, ThreadLocalRandom.current());
  }

  /**
   * Gets randomized bootstrap indices; ensures there are at least two distinct indices.
   *
   * @param length Number of random indices to generate in range [0,length)
   * @param random Source of randomness.
   * @return Randomized indices.
   */
  public static int[] getBootstrapIndices(int length, Random random) {
    return getBootstrapIndices(length, random, Math.min(2, length));
  }

  /**
   * Gets randomized bootstrap indices; ensures there are at least a few distinct indices.
   *
   * @param length      Number of random indices to generate in range [0,length)
   * @param random      Source of randomness.
   * @param minDistinct Minimum number of distinct indices.
   * @return Randomized indices.
   */
  public static int[] getBootstrapIndices(int length, Random random, int minDistinct) {
    // Handle extremely short lengths with special-case handling.
    switch (length) {
      case 0 -> {
        return new int[0];
      }
      case 1 -> {
        return new int[]{0};
      }
      case 2 -> {
        int[] indices = new int[2];
        indices[0] = random.nextBoolean() ? 0 : 1;
        indices[1] = random.nextBoolean() ? 0 : 1;
        return indices;
      }
      // Default: leave switch and handle general case.
    }

    // General case.
    int[] indices = random.ints(length, 0, length).toArray();
    long distinctVal = stream(indices).distinct().count();
    int ctr = 0;
    while (distinctVal <= minDistinct) {
      logger.info(
          format(" Regenerating array (iteration %d): only %d distinct values found for length %d.",
              ++ctr, distinctVal, length));
      indices = random.ints(length, 0, length).toArray();
      distinctVal = stream(indices).distinct().count();
    }
    return indices;
  }

  /**
   * Get bootstrap Enthalpy results for each window.
   *
   * @return Bootstrap Enthalpy results for each window.
   */
  public SummaryStatistics[] getEnthalpyResults() {
    return enthalpyResults;
  }

  /**
   * Perform bootstrap analysis.
   *
   * @param trials Number of trials.
   */
  public void bootstrap(long trials) {
    bootstrap(trials, DEFAULT_LOG_INTERVAL);
  }

  /**
   * Perform bootstrap analysis.
   *
   * @param trials      Number of trials.
   * @param logInterval Interval between logging statements.
   */
  public void bootstrap(long trials, long logInterval) {
    RunningStatistics[] windows = new RunningStatistics[nWindows];
    RunningStatistics[] enthalpyWindows = new RunningStatistics[nWindows];
    for (int i = 0; i < nWindows; i++) {
      windows[i] = new RunningStatistics();
      enthalpyWindows[i] = new RunningStatistics();
    }

    for (long i = 0; i < trials; i++) {
      if ((i + 1) % logInterval == 0) {
        logger.fine(format(" Bootstrap Trial %d", i + 1));
      }

      estimate.estimateDG(true);

      double[] fe = estimate.getFreeEnergyDifferences();
      double[] enthalpy = estimate.getEnthalpyDifferences();
      for (int j = 0; j < nWindows; j++) {
        windows[j].addValue(fe[j]);
        enthalpyWindows[j].addValue(enthalpy[j]);
      }
    }

    for (int i = 0; i < nWindows; i++) {
      freeEnergyDifferenceResults[i] = new SummaryStatistics(windows[i]);
      enthalpyResults[i] = new SummaryStatistics(enthalpyWindows[i]);
    }
  }

  /**
   * Get the total free energy difference estimate from bootstrap analysis.
   *
   * @return The total free energy difference estimate.
   */
  public double getTotalFreeEnergyDifference() {
    return getTotalFreeEnergyDifference(getFreeEnergyDifferences());
  }

  /**
   * Get the total enthalpy estimate from bootstrap analysis.
   *
   * @return The total enthalpy estimate.
   */
  public double getTotalEnthalpyChange() {
    return getTotalEnthalpyChange(getEnthalpyChanges());
  }

  /**
   * Get the total entropy change (-TdS) from bootstrap analysis.
   *
   * @return The total entropy change (-TdS).
   */
  public double getTotalEntropyChange() {
    // dG = dH - TdS
    double dG = getTotalFreeEnergyDifference();
    double dH = getTotalEnthalpyChange();
    return dG - dH;
  }

  /**
   * Get bootstrap free energy estimate for each window.
   * <p>
   * getFreeEnergyDifferences
   *
   * @return Return the bootstrap free energy difference estimate for each window.
   */
  public double[] getFreeEnergyDifferences() {
    return stream(freeEnergyDifferenceResults).mapToDouble(SummaryStatistics::getMean).toArray();
  }

  /**
   * Get bootstrap enthalpy estimate for each window.
   *
   * @return Return the bootstrap enthalpy estimate for each window.
   */
  public double[] getEnthalpyChanges() {
    return stream(enthalpyResults).mapToDouble(SummaryStatistics::getMean).toArray();
  }

  /**
   * Get bootstrap entropy estimate for each window (-TdS).
   *
   * @return Return the bootstrap entropy estimate (-TdS) for each window.
   */
  public double[] getEntropyChanges() {
    double[] dG = getFreeEnergyDifferences();
    double[] dH = getEnthalpyChanges();
    double[] dS = new double[nWindows];
    for (int i = 0; i < nWindows; i++) {
      dS[i] = dG[i] - dH[i];
    }
    return dS;
  }

  /**
   * Get the total free energy difference uncertainty estimate from bootstrap analysis.
   *
   * @return The total free energy difference uncertainty estimate.
   */
  public double getTotalFEDifferenceUncertainty() {
    return getTotalFEDifferenceUncertainty(getFEDifferenceVariances());
  }

  /**
   * Get the total enthalpy uncertainty estimate from bootstrap analysis.
   *
   * @return The total enthalpy uncertainty estimate.
   */
  public double getTotalEnthalpyUncertainty() {
    return getTotalEnthalpyUncertainty(getEnthalpyVariances());
  }

  /**
   * Get the total entropy uncertainty estimate from bootstrap analysis.
   * This is computed as the sqrt of the sum of the free energy and enthalpy variances.
   *
   * @return The total enthalpy uncertainty estimate.
   */
  public double getTotalEntropyUncertainty() {
    double dG = getTotalFEDifferenceUncertainty();
    double dH = getTotalEnthalpyUncertainty();
    return sqrt(dG * dG + dH * dH);
  }

  /**
   * Get the free energy difference uncertainties (standard deviations) from bootstrap analysis for each
   * window.
   *
   * @return Returns free energy difference standard deviation estimate for each window.
   */
  public double[] getFEDifferenceStdDevs() {
    return stream(freeEnergyDifferenceResults).mapToDouble(SummaryStatistics::getSd).toArray();
  }

  /**
   * Get the enthalpy standard deviation estimate from bootstrap analysis for each window.
   *
   * @return Returns enthalpy standard deviation estimate for each window.
   */
  public double[] getEnthalpyStdDevs() {
    return stream(enthalpyResults).mapToDouble(SummaryStatistics::getSd).toArray();
  }

  /**
   * Get the entropy standard deviation estimate from bootstrap analysis for each window.
   * <p>
   * This is computed as the sqrt of the sum of the free energy and enthalpy variances for each window.
   *
   * @return Returns enthalpy standard deviation estimate for each window.
   */
  public double[] getEntropyStdDevs() {
    double[] dG = getFEDifferenceStdDevs();
    double[] dH = getEnthalpyStdDevs();
    double[] dS = new double[nWindows];
    for (int i = 0; i < nWindows; i++) {
      dS[i] = sqrt(dG[i] * dG[i] + dH[i] * dH[i]);
    }
    return dS;
  }

  /**
   * Get the total free energy difference estimate from per window bootstrap analysis.
   *
   * @param freeEnergyDifferences The free energy difference estimate for each window.
   * @return The total free energy difference estimate from bootstrap analysis.
   */
  public double getTotalFreeEnergyDifference(double[] freeEnergyDifferences) {
    return estimate.getTotalFreeEnergyDifference(freeEnergyDifferences);
  }

  /**
   * Get the free energy difference variance estimate from bootstrap analysis for each window.
   *
   * @return Returns free energy difference variance estimate for each window.
   */
  public double[] getFEDifferenceVariances() {
    return stream(freeEnergyDifferenceResults).mapToDouble(SummaryStatistics::getVar).toArray();
  }

  /**
   * Get the total free energy difference uncertainty estimate from per window bootstrap analysis.
   *
   * @param variances The free energy difference variance estimate (not uncertainty) for each window.
   * @return The total free energy difference uncertainty from bootstrap analysis.
   */
  public double getTotalFEDifferenceUncertainty(double[] variances) {
    return estimate.getTotalFEDifferenceUncertainty(variances);
  }

  /**
   * Get the total enthalpy estimate from per window bootstrap analysis.
   *
   * @param enthalpyChanges The enthalpy estimate for each window.
   * @return The total enthalpy difference estimate from bootstrap analysis.
   */
  public double getTotalEnthalpyChange(double[] enthalpyChanges) {
    return estimate.getTotalFreeEnergyDifference(enthalpyChanges);
  }

  /**
   * Get the total enthalpy uncertainty estimate from per window bootstrap analysis.
   *
   * @param variances The enthalpy variance estimate (not uncertainty) for each window.
   * @return The total enthalpy uncertainty estimate from bootstrap analysis.
   */
  public double getTotalEnthalpyUncertainty(double[] variances) {
    return estimate.getTotalEnthalpyUncertainty(variances);
  }

  /**
   * Get the enthalpy variance estimate from bootstrap analysis for each window.
   *
   * @return Returns enthalpy variance estimate for each window.
   */
  public double[] getEnthalpyVariances() {
    return stream(enthalpyResults).mapToDouble(SummaryStatistics::getVar).toArray();
  }
}
