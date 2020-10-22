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
package ffx.numerics.estimator;

import static java.lang.String.format;
import static java.util.Arrays.stream;

import ffx.numerics.math.RunningStatistics;
import ffx.numerics.math.SummaryStatistics;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

/**
 * Bootstrap Free Energy Estimate.
 */
public class EstimateBootstrapper {

  private static final Logger logger = Logger.getLogger(EstimateBootstrapper.class.getName());
  private static final long DEFAULT_LOG_INTERVAL = 25;

  private final BootstrappableEstimator estimate;
  private final int nWindows;
  private final SummaryStatistics[] bootstrapFEResults;
  private final SummaryStatistics[] bootstrapEnthalpyResults;

  public EstimateBootstrapper(BootstrappableEstimator estimator) {
    this.estimate = estimator;
    nWindows = estimate.numberOfBins();
    bootstrapFEResults = new SummaryStatistics[nWindows];
    bootstrapEnthalpyResults = new SummaryStatistics[nWindows];
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
   * @param length Number of random indices to generate in range [0,length)
   * @param random Source of randomness.
   * @param minDistinct Minimum number of distinct indices.
   * @return Randomized indices.
   */
  public static int[] getBootstrapIndices(int length, Random random, int minDistinct) {
    // Handle extremely short lengths with special-case handling.
    switch (length) {
      case 0:
        return new int[0];
      case 1:
        return new int[] {0};
      case 2:
        int[] indices = new int[2];
        indices[0] = random.nextBoolean() ? 0 : 1;
        indices[1] = random.nextBoolean() ? 0 : 1;
        return indices;
      // Default: leave switch and handle general case.
    }

    // General case.
    int[] indices = random.ints(length, 0, length).toArray();
    long distinctVals = stream(indices).distinct().count();
    int ctr = 0;
    while (distinctVals <= minDistinct) {
      logger.info(
          format(" Regenerating array (iteration %d): only %d distinct values found for length %d.",
              ++ctr, distinctVals, length));
      indices = random.ints(length, 0, length).toArray();
      distinctVals = stream(indices).distinct().count();
    }
    return indices;
  }

  public SummaryStatistics[] getBootstrapEnthalpyResults() {
    return bootstrapEnthalpyResults;
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
   * @param trials Number of trials.
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

      double[] fe = estimate.getBinEnergies();
      double[] enthalpy = estimate.getBinEthalpies();
      for (int j = 0; j < nWindows; j++) {
        windows[j].addValue(fe[j]);
        enthalpyWindows[j].addValue(enthalpy[j]);
      }
    }

    for (int i = 0; i < nWindows; i++) {
      bootstrapFEResults[i] = new SummaryStatistics(windows[i]);
      bootstrapEnthalpyResults[i] = new SummaryStatistics(enthalpyWindows[i]);
    }
  }

  /**
   * Get bootstrap free energy estimate for each window.
   *
   * @return Return the bootstrap free energy difference estimate for each windows.
   */
  public double[] getFE() {
    return stream(bootstrapFEResults).mapToDouble(SummaryStatistics::getMean).toArray();
  }
  /**
   * Get bootstrap enthalpy estimate for each window.
   *
   * @return Return the bootstrap enthalpy estimate for each windows.
   */
  public double[] getEnthalpy() {
    return stream(bootstrapEnthalpyResults).mapToDouble(SummaryStatistics::getMean).toArray();
  }
  /**
   * Get the total free energy difference estimate from bootstrap analysis.
   *
   * @return The total free energy difference estimate.
   */
  public double getTotalFE() {
    return getTotalFE(getFE());
  }
  /**
   * Get the total enthalpy estimate from bootstrap analysis.
   *
   * @return The total enthalpy estimate.
   */

  public double getTotalEnthalpy() {
    return getTotalEnthalpy(getEnthalpy());
  }


  /**
   * Get the total free energy difference estimate from per window bootstrap analysis.
   *
   * @param fe The free energy difference estimate for each window.
   * @return The total free energy difference estimate from bootstrap analysis.
   */
  public double getTotalFE(double[] fe) {
    return estimate.sumBootstrapResults(fe);
  }
  /**
   * Get the total enthalpy estimate from per window bootstrap analysis.
   *
   * @param enthalpy The enthalpy estimate for each window.
   * @return The total enthalpy difference estimate from bootstrap analysis.
   */
  public double getTotalEnthalpy(double[] enthalpy) {
    return estimate.sumBootstrapResults(enthalpy);
  }
  /**
   * Get the total free energy difference variance estimate from bootstrap analysis.
   *
   * @return The total free energy difference variance estimate.
   */
  public double getTotalUncertainty() {
    return getTotalUncertainty(getVariance());
  }

  /**
   * Get the total free energy difference estimate from per window bootstrap analysis.
   *
   * @param var The free energy difference variance estimate (not uncertainty) for each window.
   * @return The total free energy difference variance estimate from bootstrap analysis.
   */
  public double getTotalUncertainty(double[] var) {
    return estimate.sumBootstrapUncertainty(var);
  }
  /**
   * Get the total enthalpy variance estimate from bootstrap analysis.
   *
   * @return The total enthalpy variance estimate.
   */
  public double getTotalEnthalpyUncertainty() {
    return getTotalEnthalpyUncertainty(getEnthalpyVariance());
  }
  /**
   * Get the total enthalpy estimate from per window bootstrap analysis.
   *
   * @param var The enthalpy variance estimate (not uncertainty) for each window.
   * @return The total enthalpy variance estimate from bootstrap analysis.
   */
  public double getTotalEnthalpyUncertainty(double[] var) {
    return estimate.sumBootstrapEnthalpyUncertainty(var);
  }

  /**
   * Get the free energy difference standard deviation estimate from bootstrap analysis for each
   * window.
   *
   * @return Returns free energy difference standard deviation estimate for each window.
   */
  public double[] getUncertainty() {
    return stream(bootstrapFEResults).mapToDouble(SummaryStatistics::getSd).toArray();
  }
  /**
   * Get the enthalpy standard deviation estimate from bootstrap analysis for each
   * window.
   *
   * @return Returns enthalpy standard deviation estimate for each window.
   */
  public double[] getEnthalpyUncertainty() {
    return stream(bootstrapEnthalpyResults).mapToDouble(SummaryStatistics::getSd).toArray();
  }

  /**
   * Get the free energy difference variance estimate from bootstrap analysis for each window.
   *
   * @return Returns free energy difference variance estimate for each window.
   */
  public double[] getVariance() {
    return stream(bootstrapFEResults).mapToDouble(SummaryStatistics::getVar).toArray();
  }
  /**
   * Get the enthalpy variance estimate from bootstrap analysis for each window.
   *
   * @return Returns enthalpy variance estimate for each window.
   */
  public double[] getEnthalpyVariance() {
    return stream(bootstrapEnthalpyResults).mapToDouble(SummaryStatistics::getVar).toArray();
  }
}
