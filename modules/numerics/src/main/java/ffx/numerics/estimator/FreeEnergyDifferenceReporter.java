//******************************************************************************
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
//******************************************************************************
package ffx.numerics.estimator;

import ffx.numerics.math.BootStrapStatistics;

import java.util.logging.Logger;

import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.min;

/**
 * The FreeEnergyDifferenceReporter class
 * reports free energy differences using the Bennett Acceptance Ratio (BAR) method.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class FreeEnergyDifferenceReporter {

  private static final Logger logger = Logger.getLogger(FreeEnergyDifferenceReporter.class.getName());

  /**
   * The number of states in the simulation.
   */
  private final int nStates;

  /**
   * The energiesLow array stores energy values for state L snapshots evaluated at L-dL.
   * The array is of size [nStates][nSnapshots].
   */
  private final double[][] energiesLow;

  /**
   * The energiesAt array stores energy values for state L snapshots evaluated at L.
   * The array is of size [nStates][nSnapshots].
   */
  private final double[][] energiesAt;

  /**
   * The energiesHigh array stores energy values for state L snapshots evaluated at L + dL.
   * The array is of size [nStates][nSnapshots].
   */
  private final double[][] energiesHigh;

  /**
   * The volume array stores volume values for state L at each snapshot.
   * The array is of size [nStates][nSnapshots].
   */
  private final double[][] volume;

  /**
   * The lambda values for each state.
   */
  private final double[] lambdaValues;
  /**
   * The temperature values for each state.
   */
  private final double[] temperature;
  /**
   * The convergence criterion for the BAR iteration.
   */
  private final double eps;
  /**
   * The number of iterations to be used for BAR.
   */
  private final int nIterations;

  /**
   * Maximum number of trials to be used for bootstrap.
   */
  final long MAX_BOOTSTRAP_TRIALS = 100000L;

  /**
   * Free energy difference from forward FEP.
   */
  private double forwardTotalFEDifference;
  /**
   * Enthalpy difference from forward FEP.
   */
  private double forwardTotalEnthalpyChange;
  /**
   * Entropy difference from forward FEP.
   */
  private double forwardTotalEntropyChange;
  /**
   * Free energy difference from backward FEP.
   */
  private double backwardTotalFEDifference;
  /**
   * Enthalpy difference from backward FEP.
   */
  private double backwardTotalEnthalpyChange;
  /**
   * Entropy difference from backward FEP.
   */
  private double backwardTotalEntropyChange;
  /**
   * Free energy difference from BAR.
   */
  private double barIterTotalFEDiff;
  /**
   * Free energy difference from BAR using boot-strapping.
   */
  private double barBSTotalFEDiff;
  /**
   * Enthalpy difference from BAR.
   */
  private double barBSTotalEnthalpyChange;
  /**
   * Entropy difference from BAR.
   */
  private double barBSTotalEntropyChange;

  /**
   * Report Free Energy Differences based on a series of states.
   *
   * @param nStates      The number of states.
   * @param lambdaValues The lambda value for each state.
   * @param temperature  The temperature for each state.
   * @param eps          The BAR convergence criteria.
   * @param nIterations  The BAR maximum number of iterations.
   * @param energiesLow  The energy for each snapshot from state L evaluated at L - dL.
   * @param energiesAt   The energy for each snapshot from state L evaluated at L.
   * @param energiesHigh The energy for each snapshot from state L evaluated at L + dL.
   * @param volume       The volume for each snapshot from state L.
   */
  public FreeEnergyDifferenceReporter(int nStates, double[] lambdaValues, double[] temperature, double eps, int nIterations,
                                      double[][] energiesLow, double[][] energiesAt, double[][] energiesHigh, double[][] volume) {
    this.nStates = nStates;
    this.energiesLow = energiesLow;
    this.energiesAt = energiesAt;
    this.energiesHigh = energiesHigh;
    this.volume = volume;
    this.lambdaValues = lambdaValues;
    this.temperature = temperature;
    this.eps = eps;
    this.nIterations = nIterations;
  }

  /**
   * Report the free energy differences.
   */
  public void report() {
    // Compute the mean and standard deviation of the energy for each state.
    double[] energyMean = new double[nStates];
    double[] energySD = new double[nStates];
    double[] energyVar = new double[nStates];
    for (int state = 0; state < nStates; state++) {
      BootStrapStatistics energyStats = new BootStrapStatistics(energiesAt[state]);
      energyMean[state] = energyStats.mean;
      energySD[state] = energyStats.sd;
      energyVar[state] = energyStats.var;
    }

    // Create the BAR instance.
    BennettAcceptanceRatio bar = new BennettAcceptanceRatio(lambdaValues, energiesLow,
        energiesAt, energiesHigh, temperature, eps, nIterations);

    barIterTotalFEDiff = bar.getTotalFreeEnergyDifference();
    double barIterFEDiffTotalUncertainty = bar.getTotalFEDifferenceUncertainty();
    double[] barIterFEDifferences = bar.getFreeEnergyDifferences();
    double[] barIterFEDiffUncertainties = bar.getFEDifferenceUncertainties();

    EstimateBootstrapper barBS = new EstimateBootstrapper(bar);
    EstimateBootstrapper forwardBS = new EstimateBootstrapper(bar.getInitialForwardsGuess());
    EstimateBootstrapper backwardBS = new EstimateBootstrapper(bar.getInitialBackwardsGuess());

    int numSnapshots = energiesAt[0].length;
    long bootstrap = min(MAX_BOOTSTRAP_TRIALS, numSnapshots);
    logger.info(format(" Number of bootstrap trials: %d", numSnapshots));

    long time = -System.nanoTime();
    forwardBS.bootstrap(bootstrap);
    time += System.nanoTime();
    logger.fine(format(" Forward FEP Bootstrap Complete:      %7.4f sec", time * NS2SEC));
    forwardTotalFEDifference = forwardBS.getTotalFreeEnergyDifference();
    double forwardTotalFEDiffUncertainty = forwardBS.getTotalFEDifferenceUncertainty();
    forwardTotalEnthalpyChange = forwardBS.getTotalEnthalpyChange();
    double forwardTotalEnthalpyUncertainty = forwardBS.getTotalEnthalpyUncertainty();

    time = -System.nanoTime();
    backwardBS.bootstrap(bootstrap);
    time += System.nanoTime();
    logger.fine(format(" Backward FEP Bootstrap Complete:     %7.4f sec", time * NS2SEC));
    backwardTotalFEDifference = backwardBS.getTotalFreeEnergyDifference();
    double backwardTotalFEDiffUncertainty = backwardBS.getTotalFEDifferenceUncertainty();
    backwardTotalEnthalpyChange = backwardBS.getTotalEnthalpyChange();
    double backwardTotalEnthalpyUncertainty = backwardBS.getTotalEnthalpyUncertainty();

    time = -System.nanoTime();
    barBS.bootstrap(bootstrap);
    time += System.nanoTime();
    logger.fine(format(" BAR Bootstrap Complete:              %7.4f sec", time * NS2SEC));

    barBSTotalFEDiff = barBS.getTotalFreeEnergyDifference();
    double barBSTotalFEDiffUncertainty = barBS.getTotalFEDifferenceUncertainty();
    barBSTotalEnthalpyChange = barBS.getTotalEnthalpyChange();
    double barBSTotalEnthalpyUncertainty = barBS.getTotalEnthalpyUncertainty();

    if (nStates > 2) {
      logger.info("\n Window Free Energy Differences (kcal/mol)");
      logger.info("        Forward FEP           Backward FEP          BAR Iteration         BAR Bootstrap");
      double[] forwardFEDifferences = forwardBS.getFreeEnergyDifferences();
      double[] forwardFEDUncertainty = forwardBS.getFEDifferenceStdDevs();
      double[] backwardFEDifferences = backwardBS.getFreeEnergyDifferences();
      double[] backwardFEDUncertainty = backwardBS.getFEDifferenceStdDevs();
      double[] barBSFE = barBS.getFreeEnergyDifferences();
      double[] barBSUncertaintyFE = barBS.getFEDifferenceStdDevs();
      for (int n = 0; n < nStates - 1; n++) {
        logger.info(format(" %2d %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f",
            n, forwardFEDifferences[n], forwardFEDUncertainty[n],
            backwardFEDifferences[n], backwardFEDUncertainty[n],
            barIterFEDifferences[n], barIterFEDiffUncertainties[n],
            barBSFE[n], barBSUncertaintyFE[n]));
      }
    }

    logger.info("\n Total Free Energy Difference (kcal/mol)");
    logger.info("        Forward FEP           Backward FEP          BAR Iteration         BAR Bootstrap");
    logger.info(format("    %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f",
        forwardTotalFEDifference, forwardTotalFEDiffUncertainty,
        backwardTotalFEDifference, backwardTotalFEDiffUncertainty,
        barIterTotalFEDiff, barIterFEDiffTotalUncertainty,
        barBSTotalFEDiff, barBSTotalFEDiffUncertainty));

    logger.info("\n Enthalpy from Potential Energy Averages (kcal/mol)\n");
    for (int n = 0; n < nStates; n++) {
      logger.info(format(" State %2d:     %12.4f +/- %6.4f", n, energyMean[n], energySD[n]));
    }
    double enthalpyDiff = energyMean[nStates - 1] - energyMean[0];
    double enthalpyDiffSD = Math.sqrt(energyVar[nStates - 1] + energyVar[0]);
    logger.info(format(" Enthalpy via Direct Estimate:   %12.4f +/- %6.4f", enthalpyDiff, enthalpyDiffSD));

    // Enthalpy Differences
    if (nStates > 2) {
      logger.info("\n Window Enthalpy Differences (kcal/mol)");
      logger.info("        Forward FEP           Backward FEP          BAR Bootstrap");
      double[] forwardEnthalpyDifferences = forwardBS.getEnthalpyChanges();
      double[] forwardEnthalpyUncertainty = forwardBS.getEnthalpyStdDevs();
      double[] backwardEnthalpyDifferences = backwardBS.getEnthalpyChanges();
      double[] backwardEnthalpyUncertainty = backwardBS.getEnthalpyStdDevs();
      double[] barBSenthalpy = barBS.getEnthalpyChanges();
      double[] barBSenthalpyUncertainty = barBS.getEnthalpyStdDevs();
      for (int n = 0; n < nStates - 1; n++) {
        logger.info(format(" %2d %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f",
            n, forwardEnthalpyDifferences[n], forwardEnthalpyUncertainty[n],
            backwardEnthalpyDifferences[n], backwardEnthalpyUncertainty[n],
            barBSenthalpy[n], barBSenthalpyUncertainty[n]));
      }
    }
    logger.info("\n Total Enthalpy Difference (kcal/mol)");
    logger.info("        Forward FEP           Backward FEP          BAR Bootstrap");
    logger.info(format("    %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f",
        forwardTotalEnthalpyChange, forwardTotalEnthalpyUncertainty,
        backwardTotalEnthalpyChange, backwardTotalEnthalpyUncertainty,
        barBSTotalEnthalpyChange, barBSTotalEnthalpyUncertainty));

    // Entropy Differences
    if (nStates > 2) {
      logger.info("\n Window Entropy Differences -T*dS (kcal/mol)");
      logger.info("        Forward FEP           Backward FEP          BAR Bootstrap");
      double[] forwardEntropyDifferences = forwardBS.getEntropyChanges();
      double[] forwardEntropyUncertainty = forwardBS.getEntropyStdDevs();
      double[] backwardEntropyDifferences = backwardBS.getEntropyChanges();
      double[] backwardEntropyUncertainty = backwardBS.getEntropyStdDevs();
      double[] barBSEntropy = barBS.getEntropyChanges();
      double[] barBSEntropyUncertainty = barBS.getEntropyStdDevs();
      for (int n = 0; n < nStates - 1; n++) {
        logger.info(format(" %2d %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f",
            n, forwardEntropyDifferences[n], forwardEntropyUncertainty[n],
            backwardEntropyDifferences[n], backwardEntropyUncertainty[n],
            barBSEntropy[n], barBSEntropyUncertainty[n]));
        ;
      }
    }
    forwardTotalEntropyChange = forwardBS.getTotalEntropyChange();
    double forwardTotalEntropyUncertainty = forwardBS.getTotalEntropyUncertainty();
    backwardTotalEntropyChange = backwardBS.getTotalEntropyChange();
    double backwardTotalEntropyUncertainty = backwardBS.getTotalEntropyUncertainty();
    barBSTotalEntropyChange = barBS.getTotalEntropyChange();
    double barBSTotalEntropyUncertainty = barBS.getTotalEntropyUncertainty();
    logger.info("\n Total Entropy Difference -T*dS (kcal/mol)");
    logger.info("        Forward FEP           Backward FEP          BAR Bootstrap");
    logger.info(format("    %10.4f +/- %6.4f %10.4f +/- %6.4f %10.4f +/- %6.4f",
        forwardTotalEntropyChange, forwardTotalEntropyUncertainty,
        backwardTotalEntropyChange, backwardTotalEntropyUncertainty,
        barBSTotalEntropyChange, barBSTotalEntropyUncertainty));
  }

  /**
   * Free energy difference from forward FEP.
   *
   * @return the forward total free energy difference.
   */
  public double getForwardTotalFEDifference() {
    return forwardTotalFEDifference;
  }

  /**
   * Enthalpy difference from forward FEP.
   *
   * @return the forward total enthalpy change.
   */
  public double getForwardTotalEnthalpyChange() {
    return forwardTotalEnthalpyChange;
  }

  /**
   * Entropy difference from forward FEP.
   *
   * @return the forward total entropy change.
   */
  public double getForwardTotalEntropyChange() {
    return forwardTotalEntropyChange;
  }

  /**
   * Free energy difference from backward FEP.
   *
   * @return the backward total free energy difference.
   */
  public double getBackwardTotalFEDifference() {
    return backwardTotalFEDifference;
  }

  /**
   * Enthalpy difference from backward FEP.
   *
   * @return the backward total enthalpy change.
   */
  public double getBackwardTotalEnthalpyChange() {
    return backwardTotalEnthalpyChange;
  }

  /**
   * Entropy difference from backward FEP.
   *
   * @return the backward total entropy change.
   */
  public double getBackwardTotalEntropyChange() {
    return backwardTotalEntropyChange;
  }

  /**
   * Free energy difference from BAR.
   *
   * @return the BAR iteration total free energy difference.
   */
  public double getBarIterTotalFEDiff() {
    return barIterTotalFEDiff;
  }

  /**
   * Free energy difference from BAR using boot-strapping.
   *
   * @return the BAR boot-strap total free energy difference.
   */
  public double getBarBSTotalFEDiff() {
    return barBSTotalFEDiff;
  }

  /**
   * Enthalpy difference from BAR.
   *
   * @return the BAR boot-strap total enthalpy change.
   */
  public double getBarBSTotalEnthalpyChange() {
    return barBSTotalEnthalpyChange;
  }

  /**
   * Entropy difference from BAR.
   *
   * @return the BAR boot-strap total entropy change.
   */
  public double getBarBSTotalEntropyChange() {
    return barBSTotalEntropyChange;
  }
}
