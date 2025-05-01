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

import ffx.numerics.math.SummaryStatistics;

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.estimator.EstimateBootstrapper.getBootstrapIndices;
import static ffx.numerics.estimator.Zwanzig.Directionality.BACKWARDS;
import static ffx.numerics.estimator.Zwanzig.Directionality.FORWARDS;
import static ffx.numerics.math.ScalarMath.fermiFunction;
import static ffx.utilities.Constants.R;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static java.util.Arrays.stream;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The Bennett Acceptance Ratio class implements the Bennett Acceptance Ratio (BAR) statistical
 * estimator, based on the Tinker implementation.
 *
 * <p>Literature References (from Tinker): C. H. Bennett, "Efficient Estimation of Free Energy
 * Differences from Monte Carlo Data", Journal of Computational Physics, 22, 245-268 (1976)
 *
 * <p>M. A. Wyczalkowski, A. Vitalis and R. V. Pappu, "New Estimators for Calculating Solvation
 * Entropy and Enthalpy and Comparative Assessments of Their Accuracy and Precision, Journal of
 * Physical Chemistry, 114, 8166-8180 (2010) [modified BAR algorithm, non-implemented
 * entropy/enthalpy]
 *
 * <p>K. B. Daly, J. B. Benziger, P. G. Debenedetti and A. Z. Panagiotopoulos, "Massively Parallel
 * Chemical Potential Calculation on Graphics Processing Units", Computer Physics Communications,
 * 183, 2054-2062 (2012) [non-implemented NPT modification]
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class BennettAcceptanceRatio extends SequentialEstimator implements BootstrappableEstimator {

  private static final Logger logger = Logger.getLogger(BennettAcceptanceRatio.class.getName());

  /**
   * Default BAR convergence tolerance.
   */
  private static final double DEFAULT_TOLERANCE = 1.0E-7;
  /**
   * Default maximum number of BAR iterations.
   */
  private static final int DEFAULT_MAX_BAR_ITERATIONS = 100;
  /**
   * Number of state pairs.
   */
  private final int nWindows;
  /**
   * BAR convergence tolerance.
   */
  private final double tolerance;
  /**
   * BAR maximum number of iterations.
   */
  private final int nIterations;
  /**
   * Forward Zwanzig instance.
   */
  private final Zwanzig forwardsFEP;
  /**
   * Backward Zwanzig instance.
   */
  private final Zwanzig backwardsFEP;
  /**
   * Random number generator for bootstrapping.
   */
  private final Random random;
  /**
   * Total BAR free-energy difference estimate.
   */
  private double totalFreeEnergyDifference;
  /**
   * Total BAR free-energy difference uncertainty.
   */
  private double totalFEDifferenceUncertainty;
  /**
   * BAR free-energy difference estimates.
   */
  private final double[] freeEnergyDifferences;
  /**
   * BAR free-energy difference uncertainties.
   */
  private final double[] freeEnergyDifferenceUncertainties;
  /**
   * BAR Enthalpy estimates
   */
  private final double[] enthalpyDifferences;
  /**
   * Forward Zwanzig free-energy difference estimates.
   */
  private final double[] forwardZwanzigFEDifferences;
  /**
   * Backward Zwanzig free-energy difference estimates.
   */
  private final double[] backwardZwanzigFEDifferences;

  /**
   * Constructs a BAR estimator and obtains an initial free energy estimate.
   *
   * @param lambdaValues   Value of lambda for each state.
   * @param eLambdaMinusdL Energies of state L samples at L+dL.
   * @param eLambda        Energies of state L samples at L.
   * @param eLambdaPlusdL  Energies of state L samples at L+dL.
   * @param temperature    Temperature of each state.
   */
  public BennettAcceptanceRatio(double[] lambdaValues, double[][] eLambdaMinusdL, double[][] eLambda,
                                double[][] eLambdaPlusdL, double[] temperature) {
    this(lambdaValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperature, DEFAULT_TOLERANCE);
  }

  /**
   * Constructs a BAR estimator and obtains an initial free energy estimate.
   *
   * @param lambdaValues   Value of lambda for each state.
   * @param eLambdaMinusdL Energies of state L samples at L+dL.
   * @param eLambda        Energies of state L samples at L.
   * @param eLambdaPlusdL  Energies of state L samples at L+dL.
   * @param temperature    Temperature of each state.
   * @param tolerance      Convergence criterion in kcal/mol for BAR iteration.
   */
  public BennettAcceptanceRatio(double[] lambdaValues, double[][] eLambdaMinusdL, double[][] eLambda,
                                double[][] eLambdaPlusdL, double[] temperature, double tolerance) {
    this(lambdaValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperature, tolerance, DEFAULT_MAX_BAR_ITERATIONS);
  }

  /**
   * Constructs a BAR estimator and obtains an initial free energy estimate.
   *
   * @param lambdaValues   Value of lambda for each state.
   * @param eLambdaMinusdL Energies of state L samples at L+dL.
   * @param eLambda        Energies of state L samples at L.
   * @param eLambdaPlusdL  Energies of state L samples at L+dL.
   * @param temperature    Temperature of each state.
   * @param tolerance      Convergence criterion in kcal/mol for BAR iteration.
   * @param nIterations    Maximum number of iterations for BAR.
   */
  public BennettAcceptanceRatio(double[] lambdaValues, double[][] eLambdaMinusdL, double[][] eLambda,
                                double[][] eLambdaPlusdL, double[] temperature, double tolerance, int nIterations) {

    super(lambdaValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperature);

    // Used to seed an initial guess.
    forwardsFEP = new Zwanzig(lambdaValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperature, FORWARDS);
    backwardsFEP = new Zwanzig(lambdaValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperature, BACKWARDS);

    nWindows = nStates - 1;
    forwardZwanzigFEDifferences = forwardsFEP.getFreeEnergyDifferences();
    backwardZwanzigFEDifferences = backwardsFEP.getFreeEnergyDifferences();

    freeEnergyDifferences = new double[nWindows];
    freeEnergyDifferenceUncertainties = new double[nWindows];
    enthalpyDifferences = new double[nWindows];
    this.tolerance = tolerance;
    this.nIterations = nIterations;
    random = new Random();

    estimateDG();
  }

  /**
   * Calculates the Fermi function for the differences used in estimating c.
   *
   * <p>f(x) = 1 / (1 + exp(x)) x = (e1 - e0 + c) * invRT
   *
   * @param e0         Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1         Unperturbed energy (to be subtracted; evaluated at L).
   * @param fermiDiffs Array to be filled with Fermi differences.
   * @param len        Number of energies.
   * @param c          Prior best estimate of the BAR offset/free energy.
   * @param invRT      1.0 / ideal gas constant * temperature.
   */
  private static void fermiDiffIterative(double[] e0, double[] e1, double[] fermiDiffs, int len,
                                         double c, double invRT) {
    for (int i = 0; i < len; i++) {
      fermiDiffs[i] = fermiFunction(invRT * (e0[i] - e1[i] + c));
    }
    if (stream(fermiDiffs).sum() == 0) {
      logger.warning(format(" Input Fermi with length %3d should not be permitted: c: %9.4f invRT: %9.4f Fermi output: %9.4f", len, c, invRT, stream(fermiDiffs).sum()));
    }
  }

  /**
   * Calculates forward alpha and fbsum for BAR Enthalpy estimation.
   *
   * @param e0    Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1    Unperturbed energy (to be subtracted; evaluated at L).
   * @param len   Number of energies.
   * @param c     Prior best estimate of the BAR offset/free energy.
   * @param invRT 1.0 / ideal gas constant * temperature.
   * @param ret   Return alpha and fbsum.
   */
  private void calcAlphaForward(double[] e0, double[] e1, int len, double c,
                                double invRT, double[] ret) {
    double fsum = 0;
    double fvsum = 0;
    double fbvsum = 0;
    double vsum = 0;
    double fbsum = 0;
    for (int i = 0; i < len; i++) {
      double fore = fermiFunction(invRT * (e1[i] - e0[i] - c));
      double back = fermiFunction(invRT * (e0[i] - e1[i] + c));
      fsum += fore;
      fvsum += fore * e0[i];
      fbvsum += fore * back * (e1[i] - e0[i]);
      vsum += e0[i];
      fbsum += fore * back;
    }
    double alpha = fvsum - (fsum * (vsum / len)) + fbvsum;
    ret[0] = alpha;
    ret[1] = fbsum;
  }

  /**
   * Calculates backward alpha and fbsum for BAR Enthalpy estimation.
   *
   * @param e0    Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1    Unperturbed energy (to be subtracted; evaluated at L).
   * @param len   Number of energies.
   * @param c     Prior best estimate of the BAR offset/free energy.
   * @param invRT 1.0 / ideal gas constant * temperature.
   * @param ret   Return alpha and fbsum.
   */
  private void calcAlphaBackward(double[] e0, double[] e1, int len, double c,
                                 double invRT, double[] ret) {
    double bsum = 0;
    double bvsum = 0;
    double fbvsum = 0;
    double vsum = 0;
    double fbsum = 0;
    for (int i = 0; i < len; i++) {
      double fore = fermiFunction(invRT * (e1[i] - e0[i] - c));
      double back = fermiFunction(invRT * (e0[i] - e1[i] + c));
      bsum += back;
      bvsum += back * e1[i];
      fbvsum += fore * back * (e1[i] - e0[i]);
      vsum += e1[i];
      fbsum += fore * back;
    }
    double alpha = bvsum - (bsum * (vsum / len)) - fbvsum;
    ret[0] = alpha;
    ret[1] = fbsum;
  }


  /**
   * Calculates the Fermi function for the differences used in estimating c, using bootstrap sampling
   * (choosing random indices with replacement rather than scanning through them all).
   *
   * <p>f(x) = 1 / (1 + exp(x)) x = (e1 - e0 + c) * invRT
   *
   * @param e0         Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1         Unperturbed energy (to be subtracted; evaluated at L).
   * @param fermiDiffs Array to be filled with Fermi differences.
   * @param len        Number of energies.
   * @param c          Prior best estimate of the BAR offset/free energy.
   * @param invRT      1.0 / ideal gas constant * temperature.
   */
  private static void fermiDiffBootstrap(double[] e0, double[] e1, double[] fermiDiffs,
                                         int len, double c, double invRT, int[] bootstrapSamples) {
    for (int indexI = 0; indexI < len; indexI++) {
      int i = bootstrapSamples[indexI];
      fermiDiffs[indexI] = fermiFunction(invRT * (e0[i] - e1[i] + c));
    }
  }


  /**
   * Computes one half of the BAR variance.
   *
   * @param meanFermi   Mean Fermi value for either state 0 or state 1.
   * @param meanSqFermi Mean squared Fermi value for either state 0 or state 1.
   * @param len         Number of values.
   * @return One half of BAR variance.
   */
  private static double uncertaintyCalculation(double meanFermi, double meanSqFermi, int len) {
    double sqMeanFermi = meanFermi * meanFermi;
    return ((meanSqFermi - sqMeanFermi) / len) / sqMeanFermi;
  }

  /**
   * Returns the backwards Zwanzig estimator used to seed BAR.
   *
   * @return A backwards Zwanzig estimator.
   */
  public Zwanzig getInitialBackwardsGuess() {
    return backwardsFEP;
  }

  /**
   * Returns the forwards Zwanzig estimator used to seed BAR.
   *
   * @return A forwards Zwanzig estimator.
   */
  public Zwanzig getInitialForwardsGuess() {
    return forwardsFEP;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public BennettAcceptanceRatio copyEstimator() {
    return new BennettAcceptanceRatio(lamValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperatures, tolerance, nIterations);
  }

  /**
   * Main driver for estimation of BAR free energy differences.
   * <p>
   * Based on Tinker implementation, which uses the substitution proposed in
   * Wyczalkowski, Vitalis and Pappu 2010.
   *
   * @param randomSamples Whether to use random sampling (for bootstrap analysis).
   */
  @Override
  public final void estimateDG(final boolean randomSamples) {
    double cumDG = 0;
    fill(freeEnergyDifferences, 0);
    fill(freeEnergyDifferenceUncertainties, 0);
    fill(enthalpyDifferences, 0);

    // Avoid duplicate warnings when bootstrapping.
    Level warningLevel = randomSamples ? Level.FINE : Level.WARNING;

    for (int i = 0; i < nWindows; i++) {
      // Free energy estimate/shift constant.
      if (isNaN(forwardZwanzigFEDifferences[i]) || isInfinite(forwardZwanzigFEDifferences[i])
          || isNaN(backwardZwanzigFEDifferences[i]) || isInfinite(backwardZwanzigFEDifferences[i])) {
        logger.warning(format(" Window %3d bin energies produced unreasonable value(s) for forward Zwanzig (%8.4f) and/or backward Zwanzig (%8.4f)", i, forwardZwanzigFEDifferences[i], backwardZwanzigFEDifferences[i]));
      }
      double c = 0.5 * (forwardZwanzigFEDifferences[i] + backwardZwanzigFEDifferences[i]);

      if (!randomSamples) {
        logger.fine(format(" BAR Iteration   %2d: %12.4f Kcal/mol", 0, c));
      }

      double cold = c;
      int len0 = eLambda[i].length;
      int len1 = eLambda[i + 1].length;

      if (len0 == 0 || len1 == 0) {
        freeEnergyDifferences[i] = c;
        logger.log(warningLevel, format(" Window %d has no snapshots at one end (%d, %d)!", i, len0, len1));
        continue;
      }

      // Ratio of the number of snaps: Tinker equivalent: rfrm
      double sampleRatio = ((double) len0) / ((double) len1);

      // Fermi differences.
      double[] fermi0 = new double[len0];
      double[] fermi1 = new double[len1];
      double[] ret = new double[2];

      // Ideal gas constant * temperature, or its inverse.
      double rta = R * temperatures[i];
      double rtb = R * temperatures[i + 1];
      double rtMean = 0.5 * (rta + rtb);
      double invRTA = 1.0 / rta;
      double invRTB = 1.0 / rtb;

      // Summary statistics for Fermi differences for the upper half.
      SummaryStatistics s1 = null;
      // Summary statistics for Fermi differences for the lower half.
      SummaryStatistics s0 = null;

      // Each BAR convergence cycle needs to operate on the same set of indices.
      int[] bootstrapSamples0 = null;
      int[] bootstrapSamples1 = null;

      if (randomSamples) {
        bootstrapSamples0 = getBootstrapIndices(len0, random);
        bootstrapSamples1 = getBootstrapIndices(len1, random);
      }

      int cycleCounter = 0;
      boolean converged = false;
      while (!converged) {
        if (randomSamples) {
          fermiDiffBootstrap(eLambdaPlusdL[i], eLambda[i], fermi0, len0, -c, invRTA, bootstrapSamples0);
          fermiDiffBootstrap(eLambdaMinusdL[i + 1], eLambda[i + 1], fermi1, len1, c, invRTB, bootstrapSamples1);
        } else {
          fermiDiffIterative(eLambdaPlusdL[i], eLambda[i], fermi0, len0, -c, invRTA);
          fermiDiffIterative(eLambdaMinusdL[i + 1], eLambda[i + 1], fermi1, len1, c, invRTB);
        }

        s0 = new SummaryStatistics(fermi0);
        s1 = new SummaryStatistics(fermi1);
        double ratio = stream(fermi1).sum() / stream(fermi0).sum();

        c += rtMean * log(sampleRatio * ratio);

        converged = (abs(c - cold) < tolerance);

        if (!converged && ++cycleCounter > nIterations) {
          throw new IllegalArgumentException(
              format(" BAR required too many iterations (%d) to converge! (%9.8f > %9.8f)", cycleCounter, abs(c - cold), tolerance));
        }

        if (!randomSamples) {
          logger.fine(format(" BAR Iteration   %2d: %12.4f Kcal/mol", cycleCounter, c));
        }
        cold = c;
      }

      freeEnergyDifferences[i] = c;
      cumDG += c;
      double sqFermiMean0 = new SummaryStatistics(stream(fermi0).map((double d) -> d * d).toArray()).mean;
      double sqFermiMean1 = new SummaryStatistics(stream(fermi1).map((double d) -> d * d).toArray()).mean;
      freeEnergyDifferenceUncertainties[i] = sqrt(uncertaintyCalculation(s0.mean, sqFermiMean0, len0)
          + uncertaintyCalculation(s1.mean, sqFermiMean1, len1));

      calcAlphaForward(eLambda[i], eLambdaPlusdL[i], len0, c, invRTA, ret);
      double alpha0 = ret[0];
      double fbsum0 = ret[1];

      calcAlphaBackward(eLambdaMinusdL[i + 1], eLambda[i + 1], len1, c, invRTB, ret);
      double alpha1 = ret[0];
      double fbsum1 = ret[1];

      double hBar = (alpha0 - alpha1) / (fbsum0 + fbsum1);
      enthalpyDifferences[i] = hBar;
    }

    totalFreeEnergyDifference = cumDG;
    totalFEDifferenceUncertainty = sqrt(stream(freeEnergyDifferenceUncertainties).map((double d) -> d * d).sum());
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getFreeEnergyDifferences() {
    return copyOf(freeEnergyDifferences, nWindows);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getFEDifferenceUncertainties() {
    return copyOf(freeEnergyDifferenceUncertainties, nWindows);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalFreeEnergyDifference() {
    return totalFreeEnergyDifference;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalFEDifferenceUncertainty() {
    return totalFEDifferenceUncertainty;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfBins() {
    return nWindows;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnthalpyDifference() {
    return getTotalEnthalpyDifference(enthalpyDifferences);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getEnthalpyDifferences() {
    return enthalpyDifferences;
  }
}
