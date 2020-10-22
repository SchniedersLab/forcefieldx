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

import static ffx.numerics.estimator.EstimateBootstrapper.getBootstrapIndices;
import static ffx.numerics.estimator.Zwanzig.Directionality.BACKWARDS;
import static ffx.numerics.estimator.Zwanzig.Directionality.FORWARDS;
import static ffx.numerics.math.ScalarMath.fermiFunction;
import static java.lang.String.format;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static java.util.Arrays.stream;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.math.SummaryStatistics;
import ffx.utilities.Constants;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

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
   * Maximum number of BAR iterations.
   */
  private static final int MAX_ITERS = 100;
  /**
   * Number of simulation windows.
   */
  private final int nWindows;
  /**
   * Forward Zwanzig free-energy difference estimates.
   */
  private final double[] forwardZwanzig;
  /**
   * Backward Zwanzig free-energy difference estimates.
   */
  private final double[] backwardZwanzig;
  /**
   * BAR free-energy difference estimates.
   */
  private final double[] barEstimates;
  /**
   * BAR free-energy difference uncertainties.
   */
  private final double[] barUncertainties;
  /**
   * BAR convergence tolerance.
   */
  private final double tolerance;
  /**
   * Forward Zwanzig instance.
   */
  private final Zwanzig forwardsFEP;
  /**
   * Backward Zwanzig instance.
   */
  private final Zwanzig backwardsFEP;
  private final Random random;
  /**
   * Total BAR free-energy difference estimate.
   */
  private double totalBAREstimate;
  /**
   * Total BAR free-energy difference uncertainty.
   */
  private double totalBARUncertainty;
  /**
   * BAR Enthalpy estimates
   */
  private double[] barEnthalpy;
  /**
   * Alpha for BAR Enthalpy calculations
   */
  private double alpha;
  /**
   * sum for BAR Enthalpy Calculations
   */
  private double fbsum;

  /**
   * Constructs a BAR estimator and obtains an initial free energy estimate.
   *
   * @param lambdaValues Values of lambda used.
   * @param energiesLow Energies of trajectory i at lambda (i-1).
   * @param energiesAt Energies of trajectory i at lambda i.
   * @param energiesHigh Energies of trajectory i at lambda (i+1).
   * @param temperature Temperature of each trajectory.
   */
  public BennettAcceptanceRatio(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt,
      double[][] energiesHigh, double[] temperature) {
    this(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, DEFAULT_TOLERANCE);
  }

  /**
   * Constructs a BAR estimator and obtains an initial free energy estimate.
   *
   * @param lambdaValues Values of lambda used.
   * @param energiesLow Energies of trajectory i at lambda (i-1).
   * @param energiesAt Energies of trajectory i at lambda i.
   * @param energiesHigh Energies of trajectory i at lambda (i+1).
   * @param temperature Temperature of each trajectory.
   * @param tolerance Convergence criterion in kcal/mol for BAR iteration.
   */
  public BennettAcceptanceRatio(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt,
      double[][] energiesHigh, double[] temperature, double tolerance) {

    super(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature);

    // Used to seed an initial guess.
    forwardsFEP = new Zwanzig(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature,
        FORWARDS);
    backwardsFEP = new Zwanzig(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature,
        BACKWARDS);

    nWindows = nTrajectories - 1;
    forwardZwanzig = forwardsFEP.getBinEnergies();
    backwardZwanzig = backwardsFEP.getBinEnergies();

    barEstimates = new double[nWindows];
    barUncertainties = new double[nWindows];
    barEnthalpy = new double[nWindows];
    this.tolerance = tolerance;
    random = new Random();

    estimateDG();
  }

  /**
   * Calculates the Fermi function for the differences used in estimating c.
   *
   * <p>f(x) = 1 / (1 + exp(x)) x = (e1 - e0 + c) * invRT
   *
   * @param e0 Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1 Unperturbed energy (to be subtracted; evaluated at L).
   * @param fermiDiffs Array to be filled with Fermi differences.
   * @param len Number of energies.
   * @param c Prior best estimate of the BAR offset/free energy.
   * @param invRT 1.0 / ideal gas constant * temperature.
   */
  private static void fermiDiffIterative(double[] e0, double[] e1, double[] fermiDiffs, int len,
      double c, double invRT) {
    for (int i = 0; i < len; i++) {
      fermiDiffs[i] = fermiFunction(invRT * (e0[i] - e1[i] + c));
    }
  }

  /**
   * Calculates forward alpha and fbsum for BAR Enthalpy calculations
   * @param e0 Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1 Unperturbed energy (to be subtracted; evaluated at L).
   * @param len Number of energies.
   * @param c Prior best estimate of the BAR offset/free energy.
   * @param invRT 1.0 / ideal gas constant * temperature.
   */
  private void calcAlphaForward( double[] e0, double[] e1, int len, double c, double invRT) {

    double fsum = 0 ;
    double fvsum = 0;
    double fbvsum = 0;
    double vsum = 0;
    fbsum = 0;
    for (int i = 0; i < len; i++) {
      double fore = fermiFunction(invRT * (e1[i] - e0[i] - c));
      double back = fermiFunction(invRT * (e0[i] - e1[i] + c));
      fsum += fore;
      fvsum  += fore*e0[i];
      fbvsum  += fore*back*(e1[i]-e0[i]);
      vsum  += e0[i];
      fbsum  += fore*back;

    }
    alpha = fvsum - (fsum*(vsum/len)) + fbvsum;


  }
  /**
   * Calculates backward alpha and fbsum for  BAR Enthalpy calculations
   * @param e0 Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1 Unperturbed energy (to be subtracted; evaluated at L).
   * @param len Number of energies.
   * @param c Prior best estimate of the BAR offset/free energy.
   * @param invRT 1.0 / ideal gas constant * temperature.
   */
  private void calcAlphaBackward( double[] e0, double[] e1, int len, double c, double invRT) {

    double bsum = 0;
    double bvsum = 0;
    double fbvsum = 0;
    double vsum = 0;
    fbsum = 0;
    for (int i = 0; i < len; i++) {
      double fore = fermiFunction(invRT * (e1[i] - e0[i] - c));
      double back = fermiFunction(invRT * (e0[i] - e1[i] + c));
      bsum += back;
      bvsum  += back*e1[i];
      fbvsum  += fore*back*(e1[i]-e0[i]);
      vsum  += e1[i];
      fbsum  += fore*back;

    }
    alpha = bvsum - (bsum*(vsum/len)) + fbvsum;


  }


  /**
   * Calculates the Fermi function for the differences used in estimating c, using bootstrap sampling
   * (choosing random indices with replacement rather than scanning through them all).
   *
   * <p>f(x) = 1 / (1 + exp(x)) x = (e1 - e0 + c) * invRT
   *
   * @param e0 Perturbed energy (to be added; evaluated at L +/- dL).
   * @param e1 Unperturbed energy (to be subtracted; evaluated at L).
   * @param fermiDiffs Array to be filled with Fermi differences.
   * @param len Number of energies.
   * @param c Prior best estimate of the BAR offset/free energy.
   * @param invRT 1.0 / ideal gas constant * temperature.
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
   * @param meanFermi Mean Fermi value for either state 0 or state 1.
   * @param meanSqFermi Mean squared Fermi value for either state 0 or state 1.
   * @param len Number of values.
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
   * Main driver for estimation of delta-G. Based on Tinker implementation, which uses the
   * substitution proposed in Wyczalkowski, Vitalis and Pappu 2010.
   */
  @Override
  public void estimateDG() {
    estimateDG(false);
  }

  /**
   * Returns the forwards Zwanzig estimator used to seed BAR.
   *
   * @return A forwards Zwanzig estimator.
   */
  public Zwanzig getInitialForwardsGuess() {
    return forwardsFEP;
  }

  /** {@inheritDoc} */
  @Override
  public BennettAcceptanceRatio copyEstimator() {
    return new BennettAcceptanceRatio(lamVals, eLow, eAt, eHigh, temperatures, tolerance);
  }

  /**
   * Main driver for estimation of delta-G. Based on Tinker implementation, which uses the
   * substitution proposed in Wyczalkowski, Vitalis and Pappu 2010.
   *
   * @param randomSamples Whether to use random sampling (for bootstrap analysis).
   */
  @Override
  public void estimateDG(final boolean randomSamples) {
    double cumDG = 0;
    fill(barEstimates, 0);
    fill(barUncertainties, 0);
    fill(barEnthalpy,0);

    // Avoid duplicate warnings when bootstrapping.
    Level warningLevel = randomSamples ? Level.FINE : Level.WARNING;

    for (int i = 0; i < nWindows; i++) {
      // Free energy estimate/shift constant.
      double c = 0.5 * (forwardZwanzig[i] + backwardZwanzig[i]);

      if (!randomSamples) {
        logger.fine(format(" BAR Iteration   %2d: %12.4f Kcal/mol", 0, c));
      }

      double cold = c;
      int len0 = eAt[i].length;
      int len1 = eAt[i + 1].length;

      if (len0 == 0 || len1 == 0) {
        barEstimates[i] = c;
        logger.log(warningLevel,
            format(" Window %d has no snapshots at one end (%d, %d)!", i, len0, len1));
        continue;
      }

      // Ratio of the number of samples: Tinker equivalent: rfrm
      double sampleRatio = ((double) len0) / ((double) len1);

      // Fermi differences.
      double[] fermi0 = new double[len0];
      double[] fermi1 = new double[len1];

      // Ideal gas constant * temperature, or its inverse.
      double rta = Constants.R * temperatures[i];
      double rtb = Constants.R * temperatures[i + 1];
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
          fermiDiffBootstrap(eHigh[i], eAt[i], fermi0, len0, -c, invRTA, bootstrapSamples0);
          fermiDiffBootstrap(eLow[i + 1], eAt[i + 1], fermi1, len1, c, invRTB, bootstrapSamples1);
        } else {
          fermiDiffIterative(eHigh[i], eAt[i], fermi0, len0, -c, invRTA);
          fermiDiffIterative(eLow[i + 1], eAt[i + 1], fermi1, len1, c, invRTB);
        }

        s0 = new SummaryStatistics(fermi0);
        s1 = new SummaryStatistics(fermi1);
        double ratio = stream(fermi1).sum() / stream(fermi0).sum();

        c += rtMean * log(sampleRatio * ratio);

        converged = (abs(c - cold) < tolerance);
        cold = c;

        if (++cycleCounter > MAX_ITERS) {
          throw new IllegalArgumentException(
              format(" BAR required too many iterations (%d) to converge!", cycleCounter));
        }

        if (!randomSamples) {
          logger.fine(format(" BAR Iteration   %2d: %12.4f Kcal/mol", cycleCounter, c));
        }
      }

      barEstimates[i] = c;
      cumDG += c;
      double sqFermiMean0 = new SummaryStatistics(
          stream(fermi0).map((double d) -> d * d).toArray()).mean;
      double sqFermiMean1 = new SummaryStatistics(
          stream(fermi1).map((double d) -> d * d).toArray()).mean;
      barUncertainties[i] = sqrt(uncertaintyCalculation(s0.mean, sqFermiMean0, len0)
          + uncertaintyCalculation(s1.mean, sqFermiMean1, len1));


      calcAlphaForward(eAt[i],eHigh[i],len0,c,invRTA);
      double alpha0 = alpha;
      double fbsum0 = fbsum;

      calcAlphaBackward(eLow[i+1],eAt[i+1],len1,c,invRTB);
      double alpha1 = alpha;
      double fbsum1 = fbsum;

      double hBar = 0;
      hBar = (alpha0 - alpha1) / (fbsum0+fbsum1);
      barEnthalpy[i] = hBar;

    }




    totalBAREstimate = cumDG;
    totalBARUncertainty = sqrt(stream(barUncertainties).map((double d) -> d * d).sum());
  }

  /** {@inheritDoc} */
  @Override
  public double[] getBinEnergies() {
    return copyOf(barEstimates, nWindows);
  }

  /** {@inheritDoc} */
  @Override
  public double[] getBinUncertainties() {
    return copyOf(barUncertainties, nWindows);
  }

  /** {@inheritDoc} */
  @Override
  public double getFreeEnergy() {
    return totalBAREstimate;
  }

  /** {@inheritDoc} */
  @Override
  public double getUncertainty() {
    return totalBARUncertainty;
  }

  /** {@inheritDoc} */
  @Override
  public int numberOfBins() {
    return nWindows;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getBinEthalpies() {
    return barEnthalpy;
  }
}
