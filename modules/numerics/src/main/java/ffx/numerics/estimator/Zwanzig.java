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

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

import static ffx.numerics.estimator.EstimateBootstrapper.getBootstrapIndices;
import static ffx.utilities.Constants.R;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;

/**
 * The Zwanzig class implements exponential averaging/free energy perturbation using the Zwanzig
 * relationship, in either the forwards or backwards direction (not both).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class Zwanzig extends SequentialEstimator implements BootstrappableEstimator {

  private static final Logger logger = Logger.getLogger(Zwanzig.class.getName());

  /**
   * Directionality of the Zwanzig estimation (forwards perturbation or backwards perturbation).
   */
  public final Directionality directionality;
  /**
   * Whether the Zwanzig estimation is forwards or backwards.
   */
  private final boolean forwards;
  /**
   * Number of windows.
   */
  private final int nWindows;
  /**
   * Free energy difference for each window.
   */
  private final double[] freeEnergyDifferences;
  /**
   * Enthalpy difference for each window.
   */
  private final double[] enthalpyDifferences;
  /**
   * Free energy difference uncertainty for each window.
   */
  private final double[] freeEnergyDifferenceUncertainties;
  /**
   * Random number generator for bootstrapping.
   */
  private final Random random;
  /**
   * Total free energy difference as a sum over windows.
   */
  private double totalFreeEnergyDifference;
  /**
   * Total free energy difference uncertainty:
   * totalFreeEnergyUncertainty = Sqrt [ Sum over Windows [ Window Variance ] ]
   */
  private double totalFreeEnergyDifferenceUncertainty;
  /**
   * The total enthalpy difference.
   */
  private double totalEnthalpyDifference;

  /**
   * Estimates a free energy using the Zwanzig relationship. The temperature array can be of length 1
   * if all elements are meant to be the same temperature.
   *
   * <p>The first dimension of the energies arrays corresponds to the lambda values/windows. The
   * second dimension (can be of uneven length) corresponds to potential energies of snapshots
   * sampled from that lambda value, calculated either at that lambda value, the lambda value below,
   * or the lambda value above. The arrays eLambdaMinusdL[0] and eLambdaPlusdL[n-1] is expected to be all
   * NaN.
   *
   * @param lambdaValues   Lambda values for the samples.
   * @param eLambdaMinusdL Potential energies of state L at lambda L-dL. Ignored for forwards FEP.
   * @param eLambda        Potential energies of state L at lambda L.
   * @param eLambdaPlusdL  Potential energies of state L at lambda L+dL. Ignored for backwards FEP.
   * @param temperature    Temperature each lambda window was run at (single-element indicates identical temperatures).
   * @param directionality Forwards vs. backwards FEP.
   */
  public Zwanzig(double[] lambdaValues, double[][] eLambdaMinusdL, double[][] eLambda,
                 double[][] eLambdaPlusdL, double[] temperature, Directionality directionality) {
    super(lambdaValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperature);
    this.directionality = directionality;
    nWindows = nStates - 1;

    freeEnergyDifferences = new double[nWindows];
    freeEnergyDifferenceUncertainties = new double[nWindows];
    enthalpyDifferences = new double[nWindows];

    forwards = directionality.equals(Directionality.FORWARDS);
    random = new Random();

    estimateDG();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Zwanzig copyEstimator() {
    return new Zwanzig(lamValues, eLambdaMinusdL, eLambda, eLambdaPlusdL, temperatures, directionality);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public final void estimateDG(final boolean randomSamples) {
    double cumDG = 0;

    Level warningLevel = randomSamples ? Level.FINE : Level.WARNING;

    for (int i = 0; i < nWindows; i++) {
      final int windowIndex = i + (forwards ? 0 : 1);
      final double[] referenceEnergy = eLambda[windowIndex];
      final double[] perturbedEnergy = forwards ? eLambdaPlusdL[windowIndex] : eLambdaMinusdL[windowIndex];
      final double sign = forwards ? 1.0 : -1.0;
      final double kT = temperatures[windowIndex] * R;
      final double beta = 1.0 / kT;
      final int numSamples = referenceEnergy.length;
      if (numSamples == 0) {
        logger.log(warningLevel, " Skipping frame " + i + " due to lack of snapshots!");
        continue;
      }

      // With no iteration-to-convergence, generating a fresh random index is OK.
      int[] sampleIndices = randomSamples ? getBootstrapIndices(numSamples, random) : IntStream.range(0, numSamples).toArray();

      double boltzmannFactorSum = 0.0;
      double weightedMeanEnergy = 0.0;
      double meanEnergy = 0.0;
      for (int j = 0; j < numSamples; j++) {
        int index = sampleIndices[j];
        final double dE = perturbedEnergy[index] - referenceEnergy[index];
        final double weight = exp(-beta * dE);
        boltzmannFactorSum += weight;
        if (forwards) {
          meanEnergy += referenceEnergy[index];
          weightedMeanEnergy += perturbedEnergy[index] * weight;
        } else {
          meanEnergy += perturbedEnergy[index];
          weightedMeanEnergy += referenceEnergy[index] * weight;
        }
      }
      meanEnergy = meanEnergy / numSamples;
      double dG = -sign * kT * log(boltzmannFactorSum / numSamples);

      if (isNaN(dG) || isInfinite(dG)) {
        logger.severe(format(" Change in free energy (%9.4f) for window (%2d of %2d) failed. " +
                "Sign: %9.4f, Beta: %9.4f, Temp: %9.4f, Sum: %9.4f, Len: %3d, Log: %9.4f",
            dG, i, nWindows, sign, kT, temperatures[windowIndex], boltzmannFactorSum, numSamples, log(boltzmannFactorSum / numSamples)));
      }

      freeEnergyDifferences[i] = dG;
      enthalpyDifferences[i] = sign * ((weightedMeanEnergy / boltzmannFactorSum) - meanEnergy);
      freeEnergyDifferenceUncertainties[i] = 0.0;
      totalEnthalpyDifference += enthalpyDifferences[i];
      cumDG += dG;
    }

    totalFreeEnergyDifference = cumDG;
    totalFreeEnergyDifferenceUncertainty = 0.0;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public final void estimateDG() {
    estimateDG(false);
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
  public double[] getFreeEnergyDifferences() {
    return copyOf(freeEnergyDifferences, nWindows);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalFEDifferenceUncertainty() {
    return totalFreeEnergyDifferenceUncertainty;
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
  public int getNumberOfBins() {
    return nWindows;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnthalpyDifference() {
    return totalEnthalpyDifference;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getEnthalpyDifferences() {
    return copyOf(enthalpyDifferences, nWindows);
  }

  /**
   * Directionality of the Zwanzig estimation (forwards perturbation or backwards perturbation).
   */
  public enum Directionality {
    /**
     * Forwards perturbation.
     */
    FORWARDS,
    /**
     * Backwards perturbation.
     */
    BACKWARDS
  }
}
