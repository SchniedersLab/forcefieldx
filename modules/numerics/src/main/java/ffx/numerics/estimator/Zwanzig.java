// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
import static java.util.Arrays.copyOf;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;

import ffx.utilities.Constants;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

/**
 * The Zwanzig class implements exponential averaging/free energy perturbation using the Zwanzig
 * relationship, in either the forwards or backwards direction (not both).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class Zwanzig extends SequentialEstimator implements BootstrappableEstimator {

  private static final Logger logger = Logger.getLogger(SequentialEstimator.class.getName());
  public final Directionality directionality;
  private final boolean forwards;
  /**
   * Number of windows.
   */
  private final int nWindows;
  /**
   * Free energy difference for each window.
   */
  private final double[] windowFreeEnergyDifferences;
  private final double[] windowEnthalpies;
  /**
   * Free energy difference uncertainty for each window.
   */
  private final double[] windowFreeEnergyUncertainties;
  private final Random random;
  /**
   * Total free energy difference as a sum over windows.
   */
  private double totalFreeEnergyDifference;
  /**
   * Total free energy difference uncertainty: totalFreeEnergyUncertainty = Sqrt [ Sum over Windows [
   * Window Uncertainty Squared ] ]
   */
  private double totalFreeEnergyUncertainty;

  /**
   * Total Enthalpy
   */
  private double totalEnthalpy;

  /**
   * Estimates a free energy using the Zwanzig relationship. The temperature array can be of length 1
   * if all elements are meant to be the same temperature.
   *
   * <p>The first dimension of the energies arrays corresponds to the lambda values/windows. The
   * second dimension (can be of uneven length) corresponds to potential energies of snapshots
   * sampled from that lambda value, calculated either at that lambda value, the lambda value below,
   * or the lambda value above. The arrays energiesLow[0] and energiesHigh[n-1] is expected to be all
   * NaN.
   *
   * @param lambdaValues Values of lambda dynamics was run at.
   * @param energiesLow Potential energies of trajectory L at lambda L-dL. Ignored for forwards
   *     FEP.
   * @param energiesAt Potential energies of trajectory L at lambda L.
   * @param energiesHigh Potential energies of trajectory L at lambda L+dL. Ignored for backwards
   *     FEP.
   * @param temperature Temperature each lambda window was run at (single-element indicates
   *     identical temperatures).
   * @param directionality Forwards vs. backwards FEP.
   */
  public Zwanzig(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt,
      double[][] energiesHigh,
      double[] temperature, Directionality directionality) {
    super(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature);
    this.directionality = directionality;
    nWindows = nTrajectories - 1;

    windowFreeEnergyDifferences = new double[nWindows];
    windowFreeEnergyUncertainties = new double[nWindows];
    windowEnthalpies = new double[nWindows];

    forwards = directionality.equals(Directionality.FORWARDS);
    random = new Random();

    estimateDG();
  }

  /** {@inheritDoc} */
  @Override
  public Zwanzig copyEstimator() {
    return new Zwanzig(lamVals, eLow, eAt, eHigh, temperatures, directionality);
  }

  /** {@inheritDoc} */
  @Override
  public void estimateDG(final boolean randomSamples) {
    double cumDG = 0;

    Level warningLevel = randomSamples ? Level.FINE : Level.WARNING;

    for (int i = 0; i < nWindows; i++) {
      int windowIndex = forwards ? 0 : 1;
      windowIndex += i;
      double[] e1 = eAt[windowIndex];
      double[] e2 = forwards ? eHigh[windowIndex] : eLow[windowIndex];
      double sign = forwards ? 1.0 : -1.0;
      double beta = -temperatures[windowIndex] * Constants.R;
      double invBeta = 1.0 / beta;

      int len = e1.length;
      if (len == 0) {
        logger.log(warningLevel, " Skipping frame " + i + " due to lack of snapshots!");
        continue;
      }

      // With no iteration-to-convergence, generating a fresh random index is OK.
      int[] samples =
          randomSamples ? getBootstrapIndices(len, random) : IntStream.range(0, len).toArray();

      double sum = 0.0;
      double num = 0.0;
      double term = 0.0;
      double uAve = 0.0;
      for (int indJ = 0; indJ < len; indJ++) {
        int j = samples[indJ];
        double dE = e2[j] - e1[j];
        if (sign > 0){
          uAve += e1[j];
          term = exp(invBeta*dE);
          sum += term;
          num += e2[j] * term;
        } else {
          uAve += e2[j];
          term = exp(invBeta*dE);
          sum += term;
          num += e1[j] * term;
        }


      }

      uAve = uAve/len;

      double dG = sign * beta * log(sum / len);
      windowFreeEnergyDifferences[i] = dG;
      windowEnthalpies[i] = sign*((num / sum) - uAve);
      windowFreeEnergyUncertainties[i] = 0.0;
      totalEnthalpy += windowEnthalpies[i];
      cumDG += dG;
    }

    totalFreeEnergyDifference = cumDG;
    totalFreeEnergyUncertainty = 0.0;
  }

  /** {@inheritDoc} */
  @Override
  public void estimateDG() {
    estimateDG(false);
  }

  /** {@inheritDoc} */
  @Override
  public double[] getBinEnergies() {
    return copyOf(windowFreeEnergyDifferences, nWindows);
  }

  /** {@inheritDoc} */
  @Override
  public double[] getBinUncertainties() {
    return copyOf(windowFreeEnergyUncertainties, nWindows);
  }

  /** {@inheritDoc} */
  @Override
  public double getFreeEnergy() {
    return totalFreeEnergyDifference;
  }

  /** {@inheritDoc} */
  @Override
  public double getUncertainty() {
    return totalFreeEnergyUncertainty;
  }

  /** {@inheritDoc} */
  @Override
  public int numberOfBins() {
    return nWindows;
  }

  /** {@inheritDoc} */
  @Override
  public double[] getBinEthalpies() {
    return copyOf(windowEnthalpies, nWindows);
  }

  /**
   * Directionality of the Zwanzig estimation (forwards perturbation or backwards perturbation).
   * TODO: Implement bidirectional Zwanzig with simple estimation (i.e. 0.5*(forwards + backward)).
   */
  public enum Directionality {
    FORWARDS,
    BACKWARDS
  }
}
