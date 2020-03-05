//******************************************************************************
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
//******************************************************************************
package ffx.numerics.estimator;

import ffx.numerics.math.FFXSummaryStatistics;
import ffx.numerics.math.ScalarMath;
import ffx.utilities.Constants;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.sqrt;

import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;

/**
 * The Bennett Acceptance Ratio class implements the Bennett Acceptance Ratio (BAR)
 * statistical estimator, based on the Tinker implementation.
 *
 * Literature References (from Tinker):
 * C. H. Bennett, "Efficient Estimation of Free Energy Differences
 * from Monte Carlo Data", Journal of Computational Physics, 22,
 * 245-268 (1976)
 * c
 * K. B. Daly, J. B. Benziger, P. G. Debenedetti and
 * A. Z. Panagiotopoulos, "Massively Parallel Chemical Potential
 * Calculation on Graphics Processing Units", Computer Physics
 * Communications, 183, 2054-2062 (2012)  [modification for NPT]
 * c
 * M. A. Wyczalkowski, A. Vitalis and R. V. Pappu, "New Estimators
 * for Calculating Solvation Entropy and Enthalpy and Comparative
 * Assessments of Their Accuracy and Precision, Journal of Physical
 * Chemistry, 114, 8166-8180 (2010)  [entropy and enthalpy]
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class BennettAcceptanceRatio extends SequentialEstimator implements BootstrappableEstimator {
    private static final Logger logger = Logger.getLogger(BennettAcceptanceRatio.class.getName());

    private final int nWindows;
    private final double[] estForwards;
    private final double[] estBackwards;
    private final double[] dGs;
    private final double[] uncerts;
    private final double tolerance;
    private double totDG;
    private double totUncert;
    private static final double DEFAULT_TOLERANCE = 1.0E-7;
    private static final int MAX_ITERS = 100;
    private final Random random = new Random(); // TODO: Use provided seed.

    // Hang onto these in case the end-user wants them?
    private final Zwanzig forwardsFEP;
    private final Zwanzig backwardsFEP;

    /**
     * Constructs a BAR estimator and obtains an initial free energy estimate.
     *
     * @param lambdaValues Values of lambda used.
     * @param energiesLow  Energies of trajectory i at lambda (i-1).
     * @param energiesAt   Energies of trajectory i at lambda i.
     * @param energiesHigh Energies of trajectory i at lambda (i+1).
     * @param temperature  Temperature of each trajectory.
     */
    public BennettAcceptanceRatio(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt, double[][] energiesHigh, double[] temperature) {
        this(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, DEFAULT_TOLERANCE);
    }

    /**
     * Constructs a BAR estimator and obtains an initial free energy estimate.
     *
     * @param lambdaValues Values of lambda used.
     * @param energiesLow  Energies of trajectory i at lambda (i-1).
     * @param energiesAt   Energies of trajectory i at lambda i.
     * @param energiesHigh Energies of trajectory i at lambda (i+1).
     * @param temperature  Temperature of each trajectory.
     * @param tolerance    Convergence criterion in kcal/mol for BAR iteration.
     */
    public BennettAcceptanceRatio(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt, double[][] energiesHigh, double[] temperature, double tolerance) {
        super(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature);
        // Used to seed an initial guess.
        forwardsFEP = new Zwanzig(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, Zwanzig.Directionality.FORWARDS);
        backwardsFEP = new Zwanzig(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, Zwanzig.Directionality.BACKWARDS);

        nWindows = nTrajectories - 1;
        estForwards = forwardsFEP.getBinEnergies();
        estBackwards = backwardsFEP.getBinEnergies();

        dGs = new double[nWindows];
        uncerts = new double[nWindows];
        this.tolerance = tolerance;
        estimateDG();
    }

    /**
     * Main driver for estimation of delta-G. Based on Tinker implementation, which uses the substitution
     * proposed in Wyczalkowski, Vitalis and Pappu 2010.
     */
    @Override
    public void estimateDG() {
        estimateDG(false);
    }

    /**
     * Main driver for estimation of delta-G. Based on Tinker implementation, which uses the substitution
     * proposed in Wyczalkowski, Vitalis and Pappu 2010.
     *
     * @param randomSamples Whether to use random sampling (for bootstrap analysis).
     */
    @Override
    public void estimateDG(final boolean randomSamples) {
        double cumDG = 0;
        Arrays.fill(dGs, 0);
        Arrays.fill(uncerts, 0);
        for (int i = 0; i < nWindows; i++) {
            boolean converged = false;
            double c = 0.5 * (estForwards[i] + estBackwards[i]); // Free energy estimate, closely related to the original Bennett shift constant.
            double cold = 0; // Prior value of c.
            int len0 = eAt[i].length;
            int len1 = eAt[i+1].length;

            if (len0 == 0 || len1 == 0) {
                dGs[i] = c;
                logger.warning(String.format(" Window %d has no snapshots at one end (%d, %d)!", i, len0, len1));
                continue;
            }

            // Ratio of the number of samples: Tinker equivalent: rfrm
            double sampleRatio = ((double) len0) / ((double) len1);

            // Fermi differences.
            double[] fermi1 = new double[len1];
            double[] fermi0 = new double[len0];

            // Ideal gas constant * temperature, or its inverse.
            double rta = Constants.R * temperatures[i];
            double rtb = Constants.R * temperatures[i+1];
            double rtMean = 0.5 * (rta + rtb);
            double invRTA = 1.0 / rta;
            double invRTB = 1.0 / rtb;

            FFXSummaryStatistics s1 = null; // Summary statistics for Fermi differences for the upper half.
            FFXSummaryStatistics s0 = null; // Summary statistics for Fermi differences for the lower half.

            // Each BAR convergence cycle needs to operate on the same set of indices.
            int[] bootstrapSamples0 = null;
            int[] bootstrapSamples1 = null;

            if (randomSamples) {
                bootstrapSamples0 = EstimateBootstrapper.getBootstrapIndices(len1, random);
                bootstrapSamples1 = EstimateBootstrapper.getBootstrapIndices(len0, random);
            }

            int cycleCounter = 0;
            while(!converged) {
                if (randomSamples) {
                    fermiDiffBootstrap(eLow[i + 1], eAt[i + 1], fermi1, len1, c, invRTB, bootstrapSamples1);
                    fermiDiffBootstrap(eHigh[i], eAt[i], fermi0, len0, -c, invRTA, bootstrapSamples0);
                } else {
                    fermiDiffIterative(eLow[i + 1], eAt[i + 1], fermi1, len1, c, invRTB);
                    fermiDiffIterative(eHigh[i], eAt[i], fermi0, len0, -c, invRTA);
                }

                s1 = new FFXSummaryStatistics(fermi1);
                s0 = new FFXSummaryStatistics(fermi0);

                c = rtMean * log(sampleRatio  * (s1.mean / s0.mean)) + cold;
                converged = (Math.abs(c - cold) < tolerance);
                cold = c;

                if (++cycleCounter > MAX_ITERS) {
                    throw new IllegalArgumentException(String.format(" BAR required too many iterations (%d) to converge!", cycleCounter));
                }
            }

            dGs[i] = c;
            cumDG += c;
            double sqFermiMean0 = new FFXSummaryStatistics(Arrays.stream(fermi0).map((double d) -> d*d).toArray()).mean;
            double sqFermiMean1 = new FFXSummaryStatistics(Arrays.stream(fermi1).map((double d) -> d*d).toArray()).mean;

            uncerts[i] = sqrt(uncertCalc(s0.mean, sqFermiMean0, len0) + uncertCalc(s1.mean, sqFermiMean1, len1));
        }

        totDG = cumDG;
        totUncert = sqrt(Arrays.stream(uncerts).map((double d) -> d*d).sum());
    }

    /**
     * Calculates the Fermi function for the differences used in estimating c.
     *
     * f(x) = 1 / (1 + exp(x))
     * x = (e1 - e0 + c) * invRT
     *
     * @param e0         Potential energies to be subtracted.
     * @param e1         Potential energies to be added.
     * @param fermiDiffs Array to be filled with Fermi differences.
     * @param len        Number of energies.
     * @param c          Prior best estimate of the BAR offset/free energy.
     * @param invRT      1.0 / ideal gas constant * temperature.
     */
    private static void fermiDiffIterative(double[] e0, double[] e1, double[] fermiDiffs, int len, double c, double invRT) {
        for (int i = 0; i < len; i++) {
            fermiDiffs[i] = ScalarMath.fermiFunction(invRT * (e1[i] - e0[i] + c));
        }
    }

    /**
     * Calculates the Fermi function for the differences used in estimating c, using bootstrap sampling
     * (choosing random indices w/ replacement rather than scanning through them all).
     *
     * f(x) = 1 / (1 + exp(x))
     * x = (e1 - e0 + c) * invRT
     *
     * @param e0         Potential energies to be subtracted.
     * @param e1         Potential energies to be added.
     * @param fermiDiffs Array to be filled with Fermi differences.
     * @param len        Number of energies.
     * @param c          Prior best estimate of the BAR offset/free energy.
     * @param invRT      1.0 / ideal gas constant * temperature.
     */
    private void fermiDiffBootstrap(double[] e0, double[] e1, double[] fermiDiffs, int len, double c, double invRT, int[] bootstrapSamples) {
        for (int indexI = 0; indexI < len; indexI++) {
            int i = bootstrapSamples[indexI];
            fermiDiffs[indexI] = ScalarMath.fermiFunction(invRT * (e1[i] - e0[i] + c));
        }
    }

    /**
     * Computes one half of the BAR variance.
     *
     * @param meanFermi   Mean Fermi value for either state 0 or state 1.
     * @param meanSqFermi Mean squared Fermi value for either state 0 or state 1.
     * @param len         Number of values.
     * @return            One half of BAR variance.
     */
    private static double uncertCalc(double meanFermi, double meanSqFermi, int len) {
        double sqMeanFermi = meanFermi * meanFermi;
        return ((meanSqFermi - sqMeanFermi) / len) / sqMeanFermi;
    }

    @Override
    public boolean isBidirectional() {
        return true;
    }

    @Override
    public double[] getBinEnergies() {
        return Arrays.copyOf(dGs, nWindows);
    }

    @Override
    public double[] getBinUncertainties() {
        return Arrays.copyOf(uncerts, nWindows);
    }

    @Override
    public int numberOfBins() {
        return nWindows;
    }

    @Override
    public double getFreeEnergy() {
        return totDG;
    }

    @Override
    public double getUncertainty() {
        return totUncert;
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
     * Returns the backwards Zwanzig estimator used to seed BAR.
     *
     * @return A backwards Zwanzig estimator.
     */
    public Zwanzig getInitialBackwardsGuess() {
        return backwardsFEP;
    }
}
