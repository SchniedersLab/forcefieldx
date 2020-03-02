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
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * The Bennett Acceptance Ratio class implements the Bennett Acceptance Ratio (BAR)
 * statistical estimator.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class BennettAcceptanceRatio extends SequentialEstimator {
    private static final Logger logger = Logger.getLogger(BennettAcceptanceRatio.class.getName());

    private final int nWindows;
    private final double[] dGs;
    private final double[] uncerts;
    private final double totDG;
    private final double totUncert;
    private static final double DEFAULT_TOLERANCE = 1.0E-10;
    private static final int MAX_ITERS = 100;

    // Hang onto these in case the end-user wants them?
    private final SequentialEstimator forwardsFEP;
    private final SequentialEstimator backwardsFEP;

    public BennettAcceptanceRatio(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt, double[][] energiesHigh, double[] temperature) {
        this(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, DEFAULT_TOLERANCE);
    }
    
    public BennettAcceptanceRatio(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt, double[][] energiesHigh, double[] temperature, double tolerance) {
        super(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature);
        // Used to seed an initial guess.
        forwardsFEP = new Zwanzig(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, Zwanzig.Directionality.FORWARDS);
        backwardsFEP = new Zwanzig(lambdaValues, energiesLow, energiesAt, energiesHigh, temperature, Zwanzig.Directionality.BACKWARDS);

        nWindows = nTrajectories - 1;
        double[] estForwards = forwardsFEP.getWindowEnergies();
        double[] estBackwards = backwardsFEP.getWindowEnergies();

        dGs = new double[nWindows];
        uncerts = new double[nWindows];

        double cumDG = 0;
        double cumUncert = 0;
        for (int i = 0; i < nWindows; i++) {
            boolean converged = false;
            double lastG = 0.5 * (estBackwards[i] + estForwards[i]);
            double dG = lastG;
            int len0 = eAt[i].length;
            int len1 = eAt[i+1].length;
            double sampleRatio = ((double) len1) / ((double) len0);
            double logSampleRatio = FastMath.log(sampleRatio);
            double c = estimateC(FastMath.exp(dG), sampleRatio);

            double[] diffsAbove = new double[len1];
            double[] diffsBelow = new double[len0];
            double[] fermiAbove = new double[len1];
            double[] fermiBelow = new double[len0];
            FFXSummaryStatistics numeratorStats;
            FFXSummaryStatistics denominatorStats;
            double meanFermiAbove = 1;
            double meanFermiBelow = 1;

            logger.fine(String.format(" Window %d initial guess: dG=%.8f, c=%.8f from sampleRatio %.4f fore %.6f back %.6f", i, dG, c, sampleRatio, estForwards[i], estBackwards[i]));

            int cycleCounter = 0;
            while (!converged) {
                getScalarDiffs(eLow[i+1], eAt[i+1], diffsAbove, c);
                getScalarDiffs(eHigh[i], eAt[i], diffsBelow, -c);
                fermiDiffs(diffsAbove, fermiAbove, len1);
                fermiDiffs(diffsBelow, fermiBelow, len0);

                numeratorStats = new FFXSummaryStatistics(fermiAbove);
                denominatorStats = new FFXSummaryStatistics(fermiBelow);
                meanFermiAbove = numeratorStats.mean;
                meanFermiBelow = denominatorStats.mean;

                double fermiRatio = meanFermiAbove / meanFermiBelow;
                double partitionRatio = fermiRatio * FastMath.exp(c);
                c = estimateC(partitionRatio, sampleRatio);
                dG = estimateDG(fermiRatio, c, logSampleRatio);

                converged = Math.abs(dG - lastG) < tolerance;
                logger.info(String.format(" Numerator mean: %.6f. Denominator mean: %.6f", numeratorStats.mean, denominatorStats.mean));
                logger.info(String.format(" Convergence cycle %d: dG is %.6f, c is %.6f, convergence found %b", ++cycleCounter, dG, c, converged));
                lastG = dG;
                if (cycleCounter > MAX_ITERS) {
                    throw new IllegalArgumentException(" BAR did not converge in " + MAX_ITERS + " iterations!");
                }
            }

            logger.info("\n");

            dGs[i] = dG;

            /*
             * NOTE: According to Pymbar (https://github.com/choderalab/pymbar/blob/master/pymbar/bar.py),
             * there is a typo in the original Bennett publication equation 10a, where the second denominator
             * should be n1*<f>_1^2, not n0*<f>_0^2
             */
            double[] sqFermiAbove = getSquareFermiDiffs(diffsAbove);
            double[] sqFermiBelow = getSquareFermiDiffs(diffsBelow);
            double meanSqFermiAbove = new FFXSummaryStatistics(sqFermiAbove).mean;
            double meanSqFermiBelow = new FFXSummaryStatistics(sqFermiBelow).mean;

            double sqMeanFermiBelow = meanFermiBelow * meanFermiBelow;
            double sqMeanFermiAbove = meanFermiAbove * meanFermiAbove;

            double lhs = (meanSqFermiBelow - sqMeanFermiBelow) / (len0 * sqMeanFermiBelow);
            double rhs = (meanSqFermiAbove - sqMeanFermiAbove) / (len1 * sqMeanFermiAbove);
            double var = lhs + rhs;

            uncerts[i] = Math.sqrt(var);
            cumDG += dG;
            cumUncert += uncerts[i];
        }

        totDG = cumDG;
        totUncert = cumUncert;
    }

    private static double estimateC(double partitionRatio, double sampleRatio) {
        return FastMath.log(partitionRatio, sampleRatio);
    }

    private static double estimateDG(double fermiRatio, double c, double logSampleRatio) {
        return FastMath.log(fermiRatio) + c - logSampleRatio;
    }

    private static void getScalarDiffs(double[] u0, double[] u1, double[] diffs, double c) {
        int len = u0.length;
        assert len == u1.length;
        for (int i = 0; i < len; i++) {
            diffs[i] = u0[i] - u1[i] + c;
        }
    }

    private static void fermiDiffs(double[] scalarDiffs, double[] fermiDiffs, int len) {
        for (int i = 0; i < len; i++) {
            fermiDiffs[i] = ScalarMath.fermiFunction(scalarDiffs[i]);
        }
    }

    private static double[] getSquareFermiDiffs(double[] scalarDiffs) {
        return Arrays.stream(scalarDiffs).map(ScalarMath::fermiFunction).map((double d) -> d*d).toArray();
    }

    @Override
    public boolean isBidirectional() {
        return true;
    }

    @Override
    public double[] getWindowEnergies() {
        return Arrays.copyOf(dGs, nWindows);
    }

    @Override
    public double[] getWindowUncertainties() {
        return Arrays.copyOf(uncerts, nWindows);
    }

    @Override
    public double getFreeEnergy() {
        return totDG;
    }

    @Override
    public double getUncertainty() {
        return totUncert;
    }

    public SequentialEstimator getInitialForwardsGuess() {
        return forwardsFEP;
    }

    public SequentialEstimator getInitialBackwardsGuess() {
        return backwardsFEP;
    }
}
