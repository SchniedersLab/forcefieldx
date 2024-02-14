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

import java.util.Timer;
import ffx.utilities.Constants;
import java.util.Arrays;
import java.util.Random;
import java.util.logging.Logger;
import static ffx.numerics.estimator.EstimateBootstrapper.getBootstrapIndices;
import static ffx.numerics.estimator.Zwanzig.Directionality.BACKWARDS;
import static ffx.numerics.estimator.Zwanzig.Directionality.FORWARDS;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.stream;
import static org.apache.commons.lang3.ArrayFill.fill;
import static org.apache.commons.lang3.math.NumberUtils.min;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * The MultistateBennettAcceptanceRatio class defines a statistical estimator based on a generalization
 * to the Bennett Acceptance Ratio (BAR) method for multiple lambda windows. It requires an input of
 * K X N array of energies (every window at every snap at every lambda value). Two implementations exist
 * and one is commented out. The first is slower and the second is based on pymbar code.
 *
 * This class implements the method discussed in:
 *      Shirts, M. R. and Chodera, J. D. (2008) Statistically optimal analysis of samples from multiple equilibrium
 *      states. J. Chem. Phys. 129, 124105. doi:10.1063/1.2978177
 *
 * @author Matthew J. Speranza
 * @since 1.0
 */
public class MultistateBennettAcceptanceRatio extends SequentialEstimator implements BootstrappableEstimator {
    private static final Logger logger = Logger.getLogger(MultistateBennettAcceptanceRatio.class.getName());

    /**
     * Default BAR convergence tolerance.
     */
    private static final double DEFAULT_TOLERANCE = 1.0E-7;
    /**
     * Maximum number of MBAR iterations.
     */
    private static final int MAX_ITERS = 100;
    /**
     * Number of simulation windows.
     */
    private final int nWindows;
    /**
     * Forward Zwanzig free-energy difference estimates.
     */
    private double[] forwardZwanzig;
    /**
     * Backward Zwanzig free-energy difference estimates.
     */
    private double[] backwardZwanzig;
    /**
     * MBAR free-energy difference estimates.
     */
    private final double[] mbarEstimates;
    /**
     * MBAR free-energy difference uncertainties.
     */
    private final double[] mbarUncertainties;
    /**
     * BAR convergence tolerance.
     */
    private final double tolerance;
    /**
     * Forward Zwanzig instance.
     */
    private Zwanzig forwardsFEP;
    /**
     * Backward Zwanzig instance.
     */
    private Zwanzig backwardsFEP;
    private Random random;
    /**
     * MBAR free-energy estimates at each lambda value.
     */
    private double[] mbarFreeEnergies;
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
    private final double[] mbarEnthalpy;
    private int iter;

    private SeedType seedType;

    public enum SeedType {BAR, ZWANZIG, ZEROS}

    public MultistateBennettAcceptanceRatio(double[] lambdaValues, double[][][] energiesAll, double[] temperature) {
        this(lambdaValues, energiesAll, temperature, DEFAULT_TOLERANCE, SeedType.ZWANZIG);
    }

    public MultistateBennettAcceptanceRatio(double[] lambdaValues, double[][][] energiesAll, double[] temperature,
                                            double tolerance, SeedType seedType) {
        super(lambdaValues, energiesAll, temperature);
        this.tolerance = tolerance;
        this.seedType = seedType;
        nWindows = lambdaValues.length - 1; // Between each lambda value is a window
        // MBAR calculates free energy at each lambda value (only the differences between them have physical significance)
        mbarFreeEnergies = new double[lambdaValues.length];
        mbarEstimates = new double[nWindows];
        mbarUncertainties = new double[nWindows];
        mbarEnthalpy = new double[nWindows];
        random = new Random();
        estimateDG();
    }

    private void seedEnergies() {
        switch (seedType) {
            case SeedType.BAR:
                try {
                    SequentialEstimator barEstimator = new BennettAcceptanceRatio(lamValues, eLow, eAt, eHigh, temperatures);
                    mbarFreeEnergies[0] = 0.0;
                    double[] barEstimates = barEstimator.getBinEnergies();
                    for (int i = 0; i < nWindows; i++) {
                        mbarFreeEnergies[i + 1] = mbarFreeEnergies[i] + barEstimates[i];
                    }
                    break;
                } catch (IllegalArgumentException e) {
                    logger.warning(" BAR failed to converge. Using Zwanzig seed.");
                    seedType = SeedType.ZWANZIG;
                    seedEnergies();
                    return;
                }
            case SeedType.ZWANZIG:
                forwardsFEP = new Zwanzig(lamValues, eLow, eAt, eHigh, temperatures, FORWARDS);
                backwardsFEP = new Zwanzig(lamValues, eLow, eAt, eHigh, temperatures, BACKWARDS);
                forwardZwanzig = forwardsFEP.getBinEnergies();
                backwardZwanzig = backwardsFEP.getBinEnergies();
                mbarFreeEnergies[0] = 0.0;
                for (int i = 0; i < nWindows; i++) {
                    mbarFreeEnergies[i + 1] = mbarFreeEnergies[i] + .5 * (forwardZwanzig[i] + backwardZwanzig[i]);
                }
                break;
            case SeedType.ZEROS:
                break;
            default:
                throw new IllegalArgumentException("Seed type not supported");
        }
    }

    @Override
    public void estimateDG() {
        estimateDG(false);
    }

    // Slow implementation
    /*
    @Override
    public void estimateDG(boolean randomSamples) {
        double[] prevMBAR;
        int iter = 0;

        // Precompute values for each lambda window
        double[] rtValues = new double[temperatures.length];
        double[] invRTValues = new double[temperatures.length];
        for(int i = 0; i < temperatures.length; i++){
            rtValues[i] = Constants.R * temperatures[i];
            invRTValues[i] = 1.0 / rtValues[i];
        }
        // Find the smallest length of eAllFlat array (look for INF values)
        int minSnaps = eAllFlat[0].length;

        logger.info("Using " + minSnaps + " snapshots for MBAR calculations. This is K (runs) * N (min(NumSnaps)).");
        int numSnaps = min(eAllFlat[0].length, minSnaps);
        int totalSnaps = numSnaps * nWindows;

        // Sample random snapshots from each window
        int[] indices = new int[numSnaps];
        if(randomSamples){
            indices = getBootstrapIndices(numSnaps, random);
        } else {
            for(int i = 0; i < numSnaps; i++){
                indices[i] = i;
            }
        }

        // Self-consistent iteration is used over Newton-Raphson because it is more stable & simpler to implement
        do {
            prevMBAR = copyOf(mbarFreeEnergies, mbarFreeEnergies.length);
            // Calculate new MBAR free energy estimates. Loosely based on pymbar code.
            double[] tempMBAR = new double[mbarFreeEnergies.length];
            for(int state = 0; state < mbarFreeEnergies.length; state++) { // For each lambda value
                double[] numeratorTerms = new double[numSnaps];
                double[] denominatorTerms = new double[numSnaps];
                for (int n = 0; n < numSnaps; n++) {
                    numeratorTerms[indices[n]] = -eAllFlat[state][indices[n]] * invRTValues[state];
                    double[] temp = new double[lamValues.length];
                    double[] tempB = new double[lamValues.length];
                    for (int k = 0; k < lamValues.length; k++) { // Denominator sum over K different lambda evaluations
                        temp[k] = prevMBAR[k] + (-eAllFlat[k][indices[n]] * invRTValues[k]);
                        tempB[k] = numSnaps;
                    }
                    denominatorTerms[n] = logSumExp(temp, tempB);
                }

                // Subtract all dt from all nt with a stream
                double[] temp = new double[numSnaps];
                for(int i = 0; i < numSnaps; i++){
                    temp[i] =  -denominatorTerms[i] + numeratorTerms[i];
                }
                double sum = -1.0 * logSumExp(temp);
                tempMBAR[state] = sum;
            }
            // Constrain f1=0 over the course of iterations to prevent uncontrolled growth in magnitude
            double f1 = tempMBAR[0];
            for (int i = 0; i < mbarFreeEnergies.length; i++) {
                mbarFreeEnergies[i] = tempMBAR[i] - f1;
            }
            iter++;
        } while(!converged(prevMBAR));

        for(int i = 0; i < mbarFreeEnergies.length; i++){
            mbarFreeEnergies[i] = mbarFreeEnergies[i] * rtValues[i];
        }

        logger.info(" MBAR converged after " + iter + " iterations.");
        for(int i = 0; i < nWindows; i++){
            mbarEstimates[i] = mbarFreeEnergies[i+1] - mbarFreeEnergies[i];
            logger.info(" MBAR free energy difference estimate for window " + i + " -> " + (i+1) + ": " + mbarEstimates[i]);
        }
        totalMBAREstimate = stream(mbarEstimates).sum();
        logger.info(" Total MBAR free energy difference estimate: " + totalMBAREstimate);
    }
     */

    /**
     * Implementation based on pymbar code. Precomputes values & uses the mixture distribution
     * idea for faster calculations.
     */
    @Override
    public void estimateDG(boolean randomSamples) {
        Arrays.fill(mbarFreeEnergies, 0.0); // Bootstrap needs resetting
        seedEnergies();
        // Throw error if MBAR contains NaNs or Infs
        if (stream(mbarFreeEnergies).anyMatch(Double::isInfinite) || stream(mbarFreeEnergies).anyMatch(Double::isNaN)) {
            throw new IllegalArgumentException("MBAR contains NaNs or Infs after seeding.");
        }
        double[] prevMBAR;
        iter = 0;

        // Precompute values for each lambda window
        double[] rtValues = new double[temperatures.length];
        double[] invRTValues = new double[temperatures.length];
        for (int i = 0; i < temperatures.length; i++) {
            rtValues[i] = Constants.R * temperatures[i];
            invRTValues[i] = 1.0 / rtValues[i];
        }
        // Find the smallest length of eAllFlat array (look for INF values)
        int minSnaps = Integer.MAX_VALUE;
        for (double[] trajectory : eAllFlat) {
            minSnaps = min(minSnaps, trajectory.length);
        }
        int numSnaps = min(eAllFlat[0].length, minSnaps);
        int totalSnaps = numSnaps * eAllFlat.length;

        // Sample random snapshots from each window
        int[][] indices = new int[mbarFreeEnergies.length][numSnaps];
        if (randomSamples) {
            int[] randomIndices = getBootstrapIndices(numSnaps, random);
            for (int i = 0; i < mbarFreeEnergies.length; i++) {
                indices[i] = randomIndices; // Use the same snapshots in each lambda window
            }
        } else {
            for (int i = 0; i < numSnaps; i++) {
                for (int j = 0; j < mbarFreeEnergies.length; j++) {
                    indices[j][i] = i;
                }
            }
        }

        // Precompute energies since they don't change
        double[][] u_kn = new double[mbarFreeEnergies.length][numSnaps];
        double[] N_k = new double[mbarFreeEnergies.length];
        for (int state = 0; state < mbarFreeEnergies.length; state++) { // For each lambda value
            for (int n = 0; n < numSnaps; n++) {
                u_kn[state][n] = eAllFlat[state][indices[state][n]] * invRTValues[state];
            }
            N_k[state] = ((double) numSnaps) / totalSnaps;
        }
        // Self-consistent iteration is used  because it is more stable & simpler to implement than grad descent
        double omega = 1.5;
        do {
            prevMBAR = copyOf(mbarFreeEnergies, mbarFreeEnergies.length);
            mbarFreeEnergies = selfConsistentUpdate(u_kn, N_k, mbarFreeEnergies);
            // Apply SOR
            for (int i = 0; i < mbarFreeEnergies.length; i++) {
                mbarFreeEnergies[i] = omega * mbarFreeEnergies[i] + (1-omega) * prevMBAR[i];
            }
            // Throw error if MBAR contains NaNs or Infs
            if (stream(mbarFreeEnergies).anyMatch(Double::isInfinite) || stream(mbarFreeEnergies).anyMatch(Double::isNaN)) {
                throw new IllegalArgumentException("MBAR contains NaNs or Infs after iteration " + iter);
            }
            iter++;
        } while (!converged(prevMBAR));
        logger.info(" MBAR converged after " + iter + " iterations with omega " + omega + ".");

        // Convert to kcal/mol & calculate differences/sums
        for (int i = 0; i < mbarFreeEnergies.length; i++) {
            mbarFreeEnergies[i] = mbarFreeEnergies[i] * rtValues[i];
        }

        for (int i = 0; i < nWindows; i++) {
            mbarEstimates[i] = mbarFreeEnergies[i + 1] - mbarFreeEnergies[i];
        }
        totalMBAREstimate = stream(mbarEstimates).sum();
    }

    public static double[] selfConsistentUpdate(double[][] u_kn, double[] N_k, double[] f_k) {
        int nStates = f_k.length;
        double[] updatedF_k = new double[nStates];
        double[] log_denom_n = new double[u_kn[0].length];
        double[][] logDiff = new double[u_kn.length][u_kn[0].length];
        double maxLogDiff = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < u_kn[0].length; i++) {
            double[] temp = new double[nStates];
            double maxTemp = Double.NEGATIVE_INFINITY;
            for (int j = 0; j < nStates; j++) {
                temp[j] = f_k[j] - u_kn[j][i];
                if (temp[j] > maxTemp) {
                    maxTemp = temp[j];
                }
            }
            log_denom_n[i] = logSumExp(temp, N_k, maxTemp);
            for (int j = 0; j < nStates; j++) {
                logDiff[j][i] = -log_denom_n[i] - u_kn[j][i];
                if (logDiff[j][i] > maxLogDiff) {
                    maxLogDiff = logDiff[j][i];
                }
            }
        }

        for (int i = 0; i < nStates; i++) {
            updatedF_k[i] = -1.0 * logSumExp(logDiff[i], maxLogDiff);
        }

        // Constrain f1=0 over the course of iterations to prevent uncontrolled growth in magnitude
        double norm = updatedF_k[0];
        for (int i = 0; i < nStates; i++) {
            updatedF_k[i] = updatedF_k[i] - norm;
        }

        return updatedF_k;
    }

    /**
     * Calculates the log of the sum of the exponentials of the given values.
     *
     * @param values
     * @return the sum
     */
    private static double logSumExp(double[] values, double max) {
        double[] b = fill(new double[values.length], 1.0);
        return logSumExp(values, b, max);
    }

    /**
     * Calculates the log of the sum of the exponentials of the given values.
     *
     * @param values
     * @return the sum
     */
    private static double logSumExp(double[] values, double[] b, double max) {
        // ChatGPT mostly wrote this and I tweaked it to match more closely with scipy's logsumexp implementation
        // Find the maximum value in the array.
        assert values.length == b.length : "values and b must be the same length";

        // Subtract the maximum value from each value in the array, exponentiate the result, and add up these values.
        double sum = 0.0;
        for (int i = 0; i < values.length; i++) {
            sum += b[i] * exp(values[i] - max);
        }

        // Take the natural logarithm of the sum and add the maximum value back in.
        return max + log(sum);
    }

    private boolean converged(double[] prevMBAR) {
        double[] differences = new double[prevMBAR.length];
        for (int i = 0; i < prevMBAR.length; i++) {
            differences[i] = abs(prevMBAR[i] - mbarFreeEnergies[i]);
        }
        return stream(differences).allMatch(d -> d < tolerance);
    }

    public BennettAcceptanceRatio getBAR() {
        return new BennettAcceptanceRatio(lamValues, eLow, eAt, eHigh, temperatures);
    }

    @Override
    public MultistateBennettAcceptanceRatio copyEstimator() {
        return new MultistateBennettAcceptanceRatio(lamValues, eAll, temperatures, tolerance, seedType);
    }

    @Override
    public double[] getBinEnergies() {
        return mbarEstimates;
    }

    @Override
    public double[] getBinUncertainties() {
        return mbarUncertainties;
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
        return nWindows;
    }

    @Override
    public double[] getBinEnthalpies() {
        return mbarEnthalpy;
    }


    public static class HarmonicOscillatorsTestCase {

        private double beta;
        private double[] O_k;
        private int n_states;
        private double[] K_k;

        /**
         * Constructor for HarmonicOscillatorsTestCase
         * @param O_k array of equilibrium positions
         * @param K_k array of spring constants
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

            if (observable.equals("position")) {
                return analyticalMeans();
            } else if (observable.equals("potential energy")) {
                for (int i = 0; i < n_states; i++) {
                    result[i] = 0.5 / beta;
                }
            } else if (observable.equals("position^2")) {
                for (int i = 0; i < n_states; i++) {
                    result[i] = 1.0 / (beta * K_k[i]) + Math.pow(O_k[i], 2);
                }
            } else if (observable.equals("RMS displacement")) {
                return analyticalStandardDeviations();
            }

            return result;
        }

        public double[] analyticalFreeEnergies() {
            int subtractComponentIndex = 0;
            double[] fe = new double[n_states];
            double subtract = 0.0;
            for (int i = 0; i < n_states; i++) {
                fe[i] = -0.5 * Math.log(2 * Math.PI / (beta * K_k[i]));
                if(i == 0){
                    subtract = fe[subtractComponentIndex];
                }
                fe[i] -= subtract;
            }
            return fe;
        }

        public double[] analyticalEntropies(int subtractComponent) {
            double[] entropies = new double[n_states];
            double[] potentialEnergy = analyticalObservable("potential energy");
            double[] freeEnergies = analyticalFreeEnergies();

            for (int i = 0; i < n_states; i++) {
                entropies[i] = potentialEnergy[i] - freeEnergies[i];
            }

            return entropies;
        }

        /**
         * Sample from harmonic oscillator w/ gaussian & std
         * @param N_k number of samples per state
         * @param mode only u_kn -> return K x N_tot matrix where u_kn[k,n] is reduced potential of sample n evaluated at state k
         * @return u_kn[k,n] is reduced potential of sample n evaluated at state k
         */
        public Object[] sample(int[] N_k, String mode) {
            Random random = new Random();

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

                // Number of samples
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
            }

            // Setting corrections
            if ("u_kn".equals(mode)) {
                return new Object[]{x_n, u_kn, N_k, s_n};
            } else if ("u_kln".equals(mode)) {
                return new Object[]{x_n, u_kln, N_k, s_n};
            } else {
                throw new IllegalArgumentException("Unknown mode: " + mode);
            }
        }

        public static Object[] evenlySpacedOscillators(
                int n_states, int n_samplesPerState, double lower_O_k, double upper_O_k,
                double lower_K_k, double upper_K_k, Long seed) {
            Random random = new Random(seed);

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
            Object[] result = testCase.sample(N_k, "u_kn");

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
            Object[] sampleResult = testCase.sample(N_k, setting);

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

    public static void main(String[] args) {
        // Define parameters for the harmonic oscillators
        //double[] O_k = {0, .1, .2, .3, .4, .5};
        // double[] O_k = {0, 1, 2, 3, 4};
        // double[] K_k = {1, 2, 4, 8, 16};
        // int[] N_k = {10, 20, 30, 40, 50};

        double[] O_k = {0, 1, 2, 3, 4};
        double[] K_k = {1, 5, 10, 20, 50};
        int[] N_k = {10000, 10000, 10000, 10000, 10000};
        double beta = 1.0;

        Long seed = System.currentTimeMillis();

        // Create an instance of HarmonicOscillatorsTestCase
        HarmonicOscillatorsTestCase testCase = new HarmonicOscillatorsTestCase(O_k, K_k, beta);

        // Generate sample data
        String setting = "u_kln";
        System.out.print("Generating sample data... ");
        Object[] sampleResult = testCase.sample(N_k, setting);
        System.out.println("done. \n");
        double[][][] u_kln = (double[][][]) sampleResult[1];
        double[] temps = {1 / Constants.R};

        // Create an instance of MultistateBennettAcceptanceRatio
        System.out.print("Creating MBAR instance and estimateDG() with standard tol & ZERO seeding to reduce dependancy issues... ");
        MultistateBennettAcceptanceRatio mbar = new MultistateBennettAcceptanceRatio(O_k, u_kln, temps, 1.0E-7, SeedType.ZEROS);
        double[] mbarFEEstimates = Arrays.copyOf(mbar.mbarFreeEnergies, mbar.mbarFreeEnergies.length);

        EstimateBootstrapper bootstrapper = new EstimateBootstrapper(mbar);
        bootstrapper.bootstrap(10);
        System.out.println("done. \n");

        // Get the analytical free energy differences
        double[] analyticalFreeEnergies = testCase.analyticalFreeEnergies();
        // Calculate the free energy differences from analyticalFreeEnergies & log
        double[] analyticalEstimates = new double[analyticalFreeEnergies.length - 1];
        for (int i = 0; i < analyticalEstimates.length; i++) {
            analyticalEstimates[i] = analyticalFreeEnergies[i + 1] - analyticalFreeEnergies[i];
        }
        // Calculate the error
        double[] error = new double[analyticalFreeEnergies.length];
        for (int i = 0; i < error.length; i++) {
            error[i] = - mbarFEEstimates[i] + analyticalFreeEnergies[i];
        }

        // Compare the calculated free energy differences with the analytical ones
        System.out.println("MBAR Free Energies:       " + Arrays.toString(mbarFEEstimates));
        System.out.println("Analytical Free Energies: " + Arrays.toString(analyticalFreeEnergies));
        System.out.println("Free Energy Error:        " + Arrays.toString(error));
        System.out.println();

        // Get the calculated free energy differences
        double[] mbarBootstrappedEstimates = bootstrapper.getFE();
        // Calculate the error
        double[] errors = new double[mbarBootstrappedEstimates.length];
        for (int i = 0; i < errors.length; i++) {
            errors[i] = - mbarBootstrappedEstimates[i] + analyticalEstimates[i];
        }

        System.out.println("MBAR Bootstrapped Estimates: " + Arrays.toString(mbarBootstrappedEstimates));
        System.out.println("Analytical Estimates:        " + Arrays.toString(analyticalEstimates));
        System.out.println("Free Energy Error:           " + Arrays.toString(errors));
    }
}
