package ffx.numerics.estimator;

import ffx.numerics.math.FFXRunningStatistics;
import ffx.numerics.math.FFXSummaryStatistics;

import java.util.Arrays;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;

public class EstimateBootstrapper {
    private static final Logger logger = Logger.getLogger(EstimateBootstrapper.class.getName());
    private static final long DEFAULT_LOG_INTERVAL = 25;

    private final BootstrappableEstimator estimate;
    private final int nWindows;
    private final FFXSummaryStatistics[] bootstrapResults;

    public EstimateBootstrapper(BootstrappableEstimator estimator) {
        this.estimate = estimator;
        nWindows = estimate.numberOfBins();
        bootstrapResults = new FFXSummaryStatistics[nWindows];
    }

    public void bootstrap(long trials) {
        bootstrap(trials, DEFAULT_LOG_INTERVAL);
    }

    public void bootstrap(long trials, long logInterval) {
        FFXRunningStatistics[] windows = new FFXRunningStatistics[nWindows];
        for (int i = 0; i < nWindows; i++) {
            windows[i] = new FFXRunningStatistics();
        }

        // TODO: Parallelize this loop, because so long as we construct duplicate estimators/accumulators, it should be trivially parallelizable.
        for (long i = 0; i < trials; i++) {
            if ((i+1) % logInterval == 0) {
                logger.info(String.format(" Bootstrap Trial %d", i+1));
            }

            estimate.estimateDG(true);
            double[] fe = estimate.getBinEnergies();
            for (int j = 0; j < nWindows; j++) {
                windows[j].addValue(fe[j]);
            }
        }

        for (int i = 0; i < nWindows; i++) {
            bootstrapResults[i] = new FFXSummaryStatistics(windows[i]);
        }
    }

    public double[] getFE() {
        return Arrays.stream(bootstrapResults).mapToDouble(FFXSummaryStatistics::getMean).toArray();
    }

    public double[] getUncertainty() {
        return Arrays.stream(bootstrapResults).mapToDouble(FFXSummaryStatistics::getSd).toArray();
    }

    public double[] getVariance() {
        return Arrays.stream(bootstrapResults).mapToDouble(FFXSummaryStatistics::getVar).toArray();
    }

    // The following four methods are included primarily in case of non-sequential estimators (e.g. MBAR) which aren't a simple summation.

    public double getTotalFE() {
        return getTotalFE(getFE());
    }

    public double getTotalFE(double[] fe) {
        return estimate.sumBootstrapResults(fe);
    }

    public double getTotalUncertainty() {
        return getTotalUncertainty(getVariance());
    }

    public double getTotalUncertainty(double[] var) {
        return estimate.sumBootstrapUncertainty(var);
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
        int[] indices = random.ints(length, 0, length).toArray();
        long distinctVals = Arrays.stream(indices).distinct().count();
        while (distinctVals <= minDistinct) {
            logger.info(" Regenerating array: only " + distinctVals + " distinct values found.");
            indices = random.ints(length, 0, length).toArray();
            distinctVals = Arrays.stream(indices).distinct().count();
        }
        return indices;
    }
}
