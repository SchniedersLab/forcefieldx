package ffx.numerics.estimator;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import static java.util.Arrays.stream;

import ffx.numerics.math.RunningStatistics;
import ffx.numerics.math.SummaryStatistics;

public class EstimateBootstrapper {
    private static final Logger logger = Logger.getLogger(EstimateBootstrapper.class.getName());
    private static final long DEFAULT_LOG_INTERVAL = 25;

    private final BootstrappableEstimator estimate;
    private final int nWindows;
    private final SummaryStatistics[] bootstrapResults;

    public EstimateBootstrapper(BootstrappableEstimator estimator) {
        this.estimate = estimator;
        nWindows = estimate.numberOfBins();
        bootstrapResults = new SummaryStatistics[nWindows];
    }

    public void bootstrap(long trials) {
        bootstrap(trials, DEFAULT_LOG_INTERVAL);
    }

    public void bootstrap(long trials, long logInterval) {
        RunningStatistics[] windows = new RunningStatistics[nWindows];
        for (int i = 0; i < nWindows; i++) {
            windows[i] = new RunningStatistics();
        }

        // TODO: Parallelize this loop, because so long as we construct duplicate estimators/accumulators, it should be trivially parallelizable.
        for (long i = 0; i < trials; i++) {
            if ((i + 1) % logInterval == 0) {
                logger.info(String.format(" Bootstrap Trial %d", i + 1));
            }

            estimate.estimateDG(true);
            double[] fe = estimate.getBinEnergies();
            for (int j = 0; j < nWindows; j++) {
                windows[j].addValue(fe[j]);
            }
        }

        for (int i = 0; i < nWindows; i++) {
            bootstrapResults[i] = new SummaryStatistics(windows[i]);
        }
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
        long distinctVals = stream(indices).distinct().count();
        while (distinctVals <= minDistinct) {
            logger.info(" Regenerating array: only " + distinctVals + " distinct values found.");
            indices = random.ints(length, 0, length).toArray();
            distinctVals = stream(indices).distinct().count();
        }
        return indices;
    }

    public double[] getFE() {
        return stream(bootstrapResults).mapToDouble(SummaryStatistics::getMean).toArray();
    }

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

    public double[] getUncertainty() {
        return stream(bootstrapResults).mapToDouble(SummaryStatistics::getSd).toArray();
    }

    public double[] getVariance() {
        return stream(bootstrapResults).mapToDouble(SummaryStatistics::getVar).toArray();
    }
}
