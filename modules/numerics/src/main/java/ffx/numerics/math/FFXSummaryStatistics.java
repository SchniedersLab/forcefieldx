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
package ffx.numerics.math;

import org.apache.commons.math3.distribution.TDistribution;

import java.util.Arrays;

/**
 * The FFXSummaryStatistics class uses online, stable algorithms to calculate summary
 * statistics from double arrays/lists, primarily mean, variance, standard deviation,
 * max, min, sum, and count. This is intended for accuracy and numerical stability,
 * not necessarily for performance (e.g. using Kahan summation).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class FFXSummaryStatistics {
    // Weight-sensitive values.
    public final double mean;
    public final double var;
    public final double varPopulation;
    public final double sd;
    public final double sdPopulation;
    public final double sumWeights;
    // Weight-insensitive values.
    public final double min;
    public final double max;
    public final long count;
    public final double sum;
    public final long dof;
    private final TDistribution tDist;

    private final String descString;

    /**
     * Builds a static view of a running statistic.
     *
     * @param rs Running statistic.
     */
    public FFXSummaryStatistics(FFXRunningStatistics rs) {
        mean = rs.getMean();
        var = rs.getVariance();
        varPopulation = rs.getPopulationVariance();
        sumWeights = rs.getWeight();
        min = rs.getMin();
        max = rs.getMax();
        count = rs.getCount();
        sum = rs.getSum();
        dof = rs.getDOF();
        tDist = (dof > 0) ? new TDistribution(dof) : null;
        sd = rs.getStandardDeviation();
        sdPopulation = rs.getPopulationStandardDeviation();
        descString = String.format(" Summary of %d observations: sum is %17.14g, mean is %17.14g, min is %17.14g, " +
                        "max is %17.14g, and the sum of weights is %17.14g" +
                        "\nSample standard deviation: %17.14g (dof = %d)" +
                        "\nPopulation standard deviation: %17.14g (dof = %d)",
                count, sum, mean, min, max, sumWeights, sd, dof, sdPopulation, count);
    }

    /**
     * Constructs a static summary of a statistic from provided values.
     * Assumes weights are all constant (1.0).
     * Assumes all values will be used.
     *
     * @param values  Values to summarize.
     */
    public FFXSummaryStatistics(double[] values) {
        this(values, null, 0, values.length, 1);
    }

    /**
     * Constructs a static summary of a statistic from provided values.
     * Assumes weights are all constant (1.0).
     * Assumes all values from first to end will be used.
     *
     * @param values  Values to summarize.
     * @param first   First value to use.
     */
    public FFXSummaryStatistics(double[] values, int first) {
        this(values, null, first, values.length, 1);
    }

    /**
     * Constructs a static summary of a statistic from provided values.
     * Assumes weights are all constant (1.0).
     * Assumes a stride of 1.
     *
     * @param values  Values to summarize.
     * @param first   First value to use.
     * @param last    Last value to use.
     */
    public FFXSummaryStatistics(double[] values, int first, int last) {
        this(values, null, first, last, 1);
    }

    /**
     * Constructs a static summary of a statistic from provided values.
     * Assumes weights are all constant (1.0).
     *
     * @param values  Values to summarize.
     * @param first   First value to use.
     * @param last    Last value to use.
     * @param stride  Stride between values used.
     */
    public FFXSummaryStatistics(double[] values, int first, int last, int stride) {
        this(values, null, first, last, stride);
    }

    /**
     * Constructs a static summary of a statistic from provided values.
     *
     * @param values  Values to summarize.
     * @param weights Weights for each value.
     * @param first   First value to use.
     * @param last    Last value to use.
     * @param stride  Stride between values used.
     */
    public FFXSummaryStatistics(double[] values, double[] weights, int first, int last, int stride) {
        if (values == null) {
            throw new IllegalArgumentException(" Cannot have null values!");
        }
        int nVals = values.length;

        if (first < 0 || first > (nVals - 1)) {
            throw new IllegalArgumentException(String.format(" First entry %d was not in valid range 0-%d (0 to length of values - 1)", first, nVals - 1));
        }
        if (last <= first || last > nVals) {
            throw new IllegalArgumentException(String.format(" Last entry %d was not in valid range %d-%d (first+1 to length of values", last, (first + 1), nVals));
        }

        if (weights == null) {
            weights = new double[nVals];
            Arrays.fill(weights, 1.0);
        }

        int tempCount = (last - first);
        if (tempCount % stride == 0) {
            count = tempCount / stride;
        } else {
            count = (tempCount / stride) + 1;
        }
        assert count > 0;

        if (count == 1) {
            mean = values[first];
            var = Double.NaN;
            varPopulation = 0;
            sd = Double.NaN;
            sdPopulation = 0;
            min = mean;
            max = mean;
            sum = mean;
            sumWeights = weights[first];
            dof = 0;
            tDist = null;

            descString = String.format(" Summary of single observation: value is %17.14g", mean);
        } else {
            double meanAcc = 0;
            double varAcc = 0;
            double minAcc = Double.MAX_VALUE;
            double maxAcc = Double.MIN_VALUE;
            double sumAcc = 0;
            double comp = 0;
            double weightAcc = 0;

            for (int i = 0; i < count; i++) {
                int ii = first + (i * stride);
                assert ii < last && ii < nVals;
                double val = values[ii];
                double priorMean = meanAcc;
                double weight = weights[ii];
                weightAcc += weight;

                double y = val - comp;
                double t = sumAcc + y;
                comp = (t - sumAcc) - y;
                sumAcc = t;

                minAcc = Math.min(minAcc, val);
                maxAcc = Math.max(maxAcc, val);

                double invCount = weight / weightAcc;
                meanAcc += ((val - meanAcc) * invCount);
                // TODO: Check correctness of variance accumulation when weight != 1.0.
                varAcc += ((val - priorMean) * (val - meanAcc)) * weight;
            }

            min = minAcc;
            max = maxAcc;
            mean = meanAcc;
            sum = sumAcc;
            sumWeights = weightAcc;
            // Mean via Kahan summation and online mean estimation should be pretty close to each other.
            assert Math.abs(mean - (sum / count)) < 1.0E-11;
            varPopulation = varAcc / count;
            sdPopulation = Math.sqrt(varPopulation);
            dof = count - 1;
            var = varAcc / dof;
            sd = Math.sqrt(var);
            tDist = new TDistribution(dof);

            descString = String.format(" Summary of %d observations: sum is %17.14g, mean is %17.14g, min is %17.14g, " +
                    "max is %17.14g, and the sum of weights is %17.14g" +
                    "\nSample standard deviation: %17.14g (dof = %d)" +
                    "\nPopulation standard deviation: %17.14g (dof = %d)",
                    count, sum, mean, min, max, sumWeights, sd, dof, sdPopulation, count);
        }
    }

    /**
     * Computes a 95% confidence interval based on a Student's T-distribution.
     *
     * @return 95% confidence interval.
     */
    public double confidenceInterval() {
        return confidenceInterval(0.05);
    }

    /**
     * Computes a confidence interval based on a Student's T-distribution.
     *
     * @param alpha Alpha (e.g. 0.05 for a 95% CI).
     * @return Confidence interval.
     */
    public double confidenceInterval(double alpha) {
        if (dof == 0) {
            throw new IllegalArgumentException(" Cannot calculate confidence intervals when there are no degrees of freedom!");
        }
        double critVal = tDist.inverseCumulativeProbability(0.5 * (1.0 - alpha));
        return critVal * sd / Math.sqrt(count);
    }

    public double getMean() {
        return mean;
    }

    public double getVar() {
        return var;
    }

    public double getSd() {
        return sd;
    }

    @Override
    public String toString() {
        return descString;
    }
}
