
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
package ffx.numerics.math;


import static java.lang.String.format;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import java.util.Random;
import org.apache.commons.math3.distribution.TDistribution;

/**
 * The BootStrapStatistics class uses bootstrapping to estimate statistics from a
 * given population.
 *
 * @author Michael J. Schnieders
 * @author Rose Gogal
 * @since 1.0
 */
public class BootStrapStatistics {

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
     * Constructs a static summary of a statistic from provided values. Assumes weights are all
     * constant (1.0). Assumes all values will be used.
     *
     * @param values Values to summarize.
     */
    public BootStrapStatistics(double[] values) {
        this(values, null, 0, values.length, 1);
    }

    /**
     * Constructs a static summary of a statistic from provided values. Assumes weights are all
     * constant (1.0). Assumes all values from first to end will be used.
     *
     * @param values Values to summarize.
     * @param first First value to use.
     */
    public BootStrapStatistics(double[] values, int first) {
        this(values, null, first, values.length, 1);
    }

    /**
     * Constructs a static summary of a statistic from provided values. Assumes weights are all
     * constant (1.0). Assumes a stride of 1.
     *
     * @param values Values to summarize.
     * @param first First value to use.
     * @param last Last value to use.
     */
    public BootStrapStatistics(double[] values, int first, int last) {
        this(values, null, first, last, 1);
    }

    /**
     * Constructs a static summary of a statistic from provided values. Assumes weights are all
     * constant (1.0).
     *
     * @param values Values to summarize.
     * @param first First value to use.
     * @param last Last value to use.
     * @param stride Stride between values used.
     */
    public BootStrapStatistics(double[] values, int first, int last, int stride) {
        this(values, null, first, last, stride);
    }

    /**
     * Constructs a static summary of a statistic from provided values.
     *
     * @param values Values to summarize.
     * @param weights Weights for each value.
     * @param first First value to use.
     * @param last Last value to use.
     * @param stride Stride between values used.
     */
    public BootStrapStatistics(double[] values, double[] weights, int first, int last, int stride) {
        if (values == null) {
            throw new IllegalArgumentException(" Cannot have null values!");
        }
        int nVals = values.length;

        if (first < 0 || first > (nVals - 1)) {
            throw new IllegalArgumentException(
                    format(" First entry %d was not in valid range 0-%d (0 to length of values - 1)",
                            first, nVals - 1));
        }
        if (last <= first || last > nVals) {
            throw new IllegalArgumentException(
                    format(" Last entry %d was not in valid range %d-%d (first+1 to length of values",
                            last, (first + 1), nVals));
        }

        if (weights == null) {
            weights = new double[nVals];
            fill(weights, 1.0);
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
            descString = format(" Summary of single observation: value is %17.14g", mean);
        } else {
            // Collect Bootstrap Results
            RunningStatistics bootstrapRunningStatistics = new RunningStatistics();
            Random random = new Random();

            for (int bs = 0; bs < count; bs++) {
                // Collect the mean for one Bootstrap round.
                RunningStatistics bootstrapRound = new RunningStatistics();
                for (int i=0; i < count; i++) {
                    int j = random.nextInt((int) count);
                    bootstrapRound.addValue(values[j], weights[j]);
                }
                // Add the mean from this round.
                bootstrapRunningStatistics.addValue(bootstrapRound.getMean());
            }

            min = bootstrapRunningStatistics.getMin();
            max = bootstrapRunningStatistics.getMax();
            mean = bootstrapRunningStatistics.getMean();
            sum = bootstrapRunningStatistics.getSum();
            sumWeights = bootstrapRunningStatistics.getWeight();
            varPopulation = bootstrapRunningStatistics.getPopulationVariance();
            sdPopulation = bootstrapRunningStatistics.getPopulationStandardDeviation();
            dof = bootstrapRunningStatistics.getDOF();
            var = bootstrapRunningStatistics.getVariance();
            sd = bootstrapRunningStatistics.getStandardDeviation();
            tDist = new TDistribution(dof);
            descString =
                    format(
                            " Summary of %d observations: sum is %17.14g, mean is %17.14g, min is %17.14g, "
                                    + "max is %17.14g, and the sum of weights is %17.14g"
                                    + "\nSample standard deviation: %17.14g (dof = %d)"
                                    + "\nPopulation standard deviation: %17.14g (dof = %d)",
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
            throw new IllegalArgumentException(
                    " Cannot calculate confidence intervals when there are no degrees of freedom!");
        }
        double critVal = tDist.inverseCumulativeProbability(0.5 * (1.0 - alpha));
        return critVal * sd / sqrt(count);
    }

    /**
     * The mean.
     *
     * @return Return the mean.
     */
    public double getMean() {
        return mean;
    }

    /**
     * The standard deviation.
     *
     * @return Return the standard deviation.
     */
    public double getSd() {
        return sd;
    }

    /**
     * The variance.
     *
     * @return Return the variance.
     */
    public double getVar() {
        return var;
    }

    /** ${@inheritDoc} */
    @Override
    public String toString() {
        return descString;
    }

    /**
     * Describe the Summary Statistics.
     *
     * @return Return the description.
     */
    public String describe() {
        return format(" Mean: %12.6f +/-%12.6f, Min/Max: %12.6f/%12.6f", mean, sd, min, max);
    }

}
