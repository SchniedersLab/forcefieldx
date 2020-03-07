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

/**
 * The FFXRunningStatistics class uses online, stable algorithms to calculate summary
 * statistics from a source of doubles, primarily mean, variance, standard deviation,
 * max, min, sum, and count. This is intended for accuracy and numerical stability,
 * not necessarily for performance (e.g. using Kahan summation).
 *
 * This is effectively a dynamic version of FFXSummaryStatistics.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class FFXRunningStatistics {
    private double meanAcc = 0;
    private double varAcc = 0;
    private double minAcc = Double.MAX_VALUE;
    private double maxAcc = Double.MIN_VALUE;
    private double sumAcc = 0;
    private double comp = 0;
    private double weightAcc = 0;
    private long count = 0;
    private long dof = -1;

    /**
     * Constructs new running statistics accumulator.
     */
    public FFXRunningStatistics() {
        // Empty constructor; all variables are initialized at definition.
    }

    /**
     * Add a value and update key variables.
     *
     * @param val Value to add.
     */
    public void addValue(double val) {
        addValue(val, 1.0);
    }

    /**
     * Add a value and update key variables.
     *
     * @param val    Value to add.
     * @param weight Weight to give the value.
     */
    public void addValue(double val, double weight) {
        assert Double.isFinite(val);
        assert Double.isFinite(weight);
        assert weight > 0.0;
        ++count;
        ++dof;
        double priorMean = meanAcc;
        weightAcc += weight;
        double y = val - comp;
        double t = sumAcc + y;
        comp = (t - sumAcc) - y;
        sumAcc = t;
        
        minAcc = Math.min(minAcc, val);
        maxAcc = Math.max(maxAcc, val);
        double invCount = 1.0 / weightAcc;
        meanAcc += ((val - meanAcc) * invCount);
        varAcc += ((val - priorMean) * (val - meanAcc)) * weight;
        if (Double.isNaN(varAcc)) {
            throw new IllegalArgumentException(String.format(" Val %.5f w/ wt %.3f resulted in NaN varAcc; current state %s", val, weight, new FFXSummaryStatistics(this).toString()));
        }
    }


    /**
     * Gets the mean as of the last value added.
     *
     * @return Current running mean.
     */
    public double getMean() {
        return meanAcc;
    }

    public double getVariance() {
        return varAcc / ((double) dof);
    }

    public double getPopulationVariance() {
        return varAcc / ((double) count);
    }

    public double getStandardDeviation() {
        return Math.sqrt(getVariance());
    }

    public double getPopulationStandardDeviation() {
        return Math.sqrt(getPopulationVariance());
    }

    public double getMin() {
        return minAcc;
    }

    public double getMax() {
        return maxAcc;
    }

    public double getSum() {
        return sumAcc;
    }

    public double getWeight() {
        return weightAcc;
    }

    public long getCount() {
        return count;
    }

    public long getDOF() {
        return dof;
    }
}
