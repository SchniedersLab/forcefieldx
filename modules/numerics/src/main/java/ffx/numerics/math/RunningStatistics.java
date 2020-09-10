// ******************************************************************************
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
// ******************************************************************************
package ffx.numerics.math;

import static java.lang.Double.isFinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The RunningStatistics class uses online, stable algorithms to calculate summary statistics from a
 * source of doubles, including mean, variance, standard deviation, max, min, sum, and count.
 *
 * <p>This is intended for accuracy and numerical stability, not necessarily for performance (e.g.
 * using Kahan summation).
 *
 * <p>This is effectively a dynamic version of SummaryStatistics.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class RunningStatistics {

  private double meanAcc = 0;
  private double varAcc = 0;
  private double minAcc = Double.MAX_VALUE;
  private double maxAcc = Double.MIN_VALUE;
  private double sumAcc = 0;
  private double comp = 0;
  private double weightAcc = 0;
  private long count = 0;
  private long dof = -1;

  /** Constructs new running statistics accumulator. */
  public RunningStatistics() {
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
   * @param val Value to add.
   * @param weight Weight to give the value.
   */
  public void addValue(double val, double weight) {
    assert isFinite(val);
    assert isFinite(weight);
    assert weight > 0.0;
    ++count;
    ++dof;
    double priorMean = meanAcc;
    weightAcc += weight;
    double y = val - comp;
    double t = sumAcc + y;
    comp = (t - sumAcc) - y;
    sumAcc = t;

    minAcc = min(minAcc, val);
    maxAcc = max(maxAcc, val);
    double invCount = 1.0 / weightAcc;
    meanAcc += ((val - meanAcc) * invCount);
    varAcc += ((val - priorMean) * (val - meanAcc)) * weight;
    if (isNaN(varAcc)) {
      throw new IllegalArgumentException(
          format(" Val %.5f w/ wt %.3f resulted in NaN varAcc; current state %s",
              val, weight, new SummaryStatistics(this).toString()));
    }
  }

  /**
   * Get the count.
   *
   * @return Returns the count.
   */
  public long getCount() {
    return count;
  }

  /**
   * Get the DOF.
   *
   * @return Returns DOF.
   */
  public long getDOF() {
    return dof;
  }

  /**
   * Get the max.
   *
   * @return Returns the max.
   */
  public double getMax() {
    return maxAcc;
  }

  /**
   * Gets the mean as of the last value added.
   *
   * @return Current running mean.
   */
  public double getMean() {
    return meanAcc;
  }

  /**
   * Get the min.
   *
   * @return Returns the min.
   */
  public double getMin() {
    return minAcc;
  }

  /**
   * Get the population standard deviations.
   *
   * @return The population standard deviation.
   */
  public double getPopulationStandardDeviation() {
    return sqrt(getPopulationVariance());
  }

  /**
   * Get the population variance.
   *
   * @return Returns the population variance.
   */
  public double getPopulationVariance() {
    return varAcc / ((double) count);
  }

  /**
   * Get the standard deviation.
   *
   * @return Returns the standard deviation.
   */
  public double getStandardDeviation() {
    return sqrt(getVariance());
  }

  /**
   * Get the sum.
   *
   * @return Returns the sum.
   */
  public double getSum() {
    return sumAcc;
  }

  /**
   * Get the variance.
   *
   * @return Returns the variance.
   */
  public double getVariance() {
    return varAcc / ((double) dof);
  }

  /**
   * Get the weight.
   *
   * @return Returns the weight.
   */
  public double getWeight() {
    return weightAcc;
  }
}
