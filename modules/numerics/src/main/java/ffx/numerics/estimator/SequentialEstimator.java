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

import static java.lang.System.arraycopy;
import static java.util.Arrays.copyOf;
import static java.util.Arrays.fill;
import static java.util.Arrays.stream;

/**
 * The SequentialEstimator abstract class defines a statistical estimator based on perturbative
 * potential energy differences between adjacent windows (e.g. exponential free energy perturbation,
 * Bennett Acceptance Ratio, etc).
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public abstract class SequentialEstimator implements StatisticalEstimator {

  protected final double[] lamVals;
  protected final double[][] eLow;
  protected final double[][] eAt;
  protected final double[][] eHigh;
  protected final double[] temperatures;
  protected final int nTrajectories;

  /**
   * The SequentialEstimator constructor largely just copies its parameters into local variables.
   * Most arrays are deep copied. The temperature array can be of length 1 if all elements are meant
   * to be the same temperature.
   *
   * <p>The first dimension of the energies arrays corresponds to the lambda values/windows. The
   * second dimension (can be of uneven length) corresponds to potential energies of snapshots
   * sampled from that lambda value, calculated either at that lambda value, the lambda value below,
   * or the lambda value above. The arrays energiesLow[0] and energiesHigh[n-1] is expected to be all
   * NaN.
   *
   * @param lambdaValues Values of lambda dynamics was run at.
   * @param energiesLow Potential energies of trajectory L at lambda L-dL.
   * @param energiesAt Potential energies of trajectory L at lambda L.
   * @param energiesHigh Potential energies of trajectory L at lambda L+dL.
   * @param temperature Temperature each lambda window was run at (single-element indicates
   *     identical temperatures).
   */
  public SequentialEstimator(double[] lambdaValues, double[][] energiesLow, double[][] energiesAt,
      double[][] energiesHigh, double[] temperature) {
    nTrajectories = lambdaValues.length;

    assert stream(energiesLow[0]).allMatch(Double::isNaN)
        && stream(energiesHigh[nTrajectories - 1]).allMatch(Double::isNaN);

    assert nTrajectories == energiesAt.length
        && nTrajectories == energiesLow.length
        && nTrajectories == energiesHigh.length
        : "One of the energy arrays is of the incorrect length in the first dimension!";

    this.lamVals = copyOf(lambdaValues, nTrajectories);
    temperatures = new double[nTrajectories];
    if (temperature.length == 1) {
      fill(temperatures, temperature[0]);
    } else {
      arraycopy(temperature, 0, temperatures, 0, nTrajectories);
    }

    // Just in case, copy the arrays rather than storing them as provided.
    eLow = new double[nTrajectories][];
    eAt = new double[nTrajectories][];
    eHigh = new double[nTrajectories][];
    for (int i = 0; i < nTrajectories; i++) {
      eLow[i] = copyOf(energiesLow[i], energiesLow[i].length);
      eAt[i] = copyOf(energiesAt[i], energiesAt[i].length);
      eHigh[i] = copyOf(energiesHigh[i], energiesHigh[i].length);
    }
  }
}
