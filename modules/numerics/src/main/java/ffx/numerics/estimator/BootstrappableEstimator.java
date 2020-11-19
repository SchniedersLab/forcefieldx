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
package ffx.numerics.estimator;

import static java.util.Arrays.stream;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * The BootstrappableEstimator interface describes a StatisticalEstimator which can use bootstrap
 * sampling as an additional method of calculating free energy and uncertainty. These will generally
 * perform non-bootstrap estimation on construction, with estimateDG(true) called to reset the dG and
 * uncertainty estimates using bootstrapping.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public interface BootstrappableEstimator extends StatisticalEstimator {

  /**
   * Return a copy of this Estimator. Each implementation should specify its own type as the return
   * type. Intended to make parallelization of bootstrapping easy.
   *
   * @return A copy of this Estimator.
   */
  BootstrappableEstimator copyEstimator();

  /**
   * Re-calculates free energy and enthalpy.
   *
   * @param randomSamples Whether to draw random samples with replacement (one bootstrap trial).
   */
  void estimateDG(final boolean randomSamples);

  /** Re-calculates free energy and enthalpy without bootstrapping. */
  void estimateDG();

  /**
   * Obtains bootstrap free energy. Default implementation sums by-bin free energies.
   *
   * <p>May be over-ridden by non-sequential estimators like MBAR.
   *
   * @param freeEnergyDifferences By-bin bootstrap results.
   * @return Overall free energy change.
   */
  default double sumBootstrapResults(double[] freeEnergyDifferences) {
    return stream(freeEnergyDifferences).sum();
  }

  default double sumEnthalpyBootstrapResults(double[] totalEnthalpy) {
    return stream(totalEnthalpy).sum();
  }
  /**
   * Obtains bootstrap uncertainty. Default implementation is square root of summed variances.
   *
   * <p>May be over-ridden by non-sequential estimators like MBAR.
   *
   * @param variances Variance (not uncertainty) in by-bin bootstrap results.
   * @return Overall uncertainty.
   */
  default double sumBootstrapUncertainty(double[] variances) {
    return sqrt(stream(variances).sum());
  }

  default double sumBootstrapEnthalpyUncertainty(double[] enthalpyVariances) {
    return sqrt(stream(enthalpyVariances).sum());
  }
}
