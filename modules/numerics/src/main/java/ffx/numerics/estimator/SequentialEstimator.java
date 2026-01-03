// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import java.util.ArrayList;

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
 * @author Matthew J. Speranza
 * @since 1.0
 */
public abstract class SequentialEstimator implements StatisticalEstimator {

  /**
   * The lambda values at which the samples were collected.
   */
  protected final double[] lamValues;
  /**
   * The number of states from which samples were collected.
   */
  protected final int nStates;
  /**
   * The potential energy of each snapshot at lambda - dL.
   */
  protected final double[][] eLambdaMinusdL;
  /**
   * The potential energy of each snapshot at lambda.
   */
  protected final double[][] eLambda;
  /**
   * The potential energy of each snapshot at lambda + dL.
   */
  protected final double[][] eLambdaPlusdL;
  /**
   * The potential energies of the snapshots at all other lambda values.
   * <p>
   * eAll[lambdaWindow][perturbations][energies]
   */
  protected final double[][][] eAll;
  /**
   * The potential energies of the snapshots at all other lambda values.
   * <p>
   * eAllFlat[lambda][evaluationsAtThisLambdaFromAllOtherLambda]
   */
  protected double[][] eAllFlat;
  /**
   * The number of samples for each lambda state.
   */
  protected int[] nSamples;
  /**
   * The temperatures at which the samples were collected.
   */
  protected final double[] temperatures;

  /**
   * The SequentialEstimator constructor largely just copies its parameters into local variables.
   * Most arrays are duplicated (rather a just copying their reference).
   * The temperature array can be of length 1 if all elements are meant to be the same temperature.
   *
   * <p>The first dimension of the energies arrays corresponds to the lambda values. The
   * second dimension (can be of uneven length) corresponds to potential energies of snapshots
   * sampled at that lambda value, calculated either at that lambda value, the lambda value below,
   * or the lambda value above. The arrays eLambdaMindL[0] and eLambdaPlusdL[n-1] is expected to be all
   * NaN.
   *
   * @param lambdaValues   Values of lambda for sampled states.
   * @param eLambdaMinusdL Potential energies of state L at L-dL.
   * @param eLambda        Potential energies of state L at L.
   * @param eLambdaPlusdL  Potential energies of state L at L+dL.
   * @param temperature    Temperature each state (single-element for a constant temperature).
   */
  public SequentialEstimator(double[] lambdaValues, double[][] eLambdaMinusdL, double[][] eLambda,
                             double[][] eLambdaPlusdL, double[] temperature) {
    nStates = lambdaValues.length;
    eAll = null;
    eAllFlat = null;
    nSamples = null;

    assert stream(eLambdaMinusdL[0]).allMatch(Double::isNaN)
        && stream(eLambdaPlusdL[nStates - 1]).allMatch(Double::isNaN);

    assert nStates == eLambda.length
        && nStates == eLambdaMinusdL.length
        && nStates == eLambdaPlusdL.length
        : "One of the energy arrays is of the incorrect length in the first dimension!";

    this.lamValues = copyOf(lambdaValues, nStates);
    temperatures = new double[nStates];
    if (temperature.length == 1) {
      fill(temperatures, temperature[0]);
    } else {
      arraycopy(temperature, 0, temperatures, 0, nStates);
    }

    // Just in case, copy the arrays rather than storing them as provided.
    this.eLambdaMinusdL = new double[nStates][];
    this.eLambda = new double[nStates][];
    this.eLambdaPlusdL = new double[nStates][];
    for (int i = 0; i < nStates; i++) {
      if (i != 0) {
        // There is no perturbation below first state (usually L = 0).
        this.eLambdaMinusdL[i] = copyOf(eLambdaMinusdL[i], eLambdaMinusdL[i].length);
      }
      this.eLambda[i] = copyOf(eLambda[i], eLambda[i].length);
      if (i != nStates - 1) {
        // There is no perturbation above the final state (usually L = 1).
        this.eLambdaPlusdL[i] = copyOf(eLambdaPlusdL[i], eLambdaPlusdL[i].length);
      }
    }
  }


  /**
   * The SequentialEstimator constructor largely just copies its parameters into local variables.
   * Most arrays are duplicated (rather a just copying their reference).
   * The temperature array can be of length 1 if all elements are meant to be the same temperature.
   * <p>
   * This constructor is meant for lower variance estimators such as MBAR and WHAM. These methods require energy
   * evaluations from all lambda windows at all lambda values. The energiesAll array is expected to be
   * of the form energiesAll[lambdaWindow][windowPerspective][lambdaWindowSnapshotPerspectiveEnergy].
   * As an example, at the 3rd lambda window, the energiesAll[2] array should contain the energies
   * of all the snapshots from the 3rd lambda window evaluated at all lambda values. energiesAll[2][3] is a
   * list of all snapshots from lambda 3 evaluated with the potential of lambda 4. energiesAll[2][3][4] is
   * the 5th snapshot from lambda 3 evaluated with the potential of lambda 4.
   * <p>
   * This constructor also breaks energiesAll into a flattened array (across the second dimension) such that
   * the first dimension is the lambda window where the energy was evaluated and the second dimension is the
   * snaps. energiesAll is also broken down into eAt, eLow, and eHigh arrays for convenience and so that BAR
   * calculations can be performed and compared.
   *
   * @param lambdaValues Values of lambda dynamics was run at.
   * @param energiesAll  Potential energy of each sample at all other lambdas. (Missing states are NaN)
   * @param temperature  Temperature each state (single-element for a constant temperature).
   */
  public SequentialEstimator(double[] lambdaValues, double[][][] energiesAll, double[] temperature) {
    nStates = lambdaValues.length;
    assert nStates == energiesAll.length
        : "The energy arrays is of the incorrect length in the first lambda dimension!";
    assert nStates == energiesAll[0].length
        : "The energy arrays is of the incorrect length in the second lambda dimension!";

    this.lamValues = copyOf(lambdaValues, nStates);
    temperatures = new double[nStates];
    if (temperature.length == 1) {
      fill(temperatures, temperature[0]);
    } else {
      arraycopy(temperature, 0, temperatures, 0, nStates);
    }

    // Just in case, deep copy the array rather than storing them as provided.
    eAll = new double[nStates][][];
    int maxSnaps = 0;
    for (int i = 0; i < nStates; i++) {
      eAll[i] = new double[energiesAll[i].length][];
      for (int j = 0; j < energiesAll[i].length; j++) {
        eAll[i][j] = copyOf(energiesAll[i][j], energiesAll[i][j].length);
        maxSnaps = Math.max(maxSnaps, eAll[i][j].length);
      }
    }

    // Remove jagged edges from eAll with NaN values in case it hasn't been done for you
    for (int i = 0; i < nStates; i++) {
      for (int j = 0; j < nStates; j++) {
        if (eAll[i][j].length < maxSnaps) {
          double[] temp = new double[maxSnaps];
          System.arraycopy(eAll[i][j], 0, temp, 0, eAll[i][j].length);
          for (int k = eAll[i][j].length; k < maxSnaps; k++) {
            temp[k] = Double.NaN;
          }
          eAll[i][j] = temp;
        }
      }
    }

    // Flatten the eAll array into a 2D array of [lambda][allEvaluationsAtThisLambda]
    nSamples = new int[nStates];
    int[] nanCount = new int[nStates];
    eAllFlat = new double[nStates][];
    for (int i = 0; i < nStates; i++) {
      ArrayList<Double> temp = new ArrayList<>();
      for (int j = 0; j < nStates; j++) {
        int count = 0;
        int countNaN = 0;
        for (int k = 0; k < eAll[j][i].length; k++) {
          // Don't include NaN values
          if (!Double.isNaN(eAll[j][i][k])) {
            temp.add(eAll[j][i][k]);
            count++;
          } else {
            countNaN++;
          }
        }
        nSamples[j] = count;
        nanCount[j] = countNaN;
      }
      eAllFlat[i] = temp.stream().mapToDouble(Double::doubleValue).toArray();
    }

    for (int i = 0; i < nStates; i++) {
      if (nSamples[i] + nanCount[i] != maxSnaps) {
        throw new IllegalArgumentException("Lambda window " + i
            + " is not set properly. You need to fill in the missing states with NaN.");
      }
    }

    // Assert that lengths of the energiesAll arrays are correct.
    for (int i = 0; i < nStates; i++) {
      assert eAll[i].length == nStates :
          "The energy arrays is of the incorrect length in the second lambda dimension at lambda " + i + "!";
      int nSnapshots = eAll[i][0].length;
      for (int j = 0; j < nStates; j++) {
        assert eAll[i][j].length == nSnapshots :
            "The energy arrays is of the incorrect length in numSnaps dimension at lambda " +
                i + " for evaluation at lambda " + j + "!";
      }
    }

    // Initialize the eLow, eAt, and eHigh arrays to their expected values from eAll. Don't include NaN values.
    // Handle zero sample cases for BAR
    ArrayList<Integer> nonZeroSampleStates = new ArrayList<>();
    for (int i = 0; i < nStates; i++) {
      if (nSamples[i] != 0) {
        nonZeroSampleStates.add(i);
      }
    }
    eLambdaMinusdL = new double[nonZeroSampleStates.size()][eAll[0][0].length];
    fill(eLambdaMinusdL[0], Double.NaN);
    eLambda = new double[nonZeroSampleStates.size()][];
    eLambdaPlusdL = new double[nonZeroSampleStates.size()][eAll[0][0].length];
    fill(eLambdaPlusdL[nonZeroSampleStates.size() - 1], Double.NaN);
    for (int i = 0; i < nonZeroSampleStates.size(); i++) {
      int index = nonZeroSampleStates.get(i); // Contains out of bounds index for e;
      if (i != 0) {
        int indexLow = nonZeroSampleStates.get(i - 1);
        eLambdaMinusdL[i] = copyOf(eAll[index][indexLow], nSamples[index]);
      }
      eLambda[i] = copyOf(eAll[index][index], nSamples[index]);
      if (i != nonZeroSampleStates.size() - 1) {
        int indexHigh = nonZeroSampleStates.get(i + 1);
        eLambdaPlusdL[i] = copyOf(eAll[index][indexHigh], nSamples[index]);
      }
    }
  }

  /**
   * Simpler constructor for when data provided is already flattened (although it adds uncertainty about
   * snap counts, they are all set to the same number).
   *
   * @param lambdaValues Lambda values.
   * @param nSamples     Number of samples for each state.
   * @param eAllFlat     Flattened energy evaluations at all lambda values.
   * @param temperature  Temperature each state (single-element for a constant temperature).
   */
  public SequentialEstimator(double[] lambdaValues, int[] nSamples, double[][] eAllFlat, double[] temperature) {
    nStates = lambdaValues.length;
    assert nStates == eAllFlat.length
        : "The energy arrays is of the incorrect length in the first lambda dimension!";

    this.lamValues = copyOf(lambdaValues, nStates);
    temperatures = new double[nStates];
    if (temperature.length == 1) {
      fill(temperatures, temperature[0]);
    } else {
      arraycopy(temperature, 0, temperatures, 0, nStates);
    }

    this.eAllFlat = eAllFlat;
    this.nSamples = nSamples;
    // No way of knowing the snap counts for each lambda window & therefore no way to break data into these matrices,
    // so just set these all to null.
    eLambdaMinusdL = null;
    eLambda = null;
    eLambdaPlusdL = null;
    eAll = null;
  }
}
