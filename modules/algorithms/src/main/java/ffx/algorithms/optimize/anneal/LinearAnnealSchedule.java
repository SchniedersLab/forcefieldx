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
package ffx.algorithms.optimize.anneal;

/**
 * Linear temperature schedule for simulated annealing
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class LinearAnnealSchedule implements AnnealingSchedule {
  private final int nWindows;
  private final double tHigh;
  private final double tLow;
  private final double dTemp;
  private final String description;

  /**
   * Creates an exponential annealing schedule that decays as tHigh-(n*(tHigh-tLow)).
   *
   * @param nWindows Number of windows.
   * @param tLow Low temperature bound.
   * @param tHigh High temperature bound.
   */
  public LinearAnnealSchedule(int nWindows, double tLow, double tHigh) {
    assert nWindows > 1;
    assert tLow < tHigh;
    assert tLow >= 0;
    assert Double.isFinite(tHigh);

    this.nWindows = nWindows;
    this.tHigh = tHigh;
    this.tLow = tLow;
    dTemp = (tLow - tHigh) / (nWindows - 1);

    description =
        String.format(
            "Linear annealing schedule with %d windows, "
                + "initial temperature %12.7g K, final temperature %12.7g K, dT/window %12.7g K",
            nWindows, tHigh, tLow, dTemp);
  }

  @Override
  public double getHighTemp() {
    return tHigh;
  }

  @Override
  public double getLowTemp() {
    return tLow;
  }

  @Override
  public int getNumWindows() {
    return nWindows;
  }

  @Override
  public double getTemperature(int i) {
    assert i >= 0 && i < nWindows;
    return tHigh + (dTemp * i);
  }

  @Override
  public double[] getTemperatures() {
    double[] temps = new double[nWindows];
    for (int i = 0; i < nWindows; i++) {
      temps[i] = getTemperature(i);
    }
    return temps;
  }

  @Override
  public double maxWindowLength() {
    return 1.0;
  }

  @Override
  public double minWindowLength() {
    return 1.0;
  }

  @Override
  public String toString() {
    return description;
  }

  @Override
  public double totalWindowLength() {
    return nWindows;
  }

  @Override
  public double windowLength(int window) {
    return 1.0;
  }
}
