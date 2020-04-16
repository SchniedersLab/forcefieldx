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
 * Composite annealing schedule with flat ends (i.e. spends extra time at the low and high
 * temperatures).
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class FlatEndAnnealSchedule implements AnnealingSchedule {
  private final AnnealingSchedule middle;
  private final double tLow;
  private final double tHigh;
  private final double lenBefore;
  private final double lenAfter;
  private final int totWindows;
  private final boolean useBefore;
  private final boolean useAfter;
  private final String description;

  /**
   * Creates a flat-ended annealing schedule based on a provided schedule for the middle, which is
   * flat for some number of steps at the ends.
   *
   * @param middle Annealing schedule to use in the middle.
   * @param tLow Low temperature to use at the last window.
   * @param tHigh High temperature to use at the first window.
   * @param lengthBefore Relative length of the first ("flat") tempering window at the high
   *     temperature. 0 disables this end.
   * @param lengthAfter Relative length of the last ("flat") tempering window at the low
   *     temperature. 0 disables this end.
   */
  public FlatEndAnnealSchedule(
      AnnealingSchedule middle,
      double tLow,
      double tHigh,
      double lengthBefore,
      double lengthAfter) {
    assert tLow < tHigh;
    assert tLow >= 0;
    assert Double.isFinite(tHigh);
    assert lengthBefore >= 0;
    assert lengthAfter >= 0;

    this.middle = middle;
    this.tLow = tLow;
    this.tHigh = tHigh;
    this.lenBefore = lengthBefore;
    this.lenAfter = lengthAfter;

    int nWin = middle.getNumWindows();
    useBefore = lengthBefore > 0;
    if (useBefore) {
      ++nWin;
    }
    useAfter = lengthAfter > 0;
    if (useAfter) {
      ++nWin;
    }
    totWindows = nWin;

    StringBuilder sb =
        new StringBuilder(
            String.format(
                "Flat-ended annealing schedule: starts at %10.4g K, ends at %10.4g K.",
                tHigh, tLow));
    if (useBefore) {
      sb.append(
          String.format(
              "\nMain annealing is preceded by a %9.3f-long window at %10.4g K", lenBefore, tHigh));
    }
    if (useAfter) {
      sb.append(
          String.format(
              "\nMain annealing is followed by a %9.3f-long window at %10.4g K", lenAfter, tLow));
    }
    sb.append("\nMain annealing schedule: ").append(middle.toString());
    this.description = sb.toString();
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
    return totWindows;
  }

  @Override
  public double getTemperature(int i) {
    assert i >= 0 && i < totWindows;
    if (i == 0 && useBefore) {
      return tHigh;
    } else if (i == (totWindows - 1) && useAfter) {
      return tLow;
    } else {
      i = useBefore ? i - 1 : i;
      return middle.getTemperature(i);
    }
  }

  @Override
  public double[] getTemperatures() {
    double[] temps = new double[totWindows];
    for (int i = 0; i < totWindows; i++) {
      temps[i] = getTemperature(i);
    }
    return temps;
  }

  @Override
  public double maxWindowLength() {
    return Math.max(Math.max(lenBefore, lenAfter), middle.maxWindowLength());
  }

  @Override
  public double minWindowLength() {
    double len = middle.minWindowLength();
    if (useAfter) {
      len = Math.min(len, lenAfter);
    }
    if (useBefore) {
      len = Math.min(len, lenBefore);
    }
    return len;
  }

  @Override
  public String toString() {
    return description;
  }

  @Override
  public double totalWindowLength() {
    return lenBefore + lenAfter + middle.totalWindowLength();
  }

  @Override
  public double windowLength(int window) {
    assert window >= 0 && window < totWindows;

    if (useBefore && window == 0) {
      return lenBefore;
    } else if (useAfter && window == (totWindows - 1)) {
      return lenAfter;
    } else {
      window = useBefore ? window - 1 : window;
      return middle.windowLength(window);
    }
  }
}
