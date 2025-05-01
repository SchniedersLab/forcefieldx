// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
package ffx.numerics.atomic;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDoubleArray;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * PJDoubleArray implements the AtomicDoubleArray interface using the Parallel Java class
 * SharedDoubleArray.
 *
 * <p>SharedDoubleArray is multiple thread safe and uses lock-free atomic compare-and-set.
 *
 * <p>Note: Class SharedDoubleArray is implemented using class
 * java.util.concurrent.atomic.AtomicLongArray. Each double array element is stored as a long whose
 * bit pattern is the same as the double value.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PJDoubleArray implements AtomicDoubleArray {

  private static final Logger logger = Logger.getLogger(PJDoubleArray.class.getName());

  private SharedDoubleArray array;
  private int arraySize;

  /**
   * Constructor: Initialize the size and array.
   *
   * @param arraySize the desired size of the array.
   */
  public PJDoubleArray(int arraySize) {
    this.arraySize = arraySize;
    this.array = new SharedDoubleArray(arraySize);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void add(final int threadID, final int index, final double value) {
    array.getAndAdd(index, value);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void alloc(final int newSize) {
    this.arraySize = newSize;
    if (array.length() < newSize) {
      array = new SharedDoubleArray(newSize);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double get(final int index) {
    return array.get(index);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reduce(final int lowerBound, final int upperBound) {
    // Intentionally left empty: handled atomically by SharedDoubleArray.
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reduce(final ParallelTeam parallelTeam, final int lowerBound, final int upperBound) {
    // Intentionally left empty: handled atomically by SharedDoubleArray.
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reset(final int threadID, final int lowerBound, final int upperBound) {
    resetSerially(lowerBound, upperBound);
  }

  /**
   * Reset elements serially to 0.0 in the given range.
   *
   * @param lowerBound inclusive start index.
   * @param upperBound inclusive end index.
   */
  private void resetSerially(final int lowerBound, final int upperBound) {
    for (int i = lowerBound; i <= upperBound; i++) {
      array.set(i, 0.0);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reset(final ParallelTeam parallelTeam, final int lowerBound, final int upperBound) {
    try {
      parallelTeam.execute(
          new ParallelRegion() {
            @Override
            public void run() throws Exception {
              execute(
                  lowerBound,
                  upperBound,
                  new IntegerForLoop() {
                    @Override
                    public void run(final int first, final int last) {
                      resetSerially(first, last);
                    }
                  });
            }
          });
    } catch (Exception e) {
      logger.log(Level.WARNING, "Exception resetting a PJDoubleArray", e);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void scale(final int threadID, final int index, final double value) {
    final double current = array.get(index);
    array.getAndSet(index, current * value);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void set(final int threadID, final int index, final double value) {
    array.getAndSet(index, value);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int size() {
    return arraySize;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void sub(final int threadID, final int index, final double value) {
    array.getAndAdd(index, -value);
  }
}