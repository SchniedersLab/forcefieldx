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
import java.util.concurrent.atomic.DoubleAdder;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * AdderDoubleArray implements the AtomicDoubleArray interface using an array of <code>
 * java.util.concurrent.atomic.DoubleAdder</code>.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AdderDoubleArray implements AtomicDoubleArray {

  private static final Logger logger = Logger.getLogger(AdderDoubleArray.class.getName());

  /**
   * The zero value used for resetting the array.
   */
  private static final double ZERO_VALUE = 0.0;

  /**
   * Atomic operations are handled by an array of DoubleAdder instances.
   */
  private DoubleAdder[] doubleAdders;

  /**
   * The size of the array.
   */
  private int arraySize;

  /**
   * Construct an AdderDoubleArray.
   *
   * @param size Size of the array.
   */
  public AdderDoubleArray(int size) {
    this.arraySize = size;
    this.doubleAdders = createDoubleAdders(size);
  }

  /**
   * Creates the internal array of DoubleAdder instances.
   *
   * @param size The size of the array.
   * @return An array of initialized DoubleAdder objects.
   */
  private DoubleAdder[] createDoubleAdders(int size) {
    DoubleAdder[] adders = new DoubleAdder[size];
    for (int i = 0; i < size; i++) {
      adders[i] = new DoubleAdder();
    }
    return adders;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void add(int threadID, int index, double value) {
    doubleAdders[index].add(value);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void alloc(int size) {
    this.arraySize = size;
    if (doubleAdders.length < size) {
      this.doubleAdders = createDoubleAdders(size);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double get(int index) {
    return doubleAdders[index].sum();
  }

  /**
   * The AtomicDoubleArray handles the reduction automatically, so this method does nothing.
   */
  @Override
  public void reduce(int lb, int ub) {
    // Nothing to do.
  }

  /**
   * The AtomicDoubleArray handles the reduction automatically, so this method does nothing.
   */
  @Override
  public void reduce(ParallelTeam parallelTeam, int lb, int ub) {
    // Nothing to do.
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reset(int threadID, int lb, int ub) {
    resetRange(lb, ub);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reset(ParallelTeam parallelTeam, int lb, int ub) {
    try {
      parallelTeam.execute(
          new ParallelRegion() {
            @Override
            public void run() throws Exception {
              execute(lb, ub,
                  new IntegerForLoop() {
                    @Override
                    public void run(int first, int last) {
                      resetRange(first, last);
                    }
                  });
            }
          });
    } catch (Exception e) {
      logger.log(Level.WARNING, "Exception resetting an AdderDoubleArray", e);
    }
  }

  /**
   * Resets an inclusive range of DoubleAdders in the array.
   */
  private void resetRange(int start, int end) {
    for (int i = start; i <= end; i++) {
      resetAdder(doubleAdders[i]);
    }
  }

  /**
   * Resets a single DoubleAdder to the zero value.
   *
   * @param adder The DoubleAdder to reset.
   */
  private void resetAdder(DoubleAdder adder) {
    adder.reset();
    adder.add(ZERO_VALUE);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void scale(int threadID, int index, double value) {
    double current = doubleAdders[index].sumThenReset();
    doubleAdders[index].add(current * value);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void set(int threadID, int index, double value) {
    resetAdder(doubleAdders[index]);
    doubleAdders[index].add(value);
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
  public void sub(int threadID, int index, double value) {
    doubleAdders[index].add(-value);
  }
}