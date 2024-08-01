// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import edu.rit.pj.ParallelTeam;

/**
 * This interface abstracts away the implementation of maintaining a 1D double array that is operated
 * on by multiple threads.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface AtomicDoubleArray {

  /**
   * Factory method to create an AtomicDoubleArray instance.
   *
   * @param atomicDoubleArrayImpl The implementation to use.
   * @param threads               The number of threads.
   * @param size                  The size of the array.
   * @return An AtomicDoubleArray instance.
   */
  static AtomicDoubleArray atomicDoubleArrayFactory(
      AtomicDoubleArrayImpl atomicDoubleArrayImpl, int threads, int size) {
    return switch (atomicDoubleArrayImpl) {
      case ADDER -> new AdderDoubleArray(size);
      case PJ -> new PJDoubleArray(size);
      // MULTI is the default.
      default -> new MultiDoubleArray(threads, size);
    };
  }

  /**
   * Add value to the double array at the specified index.
   *
   * @param threadID the thread ID.
   * @param index    the index.
   * @param value    the value to add.
   */
  void add(int threadID, int index, double value);

  /**
   * Ensure the AtomicDoubleArray instance is greater than or equal to size.
   *
   * @param size The size of the array.
   */
  void alloc(int size);

  /**
   * Get the value of the array at the specified index. The <code>reduce</code> method should be
   * called first when using the MULTI implementation.
   *
   * @param index the index.
   * @return the value of the array at the specified index.
   */
  double get(int index);

  /**
   * Perform reduction between the given lower bound (lb) and upper bound (up) if necessary.
   *
   * @param lb the lower bound.
   * @param ub the upper bound.
   */
  void reduce(int lb, int ub);

  /**
   * Perform reduction between the given lower bound (lb) and upper bound (up) using a ParallelTeam.
   *
   * @param parallelTeam ParallelTeam to use.
   * @param lb           the lower bound.
   * @param ub           the upper bound.
   */
  void reduce(ParallelTeam parallelTeam, int lb, int ub);

  /**
   * Reset the double array to Zero.
   *
   * @param threadID the thread ID.
   * @param lb       the lower bound.
   * @param ub       the upper bound.
   */
  void reset(int threadID, int lb, int ub);

  /**
   * Reset the double array to Zero using a ParallelTeam.
   *
   * @param parallelTeam ParallelTeam to use.
   * @param lb           the lower bound.
   * @param ub           the upper bound.
   */
  void reset(ParallelTeam parallelTeam, int lb, int ub);

  /**
   * Scale the double array at the specified index by the given value.
   *
   * @param threadID the thread ID.
   * @param index    the index.
   * @param value    the value to scale by.
   */
  void scale(int threadID, int index, double value);

  /**
   * Set the double array at the specified index to the given value.
   *
   * @param threadID the thread ID.
   * @param index    the index.
   * @param value    the value to set.
   */
  void set(int threadID, int index, double value);

  /**
   * Get the size of the array.
   *
   * @return Returns the size of the array.
   */
  int size();

  /**
   * Subtract value to the double array at the specified index.
   *
   * @param threadID the thread ID.
   * @param index    the index.
   * @param value    the value to subtract.
   */
  void sub(int threadID, int index, double value);

  /**
   * AtomicDoubleArray implementations (ADDER, MULTI, PJ).
   */
  enum AtomicDoubleArrayImpl {
    /**
     * A java.util.concurrent.atomic.DoubleAdder implementation.
     */
    ADDER,
    /**
     * Each thread has its own array, and reduction is performed by the user.
     */
    MULTI,
    /**
     * Parallel Java edu.rit.pj.reduction.SharedDoubleArray implementation.
     */
    PJ
  }
}
