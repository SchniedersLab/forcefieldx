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
   * Add a value to the double array at the specified index.
   *
   * @param threadID The thread ID.
   * @param index    The index of the array.
   * @param value    The value to add.
   */
  void add(int threadID, int index, double value);

  /**
   * Ensure the AtomicDoubleArray instance has at least the specified size.
   *
   * @param size The required size of the array.
   */
  void alloc(int size);

  /**
   * Get the value of the array at the specified index.
   * Note: The `reduce` method should be called first when using the MULTI implementation.
   *
   * @param index The index of the array.
   * @return The value at the specified index.
   */
  double get(int index);

  /**
   * Perform reduction between the given lower and upper bounds, if necessary.
   *
   * @param lowerBound The lower bound of the range.
   * @param upperBound The upper bound of the range.
   */
  void reduce(int lowerBound, int upperBound);

  /**
   * Perform reduction between the given bounds using a ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam to use.
   * @param lowerBound   The lower bound of the range.
   * @param upperBound   The upper bound of the range.
   */
  void reduce(ParallelTeam parallelTeam, int lowerBound, int upperBound);

  /**
   * Reset the double array values to zero within the specified bounds.
   *
   * @param threadID   The thread ID.
   * @param lowerBound The lower bound of the reset range.
   * @param upperBound The upper bound of the reset range.
   */
  void reset(int threadID, int lowerBound, int upperBound);

  /**
   * Reset the double array values to zero within the specified bounds using a ParallelTeam.
   *
   * @param parallelTeam The ParallelTeam to use.
   * @param lowerBound   The lower bound of the reset range.
   * @param upperBound   The upper bound of the reset range.
   */
  void reset(ParallelTeam parallelTeam, int lowerBound, int upperBound);

  /**
   * Scale the value of the double array at the specified index.
   *
   * @param threadID The thread ID.
   * @param index    The index of the array.
   * @param value    The value to scale by.
   */
  void scale(int threadID, int index, double value);

  /**
   * Set the value of the double array at the specified index.
   *
   * @param threadID The thread ID.
   * @param index    The index of the array.
   * @param value    The value to set.
   */
  void set(int threadID, int index, double value);

  /**
   * Get the current size of the array.
   *
   * @return The size of the array.
   */
  int size();

  /**
   * Subtract a value from the double array at the specified index.
   *
   * @param threadID The thread ID.
   * @param index    The index of the array.
   * @param value    The value to subtract.
   */
  void sub(int threadID, int index, double value);

  /**
   * AtomicDoubleArray implementations (ADDER, MULTI, PJ).
   */
  enum AtomicDoubleArrayImpl {
    /**
     * A java.util.concurrent.atomic.DoubleAdder implementation.
     */
    ADDER {
      @Override
      public AtomicDoubleArray createInstance(int threads, int size) {
        return new AdderDoubleArray(size);
      }
    },

    /**
     * Each thread has its own array, and reduction is performed by the user.
     */
    MULTI {
      @Override
      public AtomicDoubleArray createInstance(int threads, int size) {
        return new MultiDoubleArray(threads, size);
      }
    },

    /**
     * Parallel Java edu.rit.pj.reduction.SharedDoubleArray implementation.
     */
    PJ {
      @Override
      public AtomicDoubleArray createInstance(int threads, int size) {
        return new PJDoubleArray(size);
      }
    };

    /**
     * Factory method to create an AtomicDoubleArray instance.
     *
     * @param threads The number of threads.
     * @param size    The size of the array.
     * @return A new instance of AtomicDoubleArray.
     */
    public abstract AtomicDoubleArray createInstance(int threads, int size);
  }
}