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

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The MultiDoubleArray avoids the need for Atomic variables, but at the cost of storing a full size
 * double array for each thread.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MultiDoubleArray implements AtomicDoubleArray {

  private static final Logger logger = Logger.getLogger(MultiDoubleArray.class.getName());

  private final int threadCount;

  /**
   * Storage of the array.
   * <p>
   * First dimension is the thread. Second dimension is the value.
   */
  private final double[][] array;

  private int size;

  /**
   * Constructor for MultiDoubleArray.
   *
   * @param threadCount the number of threads.
   * @param arraySize   the size of the array.
   */
  public MultiDoubleArray(int threadCount, int arraySize) {
    this.size = arraySize;
    array = new double[threadCount][arraySize];
    this.threadCount = threadCount;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void add(int threadID, int index, double value) {
    array[threadID][index] += value;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void alloc(int arraySize) {
    this.size = arraySize;
    for (int i = 0; i < threadCount; i++) {
      if (array[i].length < arraySize) {
        array[i] = new double[arraySize];
      }
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double get(int index) {
    return array[0][index];
  }

  /**
   * Reduce contributions from each thread into the main thread's array.
   *
   * @param lb lower bound index.
   * @param ub upper bound index.
   */
  @Override
  public void reduce(int lb, int ub) {
    double[] mainArray = array[0];
    for (int t = 1; t < threadCount; t++) {
      double[] threadArray = array[t];
      reduceFromThreads(mainArray, threadArray, lb, ub);
    }
  }

  /**
   * Perform reductions from one thread's array to the main array.
   *
   * @param mainArray   the main array where reductions are stored.
   * @param threadArray the thread-specific array contributing to the reduction.
   * @param lb          the lower bound for the range to reduce.
   * @param ub          the upper bound for the range to reduce.
   */
  private void reduceFromThreads(double[] mainArray, double[] threadArray, int lb, int ub) {
    for (int i = lb; i <= ub; i++) {
      mainArray[i] += threadArray[i];
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reduce(ParallelTeam parallelTeam, int lb, int ub) {
    try {
      parallelTeam.execute(new ParallelRegion() {
        @Override
        public void run() throws Exception {
          execute(lb, ub, new IntegerForLoop() {
            @Override
            public void run(int first, int last) {
              reduce(first, last);
            }
          });
        }
      });
    } catch (Exception e) {
      logger.log(Level.WARNING, "Exception reducing a MultiDoubleArray", e);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reset(int threadID, int lb, int ub) {
    Arrays.fill(array[threadID], 0.0);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void reset(ParallelTeam parallelTeam, int lb, int ub) {
    try {
      parallelTeam.execute(new ParallelRegion() {
        @Override
        public void run() throws Exception {
          execute(0, threadCount - 1, new IntegerForLoop() {
            @Override
            public void run(int first, int last) {
              for (int i = first; i <= last; i++) {
                reset(i, lb, ub);
              }
            }
          });
        }
      });
    } catch (Exception e) {
      logger.log(Level.WARNING, "Exception resetting a MultiDoubleArray", e);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void scale(int threadID, int index, double value) {
    array[threadID][index] *= value;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void set(int threadID, int index, double value) {
    array[threadID][index] = value;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int size() {
    return size;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void sub(int threadID, int index, double value) {
    array[threadID][index] -= value;
  }
}