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
   *
   * <p>First dimension is the thread. Second dimension is the value.
   */
  private final double[][] array;

  private int size;

  /**
   * Constructor for MultiDoubleArray.
   *
   * @param nThreads a int.
   * @param size a int.
   */
  public MultiDoubleArray(int nThreads, int size) {
    this.size = size;
    array = new double[nThreads][size];
    threadCount = nThreads;
  }

  /** {@inheritDoc} */
  @Override
  public void add(int threadID, int index, double value) {
    array[threadID][index] += value;
  }

  /** {@inheritDoc} */
  @Override
  public void alloc(int size) {
    this.size = size;
    for (int i = 0; i < threadCount; i++) {
      if (array[i].length < size) {
        array[i] = new double[size];
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Return a reduced value at the given index.
   */
  @Override
  public double get(int index) {
    return array[0][index];
  }

  /**
   * {@inheritDoc}
   *
   * <p>Reduce the contributions from each thread into array[0];
   */
  @Override
  public void reduce(int lb, int ub) {
    double[] gx = array[0];
    for (int t = 1; t < threadCount; t++) {
      double[] gxt = array[t];
      for (int i = lb; i <= ub; i++) {
        gx[i] += gxt[i];
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Reduce the contributions from each thread into array[0];
   */
  @Override
  public void reduce(ParallelTeam parallelTeam, int lb, int ub) {
    try {
      parallelTeam.execute(
          new ParallelRegion() {
            @Override
            public void run() throws Exception {
              execute(
                  lb,
                  ub,
                  new IntegerForLoop() {
                    @Override
                    public void run(int first, int last) {
                      reduce(first, last);
                    }
                  });
            }
          });
    } catch (Exception e) {
      logger.log(Level.WARNING, " Exception reducing an MultiDoubleArray", e);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Initialize the storage space for the specified thread.
   */
  @Override
  public void reset(int threadID, int lb, int ub) {
    Arrays.fill(array[threadID], 0.0);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Initialize the storage space for all threads.
   */
  @Override
  public void reset(ParallelTeam parallelTeam, int lb, int ub) {
    try {
      parallelTeam.execute(
          new ParallelRegion() {
            @Override
            public void run() throws Exception {
              execute(
                  0,
                  threadCount - 1,
                  new IntegerForLoop() {
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
      logger.log(Level.WARNING, " Exception resetting an MultiDoubleArray", e);
    }
  }

  /** {@inheritDoc} */
  @Override
  public void scale(int threadID, int index, double value) {
    array[threadID][index] *= value;
  }

  /** {@inheritDoc} */
  @Override
  public void set(int threadID, int index, double value) {
    array[threadID][index] = value;
  }

  /** {@inheritDoc} */
  @Override
  public int size() {
    return size;
  }

  /** {@inheritDoc} */
  @Override
  public void sub(int threadID, int index, double value) {
    array[threadID][index] -= value;
  }
}
