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

import static java.lang.String.format;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.math.Double3;
import ffx.numerics.math.DoubleMath;
import java.util.Objects;

/**
 * Implementation of maintaining a 3D double array that is operated on by multiple threads.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AtomicDoubleArray3D {

  /**
   * Each dimension is stored in its own AtomicDoubleArray.
   */
  private final AtomicDoubleArray[] atomicDoubleArray;

  /**
   * The implementation.
   */
  private final AtomicDoubleArrayImpl atomicDoubleArrayImpl;

  /**
   * Private ParallelRegion3D for parallel reduction and resets.
   */
  private final ParallelRegion3D parallelRegion3D = new ParallelRegion3D();

  /**
   * Construct an atomic 3D double array of the specified size using the specified implementation.
   *
   * @param atomicDoubleArrayImpl Implementation.
   * @param size Size of each dimension.
   */
  public AtomicDoubleArray3D(AtomicDoubleArrayImpl atomicDoubleArrayImpl, int size) {
    this(atomicDoubleArrayImpl, size, ParallelTeam.getDefaultThreadCount());
  }

  /**
   * Construct an atomic 3D double array of the specified size, using the specified implementation,
   * and the requested number of threads.
   *
   * @param atomicDoubleArrayImpl Implementation.
   * @param size Size of each dimension.
   * @param nThreads Requested number of threads (only used by the MULTI implementation).
   */
  public AtomicDoubleArray3D(AtomicDoubleArrayImpl atomicDoubleArrayImpl, int size, int nThreads) {
    atomicDoubleArray = new AtomicDoubleArray[3];
    atomicDoubleArray[0] = atomicDoubleArrayImpl.createInstance(nThreads, size);
    atomicDoubleArray[1] = atomicDoubleArrayImpl.createInstance(nThreads, size);
    atomicDoubleArray[2] = atomicDoubleArrayImpl.createInstance(nThreads, size);
    this.atomicDoubleArrayImpl = atomicDoubleArrayImpl;
  }

  /**
   * Construct an atomic 3D double array using the specified AtomicDoubleArray instances.
   *
   * @param x AtomicDoubleArray for the X dimension.
   * @param y AtomicDoubleArray for the Y dimension.
   * @param z AtomicDoubleArray for the Z dimension.
   */
  public AtomicDoubleArray3D(AtomicDoubleArray x, AtomicDoubleArray y, AtomicDoubleArray z) {
    atomicDoubleArray = new AtomicDoubleArray[3];
    atomicDoubleArray[0] = x;
    atomicDoubleArray[1] = y;
    atomicDoubleArray[2] = z;
    if (x instanceof MultiDoubleArray) {
      this.atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
    } else if (x instanceof AdderDoubleArray) {
      this.atomicDoubleArrayImpl = AtomicDoubleArrayImpl.ADDER;
    } else {
      this.atomicDoubleArrayImpl = AtomicDoubleArrayImpl.PJ;
    }
  }

  /**
   * Add to the double arrays at the specified index the given values.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param x the x value.
   * @param y the y value.
   * @param z the z value.
   */
  public void add(int threadID, int index, double x, double y, double z) {
    atomicDoubleArray[0].add(threadID, index, x);
    atomicDoubleArray[1].add(threadID, index, y);
    atomicDoubleArray[2].add(threadID, index, z);
  }

  /**
   * Add to the double arrays at the specified index the given Double3.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param d3 the Double3 containing the values to add.
   */
  public void add(int threadID, int index, Double3 d3) {
    atomicDoubleArray[0].add(threadID, index, d3.x());
    atomicDoubleArray[1].add(threadID, index, d3.y());
    atomicDoubleArray[2].add(threadID, index, d3.z());
  }

  /**
   * Ensure the AtomicDoubleArray3D instance is greater than or equal to size.
   *
   * @param size the size of the array.
   */
  public void alloc(int size) {
    atomicDoubleArray[0].alloc(size);
    atomicDoubleArray[1].alloc(size);
    atomicDoubleArray[2].alloc(size);
  }

  /**
   * Get the value of the array at the specified index after calling the <code>reduce</code> method.
   *
   * @param dim Dimension [0, 1, 2]
   * @param index the index.
   * @return the value.
   */
  public double get(int dim, int index) {
    return atomicDoubleArray[dim].get(index);
  }

  /**
   * Get the Double3 at the specified index. This is usually after calling the <code>reduce
   * </code> method.
   *
   * @param index the index.
   * @return a new Double3 with the values at the specified index.
   */
  public Double3 get(int index) {
    return new Double3(
        atomicDoubleArray[0].get(index),
        atomicDoubleArray[1].get(index),
        atomicDoubleArray[2].get(index));
  }

  /**
   * Get the X-value of the array at the specified index after calling the <code>reduce</code>
   * method.
   *
   * @param index the index.
   * @return the X value.
   */
  public double getX(int index) {
    return atomicDoubleArray[0].get(index);
  }

  /**
   * Get the Y-value of the array at the specified index after calling the <code>reduce</code>
   * method.
   *
   * @param index the index.
   * @return the Y value.
   */
  public double getY(int index) {
    return atomicDoubleArray[1].get(index);
  }

  /**
   * Get the Z-value of the array at the specified index after calling the <code>reduce</code>
   * method.
   *
   * @param index the index.
   * @return the Z value.
   */
  public double getZ(int index) {
    return atomicDoubleArray[2].get(index);
  }

  /**
   * Perform reduction between the given lower bound (lb) and upper bound (up) if necessary.
   *
   * @param lb the lower bound.
   * @param ub the upper bound.
   */
  public void reduce(int lb, int ub) {
    // Nothing to do for PJ and Adder.
    if (Objects.requireNonNull(atomicDoubleArrayImpl) == AtomicDoubleArrayImpl.MULTI) {
      atomicDoubleArray[0].reduce(lb, ub);
      atomicDoubleArray[1].reduce(lb, ub);
      atomicDoubleArray[2].reduce(lb, ub);
    }
  }

  /**
   * Perform a reduction on the entire array.
   *
   * @param parallelTeam ParallelTeam to use.
   */
  public void reduce(ParallelTeam parallelTeam) {
    if (Objects.requireNonNull(atomicDoubleArrayImpl) == AtomicDoubleArrayImpl.MULTI) {
      parallelRegion3D.setOperation(Operation.REDUCE);
      try {
        parallelTeam.execute(parallelRegion3D);
      } catch (Exception e) {
        throw new RuntimeException(e);
      }
    }
  }

  /**
   * Reset the double array to Zero.
   *
   * @param threadID the thread ID.
   * @param lb the lower bound.
   * @param ub the upper bound.
   */
  public void reset(int threadID, int lb, int ub) {
    atomicDoubleArray[0].reset(threadID, lb, ub);
    atomicDoubleArray[1].reset(threadID, lb, ub);
    atomicDoubleArray[2].reset(threadID, lb, ub);
  }

  /**
   * Reset the double array to Zero.
   *
   * @param parallelTeam a ParallelTeam.
   */
  public void reset(ParallelTeam parallelTeam) {
    parallelRegion3D.setOperation(Operation.RESET);
    try {
      parallelTeam.execute(parallelRegion3D);
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Scale the double arrays at the specified index to the given values.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param scale The value to scale by.
   */
  public void scale(int threadID, int index, double scale) {
    atomicDoubleArray[0].scale(threadID, index, scale);
    atomicDoubleArray[1].scale(threadID, index, scale);
    atomicDoubleArray[2].scale(threadID, index, scale);
  }

  /**
   * Set the double arrays at the specified index to the given values.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param x the X value.
   * @param y the Y value.
   * @param z the Z value.
   */
  public void set(int threadID, int index, double x, double y, double z) {
    atomicDoubleArray[0].set(threadID, index, x);
    atomicDoubleArray[1].set(threadID, index, y);
    atomicDoubleArray[2].set(threadID, index, z);
  }

  /**
   * Set the double arrays at the specified index to the given Double3.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param d3 the Double3 containing the values.
   */
  public void set(int threadID, int index, Double3 d3) {
    atomicDoubleArray[0].set(threadID, index, d3.x());
    atomicDoubleArray[1].set(threadID, index, d3.y());
    atomicDoubleArray[2].set(threadID, index, d3.z());
  }

  /**
   * Subtracts from the double arrays at the specified index the given values.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param x the X value.
   * @param y the Y value.
   * @param z the Z value.
   */
  public void sub(int threadID, int index, double x, double y, double z) {
    atomicDoubleArray[0].sub(threadID, index, x);
    atomicDoubleArray[1].sub(threadID, index, y);
    atomicDoubleArray[2].sub(threadID, index, z);
  }

  /**
   * Subtracts from the double arrays at the specified index the given Double3.
   *
   * @param threadID the thread ID.
   * @param index the index.
   * @param d3 the Double3 containing the values.
   */
  public void sub(int threadID, int index, Double3 d3) {
    atomicDoubleArray[0].sub(threadID, index, d3.x());
    atomicDoubleArray[1].sub(threadID, index, d3.y());
    atomicDoubleArray[2].sub(threadID, index, d3.z());
  }

  /**
   * Return a string for given index.
   *
   * @param index Index.
   * @return Returns a String for 3D vector at the given index.
   */
  public String toString(int index) {
    String defaultLabel = " " + index + ": ";
    return toString(index, defaultLabel);
  }

  /**
   * Return a string for given index.
   *
   * @param index Index.
   * @param label Label to use for the String.
   * @return Returns a String for 3D vector at the given index.
   */
  public String toString(int index, String label) {
    var d = new double[] {getX(index), getY(index), getZ(index)};
    return DoubleMath.toString(d, label);
  }

  /**
   * Return a String for entire Array, with one 3D vector per line.
   *
   * @return Returns a String for the 3D array.
   */
  public String toString() {
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < atomicDoubleArray.length; i++) {
      sb.append(toString(i)).append("\n");
    }
    return sb.toString();
  }

  /**
   * Return a String for entire Array, with one 3D vector per line.
   *
   * @param label Label may include one "%d" format conversion to include for the index of each
   *     entry.
   * @return Returns a String for the 3D array.
   */
  public String toString(String label) {
    StringBuilder sb = new StringBuilder();
    if (label.contains("%d")) {
      for (int i = 0; i < atomicDoubleArray[0].size(); i++) {
        sb.append(toString(i, format(label, i))).append("\n");
      }
    } else {
      for (int i = 0; i < atomicDoubleArray[0].size(); i++) {
        sb.append(toString(i, label)).append("\n");
      }
    }
    return sb.toString();
  }

  /**
   * Available parallel operations.
   */
  private enum Operation {
    RESET,
    REDUCE
  }

  /**
   * Private ParallelRegion for parallel reset and reduction operations.
   */
  private class ParallelRegion3D extends ParallelRegion {

    private Operation operation;
    private IntegerForLoop3D[] integerForLoop3D;

    /**
     * Instantiates a new parallel region.
     */
    public ParallelRegion3D() {
      operation = Operation.RESET;
    }

    /**
     * Sets the operation.
     *
     * @param operation the new operation
     */
    public void setOperation(Operation operation) {
      this.operation = operation;
    }

    /** {@inheritDoc} */
    @Override
    public void start() {
      int nThreads = getThreadCount();
      if (integerForLoop3D == null) {
        integerForLoop3D = new IntegerForLoop3D[nThreads];
      }
    }

    /** {@inheritDoc} */
    @Override
    public void run() throws Exception {
      int threadID = getThreadIndex();
      if (integerForLoop3D[threadID] == null) {
        integerForLoop3D[threadID] = new IntegerForLoop3D();
      }
      int size = atomicDoubleArray[0].size();
      execute(0, size - 1, integerForLoop3D[threadID]);
    }

    /**
     * Private class for parallel reset and reduction operations.
     */
    private class IntegerForLoop3D extends IntegerForLoop {

      /** {@inheritDoc} */
      @Override
      public void run(int lb, int ub) {
        int threadID = getThreadIndex();
        switch (operation) {
          case RESET -> reset(threadID, lb, ub);
          case REDUCE -> reduce(lb, ub);
        }
      }
    }
  }


}
