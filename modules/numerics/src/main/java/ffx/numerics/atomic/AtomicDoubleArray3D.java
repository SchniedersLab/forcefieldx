//******************************************************************************
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
//******************************************************************************
package ffx.numerics.atomic;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import static ffx.numerics.atomic.AtomicDoubleArray.atomicDoubleArrayFactory;

/**
 * Implementation of maintaining a 3D double array that is operated on by multiple threads.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AtomicDoubleArray3D {

    private static final Logger logger = Logger.getLogger(AtomicDoubleArray3D.class.getName());

    /**
     * Each dimension is stored in its own AtomicDoubleArray.
     */
    private final AtomicDoubleArray[] atomicDoubleArray;

    /**
     * Construct an atomic 3D double array of the specified size using the specified implementation.
     *
     * @param atomicDoubleArrayImpl Implementation.
     * @param size                  Size of each dimension.
     */
    public AtomicDoubleArray3D(AtomicDoubleArrayImpl atomicDoubleArrayImpl, int size) {
        this(atomicDoubleArrayImpl, size, ParallelTeam.getDefaultThreadCount());
    }

    /**
     * Construct an atomic 3D double array of the specified size, using the specified implementation,
     * and the requested number of threads.
     *
     * @param atomicDoubleArrayImpl Implementation.
     * @param size                  Size of each dimension.
     * @param nThreads              Requested number of threads (only used by the MULTI implementation).
     */
    public AtomicDoubleArray3D(AtomicDoubleArrayImpl atomicDoubleArrayImpl, int size, int nThreads) {
        atomicDoubleArray = new AtomicDoubleArray[3];
        atomicDoubleArray[0] = atomicDoubleArrayFactory(atomicDoubleArrayImpl, nThreads, size);
        atomicDoubleArray[1] = atomicDoubleArrayFactory(atomicDoubleArrayImpl, nThreads, size);
        atomicDoubleArray[2] = atomicDoubleArrayFactory(atomicDoubleArrayImpl, nThreads, size);
    }

    public AtomicDoubleArray3D(AtomicDoubleArray x, AtomicDoubleArray y, AtomicDoubleArray z) {
        atomicDoubleArray = new AtomicDoubleArray[3];
        atomicDoubleArray[0] = x;
        atomicDoubleArray[1] = y;
        atomicDoubleArray[2] = z;
    }

    /**
     * Ensure the AtomicDoubleArray3D instance is greater than or equal to size.
     *
     * @param size a int.
     */
    public void alloc(int size) {
        atomicDoubleArray[0].alloc(size);
        atomicDoubleArray[1].alloc(size);
        atomicDoubleArray[2].alloc(size);
    }

    /**
     * Reset the double array to Zero.
     *
     * @param threadID a int.
     * @param lb       a int.
     * @param ub       a int.
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
     * @param lb           a int.
     * @param ub           a int.
     */
    public void reset(ParallelTeam parallelTeam, int lb, int ub) {
        try {
            parallelTeam.execute(new ParallelRegion() {
                @Override
                public void run() throws Exception {
                    int nThreads = parallelTeam.getThreadCount();
                    execute(0, nThreads - 1, new IntegerForLoop() {
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
            logger.log(Level.WARNING, " Exception resetting an AtomicDoubleArray3D", e);
        }
    }

    /**
     * Add value to the double array at the specified index.
     *
     * @param threadID a int.
     * @param index    a int.
     * @param x        a double.
     * @param y        a double.
     * @param z        a double.
     */
    public void add(int threadID, int index, double x, double y, double z) {
        atomicDoubleArray[0].add(threadID, index, x);
        atomicDoubleArray[1].add(threadID, index, y);
        atomicDoubleArray[2].add(threadID, index, z);
    }

    /**
     * Add value to the double array at the specified index.
     *
     * @param threadID a int.
     * @param index    a int.
     * @param x        a double.
     * @param y        a double.
     * @param z        a double.
     */
    public void sub(int threadID, int index, double x, double y, double z) {
        atomicDoubleArray[0].sub(threadID, index, x);
        atomicDoubleArray[1].sub(threadID, index, y);
        atomicDoubleArray[2].sub(threadID, index, z);
    }

    /**
     * Perform reduction between the given lower bound (lb) and upper bound (up)
     * if necessary.
     *
     * @param lb a int.
     * @param ub a int.
     */
    public void reduce(int lb, int ub) {
        atomicDoubleArray[0].reduce(lb, ub);
        atomicDoubleArray[1].reduce(lb, ub);
        atomicDoubleArray[2].reduce(lb, ub);
    }

    /**
     * Perform reduction between the given lower bound (lb) and upper bound (up)
     * if necessary.
     *
     * @param lb a int.
     * @param ub a int.
     */
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
            logger.log(Level.WARNING, " Exception reducing an AtomicDoubleArray3D", e);
        }
    }

    /**
     * Get the value of the array at the specified index (usually subsequent to
     * calling the <code>reduce</code> method.
     *
     * @param index a int.
     * @return a double.
     */
    public double getX(int index) {
        return atomicDoubleArray[0].get(index);
    }

    /**
     * Get the value of the array at the specified index (usually subsequent to
     * calling the <code>reduce</code> method.
     *
     * @param index a int.
     * @return a double.
     */
    public double getY(int index) {
        return atomicDoubleArray[1].get(index);
    }

    /**
     * Get the value of the array at the specified index (usually subsequent to
     * calling the <code>reduce</code> method.
     *
     * @param index a int.
     * @return a double.
     */
    public double getZ(int index) {
        return atomicDoubleArray[2].get(index);
    }

    /**
     * Get the value of the array at the specified index (usually subsequent to
     * calling the <code>reduce</code> method.
     *
     * @param dim   Dimension [0, 1, 2]
     * @param index a int.
     * @return a double.
     */
    public double get(int dim, int index) {
        return atomicDoubleArray[dim].get(index);
    }
}
