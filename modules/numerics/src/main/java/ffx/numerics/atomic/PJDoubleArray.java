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
import edu.rit.pj.reduction.SharedDoubleArray;

/**
 * PJDoubleArray implements the AtomicDoubleArray interface using the Parallel
 * Java class SharedDoubleArray.
 * <p>
 * SharedDoubleArray is multiple thread safe and uses lock-free atomic
 * compare-and-set.
 * <p>
 * Note: Class SharedDoubleArray is implemented using class
 * java.util.concurrent.atomic.AtomicLongArray. Each double array element is
 * stored as a long whose bit pattern is the same as the double value.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PJDoubleArray implements AtomicDoubleArray {

    private static final Logger logger = Logger.getLogger(PJDoubleArray.class.getName());

    private SharedDoubleArray array;

    /**
     * <p>Constructor for PJDoubleArray.</p>
     *
     * @param size a int.
     */
    public PJDoubleArray(int size) {
        array = new SharedDoubleArray(size);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void alloc(int size) {
        if (array.length() < size) {
            array = new SharedDoubleArray(size);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void reset(int threadID, int lb, int ub) {
        for (int i = lb; i <= ub; i++) {
            array.set(i, 0.0);
        }
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
                    execute(lb, ub, new IntegerForLoop() {
                        @Override
                        public void run(int first, int last) {
                            reset(getThreadIndex(), first, last);
                        }
                    });
                }
            });
        } catch (Exception e) {
            logger.log(Level.WARNING, " Exception resetting an PJDoubleArray", e);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void add(int threadID, int index, double value) {
        array.getAndAdd(index, value);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void sub(int threadID, int index, double value) {
        array.getAndAdd(index, -value);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Reduction is handled atomically by the PJ SharedDoubleArray.
     */
    @Override
    public void reduce(int lb, int ub) {
        // Nothing to do.
    }

    /**
     * {@inheritDoc}
     * Reduction is handled atomically by the PJ SharedDoubleArray.
     */
    @Override
    public void reduce(ParallelTeam parallelTeam, int lb, int ub) {
        // Nothing to do.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double get(int index) {
        return array.get(index);
    }

}
