/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.numerics.atomic;

import java.util.Arrays;

/**
 * The MultiDoubleArray avoids the need for Atomic variables, but at the cost
 * of storing a full size double array for each thread.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MultiDoubleArray implements AtomicDoubleArray {

    private final int threadCount;
    private final double array[][];

    /**
     * <p>Constructor for MultiDoubleArray.</p>
     *
     * @param nThreads a int.
     * @param size     a int.
     */
    public MultiDoubleArray(int nThreads, int size) {
        array = new double[nThreads][size];
        threadCount = nThreads;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void alloc(int size) {
        for (int i = 0; i < threadCount; i++) {
            if (array[i].length < size) {
                array[i] = new double[size];
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Initialize the storage space for the specified thread, to given value,
     * after assuring adequate storage size.
     */
    @Override
    public void reset(int threadID, int lb, int ub) {
        Arrays.fill(array[threadID], 0.0);
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
    public void sub(int threadID, int index, double value) {
        array[threadID][index] -= value;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Reduce the contributions from each thread into array[0];
     */
    @Override
    public void reduce(int lb, int ub) {
        double gx[] = array[0];
        for (int t = 1; t < threadCount; t++) {
            double gxt[] = array[t];
            for (int i = lb; i <= ub; i++) {
                gx[i] += gxt[i];
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Return a reduced value at the given index.
     */
    @Override
    public double get(int index) {
        return array[0][index];
    }

}
