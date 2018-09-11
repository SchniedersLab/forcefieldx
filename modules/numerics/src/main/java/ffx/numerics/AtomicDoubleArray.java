/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
package ffx.numerics;

/**
 * This interface abstracts away the implementation of maintaining a 1D double
 * array that is operated on by multiple threads.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface AtomicDoubleArray {

    /**
     * AtomicDoubleArray is implementations (ADDER, MULTI, PJ).
     */
    public enum AtomicDoubleArrayImpl {
        ADDER, MULTI, PJ
    };

    /**
     * Ensure the AtomicDoubleArray instance is greater than or equal to size.
     *
     * @param size a int.
     */
    public void alloc(int size);

    /**
     * Reset the double array to Zero.
     *
     * @param threadID a int.
     * @param lb a int.
     * @param ub a int.
     */
    public void reset(int threadID, int lb, int ub);

    /**
     * Add value to the double array at the specified index.
     *
     * @param threadID a int.
     * @param index a int.
     * @param value a double.
     */
    public void add(int threadID, int index, double value);

    /**
     * Subtract value to the double array at the specified index.
     *
     * @param threadID a int.
     * @param index a int.
     * @param value a double.
     */
    public void sub(int threadID, int index, double value);

    /**
     * Perform reduction between the given lower bound (lb) and upper bound (up)
     * if necessary.
     *
     * @param lb a int.
     * @param ub a int.
     */
    public void reduce(int lb, int ub);

    /**
     * Get the value of the array at the specified index (usually subsequent to
     * calling the <code>reduce</code> method.
     *
     * @param index a int.
     * @return a double.
     */
    public double get(int index);

}
