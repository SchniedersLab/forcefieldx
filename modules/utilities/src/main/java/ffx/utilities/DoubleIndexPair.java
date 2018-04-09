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
package ffx.utilities;

/**
 * <p>
 * DoubleIndexPair class.</p>
 *
 * @author Jacob M. Litman
 *
 * @since 1.0
 */
public class DoubleIndexPair implements Comparable<DoubleIndexPair> {

    private final int index;
    private final double doubleValue;

    /**
     * Allows sorting of floating-point values while retaining knowledge of
     * where that value was in some original list or array.
     *
     * @param index the original index of this double.
     * @param doubleValue the double value.
     */
    public DoubleIndexPair(int index, double doubleValue) {
        this.index = index;
        this.doubleValue = doubleValue;
    }

    public int getIndex() {
        return index;
    }

    public double getDoubleValue() {
        return doubleValue;
    }

    @Override
    public int compareTo(DoubleIndexPair o) {
        if (o == null) {
            return 0;
        }
        return Double.compare(doubleValue, o.doubleValue);
    }
}
