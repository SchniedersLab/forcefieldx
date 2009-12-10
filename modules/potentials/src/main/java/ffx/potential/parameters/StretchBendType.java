/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.parameters;

import static java.lang.Math.PI;

/**
 * The StretchBendType class defines one out-of-plane angle bending energy type.
 */
public final class StretchBendType extends BaseType {

    /**
     * Atom class for this out-of-plane angle bending type.
     */
    public final int atomClasses[];
    /**
     * Force constants (Kcal/mole/Angstrom-Degrees).
     */
    public final double forceConstants[];

    /**
     * StretchBendType Constructor.
     *
     * @param atomClasses
     *            int[]
     * @param forceConstants
     *            double[]
     */
    public StretchBendType(int atomClasses[], double forceConstants[]) {
        super(ForceField.ForceFieldType.STRBND, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstants = forceConstants;
    }

    /**
     * This method sorts the atom classes as: min, c[1], max
     *
     * @param c
     *            atomClasses
     * @return lookup key
     */
    public static String sortKey(int c[]) {
        if (c == null || c.length != 3) {
            return null;
        }
        if (c[0] > c[2]) {
            int temp = c[0];
            c[0] = c[2];
            c[2] = temp;
        }
        String key = c[0] + " " + c[1] + " " + c[2];
        return key;
    }

    /**
     * Nicely formatted Stretch-Bend string.
     *
     * @return String
     */
    @Override
    public String toString() {
        return String.format("strbnd  %5d  %5d  %5d  %6.2f  %6.2f",
                atomClasses[0], atomClasses[1], atomClasses[2],
                forceConstants[0], forceConstants[1]);
    }
    
    public static final double units = PI / 180.0;
}
