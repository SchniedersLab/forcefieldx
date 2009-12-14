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
import static java.lang.Math.pow;

/**
 * The OutOfPlaneBendType class defines one Allinger style
 * out-of-plane angle bending energy type.
 */
public final class OutOfPlaneBendType extends BaseType {

    /**
     * Atom classes for this out-of-plane angle bending type.
     */
    public final int atomClasses[];
    /**
     * Force constant (Kcal/mol).
     */
    public final double forceConstant;

    /**
     * OutOfPlaneBendType Constructor.
     *
     * @param atomClasses
     *            int[]
     * @param forceConstant
     *            double
     */
    public OutOfPlaneBendType(int atomClasses[], double forceConstant) {
        super(ForceField.ForceFieldType.OPBEND, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
    }

    /**
     * This method sorts the atom classes for the out-of-plane angle bending type.
     *
     * @param c
     *            atomClasses
     * @return lookup key
     */
    public static String sortKey(int c[]) {
        if (c == null || c.length != 4) {
            return null;
        }
        String key = c[0] + " " + c[1] + " " + c[2] + " " + c[3];
        return key;
    }

    /**
     * Nicely formatted out-of-plane angle bending string.
     *
     * @return String
     */
    @Override
    public String toString() {
        return String.format("opbend  %5d  %5d  %5d  %5d  %4.2f", atomClasses[0],
                atomClasses[1], atomClasses[2],
                atomClasses[3], forceConstant);
    }
    /**
     * Cubic coefficient in out-of-plane angle bending potential.
     */
    public static final double cubic = -0.014;
    /**
     * Quartic coefficient in out-of-plane angle bending potential.
     */
    public static final double quartic = 0.000056;
    /**
     * Quintic coefficient in out-of-plane angle bending potential.
     */
    public static final double quintic = -0.0000007;
    /**
     * Sextic coefficient in out-of-plane angle bending potential.
     */
    public static final double sextic = 0.000000022;
    /**
     * Convert Out-of-Plane bending energy to kcal/mole.
     */
    public static final double units = 1.0 / pow(180.0 / PI, 2);
}
