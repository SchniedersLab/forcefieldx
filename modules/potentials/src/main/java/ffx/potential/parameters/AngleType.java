/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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

import java.util.Arrays;
import java.util.Comparator;

/**
 * The AngleType class defines one harmonic angle bend energy term.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public final class AngleType extends BaseType implements Comparator<String> {

    /**
     * Atom classes that for this Angle type.
     */
    public final int atomClasses[];
    /**
     * Force constant (Kcal/mole/radian^2).
     */
    public final double forceConstant;
    /**
     * Equilibrium angle (degrees). There can be up to three equilibrium angles,
     * depending on the number of attached hydrogens (0, 1, or 2).
     */
    public final double angle[];

    /**
     * @param atomClasses
     * @param forceConstant
     * @param angle
     */
    public AngleType(int atomClasses[], double forceConstant, double angle[]) {
        super(ForceField.ForceFieldType.ANGLE, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
        this.angle = angle;
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
     * Nicely formatted Angle bending string.
     *
     * @return String
     */
    @Override
    public String toString() {
        StringBuilder angleString = new StringBuilder(String.format(
                "angle  %5d  %5d  %5d  %6.2f", atomClasses[0], atomClasses[1],
                atomClasses[2], forceConstant));
        for (double eq : angle) {
            angleString.append(String.format("  %6.2f", eq));
        }
        return angleString.toString();
    }
    /**
     * Cubic coefficient in angle bending potential.
     */
    public static final double cubic = -0.014;
    /**
     * Quartic coefficient in angle bending potential.
     */
    public static final double quartic = 0.000056;
    /**
     * Quintic coefficient in angle bending potential.
     */
    public static final double quintic = -0.0000007;
    /**
     * Sextic coefficient in angle bending potential.
     */
    public static final double sextic = 0.000000022;
    /**
     * Convert angle bending energy to kcal/mole.
     */
    public static double units = 1.0 / pow(180.0 / PI, 2);

    @Override
    public int compare(String key1, String key2) {
        String keys1[] = key1.split(" ");
        String keys2[] = key2.split(" ");
        int c1[] = new int[3];
        int c2[] = new int[3];
        for (int i = 0; i < 3; i++) {
            c1[i] = Integer.parseInt(keys1[i]);
            c2[i] = Integer.parseInt(keys2[i]);
        }

        if (c1[1] < c2[1]) {
            return -1;
        } else if (c1[1] > c2[1]) {
            return 1;
        } else if (c1[0] < c2[0]) {
            return -1;
        } else if (c1[0] > c2[0]) {
            return 1;
        } else if (c1[2] < c2[2]) {
            return -1;
        } else if (c1[2] > c2[2]) {
            return 1;
        }

        return 0;
    }

    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (other == null || !(other instanceof AngleType)) {
            return false;
        }
        AngleType angleType = (AngleType) other;
        int c[] = angleType.atomClasses;
        if (c[0] == atomClasses[0] && c[1] == atomClasses[1] && c[2] == atomClasses[2]) {
            return true;
        }
        return false;
    }

    @Override
    public int hashCode() {
        int hash = 5;
        hash = 37 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }
}
