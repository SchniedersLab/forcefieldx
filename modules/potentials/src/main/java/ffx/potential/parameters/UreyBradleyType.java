/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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

import java.util.Arrays;
import java.util.Comparator;

/**
 * The UreyBradleyType class defines one harmonic UreyBradley cross term.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public final class UreyBradleyType extends BaseType implements Comparator<String> {

    /**
     * Atom classes that form this Urey-Bradley cross term.
     */
    public final int atomClasses[];
    /**
     * Force constant (Kcal/mole/angstroms^2).
     */
    public final double forceConstant;
    /**
     * Equilibrium 1-3 separation (Angstroms).
     */
    public final double distance;

    /**
     * UreyBradleyType constructor.
     *
     * @param atomClasses
     *            int[]
     * @param forceConstant
     *            double
     * @param distance
     *            double
     */
    public UreyBradleyType(int atomClasses[], double forceConstant,
            double distance) {
        super(ForceField.ForceFieldType.UREYBRAD, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
        this.distance = distance;
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
     * Nicely formatted Urey-Bradley string.
     *
     * @return String
     */
    @Override
    public String toString() {
        return String.format("ureybrad  %5d  %5d  %5d  %6.2f  %7.4f",
                atomClasses[0], atomClasses[1], atomClasses[2], forceConstant,
                distance);
    }
    /**
     * Convert bond stretch energy to kcal/mole.
     */
    public static final double units = 1.0;
    /**
     * Cubic coefficient in bond stretch potential.
     */
    public static final double cubic = 0.0;
    /**
     * Quartic coefficient in bond stretch potential.
     */
    public static final double quartic = 0.0;

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
        if (other == null || !(other instanceof UreyBradleyType)) {
            return false;
        }
        UreyBradleyType ureyBradleyType = (UreyBradleyType) other;
        int c[] = ureyBradleyType.atomClasses;
        if (c[0] == atomClasses[0] && c[1] == atomClasses[1] && c[2] == atomClasses[2]) {
            return true;
        }
        return false;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 79 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }
}
