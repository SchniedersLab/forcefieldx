/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
 */
package ffx.potential.parameters;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

/**
 * The ImproperType class defines an improper torsion.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class ImproperTorsionType extends BaseType implements Comparator<String> {

    /**
     * Atom classes that for this Improper Torsion angle.
     */
    public final int atomClasses[];
    /**
     * Phases in degrees.
     */
    public final double phase;
    /**
     * Force constant in kcal/mol.
     */
    public final double k;
    /**
     * Periodicity.
     */
    public final int periodicity;

    /**
     * TorsionType Constructor.
     *
     * @param atomClasses int[]
     * @param k double
     * @param phase double
     * @param periodicity int
     */
    public ImproperTorsionType(int atomClasses[], double k, double phase,
            int periodicity) {
        super(ForceField.ForceFieldType.IMPTORS, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.k = k;
        this.phase = phase;
        this.periodicity = periodicity;
    }

    /**
     * <p>
     * incrementClasses</p>
     *
     * @param increment a int.
     */
    public void incrementClasses(int increment) {
        for (int i = 0; i < atomClasses.length; i++) {
            atomClasses[i] += increment;
        }
        setKey(sortKey(atomClasses));
    }

    /**
     * Remap new atom classes to known internal ones.
     *
     * @param typeMap a lookup between new atom types and known atom types.
     */
    public void patchClasses(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        for (AtomType newType : typeMap.keySet()) {
            for (int i = 0; i < atomClasses.length; i++) {
                if (atomClasses[i] == newType.atomClass) {
                    count++;
                }
            }
        }
        if (count > 0 && count < atomClasses.length) {
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < atomClasses.length; i++) {
                    if (atomClasses[i] == newType.atomClass) {
                        AtomType knownType = typeMap.get(newType);
                        atomClasses[i] = knownType.atomClass;
                    }
                }

            }
            setKey(sortKey(atomClasses));
        }
    }

    /**
     * {@inheritDoc}
     *
     * Nicely formatted Torsion angle.
     *
     * @since 1.0
     */
    @Override
    public String toString() {
        StringBuilder imptorsBuffer = new StringBuilder("torsion");
        for (int i : atomClasses) {
            imptorsBuffer.append(String.format(" %5d", i));
        }
        imptorsBuffer.append(String.format(" %7.3f %7.3f %1d",
                k, phase, periodicity));

        return imptorsBuffer.toString();
    }

    /**
     * This method sorts the atom classes for the torsion.
     *
     * @param c atomClasses
     * @return lookup key
     * @since 1.0
     */
    public static String sortKey(int c[]) {
        if (c == null || c.length != 4) {
            return null;
        }
        if (c[1] < c[2]) {
            // Do nothing.
        } else if (c[2] < c[1]) {
            // Reverse the order.
            int temp = c[0];
            c[0] = c[3];
            c[3] = temp;
            temp = c[1];
            c[1] = c[2];
            c[2] = temp;
        } else if (c[0] <= c[3]) {
            // Do nothing.
        } else {
            // Reverse the order.
            int temp = c[0];
            c[0] = c[3];
            c[3] = temp;
            temp = c[1];
            c[1] = c[2];
            c[2] = temp;
        }
        String key = c[0] + " " + c[1] + " " + c[2] + " " + c[3];
        return key;
    }
    /**
     * Convert Torsional Angle energy to kcal/mole.
     *
     * @since 1.0
     */
    public static final double units = 0.5;

    /**
     * {@inheritDoc}
     *
     * Implements the Comparator<String> interface.
     *
     * @since 1.0
     */
    @Override
    public int compare(String s1, String s2) {
        String keys1[] = s1.split(" ");
        String keys2[] = s2.split(" ");
        int c1[] = new int[4];
        int c2[] = new int[4];

        for (int i = 0; i < 4; i++) {
            c1[i] = Integer.parseInt(keys1[i]);
            c2[i] = Integer.parseInt(keys2[i]);
        }

        if (c1[1] < c2[1]) {
            return -1;
        } else if (c1[1] > c2[1]) {
            return 1;
        } else if (c1[2] < c2[2]) {
            return -1;
        } else if (c1[2] > c2[2]) {
            return 1;
        } else if (c1[0] < c2[0]) {
            return -1;
        } else if (c1[0] > c2[0]) {
            return 1;
        } else if (c1[3] < c2[3]) {
            return -1;
        } else if (c1[3] > c2[3]) {
            return 1;
        }

        return 0;
    }

    /**
     * {@inheritDoc}
     *
     * Override the default <code>equals</code> method.
     *
     * @since 1.0
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (other == null || !(other instanceof ImproperTorsionType)) {
            return false;
        }
        ImproperTorsionType torsionType = (ImproperTorsionType) other;
        for (int i = 0; i < 4; i++) {
            if (torsionType.atomClasses[i] != atomClasses[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of the <code>hashCode</code> method.
     *
     * @since 1.0
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 89 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }
}
