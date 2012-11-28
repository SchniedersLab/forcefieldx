/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 * The PiTorsionType class defines a Pi-Torsion energy term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class PiTorsionType extends BaseType implements Comparator<String> {

    /**
     * Atom classes that form this Pi-Torsion.
     */
    public final int atomClasses[];
    /**
     * Force constant.
     */
    public final double forceConstant;

    /**
     * PiTorsionType Constructor.
     *
     * @param atomClasses int[]
     * @param forceConstant double
     */
    public PiTorsionType(int atomClasses[], double forceConstant) {
        super(ForceField.ForceFieldType.PITORS, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
    }

    /**
     * <p>incrementClasses</p>
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
     * This method sorts the atom classes as: min, max
     *
     * @param c atomClasses
     * @return lookup key
     */
    public static String sortKey(int c[]) {
        if (c == null || c.length != 2) {
            return null;
        }
        String key = null;
        int temp;
        if (c[1] <= c[0]) {
            temp = c[1];
            c[1] = c[0];
            c[0] = temp;
        }
        key = c[0] + " " + c[1];
        return key;
    }

    /**
     * {@inheritDoc}
     *
     * Nicely formatted Pi-Torsion type.
     */
    @Override
    public String toString() {
        return String.format("pitors  %5d  %5d  %4.2f", atomClasses[0],
                atomClasses[1], forceConstant);
    }
    /**
     * Convert Pi-Torsion energy to kcal/mole.
     */
    public static double units = 1.0;

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String s1, String s2) {
        String keys1[] = s1.split(" ");
        String keys2[] = s2.split(" ");

        for (int i = 0; i < 2; i++) {
            int c1 = Integer.parseInt(keys1[i]);
            int c2 = Integer.parseInt(keys2[i]);
            if (c1 < c2) {
                return -1;
            } else if (c1 > c2) {
                return 1;
            }
        }
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (other == null || !(other instanceof PiTorsionType)) {
            return false;
        }
        PiTorsionType piTorsionType = (PiTorsionType) other;
        for (int i = 0; i < 2; i++) {
            if (piTorsionType.atomClasses[i] != atomClasses[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 97 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }
}
