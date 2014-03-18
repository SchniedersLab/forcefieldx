/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 * The BondType class defines one harmonic bond stretch energy term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class BondType extends BaseType implements Comparator<String> {

    public enum BondFunction {

        HARMONIC, QUARTIC
    }

    /**
     * Atom classes that form this bond stretch.
     */
    public final int atomClasses[];
    /**
     * Force constant (Kcal/mol).
     */
    public final double forceConstant;
    /**
     * Equilibrium separation (Angstroms).
     */
    public final double distance;

    public final BondFunction bondFunction;

    /**
     * BondType constructor.
     *
     * @param atomClasses int[]
     * @param forceConstant double
     * @param distance double
     * @param bondFunction
     */
    public BondType(int atomClasses[], double forceConstant, double distance, BondFunction bondFunction) {
        super(ForceField.ForceFieldType.BOND, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
        this.distance = distance;
        this.bondFunction = bondFunction;
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
     * Nicely formatted bond stretch string.
     */
    @Override
    public String toString() {
        return String.format("bond  %5d  %5d  %6.1f  %7.4f", atomClasses[0],
                atomClasses[1], forceConstant, distance);
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
     * Convert bond stretch energy to kcal/mole.
     */
    public static final double units = 1.0;
    /**
     * Cubic coefficient in bond stretch potential.
     */
    public static final double cubic = -2.55;
    /**
     * Quartic coefficient in bond stretch potential.
     */
    public static final double quartic = 3.793125;

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String key1, String key2) {
        String keys1[] = key1.split(" ");
        String keys2[] = key2.split(" ");
        int c1[] = new int[2];
        int c2[] = new int[2];
        for (int i = 0; i < 2; i++) {
            c1[i] = Integer.parseInt(keys1[i]);
            c2[i] = Integer.parseInt(keys2[i]);
        }

        if (c1[0] < c2[0]) {
            return -1;
        } else if (c1[0] > c2[0]) {
            return 1;
        } else if (c1[1] < c2[1]) {
            return -1;
        } else if (c1[1] > c2[1]) {
            return 1;
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
        if (other == null || !(other instanceof BondType)) {
            return false;
        }
        BondType bondType = (BondType) other;
        int c[] = bondType.atomClasses;
        if (c[0] == atomClasses[0] && c[1] == atomClasses[1]) {
            return true;
        }
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }
}
