/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.parameters;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import static java.lang.String.format;

/**
 * The UreyBradleyType class defines one harmonic UreyBradley cross term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class UreyBradleyType extends BaseType implements Comparator<String> {

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
    /**
     * Atom classes that form this Urey-Bradley cross term.
     */
    public final int[] atomClasses;
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
     * @param atomClasses   Atom classes.
     * @param forceConstant Force constant (Kcal/mole/angstroms^2).
     * @param distance      Equilibrium 1-3 separation (Angstroms).
     */
    public UreyBradleyType(int[] atomClasses, double forceConstant,
                           double distance) {
        super(ForceField.ForceFieldType.UREYBRAD, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
        this.distance = distance;
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
            for (int atomClass : atomClasses) {
                if (atomClass == newType.atomClass) {
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
     * This method sorts the atom classes as: min, c[1], max
     *
     * @param c atomClasses
     * @return lookup key
     */
    public static String sortKey(int[] c) {
        if (c == null || c.length != 3) {
            return null;
        }
        if (c[0] > c[2]) {
            int temp = c[0];
            c[0] = c[2];
            c[2] = temp;
        }

        return c[0] + " " + c[1] + " " + c[2];
    }

    /**
     * <p>average.</p>
     *
     * @param ureyBradleyType1 a {@link ffx.potential.parameters.UreyBradleyType} object.
     * @param ureyBradleyType2 a {@link ffx.potential.parameters.UreyBradleyType} object.
     * @param atomClasses      an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.UreyBradleyType} object.
     */
    public static UreyBradleyType average(UreyBradleyType ureyBradleyType1,
                                          UreyBradleyType ureyBradleyType2, int[] atomClasses) {
        if (ureyBradleyType1 == null || ureyBradleyType2 == null || atomClasses == null) {
            return null;
        }

        double forceConstant = (ureyBradleyType1.forceConstant + ureyBradleyType2.forceConstant) / 2.0;
        double distance = (ureyBradleyType1.distance + ureyBradleyType2.distance) / 2.0;

        return new UreyBradleyType(atomClasses, forceConstant, distance);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted Urey-Bradley string.
     */
    @Override
    public String toString() {
        return format("ureybrad  %5d  %5d  %5d  %6.2f  %7.4f",
                atomClasses[0], atomClasses[1], atomClasses[2], forceConstant, distance);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String key1, String key2) {
        String[] keys1 = key1.split(" ");
        String[] keys2 = key2.split(" ");
        int[] c1 = new int[3];
        int[] c2 = new int[3];
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

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof UreyBradleyType)) {
            return false;
        }
        UreyBradleyType ureyBradleyType = (UreyBradleyType) other;
        int[] c = ureyBradleyType.atomClasses;

        return (c[0] == atomClasses[0] && c[1] == atomClasses[1] && c[2] == atomClasses[2]);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 3;
        hash = 79 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }

}
