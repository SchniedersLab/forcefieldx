//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.parameters;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Objects;
import static java.util.Arrays.copyOf;


/**
 * The AngleTorsionType class defines one angle-torsion energy type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class AngleTorsionType extends BaseType implements Comparator<String> {

    /**
     * Atom classes for this stretch-torsion type.
     */
    public final int[] atomClasses;
    /**
     * Force constants.
     */
    public final double[] forceConstants;
    /**
     * Convert angle-torsion to kcal/mole.
     */
    public static final double units = 1.0;

    /**
     * AngleTorsionType Constructor.
     *
     * @param atomClasses    Atomic classes.
     * @param forceConstants Force constants.
     */
    public AngleTorsionType(int[] atomClasses, double[] forceConstants) {

        // Pass the key from sorted classes to the super constructor.
        super(ForceField.ForceFieldType.ANGTORS, sortKey(atomClasses));

        this.atomClasses = atomClasses;
        this.forceConstants = forceConstants;
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
     * @return a {@link ffx.potential.parameters.AngleTorsionType} object.
     */
    public AngleTorsionType patchClasses(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        int len = atomClasses.length;

        // Check if this Type contain 1 or 2 mapped atom classes.
        for (AtomType newType : typeMap.keySet()) {

            for (int atomClass : atomClasses) {
                if (atomClass == newType.atomClass) {
                    count++;
                }
            }
        }

        // If found, create a new AngleTorsionType that bridges to known classes.
        if (count == 1 || count == 2) {
            int[] newClasses = copyOf(atomClasses, len);
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < len; i++) {
                    if (atomClasses[i] == newType.atomClass) {
                        AtomType knownType = typeMap.get(newType);
                        newClasses[i] = knownType.atomClass;
                    }
                }
            }
            return new AngleTorsionType(newClasses, forceConstants);
        }
        return null;
    }

    /**
     * This method sorts the atom classes for the angle-torsion.
     *
     * @param c atomClasses
     * @return lookup key
     * @since 1.0
     */
    public static String sortKey(int[] c) {
        return c[0] + " " + c[1] + " " + c[2] + " " + c[3];
    }

    /**
     * <p>average.</p>
     *
     * @param angleTorsionType1 a {@link ffx.potential.parameters.AngleTorsionType} object.
     * @param angleTorsionType2 a {@link ffx.potential.parameters.AngleTorsionType} object.
     * @param atomClasses       an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.AngleTorsionType} object.
     */
    public static AngleTorsionType average(AngleTorsionType angleTorsionType1,
                                           AngleTorsionType angleTorsionType2, int[] atomClasses) {
        if (angleTorsionType1 == null || angleTorsionType2 == null || atomClasses == null) {
            return null;
        }
        int len = angleTorsionType1.forceConstants.length;
        if (len != angleTorsionType2.forceConstants.length) {
            return null;
        }
        double[] forceConstants = new double[len];
        for (int i = 0; i < len; i++) {
            forceConstants[i] = (angleTorsionType1.forceConstants[i]
                    + angleTorsionType2.forceConstants[i]) / 2.0;
        }
        return new AngleTorsionType(atomClasses, forceConstants);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted Angle-Torsion string.
     */
    @Override
    public String toString() {
        return String.format("angtors  %5d  %5d  %5d  %5d  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f",
                atomClasses[0], atomClasses[1], atomClasses[2], atomClasses[3],
                forceConstants[0], forceConstants[1], forceConstants[2], forceConstants[3], forceConstants[4], forceConstants[5]);
    }

    /**
     * {@inheritDoc}
     *
     * @since 1.0
     */
    @Override
    public int compare(String s1, String s2) {
        String[] keys1 = s1.split(" ");
        String[] keys2 = s2.split(" ");
        int[] c1 = new int[4];
        int[] c2 = new int[4];

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
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof AngleTorsionType)) {
            return false;
        }
        AngleTorsionType stretchTorsionType = (AngleTorsionType) other;
        int[] c = stretchTorsionType.atomClasses;

        return (c[0] == atomClasses[0] && c[1] == atomClasses[1]
                && c[2] == atomClasses[2] && c[3] == atomClasses[3]);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return Objects.hash(atomClasses);
    }

}
