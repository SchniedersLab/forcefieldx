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

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import static java.lang.Integer.parseInt;

import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.toRadians;

/**
 * The ImproperTorsionType class defines an improper torsion.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class ImproperTorsionType extends BaseType implements Comparator<String> {

    /**
     * Atom classes that for this Improper Torsion angle.
     */
    public final int[] atomClasses;
    /**
     * Force constant in kcal/mol.
     */
    public final double k;
    /**
     * Phases in degrees.
     */
    public final double phase;
    /**
     * Periodicity (should be 2 for an Improper Torsion).
     */
    public final int periodicity;
    /**
     * Value of cos(toRadians(phase)).
     */
    public final double cos;
    /**
     * Value of sin(toRadians(phase)).
     */
    public final double sin;

    /**
     * TorsionType Constructor.
     *
     * @param atomClasses Atom classes.
     * @param k           Force constant.
     * @param phase       The phase.
     * @param periodicity The periodicity.
     */
    public ImproperTorsionType(int[] atomClasses, double k, double phase, int periodicity) {
        super(ForceField.ForceFieldType.IMPTORS, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        double symm = 1.0;
        this.periodicity = periodicity;
        this.k = k / symm;
        this.phase = phase;
        cos = cos(toRadians(phase));
        sin = sin(toRadians(phase));

        assert (periodicity == 2);
    }

    /**
     * Returns true if the atoms can be assigned this improperTorsionType.
     *
     * @param inputClasses The atom classes will be re-ordered if its member
     *                     atoms match this ImproperTorsionType. The trigonal atom will not change
     *                     position.
     * @param allowInitialWildCards Allow wildcard match to first two classes.
     * @param allowFinalWildCard Allow wildcard match for final class.
     *
     * @return True if this torsionType is assignable to the atom array.
     */
    public boolean assigned(int[] inputClasses, boolean allowInitialWildCards, boolean allowFinalWildCard) {
        // Assign the trigonal atom.
        if (inputClasses[2] != atomClasses[2]) {
            return false;
        }

        // Assign the final atom.
        if (inputClasses[3] == atomClasses[3] || (atomClasses[3] == 0 && allowFinalWildCard)) {
            // do nothing.
        } else if (inputClasses[1] == atomClasses[3]) {
            int temp = inputClasses[3];
            inputClasses[3] = inputClasses[1];
            inputClasses[1] = temp;
        } else if (inputClasses[0] == atomClasses[3]) {
            int temp = inputClasses[3];
            inputClasses[3] = inputClasses[0];
            inputClasses[0] = temp;
        } else {
            return false;
        }

        // Assign the second atom.
        if (inputClasses[1] == atomClasses[1] || (atomClasses[1] == 0 && allowInitialWildCards)) {
            // Do nothing.
        } else if (inputClasses[0] == atomClasses[1]) {
            int temp = inputClasses[1];
            inputClasses[1] = inputClasses[0];
            inputClasses[0] = temp;
        } else {
            return false;
        }

        // Assign the first atom.
        return (inputClasses[0] == atomClasses[0] || (atomClasses[0] == 0 && allowInitialWildCards));
    }

    /**
     * Check if this Improper Torsion Type is defined by 1 or more atom classes equal to zero.
     *
     * @return True if there are no zero "wildcard" atom classes for this type.
     */
    public boolean noZeroClasses() {
        if (atomClasses[0] != 0 && atomClasses[1] != 0 && atomClasses[3] != 0) {
            return true;
        }

        return false;
    }

    /**
     * <p>
     * incrementClasses</p>
     *
     * @param increment a int.
     */
    public void incrementClasses(int increment) {
        for (int i = 0; i < atomClasses.length; i++) {
            if (atomClasses[i] != 0) {
                atomClasses[i] += increment;
            }
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
     * This method sorts the atom classes for the improper torsion.
     *
     * @param c atomClasses
     * @return lookup key
     * @since 1.0
     */
    public static String sortKey(int[] c) {
        if (c == null || c.length != 4) {
            return null;
        }
        return c[0] + " " + c[1] + " " + c[2] + " " + c[3];
    }

    /**
     * Average two ImproperTorsionType instances. The atom classes that define the
     * new type must be supplied.
     *
     * @param improperTorsionType1 a {@link ffx.potential.parameters.ImproperTorsionType} object.
     * @param improperTorsionType2 a {@link ffx.potential.parameters.ImproperTorsionType} object.
     * @param atomClasses          an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.ImproperTorsionType} object.
     */
    public static ImproperTorsionType average(ImproperTorsionType improperTorsionType1,
                                              ImproperTorsionType improperTorsionType2, int[] atomClasses) {

        if (improperTorsionType1 == null || improperTorsionType2 == null || atomClasses == null) {
            return null;
        }

        int periodicity = improperTorsionType1.periodicity;
        if (periodicity != improperTorsionType2.periodicity) {
            return null;
        }

        double forceConstant = (improperTorsionType1.k + improperTorsionType2.k) / 2.0;
        double phase = (improperTorsionType1.phase + improperTorsionType2.phase) / 2.0;

        return new ImproperTorsionType(atomClasses, forceConstant, phase, periodicity);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted Torsion angle.
     *
     * @since 1.0
     */
    @Override
    public String toString() {
        StringBuilder imptorsBuffer = new StringBuilder("imptors");
        for (int i : atomClasses) {
            imptorsBuffer.append(String.format(" %5d", i));
        }
        imptorsBuffer.append(String.format(" %7.3f %7.3f %1d",
                k, phase, periodicity));

        return imptorsBuffer.toString();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Implements the Comparator interface.
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
            c1[i] = parseInt(keys1[i]);
            c2[i] = parseInt(keys2[i]);
        }

        if (c1[2] < c2[2]) {
            return -1;
        } else if (c1[2] > c2[2]) {
            return 1;
        } else if (c1[0] < c2[0]) {
            return -1;
        } else if (c1[0] > c2[0]) {
            return 1;
        } else if (c1[1] < c2[1]) {
            return -1;
        } else if (c1[1] > c2[1]) {
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
     * <p>
     * Override the default <code>equals</code> method.
     *
     * @since 1.0
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof ImproperTorsionType)) {
            return false;
        }
        ImproperTorsionType improperTorsionType = (ImproperTorsionType) other;
        for (int i = 0; i < atomClasses.length; i++) {
            if (improperTorsionType.atomClasses[i] != atomClasses[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * {@inheritDoc}
     * <p>
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
