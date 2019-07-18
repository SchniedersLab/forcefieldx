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
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.pow;

import ffx.potential.parameters.ForceField.ForceFieldType;

/**
 * The AngleType class defines one harmonic angle bend energy term.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class AngleType extends BaseType implements Comparator<String> {

    /**
     * Angle function types include harmonic or sextic.
     */
    public enum AngleFunction {

        HARMONIC, SEXTIC
    }

    /**
     * Angle modes include Normal or In-Plane
     */
    public enum AngleMode {

        NORMAL, IN_PLANE
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
    public static final double units = 1.0 / pow(180.0 / PI, 2.0);
    /**
     * Atom classes that for this Angle type.
     */
    public final int[] atomClasses;
    /**
     * Force constant (Kcal/mole/radian^2).
     */
    public final double forceConstant;
    /**
     * Equilibrium angle (degrees). There can be up to three equilibrium angles,
     * depending on the number of attached hydrogens (0, 1, or 2).
     */
    public final double[] angle;
    /**
     * The angle function in use.
     */
    public final AngleFunction angleFunction;
    /**
     * The angle mode in use.
     */
    public AngleMode angleMode;

    /**
     * <p>
     * Constructor for normal AngleType.</p>
     *
     * @param atomClasses   an array of int.
     * @param forceConstant a double.
     * @param angle         an array of double.
     * @param angleFunction the AngleFunction to apply.
     */
    public AngleType(int[] atomClasses, double forceConstant, double[] angle, AngleFunction angleFunction) {
        super(ForceFieldType.ANGLE, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
        this.angle = angle;
        this.angleFunction = angleFunction;
        this.angleMode = AngleMode.NORMAL;
    }

    /**
     * <p>
     * Constructor for In-Plane AngleType.</p>
     *
     * @param atomClasses   an array of int.
     * @param forceConstant a double.
     * @param angle         an array of double.
     * @param angleFunction the AngleFunction to apply.
     */
    public AngleType(int[] atomClasses, double forceConstant, double[] angle,
                     AngleFunction angleFunction, AngleMode angleMode) {
        super(ForceFieldType.ANGLEP, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
        this.angle = angle;
        this.angleFunction = angleFunction;
        this.angleMode = angleMode;
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
     * @return a {@link ffx.potential.parameters.AngleType} object.
     */
    public AngleType patchClasses(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        int len = atomClasses.length;

        // Look for new AngleTypes that contain 1 or 2 mapped atom classes.
        for (AtomType newType : typeMap.keySet()) {
            for (int atomClass : atomClasses) {
                if (atomClass == newType.atomClass) {
                    count++;
                }
            }
        }

        // If found, create a new AngleType that bridges to known classes.
        if (count == 1 || count == 2) {
            int[] newClasses = Arrays.copyOf(atomClasses, len);
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < len; i++) {
                    if (atomClasses[i] == newType.atomClass) {
                        AtomType knownType = typeMap.get(newType);
                        newClasses[i] = knownType.atomClass;
                    }
                }
            }
            return new AngleType(newClasses, forceConstant, angle, angleFunction);
        }

        return null;
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
     * Average two AngleType instances. The atom classes that define the new
     * type must be supplied.
     *
     * @param angleType1  a {@link ffx.potential.parameters.AngleType} object.
     * @param angleType2  a {@link ffx.potential.parameters.AngleType} object.
     * @param atomClasses an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.AngleType} object.
     */
    public static AngleType average(AngleType angleType1, AngleType angleType2, int[] atomClasses) {
        if (angleType1 == null || angleType2 == null || atomClasses == null) {
            return null;
        }
        AngleFunction angleFunction = angleType1.angleFunction;
        if (angleFunction != angleType2.angleFunction) {
            return null;
        }
        int len = angleType1.angle.length;
        if (len != angleType2.angle.length) {
            return null;
        }
        double forceConstant = (angleType1.forceConstant + angleType2.forceConstant) / 2.0;
        double[] angle = new double[len];
        for (int i = 0; i < len; i++) {
            angle[i] = (angleType1.angle[i] + angleType2.angle[i]) / 2.0;
        }

        return new AngleType(atomClasses, forceConstant, angle, angleFunction);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted Angle bending string.
     */
    @Override
    public String toString() {
        StringBuilder angleString;
        if (angleMode == AngleMode.NORMAL) {
            angleString = new StringBuilder(format(
                    "angle  %5d  %5d  %5d  %6.2f", atomClasses[0], atomClasses[1],
                    atomClasses[2], forceConstant));
        } else {
            angleString = new StringBuilder(format(
                    "anglep  %5d  %5d  %5d  %6.2f", atomClasses[0], atomClasses[1],
                    atomClasses[2], forceConstant));
        }
        for (double eq : angle) {
            angleString.append(format("  %6.2f", eq));
        }
        return angleString.toString();
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
        if (!(other instanceof AngleType)) {
            return false;
        }
        AngleType angleType = (AngleType) other;
        int[] c = angleType.atomClasses;

        return (c[0] == atomClasses[0] && c[1] == atomClasses[1] && c[2] == atomClasses[2]);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 37 * hash + Arrays.hashCode(atomClasses);
        return hash;
    }
}
