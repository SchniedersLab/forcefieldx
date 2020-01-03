//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.util.Arrays.copyOf;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.pow;

import static ffx.potential.parameters.ForceField.ForceFieldType.OPBEND;

/**
 * The OutOfPlaneBendType class defines one Allinger style out-of-plane angle
 * bending energy type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class OutOfPlaneBendType extends BaseType implements Comparator<String> {

    /**
     * A Logger for the OutOfPlaneBendType class.
     */
    private static final Logger logger = Logger.getLogger(OutOfPlaneBendType.class.getName());

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
     * <p>
     * TINKER v.5 and v.6 Units: 1.0 / (180.0/PI)^2 = 0.00030461741979 TINKER
     * v.4 Units: 0.02191418
     * <p>
     * Ratio of v.4 to v.5/6 = 0.02191418 / 1.0 / (180.0/PI)^2 = 71.94
     */
    public static final double units = 1.0 / pow(180.0 / PI, 2);
    /**
     * Atom classes for this out-of-plane angle bending type.
     */
    public final int[] atomClasses;
    /**
     * Force constant (Kcal/mol/Angstrom).
     */
    public final double forceConstant;

    /**
     * OutOfPlaneBendType Constructor.
     *
     * @param atomClasses   int[]
     * @param forceConstant double
     */
    public OutOfPlaneBendType(int[] atomClasses, double forceConstant) {
        super(OPBEND, sortKey(atomClasses));
        this.atomClasses = atomClasses;
        this.forceConstant = forceConstant;
    }

    /**
     * <p>
     * incrementClasses</p>
     *
     * @param increment a int.
     */
    public void incrementClasses(int increment) {
        for (int i = 0; i < atomClasses.length; i++) {
            if (atomClasses[i] > 0) {
                atomClasses[i] += increment;
            }
        }
        setKey(sortKey(atomClasses));
    }

    /**
     * Remap new atom classes to known internal ones.
     *
     * @param typeMap a lookup between new atom types and known atom types.
     * @return a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     */
    public OutOfPlaneBendType patchClasses(HashMap<AtomType, AtomType> typeMap) {
        int count = 0;
        int len = atomClasses.length;

        // Look for new OutOfPlaneBends that contain 1 mapped atom classes.
        for (AtomType newType : typeMap.keySet()) {

            for (int atomClass : atomClasses) {
                if (atomClass == newType.atomClass) {
                    count++;
                }
            }
        }

        // If found, create a new OutOfPlaneBend that bridges to known classes.
        if (count == 1) {
            int[] newClasses = copyOf(atomClasses, len);
            for (AtomType newType : typeMap.keySet()) {
                for (int i = 0; i < len; i++) {
                    if (atomClasses[i] == newType.atomClass) {
                        AtomType knownType = typeMap.get(newType);
                        newClasses[i] = knownType.atomClass;
                    }
                }
            }
            return new OutOfPlaneBendType(newClasses, forceConstant);
        }
        return null;
    }

    /**
     * This method sorts the atom classes for the out-of-plane angle bending
     * type.
     *
     * @param c atomClasses
     * @return lookup key
     */
    public static String sortKey(int[] c) {
        if (c == null || c.length != 4) {
            return null;
        }
        return c[0] + " " + c[1] + " " + c[2] + " " + c[3];
    }

    /**
     * Average two OutOfPlaneBendType instances. The atom classes that define
     * the new type must be supplied.
     *
     * @param outOfPlaneBendType1 a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     * @param outOfPlaneBendType2 a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     * @param atomClasses         an array of {@link int} objects.
     * @return a {@link ffx.potential.parameters.OutOfPlaneBendType} object.
     */
    public static OutOfPlaneBendType average(OutOfPlaneBendType outOfPlaneBendType1,
                                             OutOfPlaneBendType outOfPlaneBendType2, int[] atomClasses) {
        if (outOfPlaneBendType1 == null || outOfPlaneBendType2 == null || atomClasses == null) {
            return null;
        }

        double forceConstant = (outOfPlaneBendType1.forceConstant + outOfPlaneBendType2.forceConstant) / 2.0;

        return new OutOfPlaneBendType(atomClasses, forceConstant);
    }

    /**
     * Construct an OutOfPlaneBendType from an input string.
     *
     * @param input  The overall input String.
     * @param tokens The input String tokenized.
     * @return an OutOfPlaneBendType instance.
     */
    public static OutOfPlaneBendType parse(String input, String[] tokens) {
        if (tokens.length < 6) {
            logger.log(Level.WARNING, "Invalid OPBEND type:\n{0}", input);
        } else {
            try {
                int[] atomClasses = new int[4];
                atomClasses[0] = parseInt(tokens[1]);
                atomClasses[1] = parseInt(tokens[2]);
                atomClasses[2] = parseInt(tokens[3]);
                atomClasses[3] = parseInt(tokens[4]);
                double forceConstant = parseDouble(tokens[5]);
                return new OutOfPlaneBendType(atomClasses, forceConstant);
            } catch (NumberFormatException e) {
                String message = "Exception parsing OPBEND type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted out-of-plane angle bending string.
     */
    @Override
    public String toString() {
        return String.format("opbend  %5d  %5d  %5d  %5d  %6.2f", atomClasses[0],
                atomClasses[1], atomClasses[2],
                atomClasses[3], forceConstant);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String s1, String s2) {
        String[] keys1 = s1.split(" ");
        String[] keys2 = s2.split(" ");

        for (int i = 0; i < 4; i++) {
            int c1 = parseInt(keys1[i]);
            int c2 = parseInt(keys2[i]);
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
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        OutOfPlaneBendType outOfPlaneBendType = (OutOfPlaneBendType) o;
        return Arrays.equals(atomClasses, outOfPlaneBendType.atomClasses);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(atomClasses);
    }


}
