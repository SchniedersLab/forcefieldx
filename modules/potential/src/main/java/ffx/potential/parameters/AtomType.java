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

import java.util.Comparator;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import static ffx.potential.parameters.ForceField.ForceFieldType.ATOM;

/**
 * The AtomType class represents one molecular mechanics atom type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class AtomType extends BaseType implements Comparator<String> {

    /**
     * A Logger for the AngleType class.
     */
    private static final Logger logger = Logger.getLogger(AtomType.class.getName());
    /**
     * Short name (ie CH3/CH2 etc).
     */
    public final String name;
    /**
     * Description of the atom's bonding environment.
     */
    public final String environment;
    /**
     * Atomic Number.
     */
    public final int atomicNumber;
    /**
     * Atomic weight. "An atomic weight (relative atomic atomicWeight) of an
     * element from a specified source is the ratio of the average atomicWeight
     * per atom of the element to 1/12 of the atomicWeight of an atom of 12C"
     */
    public final double atomicWeight;
    /**
     * Valence number for this type.
     */
    public final int valence;
    /**
     * Atom type.
     */
    public int type;
    /**
     * Atom class.
     */
    public int atomClass;

    /**
     * AtomType Constructor.
     *
     * @param type         int
     * @param atomClass    int
     * @param name         String
     * @param environment  String
     * @param atomicNumber int
     * @param atomicWeight double
     * @param valence      int
     */
    public AtomType(int type, int atomClass, String name, String environment,
                    int atomicNumber, double atomicWeight, int valence) {
        super(ATOM, Integer.toString(type));
        this.type = type;
        this.atomClass = atomClass;
        this.name = name;
        this.environment = environment;
        this.atomicNumber = atomicNumber;
        this.atomicWeight = atomicWeight;
        this.valence = valence;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String s1, String s2) {
        int t1 = parseInt(s1);
        int t2 = parseInt(s2);
        return Integer.compare(t1, t2);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        AtomType atomType = (AtomType) o;
        return atomType.type == this.type;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return Objects.hash(type);
    }

    /**
     * Construct an AtomType from an input string.
     *
     * @param input  The overall input String.
     * @param tokens The input String tokenized.
     * @return an AtomType instance.
     */
    public static AtomType parse(String input, String[] tokens) {
        if (tokens.length < 7) {
            logger.log(Level.WARNING, "Invalid ATOM type:\n{0}", input);
        } else {
            try {
                int index = 1;
                // Atom Type
                int type = parseInt(tokens[index++]);
                // Atom Class
                int atomClass;
                // The following try/catch is a nasty hack to check for one of the
                // the following two cases:
                //
                // NUMBER TYPE CLASS IDENTIFIER ... (example is OPLSAA)
                // vs.
                // NUMBER TYPE IDENTIFIER ... (example is OPLSUA)
                //
                // If there is no atom class, a harmless exception will be caught
                // and the atomClass field will remain equal to null.
                try {
                    atomClass = parseInt(tokens[index]);
                    // If the parseInt succeeds, this force field has atom classes.
                    index++;
                } catch (NumberFormatException e) {
                    // Some force fields do not use atom classes.
                    atomClass = -1;
                }
                // Name
                String name = tokens[index].intern();
                // The "environment" string may contain spaces,
                // and is therefore surrounded in quotes located at "first" and
                // "last".
                int first = input.indexOf("\"");
                int last = input.lastIndexOf("\"");
                if (first >= last) {
                    logger.log(Level.WARNING, "Invalid ATOM type:\n{0}", input);
                    return null;
                }
                // Environment
                String environment = input.substring(first, last + 1).intern();
                // Shrink the tokens array to only include entries
                // after the environment field.
                tokens = input.substring(last + 1).trim().split(" +");
                index = 0;
                // Atomic Number
                int atomicNumber = parseInt(tokens[index++]);
                // Atomic Mass
                double mass = parseDouble(tokens[index++]);
                // Hybridization
                int hybridization = parseInt(tokens[index]);
                return new AtomType(type, atomClass, name, environment, atomicNumber, mass, hybridization);
            } catch (NumberFormatException e) {
                String message = "Exception parsing CHARGE type:\n" + input + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
        return null;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted atom type string.
     */
    @Override
    public String toString() {
        String s;
        if (atomClass >= 0) {
            s = format("atom  %5d  %5d  %-4s  %-25s  %3d  %8.4f  %d",
                    type, atomClass, name, environment, atomicNumber, atomicWeight, valence);
        } else {
            s = format("atom  %5d  %-4s  %-25s  %3d  %8.4f  %d", type,
                    name, environment, atomicNumber, atomicWeight, valence);
        }
        return s;
    }

    /**
     * <p>
     * incrementClassAndType</p>
     *
     * @param classIncrement a int.
     * @param typeIncrement  a int.
     */
    void incrementClassAndType(int classIncrement, int typeIncrement) {
        atomClass += classIncrement;
        type += typeIncrement;
        setKey(Integer.toString(type));
    }
}
