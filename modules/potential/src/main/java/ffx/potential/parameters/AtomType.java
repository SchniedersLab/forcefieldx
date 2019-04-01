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

import java.util.Comparator;

/**
 * The AtomType class represents one molecular mechanics atom type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class AtomType extends BaseType implements Comparator<String> {

    /**
     * Atom type.
     */
    public int type;
    /**
     * Atom class.
     */
    public int atomClass;
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
        super(ForceField.ForceFieldType.ATOM, Integer.toString(type));
        this.type = type;
        this.atomClass = atomClass;
        this.name = name;
        this.environment = environment;
        this.atomicNumber = atomicNumber;
        this.atomicWeight = atomicWeight;
        this.valence = valence;
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

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted atom type string.
     */
    @Override
    public String toString() {
        String s;
        if (atomClass >= 0) {
            s = String.format("atom  %5d  %5d  %-4s  %-25s  %3d  %8.4f  %d",
                    type, atomClass, name, environment, atomicNumber, atomicWeight, valence);
        } else {
            s = String.format("atom  %5d  %-4s  %-25s  %3d  %8.4f  %d", type,
                    name, environment, atomicNumber, atomicWeight, valence);
        }
        return s;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String s1, String s2) {

        int t1 = Integer.parseInt(s1);
        int t2 = Integer.parseInt(s2);

        return Integer.compare(t1, t2);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof AtomType)) {
            return false;
        }
        AtomType atomType = (AtomType) other;

        return atomType.type == this.type;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 71 * hash + this.type;
        return hash;
    }
}
