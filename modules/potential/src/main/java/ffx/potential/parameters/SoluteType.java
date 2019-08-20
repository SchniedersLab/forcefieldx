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

import ffx.potential.parameters.ForceField.ForceFieldType;

/**
 * The SoluteType class defines one implicit solvent radius.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class SoluteType extends BaseType implements Comparator<String> {

    /**
     * Atom class for this solute type.
     */
    private int atomClass;

    /**
     * Solute atomic diameter.
     */
    public final double diameter;

    /**
     * Solute atom chemical description.
     */
    public final String description;

    /**
     * <p>Constructor for SoluteType.</p>
     *
     * @param atomClass    a int.
     * @param diameter a double.
     */
    public SoluteType(int atomClass, String description, double diameter) {
        super(ForceFieldType.SOLUTE, Integer.toString(atomClass));
        this.atomClass = atomClass;
        this.diameter = diameter;
        this.description = description;
    }

    /**
     * <p>
     * incrementClasses</p>
     *
     * @param increment a int.
     */
    void incrementType(int increment) {
        atomClass += increment;
        setKey(Integer.toString(atomClass));
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted bond stretch string.
     */
    @Override
    public String toString() {
        return String.format("solute  %4d  %30s  %7.5f", atomClass, description, diameter);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String key1, String key2) {
        int type1 = Integer.parseInt(key1);
        int type2 = Integer.parseInt(key2);

        return Integer.compare(type1, type2);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object other) {
        if (other == this) {
            return true;
        }
        if (!(other instanceof SoluteType)) {
            return false;
        }
        SoluteType otherType = (SoluteType) other;
        return (otherType.atomClass == this.atomClass);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        int hash = 7;
        hash = 97 * hash + this.atomClass;
        return hash;
    }
}
