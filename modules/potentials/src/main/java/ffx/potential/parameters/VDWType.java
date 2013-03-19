/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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

import java.util.Comparator;

/**
 * The VDWType class defines a van der Waals type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public final class VDWType extends BaseType implements Comparator<String> {

    public enum RadiusSize {

        RADIUS, DIAMETER
    }

    public enum RadiusType {

        RMIN, SIGMA
    }
    /**
     * The atom class that uses this van der Waals parameter.
     */
    public int atomClass;
    /**
     * The radius of the minimum well depth energy (angstroms).
     */
    public final double radius;
    /**
     * The minimum energy of the vdw function (kcal/mol).
     */
    public final double wellDepth;
    /**
     * Reduction factor for evaluating van der Waals pairs. Valid range: 0.0 >
     * reduction <= 1.0 Usually only hydrogen atom have a reduction factor.
     * Setting the reduction to < 0.0 indicates it is not being used.
     */
    public final double reductionFactor;

    /**
     * van der Waals constructor. If the reduction factor is <= 0.0, no
     * reduction is used for this atom type.
     *

     *
     * @param atomClass int
     * @param radius double
     * @param wellDepth double
     * @param reductionFactor double
     */
    public VDWType(int atomClass, double radius, double wellDepth,
            double reductionFactor) {
        super(ForceField.ForceFieldType.VDW, Integer.toString(atomClass));
        this.atomClass = atomClass;
        this.radius = radius;
        this.wellDepth = wellDepth;
        this.reductionFactor = reductionFactor;
    }

    /**
     * <p>incrementClass</p>
     *
     * @param increment a int.
     */
    public void incrementClass(int increment) {
        atomClass += increment;
        setKey(Integer.toString(atomClass));
    }

    /**
     * {@inheritDoc}
     *
     * Nicely formatted van der Waals type.
     */
    @Override
    public String toString() {
        String vdwString;
        // No reduction factor.
        if (reductionFactor <= 0e0) {
            vdwString = String.format("vdw  %5d  %11.9f  %11.9f", atomClass,
                    radius, wellDepth);
        } else {
            vdwString = String.format("vdw  %5d  %11.9f  %11.9f  %5.3f",
                    atomClass, radius, wellDepth, reductionFactor);
        }
        return vdwString;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int compare(String s1, String s2) {

        int t1 = Integer.parseInt(s1);
        int t2 = Integer.parseInt(s2);

        if (t1 < t2) {
            return -1;
        }
        if (t1 > t2) {
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
        if (other == null || !(other instanceof VDWType)) {
            return false;
        }
        VDWType vdwType = (VDWType) other;
        if (vdwType.atomClass == this.atomClass) {
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
        hash = 97 * hash + this.atomClass;
        return hash;
    }
}
