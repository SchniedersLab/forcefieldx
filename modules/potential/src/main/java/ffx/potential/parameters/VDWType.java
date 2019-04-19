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

/**
 * The VDWType class defines a van der Waals type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
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
     * Reduction factor for evaluating van der Waals pairs. Valid range: 0.0
     * .GT. reduction .LE. 1.0 Usually only hydrogen atoms have a reduction
     * factor. Setting the reduction to .LT. 0.0 indicates it is not being used.
     */
    public final double reductionFactor;

    /**
     * van der Waals constructor. If the reduction factor is .LE. 0.0, no
     * reduction is used for this atom type.
     *
     * @param atomClass       The atom class that uses this van der Waals parameter.
     * @param radius          The radius of the minimum well depth energy (angstroms).
     * @param wellDepth       The minimum energy of the vdw function (kcal/mol).
     * @param reductionFactor Reduction factor for evaluating van der Waals pairs.
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
     * <p>
     * incrementClass</p>
     *
     * @param increment a int.
     */
    void incrementClass(int increment) {
        atomClass += increment;
        setKey(Integer.toString(atomClass));
    }

    /**
     * <p>average.</p>
     *
     * @param vdwType1  a {@link ffx.potential.parameters.VDWType} object.
     * @param vdwType2  a {@link ffx.potential.parameters.VDWType} object.
     * @param atomClass a int.
     * @return a {@link ffx.potential.parameters.VDWType} object.
     */
    public static VDWType average(VDWType vdwType1, VDWType vdwType2, int atomClass) {
        if (vdwType1 == null || vdwType2 == null) {
            return null;
        }

        double radius = (vdwType1.radius + vdwType2.radius) / 2.0;
        double wellDepth = (vdwType1.wellDepth + vdwType2.wellDepth) / 2.0;
        double reductionFactor = (vdwType1.reductionFactor + vdwType2.reductionFactor) / 2.0;

        return new VDWType(atomClass, radius, wellDepth, reductionFactor);
    }

    /**
     * {@inheritDoc}
     * <p>
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
        if (!(other instanceof VDWType)) {
            return false;
        }
        VDWType vdwType = (VDWType) other;
        return (vdwType.atomClass == this.atomClass);
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
