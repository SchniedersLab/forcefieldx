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
import static java.lang.String.format;

/**
 * The ChargeType class defines a partial atomic charge type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class ChargeType extends BaseType implements Comparator<String> {

    /**
     * The atom type that uses this charge parameter.
     */
    public int atomType;
    /**
     * Partial atomic charge in units of electrons.
     */
    public final double charge;

    /**
     * ChargeType constructor.
     *
     * @param atomType int
     * @param charge   double
     */
    public ChargeType(int atomType, double charge) {
        super(ForceField.ForceFieldType.CHARGE, "" + atomType);
        this.atomType = atomType;
        this.charge = charge;
    }

    /**
     * <p>
     * incrementType</p>
     *
     * @param increment a int.
     */
    public void incrementType(int increment) {
        this.atomType += increment;
    }

    /**
     * Average two ChargeType instances. The atom type that defines the
     * new type must be supplied.
     *
     * @param chargeType1 a {@link ffx.potential.parameters.ChargeType} object.
     * @param chargeType2 a {@link ffx.potential.parameters.ChargeType} object.
     * @param atomType    a int.
     * @return a {@link ffx.potential.parameters.ChargeType} object.
     */
    public static ChargeType average(ChargeType chargeType1, ChargeType chargeType2, int atomType) {
        if (chargeType1 == null || chargeType2 == null) {
            return null;
        }
        double charge = (chargeType1.charge + chargeType2.charge) / 2.0;

        return new ChargeType(atomType, charge);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Nicely formatted Charge type.
     */
    @Override
    public String toString() {
        return format("charge  %5d  % 7.5f", atomType, charge);
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

}
