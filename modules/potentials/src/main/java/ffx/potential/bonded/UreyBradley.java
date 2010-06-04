/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.bonded;

import static ffx.numerics.VectorMath.*;
import static ffx.potential.parameters.UreyBradleyType.*;

import java.util.logging.Logger;

import ffx.potential.parameters.UreyBradleyType;

/**
 * The UreyBradley class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class UreyBradley extends BondedTerm implements Comparable<UreyBradley> {

    private static final Logger logger = Logger.getLogger(UreyBradley.class.getName());
    private static final long serialVersionUID = 1L;
    /**
     * Force field parameters to compute the Stretch-Bend energy.
     */
    public UreyBradleyType ureyBradleyType = null;
    /**
     * The Angle this UreyBradley term is based on.
     */
    protected Angle angle = null;

    public void setUreyBradleyType(UreyBradleyType a) {
        ureyBradleyType = a;
    }

    /**
     * Constructor for the UreyBradley class.
     */
    public UreyBradley(Angle a) {
        super();
        angle = a;
        bonds = a.bonds;
        atoms = a.atoms;
        setID_Key(false);
    }

    /**
     * Update recomputes the UreyBradley's value and energy.
     */
    @Override
    public void update() {
        energy(false);
    }
    /**
     * The vector from Atom 2 to Atom 0.
     */
    protected static final double v20[] = new double[3];
    /**
     * Gradient on Atom 0.
     */
    protected static final double g0[] = new double[3];
    /**
     * Gradient on Atom 2.
     */
    protected static final double g2[] = new double[3];

    /**
     * Evaluate the Urey-Bradley energy.
     *
     * @param gradient
     *            Evaluate the gradient.
     * @return Returns the energy.
     */
    public double energy(boolean gradient) {
        diff(atoms[0].getXYZ(), atoms[2].getXYZ(), v20);
        value = r(v20);
        double dv = value - ureyBradleyType.distance;
        double dv2 = dv * dv;
        energy = units * ureyBradleyType.forceConstant * dv2 * (1.0 + cubic * dv + quartic * dv2);
        if (gradient) {
            double deddt = 2.0 * units * ureyBradleyType.forceConstant * dv * (1.0 + 1.5 * cubic * dv + 2.0 * quartic * dv2);
            double de = 0.0;
            if (value > 0.0) {
                de = deddt / value;
            }
            scalar(v20, de, g0);
            scalar(v20, -de, g2);
            atoms[0].addToXYZGradient(g0[0], g0[1], g0[2]);
            atoms[2].addToXYZGradient(g2[0], g2[1], g2[2]);
        }
        return energy;
    }

    /*
     * Log details for this Angle energy term.
     */
    public void log() {
        logger.info(String.format(" %s %6d-%s %6d-%s %6.4f  %6.4f  %10.4f",
                "Urey-Bradley", atoms[0].getXYZIndex(), atoms[0].getAtomType().name,
                atoms[2].getXYZIndex(), atoms[2].getAtomType().name, ureyBradleyType.distance, value,
                energy));
    }

    @Override
    public int compareTo(UreyBradley ub) {
        return angle.compareTo(ub.angle);
    }
}
