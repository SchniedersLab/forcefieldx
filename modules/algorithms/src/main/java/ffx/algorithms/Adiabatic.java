/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
package ffx.algorithms;

import ffx.numerics.Potential.VARIABLE_TYPE;

/**
 * The Adiabatic thermostat is for NVE simulations and does not alter
 * particle velocities.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class Adiabatic extends Thermostat {

    /**
     * <p>Constructor for Adiabatic.</p>
     *
     * @param n a int.
     * @param x an array of double.
     * @param v an array of double.
     * @param mass an array of double.
     */
    public Adiabatic(int n, double x[], double v[], double mass[], VARIABLE_TYPE type[]) {
        super(n, x, v, mass, type, 0.0);
        this.name = Thermostats.ADIABATIC;
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return String.format("%s thermostat", name);
    }

    /**
     * {@inheritDoc}
     *
     * No half-step velocity modifications are made.
     */
    @Override
    public void halfStep(double dt) {
        return;
    }

    /**
     * {@inheritDoc}
     *
     * No full-step velocity modifications are made.
     */
    @Override
    public void fullStep(double dt) {
        return;
    }
}
