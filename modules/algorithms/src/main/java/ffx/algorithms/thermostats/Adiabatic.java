/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
package ffx.algorithms.thermostats;

import ffx.numerics.Potential.VARIABLE_TYPE;

/**
 * The Adiabatic thermostat is for NVE simulations and does not alter particle
 * velocities.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class Adiabatic extends Thermostat {

    /**
     * <p>
     * Constructor for Adiabatic.</p>
     *
     * @param n a int.
     * @param x an array of double.
     * @param v an array of double.
     * @param mass an array of double.
     * @param type the VARIABLE_TYPE of each variable.
     */
    public Adiabatic(int n, double x[], double v[], double mass[], VARIABLE_TYPE type[]) {
        super(n, x, v, mass, type, 0.0);
        this.name = ThermostatEnum.ADIABATIC;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return " Adiabatic thermostat";
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
