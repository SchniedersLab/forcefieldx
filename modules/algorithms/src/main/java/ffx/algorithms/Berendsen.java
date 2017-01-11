/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.algorithms;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.Potential.VARIABLE_TYPE;

/**
 * Thermostat a molecular dynamics trajectory to an external bath using the
 * Berendsen weak-coupling thermostat.
 *
 * @author Michael J. Schnieders derived from TINKER temperature control by Alan
 * Grossfield and Jay Ponder
 *
 * @see <a href="http://link.aip.org/link/?JCP/81/3684"> H. J. C. Berendsen, J.
 * P. M. Postma, W. F. van Gunsteren, A. DiNola and J. R. Hauk, "Molecular
 * Dynamics with Coupling to an External Bath", Journal of Chemical Physics, 81,
 * 3684-3690 (1984)</a>
 *
 */
public class Berendsen extends Thermostat {

    /**
     * Berendsen time constant (psec).
     */
    private double tau;

    /**
     * <p>
     * Constructor for Berendsen.</p>
     *
     * @param n a int.
     * @param x an array of double.
     * @param v an array of double.
     * @param mass an array of double.
     * @param type the VARIABLE_TYPE of each variable.
     * @param targetTemperature a double.
     * @param tau a double.
     */
    public Berendsen(int n, double x[], double v[], double mass[],
            VARIABLE_TYPE type[], double targetTemperature, double tau) {
        super(n, x, v, mass, type, targetTemperature);
        this.name = Thermostats.BERENDSEN;
        this.tau = tau;
    }

    /**
     * <p>
     * Constructor for Berendsen.</p>
     *
     * @param n a int.
     * @param x an array of double.
     * @param v an array of double.
     * @param mass an array of double.
     * @param type the VARIABLE_TYPE of each variable.
     * @param targetTemperature a double.
     */
    public Berendsen(int n, double x[], double v[], double mass[],
            VARIABLE_TYPE type[], double targetTemperature) {
        this(n, x, v, mass, type, targetTemperature, 0.2e0);
    }

    /**
     * <p>
     * Setter for the field <code>tau</code>.</p>
     *
     * @param tau a double.
     */
    public void setTau(double tau) {
        this.tau = tau;
    }

    /**
     * <p>
     * Getter for the field <code>tau</code>.</p>
     *
     * @return a double.
     */
    public double getTau() {
        return tau;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return String.format(" Berendsen Thermostat (tau = %8.3f psec)", tau);
    }

    /**
     * {@inheritDoc}
     *
     * No velocity modifications are made by the Berendsen method at the
     * half-step.
     */
    @Override
    public void halfStep(double dt) {
        return;
    }

    /**
     * {@inheritDoc}
     *
     * Full step velocity modification.
     */
    @Override
    public void fullStep(double dt) {
        double ratio = targetTemperature / currentTemperature;
        double scale = sqrt(1.0 + (dt / tau) * (ratio - 1.0));
        for (int i = 0; i < nVariables; i++) {
            v[i] *= scale;
        }
    }
}
