/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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

import static java.lang.Math.sqrt;

/**
 * Thermostat a molecular dynamics trajectory to an external bath
 * using the Berendensen weak-coupling thermostat.
 *
 * @author Michael J. Schnieders
 *         derived from TINKER temperature control by Alan Grossfield
 *         and Jay Ponder
 *
 * @see H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
 *      A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
 *      to an External Bath", Journal of Chemical Physics, 81, 3684-3690 (1984)
 */
public class Berendsen extends Thermostat {

    private double tau;

    public Berendsen(int n, double v[], double mass[], double targetTemperature,
            double tau) {
        super(n, v, mass, targetTemperature);
        this.tau = tau;
    }

    public Berendsen(int n, double v[], double mass[], double targetTemperature) {
        this(n, v, mass, targetTemperature, 0.2e0);
    }

    public void setTau(double tau) {
        this.tau = tau;
    }

    public double getTau() {
        return tau;
    }

    /**
     * No velocity modifications are made by the Berendesen method at
     * the half-step.
     */
    @Override
    public void halfStep(double dt) {
        return;
    }

    /**
     * Full step velocity modification.
     */
    @Override
    public void fullStep(double dt) {
        kineticEnergy();
        double ratio = targetTemperature/currentTemperature;
        double scale = sqrt(1.0 + (dt/tau)*(ratio-1.0));
        for (int i=0; i<3*n; i++) {
            v[i] *= scale;
        }
    }

}
