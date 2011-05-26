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

import static java.lang.Math.exp;
import static java.lang.Math.sqrt;

import java.util.Random;

/**
 * Thermostat a molecular dynamics trajectory to an external bath
 * using the Bussi, Donadio, and Parrinello method. This method is similar to
 * Berendsen thermostat, but generates a canonical distribution.
 *
 * @author Michael J. Schnieders
 *
 *         Derived from TINKER temperature control by Alan Grossfield
 *         and Jay Ponder.
 *
 * @see <a href="http://dx.doi.org/10.1016/j.cpc.2008.01.006">
 *      G. Bussi and M. Parrinello, "Stochastic Thermostats: Comparison
 *      of Local and Global Schemes", Computer Physics Communications,
 *      179, 26-29 (2008)</a>
 *
 * @since 1.0
 */
public class Bussi extends Thermostat {

    private double tau;
    private final Random random;

    public Bussi(int dof, double x[], double v[], double mass[], double targetTemperature,
            double tau) {
        super(dof, x, v, mass, targetTemperature);
        this.name = Thermostats.BUSSI;
        this.tau = tau;
        this.random = new Random(0);
    }

    public Bussi(int dof, double x[], double v[], double mass[], double targetTemperature) {
        this(dof, x, v, mass, targetTemperature, 0.2e0);
    }

    public void setTau(double tau) {
        this.tau = tau;
    }

    public double getTau() {
        return tau;
    }

    @Override
    public String toString() {
        return String.format("%s thermostat (tau = %8.3f)", name, tau);
    }

    /**
     * No velocity modifications are made by the Bussi method at
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
        double exptau = exp(-dt / tau);
        double ratio = targetTemperature / currentTemperature;
        double rate = (1.0 - exptau) * ratio / dof;
        double r = random.nextGaussian();
        double s = 0.0;
        for (int i = 0; i < dof - 1; i++) {
            double si = random.nextGaussian();
            s += si * si;
        }
        double scale = sqrt(exptau + (s + r * r) * rate + 2.0 * r * sqrt(exptau * rate));
        if (r + sqrt(exptau / rate) < 0.0) {
            scale = -scale;
        }
        for (int i = 0; i < dof; i++) {
            v[i] *= scale;
        }
    }
}
