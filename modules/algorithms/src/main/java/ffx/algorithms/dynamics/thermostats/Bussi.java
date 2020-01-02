//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.algorithms.dynamics.thermostats;

import java.util.Collections;
import java.util.List;
import java.util.Random;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.Constraint;
import ffx.numerics.Potential.VARIABLE_TYPE;

/**
 * Thermostat a molecular dynamics trajectory to an external bath using the
 * Bussi, Donadio, and Parrinello method. This method is similar to Berendsen
 * thermostat, but generates a canonical distribution.
 *
 * @author Michael J. Schnieders
 * <p>
 * Derived from TINKER temperature control by Alan Grossfield and Jay Ponder.
 * @see <a href="http://dx.doi.org/10.1016/j.cpc.2008.01.006"> G. Bussi and M.
 * Parrinello, "Stochastic Thermostats: Comparison of Local and Global Schemes",
 * Computer Physics Communications, 179, 26-29 (2008)</a>
 * @since 1.0
 */
public class Bussi extends Thermostat {

    /**
     * Bussi thermostat time constant (psec).
     */
    private double tau;
    /**
     * The random number generator used to perturb velocities.
     */
    private final Random bussiRandom;

    /**
     * <p>
     * Constructor for Bussi.</p>
     *
     * @param n                 Number of degrees of freedom.
     * @param x                 Atomic coordinates.
     * @param v                 Velocities.
     * @param mass              Mass of each degrees of freedom.
     * @param type              the VARIABLE_TYPE of each variable.
     * @param targetTemperature The target temperature.
     * @param tau               Bussi thermostat time constant (psec).
     */
    public Bussi(int n, double[] x, double[] v, double[] mass,
                 VARIABLE_TYPE[] type, double targetTemperature,
                 double tau) {
        this(n, x, v, mass, type, targetTemperature, tau, Collections.emptyList());
    }

    public Bussi(int n, double[] x, double[] v, double[] mass,
                 VARIABLE_TYPE[] type, double targetTemperature,
                 double tau, List<Constraint> constraints) {
        super(n, x, v, mass, type, targetTemperature, constraints);
        this.name = ThermostatEnum.BUSSI;
        this.tau = tau;
        this.bussiRandom = new Random();
    }

    /**
     * <p>
     * Constructor for Bussi.</p>
     *
     * @param n                 Number of degrees of freedom.
     * @param x                 Atomic coordinates.
     * @param v                 Velocities.
     * @param mass              Mass of each degrees of freedom.
     * @param type              the VARIABLE_TYPE of each variable.
     * @param targetTemperature a double.
     */
    public Bussi(int n, double[] x, double[] v, double[] mass,
                 VARIABLE_TYPE[] type, double targetTemperature) {
        this(n, x, v, mass, type, targetTemperature, 0.2e0);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Initialize the Random number generator used to apply random forces to the
     * particles.
     */
    public void setRandomSeed(long seed) {
        bussiRandom.setSeed(seed);
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
        return format(" Bussi Thermostat (tau = %8.3f psec)", tau);
    }

    /**
     * {@inheritDoc}
     * <p>
     * No velocity modifications are made by the Bussi method at the half-step.
     */
    @Override
    public void halfStep(double dt) {
    }

    /**
     * {@inheritDoc}
     * <p>
     * Full step velocity modification.
     */
    @Override
    public void fullStep(double dt) {
        double expTau = exp(-dt / tau);
        double tempRatio = targetTemperature / currentTemperature;
        double rate = (1.0 - expTau) * tempRatio / dof;
        double r = bussiRandom.nextGaussian();
        double s = 0.0;
        for (int i = 0; i < dof - 1; i++) {
            double si = bussiRandom.nextGaussian();
            s += si * si;
        }
        double scale = expTau + (s + r * r) * rate + 2.0 * r * sqrt(expTau * rate);
        scale = sqrt(scale);
        if (r + sqrt(expTau / rate) < 0.0) {
            scale = -scale;
        }
        for (int i = 0; i < nVariables; i++) {
            v[i] *= scale;
        }
    }
}
