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

import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The abstract Thermostat class implements methods common to all thermostats
 * for initalizing velocities from a Maxwell-Boltzmann distribution and
 * computing the instantaneous temperature. Abstract methods are declared for
 * half-step and full-step modification of velocities thermostat implementions.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public abstract class Thermostat {

    private static final Logger logger = Logger.getLogger(ffx.algorithms.Thermostat.class.getName());
    /**
     * Boltzmann constant in units of g*Ang**2/ps**2/mole/K
     */
    public static final double kB = 0.83144725;
    /**
     * Conversion from kcal/mole to g*Ang**2/ps**2.
     */
    public static final double convert = 4.1840e2;

    public enum Thermostats { BERENDSEN, BUSSI, ISOTHERMAL };

    protected double targetTemperature;
    protected double currentTemperature;
    protected double kT;
    protected double kineticEnergy;
    protected int n;
    protected int dof;
    protected double v[];
    protected double mass[];

    public Thermostat(int n, double v[], double mass[], double t) {
        this.n = n;
        this.dof = n * 3;
        this.v = v;
        this.mass = mass;
        assert (v.length == dof);
        assert (mass.length == n);
        setTargetTemperature(t);
    }

    protected void log(Level level) {
        if (logger.isLoggable(level)) {
            logger.log(level, String.format("\nThermostat target temperature %6.2.", targetTemperature));
            logger.log(level, String.format("\nkT per degree of freedom:  %7.3f\n", convert * kineticEnergy / (dof * kT)));
        }
    }

    public double getCurrentTemperture() {
        return currentTemperature;
    }

    public double getKineticEnergy() {
        return kineticEnergy;
    }

    public double getTargetTemperature() {
        return targetTemperature;
    }

    /**
     * Set the target temperature.
     *
     * @param t Target temperature must be greater than absolute zero.
     *
     * @since 1.0
     */
    public void setTargetTemperature(double t) {
        /**
         * Obey the Third Law of Thermodynamics.
         */
        assert(t > 0.0);
        targetTemperature = t;
        kT = t * kB;
    }

    /**
     * Rest velocities from a Maxwell-Boltzmann distribution of momenta.
     * The varience of each independent momentum component is kT * mass.
     */
    public void maxwell() {
        Random random = new Random();
        for (int index = 0, i = 0; i < n; i++) {
            double m = mass[i];
            double variance = m * kT;
            double sd = Math.sqrt(variance);
            v[index++] = random.nextGaussian() * sd / m;
            v[index++] = random.nextGaussian() * sd / m;
            v[index++] = random.nextGaussian() * sd / m;
        }
        kineticEnergy();
    }

    /**
     * Compute the current temperature and kinetic energy of the system.
     */
    protected void kineticEnergy() {
        double e = 0.0;
        for (int index = 0, i = 0; i < n; i++) {
            double velocity = v[index++];
            double v2 = velocity * velocity;
            velocity = v[index++];
            v2 += velocity * velocity;
            velocity = v[index++];
            v2 += velocity * velocity;
            e += mass[i] * v2;
        }
        e *= 0.5 / convert;
        currentTemperature = convert * e / (0.5 * kB * dof);
        kineticEnergy = e;
    }

    public abstract void halfStep(double dt);
    public abstract void fullStep(double dt);

}
