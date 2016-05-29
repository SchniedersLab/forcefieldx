/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import java.util.Random;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.Potential;

/**
 * Stochastic dynamics time step via a velocity Verlet integration algorithm.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class Stochastic extends Integrator {

    private double vfric[];
    private double vrand[];
    private final double friction;
    private double inverseFriction;
    private double fdt;
    private double efdt;
    private double temperature;
    private final Random random;

    /**
     * Constructor for Stochastic Dynamics.
     *
     * @param friction Friction coefficient.
     * @param nVariables Number of variables.
     * @param x Variables current value.
     * @param v Current velocities.
     * @param a Current accelerations.
     * @param mass Mass of the variables.
     */
    public Stochastic(double friction, int nVariables, double x[],
            double v[], double a[], double mass[]) {
        super(nVariables, x, v, a, mass);
        this.friction = friction;
        if (friction >= 0) {
            inverseFriction = 1.0 / friction;
        } else {
            inverseFriction = Double.POSITIVE_INFINITY;
        }
        vfric = new double[nVariables];
        vrand = new double[nVariables];
        fdt = friction * dt;
        efdt = exp(-fdt);
        temperature = 298.15;
        random = new Random();
    }
    
    /**
     * Set the stochastic dynamics time-step.
     *
     * @param dt the time step.
     */
    @Override
    public void setTimeStep(double dt) {
        this.dt = dt;
        fdt = friction * dt;
        efdt = exp(-fdt);
        if (friction >= 0) {
            inverseFriction = 1.0 / friction;
        } else {
            inverseFriction = Double.POSITIVE_INFINITY;
        }
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    /**
     * Initialize the Random number generator used to apply random forces to the
     * particles.
     *
     * @param seed Random number generator seed.
     */
    public void setRandomSeed(long seed) {
        random.setSeed(seed);
    }

    /**
     * Set the frictional and random coefficients, store the current atom
     * positions, then find new atom positions and half-step velocities via
     * Verlet recursion.
     */
    @Override
    public void preForce(Potential potential) {
        for (int i = 0; i < nVariables; i++) {
            double m = mass[i];
            double pfric;
            double afric;
            double prand;
            if (fdt <= 0.0) {
                /**
                 * In the limit of no friction, SD recovers normal molecular
                 * dynamics.
                 */
                pfric = 1.0;
                vfric[i] = dt;
                afric = 0.5 * dt * dt;
                prand = 0.0;
                vrand[i] = 0.0;
            } else {
                double pterm;
                double vterm;
                double rho;
                if (fdt >= 0.05) {
                    /**
                     * Analytical expressions when the friction coefficient is
                     * large.
                     */
                    pfric = efdt;
                    vfric[i] = (1.0 - efdt) * inverseFriction;
                    afric = (dt - vfric[i]) * inverseFriction;
                    pterm = 2.0 * fdt - 3.0 + (4.0 - efdt) * efdt;
                    vterm = 1.0 - efdt * efdt;
                    rho = (1.0 - efdt) * (1.0 - efdt) / sqrt(pterm * vterm);
                } else {
                    /**
                     * Use a series expansions when friction coefficient is
                     * small.
                     */
                    double fdt2 = fdt * fdt;
                    double fdt3 = fdt * fdt2;
                    double fdt4 = fdt2 * fdt2;
                    double fdt5 = fdt2 * fdt3;
                    double fdt6 = fdt3 * fdt3;
                    double fdt7 = fdt3 * fdt4;
                    double fdt8 = fdt4 * fdt4;
                    double fdt9 = fdt4 * fdt5;
                    afric = (fdt2 / 2.0 - fdt3 / 6.0 + fdt4 / 24.0
                            - fdt5 / 120.0 + fdt6 / 720.0
                            - fdt7 / 5040.0 + fdt8 / 40320.0
                            - fdt9 / 362880.0) / (friction * friction);
                    vfric[i] = dt - friction * afric;
                    pfric = 1.0 - friction * vfric[i];
                    pterm = 2.0 * fdt3 / 3.0 - fdt4 / 2.0
                            + 7.0 * fdt5 / 30.0 - fdt6 / 12.0
                            + 31.0 * fdt7 / 1260.0 - fdt8 / 160.0
                            + 127.0 * fdt9 / 90720.0;
                    vterm = 2.0 * fdt - 2.0 * fdt2 + 4.0 * fdt3 / 3.0
                            - 2.0 * fdt4 / 3.0 + 4.0 * fdt5 / 15.0
                            - 4.0 * fdt6 / 45.0 + 8.0 * fdt7 / 315.0
                            - 2.0 * fdt8 / 315.0 + 4.0 * fdt9 / 2835.0;
                    rho = sqrt(3.0) * (0.5 - fdt / 16.0
                            - 17.0 * fdt2 / 1280.0
                            + 17.0 * fdt3 / 6144.0
                            + 40967.0 * fdt4 / 34406400.0
                            - 57203.0 * fdt5 / 275251200.0
                            - 1429487.0 * fdt6 / 13212057600.0
                            + 1877509.0 * fdt7 / 105696460800.0);
                }
                /**
                 * Compute random terms to thermostat the nonzero friction case.
                 */
                double ktm = Thermostat.kB * temperature / m;
                double psig = sqrt(ktm * pterm) / friction;
                double vsig = sqrt(ktm * vterm);
                double rhoc = sqrt(1.0 - rho * rho);
                double pnorm = random.nextGaussian();
                double vnorm = random.nextGaussian();
                prand = psig * pnorm;
                vrand[i] = vsig * (rho * pnorm + rhoc * vnorm);
            }

            /**
             * Store the current atom positions, then find new atom positions
             * and half-step velocities via Verlet recursion.
             */
            x[i] += (v[i] * vfric[i] + a[i] * afric + prand);
            v[i] = v[i] * pfric + 0.5 * a[i] * vfric[i];
        }
    }

    /**
     * Use Newton's second law to get the next acceleration and find the
     * full-step velocities using the Verlet recursion.
     */
    @Override
    public void postForce(double[] gradient) {
        for (int i = 0; i < nVariables; i++) {
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] += (0.5 * a[i] * vfric[i] + vrand[i]);
        }
    }

    /**
     * Update the integrator to be consistent with chemical perturbations. This
     * overrides the default implementation so that the vfric and vrand arrays
     * can be resized.
     *
     * @param nVariables the number of variables being integrated.
     * @param x the current value of each variable.
     * @param v the current velocity of each variable.
     * @param a the current acceleration of each variable.
     * @param aPrevious the previous acceleration of each variable.
     * @param mass the mass for each variable.
     */
    @Override
    public void setNumberOfVariables(int nVariables, double x[], double v[],
            double a[], double aPrevious[], double mass[]) {
        super.setNumberOfVariables(nVariables, x, v, a, aPrevious, mass);
        if (nVariables > vfric.length) {
            vfric = new double[nVariables];
            vrand = new double[nVariables];
        }
    }

}
