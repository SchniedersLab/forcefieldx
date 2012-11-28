/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */
package ffx.algorithms;

import ffx.numerics.Potential;

/**
 * Respa performs multiple time step molecular dynamics using the reversible
 * reference system propagation algorithm (r-RESPA) via a Verlet core with the
 * potential split into fast- and slow-evolving portions.
 *
 * The inner RESPA loop is position Verlet.
 *
 * @author Gaurav Chattree
 *
 * @since 1.0
 */
public class Respa extends Integrator {

    private double x[];
    private double v[];
    private double a[];
    private double aAlt[];
    private double mass[];
    private int nVariables;
    private double dt;
    private double dt_2;
    private double dalt;
    private double dta;
    private double dta_2;
    private int nalt;
    private double eps = .00000001;
    private double halfStepEnergy = 0;

    /**
     * Initialize Respa multiple time step molecular dynamics.
     *
     * @param nVariables
     * @param x
     * @param v
     * @param a
     * @param aPrevious
     * @param mass
     */
    public Respa(int nVariables, double x[], double v[], double a[],
            double aPrevious[], double mass[]) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.aAlt = aPrevious;
        this.mass = mass;
        dt = 1.0;
        dt_2 = 0.5 * dt;
        dalt = 0.00025;
        nalt = ((int) (dt / (dalt + eps))) + 1;
        dalt = (double) nalt;
        dta = dt / dalt;
        dta_2 = 0.5 * dta;
    }

    /**
     * Get the potential energy of the fast degrees of freedom.
     *
     * @return The potential energy of the fast degrees of freedom.
     */
    public double getHalfStepEnergy() {
        return halfStepEnergy;
    }

    /**
     * Performs the inner RESPA loop via position Verlet.
     *
     * @param potential
     */
    @Override
    public void halfStep(Potential potential) {
        double gradient[] = new double[nVariables];
        for (int i = 0; i < nVariables; i++) {
            v[i] += a[i] * dt_2;
        }

        /**
         * The inner RESPA loop.
         */
        for (int j = 0; j < nalt; j++) {
            for (int i = 0; i < nVariables; i++) {
                x[i] += v[i] * dta_2;
            }
            potential.setEnergyTermState(Potential.STATE.FAST);
            halfStepEnergy = potential.energyAndGradient(x, gradient);
            for (int i = 0; i < nVariables; i++) {
                aAlt[i] = -Thermostat.convert * gradient[i] / mass[i];
                v[i] += aAlt[i] * dta;
                x[i] += v[i] * dta_2;
            }
        }
        potential.setEnergyTermState(Potential.STATE.SLOW);
    }

    /**
     * The Respa full-step integration operation.
     *
     * @param gradient
     */
    @Override
    public void fullStep(double[] gradient) {
        for (int i = 0; i < nVariables; i++) {
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] += a[i] * dt_2;
        }
    }

    /**
     * Set outer Respa time step. The inner time step is fixed at 0.25 fsec.
     *
     * @param dt Outer time step (dt must be >= 0.5 fsec).
     */
    @Override
    public void setTimeStep(double dt) {
        if (dt < 0.0005) {
            dt = 0.0005;
        }

        this.dt = dt;
        dt_2 = 0.5 * dt;
        dalt = 0.00025;
        nalt = ((int) (dt / (dalt + eps))) + 1;
        dalt = (double) nalt;
        dta = dt / dalt;
        dta_2 = 0.5 * dta;
    }
}
