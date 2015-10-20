/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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

    private double dalt;
    private double dta;
    private double dta_2;
    private int nalt;
    private final double eps = .00000001;
    private double halfStepEnergy = 0;

    /**
     * Initialize Respa multiple time step molecular dynamics.
     *
     * @param nVariables Number of variables.
     * @param x Variables current value.
     * @param v Current velocities.
     * @param a Current accelerations.
     * @param aPrevious Previous accelerations.
     * @param mass Mass of the variables.
     */
    public Respa(int nVariables, double x[], double v[], double a[],
            double aPrevious[], double mass[]) {
        super(nVariables, x, v, a, aPrevious, mass);

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
     * @param potential the Potential for RESPA.
     */
    @Override
    public void preForce(Potential potential) {
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
                aPrevious[i] = -Thermostat.convert * gradient[i] / mass[i];
                v[i] += aPrevious[i] * dta;
                x[i] += v[i] * dta_2;
            }
        }
        potential.setEnergyTermState(Potential.STATE.SLOW);
    }

    /**
     * The Respa full-step integration operation.
     *
     * @param gradient the Gradient over all parameters.
     */
    @Override
    public void postForce(double[] gradient) {
        for (int i = 0; i < nVariables; i++) {
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] += a[i] * dt_2;
        }
    }

    /**
     * Set outer Respa time step. The inner time step is fixed at 0.25 fsec.
     *
     * @param dt Outer time step (dt must be .GE. 0.5 fsec).
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
