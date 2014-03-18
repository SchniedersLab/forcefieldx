/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 * Integrate Newton's equations of motion using a Beeman multistep recursion
 * formula; the actual coefficients are Brooks' "Better Beeman" values.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class BetterBeeman extends Integrator {

    private double x[];
    private double v[];
    private double a[];
    private double aPrevious[];
    private double mass[];
    private int nVariables;
    private double dt;
    private double dt2_8;
    private double dt_8;

    /**
     * Constructor for BetterBeeman.
     *
     * @param nVariables number of Variables.
     * @param x Cartesian coordinates (Angstroms).
     * @param v Velocities.
     * @param a Accelerations.
     * @param aPrevious Previous Accelerations.
     * @param mass Mass.
     */
    public BetterBeeman(int nVariables, double x[], double v[], double a[],
            double aPrevious[], double mass[]) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.aPrevious = aPrevious;
        this.mass = mass;

        dt = 1.0;
        dt_8 = 0.125 * dt;
        dt2_8 = dt * dt_8;
    }

    /**
     * Store the current atom positions, then find new atom positions and
     * half-step velocities via Beeman recursion.
     */
    @Override
    public void halfStep(Potential potential) {
        for (int i = 0; i < nVariables; i++) {
            double temp = 5.0 * a[i] - aPrevious[i];
            x[i] += v[i] * dt + temp * dt2_8;
            v[i] += temp * dt_8;
        }
    }

    /**
     * Use Newton's second law to get the next acceleration and find the
     * full-step velocities using the Beeman recusion.
     */
    @Override
    public void fullStep(double gradient[]) {
        for (int i = 0; i < nVariables; i++) {
            aPrevious[i] = a[i];
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] += (3.0 * a[i] + aPrevious[i]) * dt_8;
        }
    }

    @Override
    public void setTimeStep(double dt) {
        this.dt = dt;
        dt_8 = 0.125 * dt;
        dt2_8 = dt * dt_8;
    }
}
