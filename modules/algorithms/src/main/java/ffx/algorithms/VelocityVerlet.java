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
 * Integrate Newton's equations of motion using a Velocity Verlet multistep
 * recursion formula.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class VelocityVerlet extends Integrator {

    private double x[];
    private double v[];
    private double a[];
    private double mass[];
    private int nVariables;
    private double dt;
    private double dt_2;

    /**
     * Constructor for VelocityVerlet.
     *
     * @param nVariables number of Variables.
     * @param x Cartesian coordinates (Angstroms).
     * @param v Velocities.
     * @param a Accelerations.
     * @param mass Mass.
     */
    public VelocityVerlet(int nVariables, double x[], double v[], double a[],
            double mass[]) {
        this.nVariables = nVariables;
        this.x = x;
        this.v = v;
        this.a = a;
        this.mass = mass;
        dt = 1.0;
        dt_2 = dt * .5;
    }

    /**
     * Find pre-gradient velocities using and pre-gradient positions
     */
    @Override
    public void preForce(Potential potential) {
        for (int i = 0; i < nVariables; i++) {
            v[i] = v[i] + a[i] * dt_2;
        }
        for (int i = 0; i < nVariables; i++) {
            x[i] = x[i] + v[i] * dt;
        }
    }

    /**
     * Use Newton's second law to find acceleration from full-step positions and
     * find the full-step velocities using the new acceleration.
     */
    @Override
    public void postForce(double gradient[]) {

        for (int i = 0; i < nVariables; i++) {
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] = v[i] + a[i] * dt_2;
        }
    }

    @Override
    public void setTimeStep(double dt) {
        this.dt = dt;
        dt_2 = dt * .5;
    }
}
