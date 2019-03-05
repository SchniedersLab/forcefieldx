/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms.integrators;

import java.util.logging.Logger;
import static java.lang.String.format;

import ffx.algorithms.thermostats.Thermostat;
import ffx.numerics.Potential;

/**
 * Respa performs multiple time step molecular dynamics using the reversible
 * reference system propagation algorithm (r-RESPA) via a Verlet core with the
 * potential split into fast- and slow-evolving portions.
 * <p>
 * The inner RESPA loop is position Verlet.
 * <p>
 * D. D. Humphreys, R. A. Friesner and B. J. Berne, "A Multiple-Time-Step
 * Molecular Dynamics Algorithm for Macromolecules", Journal of Physical
 * Chemistry, 98, 6885-6892 (1994)
 * <p>
 * X. Qian and T. Schlick, "Efficient Multiple-Time-Step Integrators with
 * Distance-Based Force Splitting for Particle-Mesh-Ewald Molecular Dynamics
 * Simulations", Journal of Chemical Physics, 115, 4019-4029 (2001)
 *
 * @author Gaurav Chattree
 * @since 1.0
 */
public class Respa extends Integrator {

    private static final Logger logger = Logger.getLogger(Respa.class.getName());

    /**
     * Number of inner time steps.
     */
    private int innerSteps;

    /**
     * Inner time step in psec.
     */
    private double innerTimeStep;

    /**
     * Half the inner time step.
     */
    private double halfInnerTimeStep;

    private double halfStepEnergy = 0;

    /**
     * Initialize Respa multiple time step molecular dynamics.
     *
     * @param nVariables Number of variables.
     * @param x          Variables current value.
     * @param v          Current velocities.
     * @param a          Current accelerations.
     * @param aPrevious  Previous accelerations.
     * @param aPrevious  Previous accelerations.
     * @param mass       Mass of the variables.
     */
    public Respa(int nVariables, double[] x, double[] v, double[] a,
                 double[] aPrevious, double[] mass) {
        super(nVariables, x, v, a, aPrevious, mass);

        innerSteps = 4;
        innerTimeStep = dt / innerSteps;
        halfInnerTimeStep = 0.5 * innerTimeStep;
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
     * {@inheritDoc}
     * <p>
     * Performs the inner RESPA loop via position Verlet.
     */
    @Override
    public void preForce(Potential potential) {
        double[] gradient = new double[nVariables];

        /**
         * Find half-step velocities via velocity Verlet recursion
         */
        for (int i = 0; i < nVariables; i++) {
            v[i] += a[i] * dt_2;
        }

        /**
         * Initialize accelerations due to fast-evolving forces.
         */
        potential.setEnergyTermState(Potential.STATE.FAST);
        halfStepEnergy = potential.energyAndGradient(x, gradient);
        for (int i = 0; i < nVariables; i++) {
            aPrevious[i] = -Thermostat.convert * gradient[i] / mass[i];
        }

        /**
         * Complete the inner RESPA loop.
         */
        for (int j = 0; j < innerSteps; j++) {

            /**
             * Find fast-evolving velocities and positions via Verlet recursion.
             */
            for (int i = 0; i < nVariables; i++) {
                v[i] += aPrevious[i] * halfInnerTimeStep;
                x[i] += v[i] * innerTimeStep;
            }

            /**
             * Update accelerations from fast varying forces.
             */
            halfStepEnergy = potential.energyAndGradient(x, gradient);
            for (int i = 0; i < nVariables; i++) {

                /**
                 * Use Newton's second law to get fast-evolving accelerations.
                 * Update fast-evolving velocities using the Verlet recursion.
                 */
                aPrevious[i] = -Thermostat.convert * gradient[i] / mass[i];
                v[i] += aPrevious[i] * halfInnerTimeStep;
            }
        }

        /**
         * Revert to computing slowly varying forces.
         */
        potential.setEnergyTermState(Potential.STATE.SLOW);
    }

    /**
     * {@inheritDoc}
     * <p>
     * The Respa full-step integration operation.
     */
    @Override
    public void postForce(double[] gradient) {
        for (int i = 0; i < nVariables; i++) {
            a[i] = -Thermostat.convert * gradient[i] / mass[i];
            v[i] += a[i] * dt_2;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Set outer Respa time step.
     */
    @Override
    public void setTimeStep(double dt) {
        if (dt < 0.0005) {
            dt = 0.0005;
        }

        this.dt = dt;
        dt_2 = 0.5 * dt;
        innerTimeStep = dt / innerSteps;
        halfInnerTimeStep = 0.5 * innerTimeStep;

        logger.info(format(" Time step set at %f (psec) and inner time step set at %f (psec) \n", this.dt, innerTimeStep));
    }

    /**
     * Set inner Respa number of time steps.
     *
     * @param n Number of inner time steps (must be greater than or equal to 2).
     */
    public void setInnerTimeSteps(int n) {
        if (n < 2) {
            n = 2;
        }

        innerSteps = n;

        // Update inner time step
        setTimeStep(dt);
    }
}
