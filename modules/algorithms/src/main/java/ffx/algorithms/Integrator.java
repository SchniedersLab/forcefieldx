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
 * The Integrator class is responsible for propagation of degrees of freedom
 * through time. Implementations must define their behavior at half-step and
 * full-step time points.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public abstract class Integrator {

    public enum Integrators {

        BEEMAN, RESPA, STOCHASTIC
    };

    abstract public void setTimeStep(double dt);

    /**
     * Integrator halfStep operation.
     */
    abstract public void preForce(Potential potential);

    /**
     * Integrator fullStep operation.
     */
    abstract public void postForce(double gradient[]);
}
