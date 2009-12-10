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
package ffx.numerics;

/**
 * The OptimizationInterface defines methods required by optimizers.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Optimizable {

    /**
     * This method is called repeatedly by the optimizer to compute the
     * function energy and gradient.
     * 
     * @param x Input parameters.
     * @param g Output gradients with respect to each parameter.
     * 
     * @return Function value at <code>x</code>.
     * 
     * @since 1.0
     */
    public abstract double energyAndGradient(double x[], double g[]);

    /**
     * Scale the optimization problem. A good choice is the square root of the
     * median eigenvalue of a typical Hessian.
     * 
     * @param scaling The scaling value to use for each variable.
     * 
     * @since 1.0
     */
    public abstract void setOptimizationScaling(double scaling[]);

    /**
     * Get the scaling for the optimization problem.
     *
     * @return The scaling value used for each variable.
     * 
     * @since 1.0
     */
    public abstract double[] getOptimizationScaling();

}
