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
package ffx.numerics;

/**
 * The Potential interface defines methods required by an optimizer or molecular
 * dynamics.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public interface Potential {

    /**
     * Recognized variables currently include Cartesian coordinates and OTHER.
     */
    public enum VARIABLE_TYPE {

        X, Y, Z, OTHER
    };

    /**
     * This method is called repeatedly to compute the function energy.
     *
     * @param x Input parameters.
     * @return Function value at <code>x</code>.
     * @since 1.0
     */
    public abstract double energy(double x[]);

    /**
     * This method is called repeatedly to compute the function energy and
     * gradient.
     *
     * @param x Input parameters.
     * @param g Output gradients with respect to each parameter.
     * @return Function value at <code>x</code>.
     * @since 1.0
     */
    public abstract double energyAndGradient(double x[], double g[]);

    /**
     * Scale the problem. A good choice for optimization is the square root of
     * the median eigenvalue of a typical Hessian.
     *
     * @param scaling The scaling value to use for each variable.
     * @since 1.0
     */
    public abstract void setScaling(double scaling[]);

    /**
     * Get the problem scaling.
     *
     * @return The scaling value used for each variable.
     * @since 1.0
     */
    public abstract double[] getScaling();

    /**
     * Load the current value of the parameters. If the supplied array is null
     * or not large enough, a new one should be created. The filled array is
     * returned.
     *
     * @param parameters Supplied array.
     *
     * @return The array filled with parameter values.
     */
    public abstract double[] getCoordinates(double[] parameters);

    /**
     * Get the mass of each degree of freedom. This is required for molecular
     * dynamics.
     *
     * @return The mass of each degree of freedom.
     */
    public abstract double[] getMass();

    /**
     * Get the total energy of the system
     *
     * @return the total energy
     */
    public abstract double getTotalEnergy();

    /**
     * Get the number of variables being operated on.
     *
     * @return Number of variables.
     */
    public abstract int getNumberOfVariables();

    /**
     * Get the type of all variables.
     *
     * @return The VARIABLE_TYPE of each variable.
     */
    public abstract VARIABLE_TYPE[] getVariableTypes();

    /**
     * Set the state of the Potential to include FAST varying energy terms, SLOW
     * varying energy terms or BOTH.
     */
    public enum STATE {

        FAST, SLOW, BOTH
    };

    /**
     * Set the Potential Energy terms that should be active.
     *
     * @param state include FAST varying energy terms, SLOW
     * varying energy terms or BOTH.
     */
    public abstract void setEnergyTermState(STATE state);
    
    /**
     * Re-initializes the Potential and performs all necessary operations upon 
     * any change in the underlying variables.
     */
    public abstract void reInit();
}
