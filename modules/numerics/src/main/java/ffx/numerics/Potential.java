/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.numerics;

/**
 * The Potential interface defines methods required by an optimizer or molecular
 * dynamics.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface Potential {

    /**
     * Recognized variables currently include Cartesian coordinates and OTHER.
     */
    enum VARIABLE_TYPE {

        X, Y, Z, OTHER
    }

    /**
     * This method is called repeatedly to compute the function energy.
     *
     * @param x Input parameters.
     * @return Function value at <code>x</code>.
     * @since 1.0
     */
    double energy(double x[]);

    /**
     * This method is called repeatedly to compute the function energy and
     * gradient.
     *
     * @param x Input parameters.
     * @param g Output gradients with respect to each parameter.
     * @return Function value at <code>x</code>.
     * @since 1.0
     */
    double energyAndGradient(double x[], double g[]);

    /**
     * This method is called repeatedly to compute the function energy. The
     * verbose flag may not be used by all implementations.
     *
     * @param x Input parameters.
     * @param verbose Display extra information.
     * @return Function value at <code>x</code>
     */
    default double energy(double[] x, boolean verbose) {
        return energy(x);
    }

    /**
     * This method is called repeatedly to compute the function energy and
     * gradient. The verbose flag may not be used by all implementations.
     *
     * @param x Input parameters.
     * @param g Output gradients with respect to each parameter.
     * @param verbose Display extra information.
     * @return Function value at <code>x</code>.
     * @since 1.0
     */
    default double energyAndGradient(double[] x, double[] g, boolean verbose) {
        return energyAndGradient(x, g);
    }

    /**
     * Scale the problem. A good choice for optimization is the square root of
     * the median eigenvalue of a typical Hessian.
     *
     * @param scaling The scaling value to use for each variable.
     * @since 1.0
     */
    void setScaling(double scaling[]);

    /**
     * Get the problem scaling.
     *
     * @return The scaling value used for each variable.
     * @since 1.0
     */
    double[] getScaling();

    /**
     * Load the current value of the parameters. If the supplied array is null
     * or not large enough, a new one should be created. The filled array is
     * returned.
     *
     * @param parameters Supplied array.
     * @return The array filled with parameter values.
     */
    double[] getCoordinates(double[] parameters);

    /**
     * Get the mass of each degree of freedom. This is required for molecular
     * dynamics.
     *
     * @return The mass of each degree of freedom.
     */
    double[] getMass();

    /**
     * Get the total energy of the system
     *
     * @return the total energy
     */
    double getTotalEnergy();

    /**
     * Get the number of variables being operated on.
     *
     * @return Number of variables.
     */
    int getNumberOfVariables();

    /**
     * Get the type of all variables.
     *
     * @return The VARIABLE_TYPE of each variable.
     */
    VARIABLE_TYPE[] getVariableTypes();

    /**
     * <p>setVelocity.</p>
     *
     * @param velocity an array of {@link double} objects.
     */
    void setVelocity(double velocity[]);

    /**
     * <p>setAcceleration.</p>
     *
     * @param acceleration an array of {@link double} objects.
     */
    void setAcceleration(double acceleration[]);

    /**
     * <p>setPreviousAcceleration.</p>
     *
     * @param previousAcceleration an array of {@link double} objects.
     */
    void setPreviousAcceleration(double previousAcceleration[]);

    /**
     * <p>getVelocity.</p>
     *
     * @param velocity an array of {@link double} objects.
     * @return an array of {@link double} objects.
     */
    double[] getVelocity(double velocity[]);

    /**
     * <p>getAcceleration.</p>
     *
     * @param acceleration an array of {@link double} objects.
     * @return an array of {@link double} objects.
     */
    double[] getAcceleration(double acceleration[]);

    /**
     * <p>getPreviousAcceleration.</p>
     *
     * @param previousAcceleration an array of {@link double} objects.
     * @return an array of {@link double} objects.
     */
    double[] getPreviousAcceleration(double previousAcceleration[]);

    /**
     * Set the state of the Potential to include FAST varying energy terms, SLOW
     * varying energy terms or BOTH.
     */
    enum STATE {

        FAST, SLOW, BOTH
    }

    /**
     * Set the Potential Energy terms that should be active.
     *
     * @param state include FAST varying energy terms, SLOW varying energy terms
     * or BOTH.
     */
    void setEnergyTermState(STATE state);

    /**
     * Get the Potential Energy terms that is active.
     *
     * @return the STATE
     */
    STATE getEnergyTermState();

    /**
     * Destroys this Potential and frees up any associated resources, particularly worker Threads.
     * Default implementation is to return true (assume destruction successful).
     *
     * @return If resource reclamation successful, or resources already reclaimed.
     */
    default boolean destroy() {
        return true;
    }

}
