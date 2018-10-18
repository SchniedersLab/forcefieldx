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
package ffx.potential.bonded;

/**
 * The LambdaInterface should be implemented by potential energy terms that can
 * accept a lambda value from [0 .. 1] that defines a twice differentiable path
 * between states 0 and 1.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public interface LambdaInterface {

    /**
     * Set the current value of the state variable. May be ignored if lambda is not being applied.
     *
     * @param lambda a double.
     * @since 1.0
     */
    void setLambda(double lambda);

    /**
     * Get the current value of the state variable.
     *
     * @return state
     * @since 1.0
     */
    double getLambda();

    /**
     * Get the partial derivative of the energy with respect to lambda.
     *
     * @return dEdL
     * @since 1.0
     */
    double getdEdL();

    /**
     * Get the 2nd partial derivative of the energy with respect to lambda.
     *
     * @return d2EdL2
     * @since 1.0
     */
    double getd2EdL2();

    /**
     * Get the gradient of dEdL with respect to each parameter.
     *
     * @param gradient - A double array of length the number of parameters in
     *                 the model (commonly 3 * number of atoms).
     * @since 1.0
     */
    void getdEdXdL(double gradient[]);

    /**
     * Returns true if dUdL is guaranteed to be zero at 0 and 1. Default
     * implementation is to return false.
     *
     * @return True if dUdL is guaranteed 0 at endpoints.
     */
    default boolean dEdLZeroAtEnds() {
        return false;
    }
}
