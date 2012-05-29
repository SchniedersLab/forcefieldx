/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012
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

package ffx.potential;

/**
 * The LambdaInterface should be implemented by potential energy terms that
 * can accept a lambda value from [0 .. 1] that defines a twice differentiable
 * path between states 0 and 1.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public interface LambdaInterface {

    /**
     * Set the current value of the state variable.
     *
     * @param lambda a double.
     * @since 1.0
     */
    public void setLambda(double lambda);

    /**
     * Get the current value of the state variable.
     *
     * @return state
     * @since 1.0
     */
    public double getLambda();
    
    /**
     * Get the partial derivative of the energy with respect to lambda.
     *
     * @return dEdL
     * @since  1.0
     */
    public double getdEdL();
    
    /**
     * Get the 2nd partial derivative of the energy with respect to lambda.
     *
     * @return d2EdL2
     * @since  1.0
     */
    public double getd2EdL2();
    
    /**
     * Get the gradient of dEdL with respect to each parameter.
     *
     * @param gradient - A double array of length the number of parameters
     *        in the model (commonly 3 * number of atoms).
     * @since 1.0
     */
    public void getdEdXdL(double gradient[]);
    
}
