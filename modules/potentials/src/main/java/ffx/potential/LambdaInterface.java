/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
 * can accept a lambda value from [0 .. 1] that defines a smooth path between
 * state 0 to state 1.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public interface LambdaInterface {

    public void setLambda(double lambda);

    public double getLambda();
    
    public void lambdaGradient(boolean lambdaGradient);
    
    public double getdEdLambda();
    
    public double getd2EdLambda2();
    
    public void getdEdLambdaGradient(double gradient[]);
    
}
