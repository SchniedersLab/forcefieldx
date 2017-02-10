/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
 * A UnivariateSwitchingFunction describes a function of a single value (often
 * lambda), where f(lb) = 0, f(ub) = 1, and df(x)/dx >= 0 for all x lb-ub. 
 * Additionally, good switching functions have zero first and second derivatives
 * at the lower and upper bound.
 * 
 * A number of methods exist to check for various properties of a switching
 * function; these will often be implemented as simple return-boolean methods.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public interface UnivariateSwitchingFunction {
    
    /**
     * Gets the lower bound, where f(x) becomes zero.
     * @return Lower bound
     */
    public double getLowerBound();
    
    /**
     * Gets the upper bound, where f(x) becomes one.
     * @return Upper bound
     */
    public double getUpperBound();
    
    /**
     * Remains 0 below the lower bound, and 1 above the upper bound (i.e. 
     * a multiplicative switch).
     * @return df(x)/dx is zero outside range lb-ub.
     */
    public boolean constantOutsideBounds();
    
    /**
     * Remains in the range 0-1 outside the bounds. Implied to be true if 
     * constantOutsideBounds is true.
     * @return min(f(x)) = 0 and max(f(x)) = 1.
     */
    public boolean validOutsideBounds();
    
    /**
     * The highest-order derivative that is zero at the bounds.
     * @return Maximum zero derivative at bounds.
     */
    public int getHighestOrderZeroDerivative();
    
    /**
     * True if f(lb + delta) + f(ub - delta) = 1 for all delta between 0 and
     * (ub - lb). For example, a power switch with beta 1 is symmetric to unity,
     * as f(l) + f(1-l) = 1, but beta 2 produces a non-unity result, where 
     * f(0.5) + f(0.5) = 0.5.
     * 
     * @return If symmetry produces unity result.
     */
    public boolean symmetricToUnity();
    
    /**
     * Value at a point
     * @param x
     * @return f(x)
     * @throws IllegalArgumentException If f(x) is undefined at x.
     */
    public double valueAt(double x) throws IllegalArgumentException;
    
    /**
     * First derivative at a point.
     * @param x
     * @return f'(x)
     * @throws IllegalArgumentException If f'(x) is undefined at x.
     */
    public double firstDerivative(double x) throws IllegalArgumentException;
    
    /**
     * Second derivative at a point.
     * @param x
     * @return f''(x)
     * @throws IllegalArgumentException If f''(x) is undefined at x.
     */
    public double secondDerivative(double x) throws IllegalArgumentException;
    
    /**
     * N'th order derivative at a point. Should be relatively optional for any
     * order above 2.
     * @param x
     * @param order Derivative order (>= 1)
     * @return d^nf(x)/dx^n
     * @throws IllegalArgumentException If derivative undefined at x.
     */
    public double nthDerivative(double x, int order) throws IllegalArgumentException;
}
