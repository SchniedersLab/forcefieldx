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
package ffx.numerics.switching;

/**
 * A LinearDerivativeSwitch interpolates between 0 and 1 vi f(x) = 2*x - x^2.
 * <p>
 * The derivative is then linear in x: f'(x) = 2 - 2*x
 * <p>
 * Limiting behavior is given by: f(0) = 0, f(1) = 1, f'(0) = 2, f'(1) = 0.
 *
 * @author Michael J. Schnieders
 */
public class LinearDerivativeSwitch implements UnivariateSwitchingFunction {

    /**
     * Constructor for the LinearDerivativeSwitch.
     */
    public LinearDerivativeSwitch() {
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getZeroBound() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getOneBound() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean constantOutsideBounds() {
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean validOutsideBounds() {
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getHighestOrderZeroDerivative() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean symmetricToUnity() {
        return false;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        return 2 * x - (x * x);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double firstDerivative(double x) throws IllegalArgumentException {
        return 2.0 - 2.0 * x;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double secondDerivative(double x) throws IllegalArgumentException {
        return -2.0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        if (order < 1) {
            throw new IllegalArgumentException("Order must be >= 1");
        }
        switch (order) {
            case 1:
                return firstDerivative(x);
            case 2:
                return secondDerivative(x);
            default:
                return 0;
        }
    }
}
