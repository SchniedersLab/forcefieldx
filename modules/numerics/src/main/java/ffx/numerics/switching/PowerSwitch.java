//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.numerics.switching;

import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.pow;

/**
 * A PowerSwitch interpolates between 0 and 1 vi f(x) = (ax)^beta, where x must
 * be between 0 and 1/a.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class PowerSwitch implements UnivariateSwitchingFunction {

    /**
     * The multiplier a.
     */
    private final double a;
    /**
     * The power of the switch.
     */
    private final double beta;
    /**
     * The upper bound of the switch.
     */
    private final double ub;

    /**
     * Default Constructor of the PowerSwitch: constructs a linear switch.
     */
    public PowerSwitch() {
        this(1.0, 1.0);
    }

    /**
     * Constructor of the PowerSwitch.
     *
     * @param a    The upper bound of the switch is 1.0 / a.
     * @param beta The power of the function f(x) = (ax)^beta,
     */
    public PowerSwitch(double a, double beta) {
        if (a <= 0) {
            throw new IllegalArgumentException(" The constant a must be > 0");
        }
        if (beta == 0) {
            throw new IllegalArgumentException(" The exponent must be > 0 (preferably >= 1)");
        }
        this.a = a;
        this.beta = beta;
        ub = 1.0 / a;
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
    public int getHighestOrderZeroDerivative() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getOneBound() {
        return ub;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getZeroBound() {
        return 0;
    }

    /**
     * Power switch derivatives can be zero at the zero bound if the exponent
     * is greater than the derivative order.
     *
     * @return the highest order zero derivative at zero bound
     */
    @Override
    public int highestOrderZeroDerivativeAtZeroBound() {
        return beta >= 1 ? ((int) beta) - 1 : 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean symmetricToUnity() {
        return (beta == 1.0);
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
    public double firstDerivative(double x) throws IllegalArgumentException {
        x *= a;
        return beta * a * pow(x, beta - 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        x *= a;
        if (order < 1) {
            throw new IllegalArgumentException("Order must be >= 1");
        }
        switch (order) {
            case 1:
                return firstDerivative(x);
            case 2:
                return secondDerivative(x);
            default:
                double orderDiff = order - beta;
                if (orderDiff % 1.0 == 0 && orderDiff >= 1.0) {
                    return 0.0;
                } else {
                    double val = pow(x, beta - order);
                    for (int i = 0; i < order; i++) {
                        val *= (beta - i) * a;
                    }
                    return val;
                }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double secondDerivative(double x) throws IllegalArgumentException {
        x *= a;
        return beta == 1.0 ? 0.0 : beta * (beta - 1) * a * a * pow(x, beta - 2);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        x *= a;
        return pow(x, beta);
    }

    /**
     * Gets the value of beta in f(x) = (a*x)^beta
     *
     * @return Exponent of input
     */
    public double getExponent() {
        return beta;
    }

    /**
     * Gets the value of a in f(x) = (a*x)^beta.
     *
     * @return Multiplier of input
     */
    public double getMultiplier() {
        return a;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return format("Power switching function f(x) = (%8.4g * x)^%8.4g", a, beta);
    }
}
