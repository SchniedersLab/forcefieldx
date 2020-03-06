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

import org.apache.commons.math3.util.FastMath;

/**
 * Implements a bell-shaped switching function by stitching together
 * a pair of MultiplicativeSwitches. f(midpoint - 0.5*width) = 0,
 * f(midpoint) = 1, f(midpoint + 0.5*width) = 0.
 *
 * @author Jacob M. Litman
 * @author Rae Corrigan
 */
public class BellCurveSwitch implements UnivariateSwitchingFunction {
    private final double midpoint;
    private final double halfWidth;
    private final double invWidth;
    private final UnivariateSwitchingFunction switchingFunction;
    private final UnivariateSwitchingFunction secondSwitchingFunction;

    /**
     * Construct a bell curve (spliced 5-'th order Hermite splines) of width 1.0, midpoint 0.5.
     */
    public BellCurveSwitch() {
        this(0.5);
    }

    /**
     * Construct a bell curve (spliced 5-'th order Hermite splines) of width 1.0.
     *
     * @param midpoint Midpoint of the curve.
     */
    public BellCurveSwitch(double midpoint) {
        this(midpoint, 1.0);
    }

    /**
     * Construct a bell curve (spliced 5-'th order Hermite splines).
     * @param midpoint Midpoint of the curve.
     * @param width Width of the curve, between the two zero points.
     */
    public BellCurveSwitch(double midpoint, double width) {
        this.midpoint = midpoint;
        invWidth = 1.0 / width;

        halfWidth = 0.5 * width;
        switchingFunction = new MultiplicativeSwitch(midpoint, midpoint - halfWidth);
        secondSwitchingFunction = new MultiplicativeSwitch(midpoint, midpoint + halfWidth);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getZeroBound() {
        return midpoint - halfWidth;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getOneBound() {
        return midpoint + halfWidth;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean constantOutsideBounds() {
        return switchingFunction.constantOutsideBounds() && secondSwitchingFunction.constantOutsideBounds();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean validOutsideBounds() {
        return switchingFunction.validOutsideBounds() && secondSwitchingFunction.validOutsideBounds();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getHighestOrderZeroDerivative() {
        return Math.min(switchingFunction.getHighestOrderZeroDerivative(), secondSwitchingFunction.getHighestOrderZeroDerivative());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean symmetricToUnity() {
        return switchingFunction.symmetricToUnity() && secondSwitchingFunction.symmetricToUnity();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        if (x > midpoint) {
            return secondSwitchingFunction.valueAt(x);
        } else {
            return switchingFunction.valueAt(x);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double firstDerivative(double x) {
        if (x > midpoint) {
            return invWidth * secondSwitchingFunction.firstDerivative(x);
        } else {
            return invWidth * switchingFunction.firstDerivative(x);
        }
    }

    /**
     * {@inheritDoc}return max(cut, off);
     */
    @Override
    public double secondDerivative(double x) {
        if (x > midpoint) {
            return invWidth * invWidth * secondSwitchingFunction.secondDerivative(x);
        } else {
            return invWidth * invWidth * switchingFunction.secondDerivative(x);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        double mult = FastMath.pow(invWidth, order);
        if (x > midpoint) {
            return mult * secondSwitchingFunction.nthDerivative(x, order);
        } else {
            return mult * switchingFunction.nthDerivative(x, order);
        }
    }

    @Override
    public String toString() {
        return String.format(" Spliced 5'th order Hermite splines with midpoint " +
                "%11.5g, width %11.5g", midpoint, 2.0 * halfWidth);
    }
}
