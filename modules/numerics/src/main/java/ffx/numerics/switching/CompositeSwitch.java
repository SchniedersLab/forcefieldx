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

import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.max;

/**
 * The CompositeSwitch uses a primary switch in the middle, and then two secondary
 * switches at the ends of the path to smoothly switch to the primary switch. For
 * example, one can smoothly interpolate from 0/0/0 value/derivative/second derivative
 * to a linear switch by multiplying the linear switch by a MultiplicativeSwitch in the
 * range 0-0.1.
 * <p>
 * At present, there is an assumption that x gets linearly scaled when passed to the switch;
 * at the ends, s(x) = f(g(x))*h(x), where h(x) is the primary switch, g(x) = x / (switching range),
 * and f(g(x)) is the secondary switch.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class CompositeSwitch implements UnivariateSwitchingFunction {
    private static final Logger logger = Logger.getLogger(CompositeSwitch.class.getName());

    /**
     * Primary switching function to be obeyed exactly in the middle.
     */
    private final UnivariateSwitchingFunction primaryFunction;
    /**
     * Secondary switching function used to interpolate smoothly
     * (usually from 0/0/0) into primaryFunction over lb to lbPrmary.
     */
    private final UnivariateSwitchingFunction startSwitch;
    /**
     * Secondary switching function used to interpolate smoothly
     * from primaryFunction (usually to 1/0/0) over lb to lbPrmary.
     */
    private final UnivariateSwitchingFunction endSwitch;
    /**
     * Lower bound of the primary switch/upper bound of startSwitch * primarySwitch.
     */
    private final double lbPrimary;
    /**
     * Upper bound of the primary switch/lower bound of endSwitch * primarySwitch.
     */
    private final double ubPrimary;
    /**
     * Overall lower bound of the CompositeSwitch.
     */
    private final double lb;
    /**
     * Overall upper bound of the CompositeSwitch.
     */
    private final double ub;

    // Below six constants are only valid if g(x) = x / range
    // Precomputed constants used in evaluating values/derivatives from lb to lbPrimary.
    private final double multLB;
    private final double fdLB;
    private final double fdLB2;
    // Precomputed constants used in evaluating values/derivatives from ubPrimary to ub.
    private final double multUB;
    private final double fdUB;
    private final double fdUB2;

    /**
     * Builds a switch that uses MultiplicativeSwitches at the ends (0-0.1, 0.9-1.0) to smoothly interpolate a linear
     * switch between 0 and 1 with smooth 2'nd and 3'rd derivatives. Will not be quite linear from 0-0.1 and 0.9-1.0.
     */
    public CompositeSwitch() {
        this(new PowerSwitch());
    }

    /**
     * Builds a switch that uses MultiplicativeSwitches at the ends (0-0.1, 0.9-1.0) to smoothly interpolate a provided
     * switch between 0 and 1 with smooth 2'nd and 3'rd derivatives.
     *
     * @param primary Primary switch to obey exactly from 0.1-0.9.
     */
    public CompositeSwitch(UnivariateSwitchingFunction primary) {
        this(primary, new MultiplicativeSwitch(), new MultiplicativeSwitch(), 0.1, 0.9);
    }

    /**
     * Builds a composite switch in .
     *
     * @param primary   Primary switch to obey exactly from lbPrimary to ubPrimary.
     * @param start     Switch to interpolate from 0 to primary between 0 and lbPrimary; assumed to internally function from 0-1.
     * @param end       Switch to interpolate from primary to 1.0 between ubPrimary and 1; assumed to internally function from 0-1.
     * @param lbPrimary Value at which primary should begin to be obeyed exactly.
     * @param ubPrimary Value at which primary should stop being obeyed exactly.
     */
    public CompositeSwitch(UnivariateSwitchingFunction primary, UnivariateSwitchingFunction start, UnivariateSwitchingFunction end, double lbPrimary, double ubPrimary) {
        this(primary, start, end, lbPrimary, ubPrimary, 0, 1);
    }

    /**
     * Builds a composite switch in .
     *
     * @param primary   Primary switch to obey exactly from lbPrimary to ubPrimary.
     * @param start     Switch to interpolate from 0 to primary between 0 and lbPrimary; assumed to internally function from 0-1.
     * @param end       Switch to interpolate from primary to 1.0 between ubPrimary and 1; assumed to internally function from 0-1.
     * @param lbPrimary Value at which primary should begin to be obeyed exactly.
     * @param ubPrimary Value at which primary should stop being obeyed exactly.
     * @param ub        Overall upper bound of the CompositeSwitch.
     * @oaram lb        Overall lower bound of the CompositeSwitch.
     */
    public CompositeSwitch(UnivariateSwitchingFunction primary, UnivariateSwitchingFunction start, UnivariateSwitchingFunction end, double lbPrimary, double ubPrimary, double lb, double ub) {
        if (lbPrimary > ubPrimary) {
            throw new IllegalArgumentException(String.format(" Lower primary bound %10.4g was greater than upper primary bound %10.4g", lbPrimary, ubPrimary));
        }
        if (lb > ub) {
            throw new IllegalArgumentException(String.format(" Lower bound %10.4g was greater than upper bound %10.4g", lb, ub));
        }
        assert lb < lbPrimary && ub > ubPrimary;

        primaryFunction = primary;
        startSwitch = start;
        endSwitch = end;
        this.lbPrimary = lbPrimary;
        this.ubPrimary = ubPrimary;
        this.lb = lb;
        this.ub = ub;

        fdLB = lbPrimary - lb;
        multLB = 1.0 / fdLB;
        fdLB2 = fdLB * fdLB;

        fdUB = ub - ubPrimary;
        multUB = 1.0 / fdUB;
        fdUB2 = fdUB * fdUB;

        if (!testJoints()) {
            logger.warning(String.format(" Switch %s is not smooth at one of its joints!", toString()));
        }
    }

    @Override
    public boolean constantOutsideBounds() {
        return startSwitch.constantOutsideBounds() && endSwitch.constantOutsideBounds();
    }

    @Override
    public int getHighestOrderZeroDerivative() {
        return Math.max(Math.max(startSwitch.getHighestOrderZeroDerivative(), endSwitch.getHighestOrderZeroDerivative()), primaryFunction.getHighestOrderZeroDerivative());
    }

    @Override
    public double getOneBound() {
        return ub;
    }

    @Override
    public double getZeroBound() {
        return lb;
    }

    @Override
    public boolean symmetricToUnity() {
        return primaryFunction.symmetricToUnity() && startSwitch.equals(endSwitch) && (lbPrimary - lb == ub - ubPrimary);
    }

    @Override
    public boolean validOutsideBounds() {
        return startSwitch.constantOutsideBounds() && endSwitch.constantOutsideBounds();
    }

    @Override
    public double firstDerivative(double x) throws IllegalArgumentException {
        if (x < lbPrimary) {
            return fdLower(x);
        } else if (x > ubPrimary) {
            return fdUpper(x);
        } else {
            return primaryFunction.firstDerivative(x);
        }
    }

    @Override
    public double nthDerivative(double x, int order) throws IllegalArgumentException {
        switch (order) {
            case 0:
                return valueAt(x);
            case 1:
                return firstDerivative(x);
            case 2:
                return secondDerivative(x);
            default:
                throw new IllegalArgumentException(" Composite switches do not yet have support for arbitrary derivatives");
        }
    }

    @Override
    public double secondDerivative(double x) throws IllegalArgumentException {
        if (x < lbPrimary) {
            return sdLower(x);
        } else if (x > ubPrimary) {
            return sdUpper(x);
        } else {
            return primaryFunction.secondDerivative(x);
        }
    }

    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        if (x < lbPrimary) {
            return valLower(x);
        } else if (x > ubPrimary) {
            return valUpper(x);
        } else {
            return primaryFunction.valueAt(x);
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(String.format(" Composite switch with overall range %12.5g-%12.5g, " +
                "with an inner range %12.5g-%12.5g", lb, ub, lbPrimary, ubPrimary));
        sb.append("\n Primary switch: ").append(primaryFunction.toString());
        sb.append("\n Start switch:   ").append(startSwitch.toString());
        sb.append("\n End switch:     ").append(endSwitch.toString());
        return sb.toString();
    }

    private boolean approxEquals(double x1, double x2) {
        return approxEquals(x1, x2, 1E-11);
    }

    /**
     * Tests double equality to within a reasonable tolerance.
     *
     * @param x1
     * @param x2
     * @param tol
     * @return Fuzzy equality.
     */
    private boolean approxEquals(double x1, double x2, double tol) {
        double largerVal = max(abs(x1), abs(x2));
        if (largerVal == 0) {
            // If both are zero, they're both equal.
            return true;
        } else if (largerVal > 1E-6) {
            // Use normalized approximate-equals.
            return abs((x1 - x2) / largerVal) < tol;
        } else {
            // Use absolute approximate-equals to avoid numerical issues.
            return abs(x1 - x2) < tol;
        }
    }

    /**
     * Test that values, first derivatives, and second derivatives are smooth at lbPrimary and ubPrimary,
     * such that the primary switch matches the composite switch at these values.
     *
     * @return If joints are smooth.
     */
    private boolean testJoints() {
        if (!approxEquals(valLower(lbPrimary), primaryFunction.valueAt(lbPrimary))) {
            return false;
        }
        if (!approxEquals(valUpper(ubPrimary), primaryFunction.valueAt(ubPrimary))) {
            return false;
        }

        if (!approxEquals(fdLower(lbPrimary), primaryFunction.firstDerivative(lbPrimary))) {
            return false;
        }
        if (!approxEquals(fdUpper(ubPrimary), primaryFunction.firstDerivative(ubPrimary))) {
            return false;
        }

        if (!approxEquals(sdLower(lbPrimary), primaryFunction.secondDerivative(lbPrimary))) {
            return false;
        }
        return approxEquals(sdUpper(ubPrimary), primaryFunction.secondDerivative(ubPrimary));
    }

    private double lbX(double x) {
        return (x - lb) * multLB;
    }

    private double ubX(double x) {
        return (ub - x) * multUB;
    }

    /**
     * Broken out from valueAt to make testing joints easier.
     *
     * @param x x value to evalute.
     * @return f(x) using both startSwitch and primary function.
     */
    private double valLower(double x) {
        return startSwitch.valueAt(lbX(x)) * primaryFunction.valueAt(x);
    }

    /**
     * Broken out from valueAt to make testing joints easier.
     *
     * @param x x value to evalute.
     * @return f(x) using both endSwitch and primary function.
     */
    private double valUpper(double x) {
        return endSwitch.valueAt(ubX(x)) * primaryFunction.valueAt(x);
    }

    private double fdLower(double x) {
        double swX = lbX(x);
        double val = primaryFunction.firstDerivative(x) * startSwitch.valueAt(swX);
        val += primaryFunction.valueAt(x) * startSwitch.firstDerivative(swX) * fdLB;
        return val;
    }

    private double fdUpper(double x) {
        double swX = ubX(x);
        double val = primaryFunction.firstDerivative(x) * endSwitch.valueAt(swX);
        val += primaryFunction.valueAt(x) * endSwitch.firstDerivative(swX) * fdUB;
        return val;
    }

    private double sdLower(double x) {
        double swX = lbX(x);
        double val = primaryFunction.secondDerivative(x) * startSwitch.valueAt(swX);
        val += (2 * primaryFunction.firstDerivative(x) * startSwitch.firstDerivative(swX) * fdLB);
        val += (primaryFunction.valueAt(x) * startSwitch.secondDerivative(swX) * fdLB2);
        return val;
    }

    private double sdUpper(double x) {
        double swX = ubX(x);
        double val = primaryFunction.secondDerivative(x) * endSwitch.valueAt(swX);
        val += (2 * primaryFunction.firstDerivative(x) * endSwitch.firstDerivative(swX) * fdUB);
        val += (primaryFunction.valueAt(x) * endSwitch.secondDerivative(swX) * fdUB2);
        return val;
    }
}
