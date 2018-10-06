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
package ffx.potential.nonbonded;

import static org.apache.commons.math3.util.FastMath.pow;

import ffx.numerics.switching.UnivariateSwitchingFunction;

/**
 * The 6 coefficients of the multiplicative polynomial switch are unique given
 * the distances "off" and "cut". They are found by solving a system of 6
 * equations, which define the boundary conditions of the switch.
 * <br>
 * f(cut) = 1
 * <br>
 * f'(cut) = f"(cut) = 0
 * <br>
 * f(off) = f'(off) = f"(off) = 0
 *
 * @author Michael J. Schnieders
 */
public class MultiplicativeSwitch implements UnivariateSwitchingFunction {

    private final double off;
    private final double cut;
    private final double c0;
    private final double c1;
    private final double c2;
    private final double c3;
    private final double c4;
    private final double c5;
    private final double twoC2;
    private final double threeC3;
    private final double fourC4;
    private final double fiveC5;

    /**
     * Constructs a MultiplicativeSwitch that starts at 0.0, ends at 1.0,
     * and smoothly interpolates between them via a sinusoid with zero first and
     * second derivatives at 0 and 1.
     */
    public MultiplicativeSwitch() {
        this(0.0, 1.0);
    }

    /**
     * Constructs a multiplicative switch which starts at off and ends at cut,
     * which smoothly interpolates between 0-1 across that range, with
     * zero first and second derivatives at off and cut.
     *
     * @param off Zero point of the switch
     * @param cut End point of the switch
     */
    public MultiplicativeSwitch(double off, double cut) {

        this.off = off;
        this.cut = cut;

        double off2 = off * off;
        double cut2 = cut * cut;

        double denom = pow(off - cut, 5.0);
        c0 = off * off2 * (off2 - 5.0 * off * cut + 10.0 * cut2) / denom;
        c1 = -30.0 * off2 * cut2 / denom;
        c2 = 30.0 * (off2 * cut + off * cut2) / denom;
        c3 = -10.0 * (off2 + 4.0 * off * cut + cut2) / denom;
        c4 = 15.0 * (off + cut) / denom;
        c5 = -6.0 / denom;
        twoC2 = 2.0 * c2;
        threeC3 = 3.0 * c3;
        fourC4 = 4.0 * c4;
        fiveC5 = 5.0 * c5;
    }

    /**
     * Value of the switching function at r.
     *
     * @param r r
     * @return Value of switch at r
     */
    public double taper(double r) {
        // Minimize number of multiply operations by storing r^2, r^3.
        double r2 = r * r;
        double r3 = r2 * r;
        return taper(r, r2, r3, r2 * r2, r3 * r2);
    }

    /**
     * First derivative of the switching function at r.
     *
     * @param r r
     * @return First derivative of switch at r
     */
    public double dtaper(double r) {
        // Minimize number of multiply operations by storing r^2.
        double r2 = r * r;
        return dtaper(r, r2, r2 * r, r2 * r2);
    }

    /**
     * Value of the switching function at r.
     *
     * @param r  r
     * @param r2 r^2
     * @param r3 r^3
     * @param r4 r^4
     * @param r5 r^5
     * @param r2 r^2
     * @param r3 r^3
     * @param r4 r^4
     * @param r5 r^5
     * @param r2 r^2
     * @param r2 r^2
     * @param r3 r^3
     * @param r3 r^3
     * @param r4 r^4
     * @param r4 r^4
     * @param r5 r^5
     * @param r5 r^5
     * @param r2 r^2
     * @param r2 r^2
     * @param r3 r^3
     * @param r3 r^3
     * @param r4 r^4
     * @param r4 r^4
     * @param r5 r^5
     * @param r5 r^5
     * @return Value of switch at r
     */
    public double taper(double r, double r2, double r3, double r4, double r5) {
        return c5 * r5 + c4 * r4 + c3 * r3 + c2 * r2 + c1 * r + c0;
    }

    /**
     * First derivative of the switching function at r.
     *
     * @param r  r
     * @param r2 r^2
     * @param r3 r^3
     * @param r4 r^4
     * @param r2 r^2
     * @param r3 r^3
     * @param r4 r^4
     * @param r2 r^2
     * @param r2 r^2
     * @param r3 r^3
     * @param r3 r^3
     * @param r4 r^4
     * @param r4 r^4
     * @param r2 r^2
     * @param r2 r^2
     * @param r3 r^3
     * @param r3 r^3
     * @param r4 r^4
     * @param r4 r^4
     * @return First derivative of switch at r
     */
    public double dtaper(double r, double r2, double r3, double r4) {
        return fiveC5 * r4 + fourC4 * r3 + threeC3 * r2 + twoC2 * r + c1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getZeroBound() {
        return off < cut ? off : cut;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getOneBound() {
        return cut > off ? cut : off;
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
        return 2;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean symmetricToUnity() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double valueAt(double x) throws IllegalArgumentException {
        return taper(x);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double firstDerivative(double x) {
        return dtaper(x);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double secondDerivative(double x) {
        /*double val = 20.0 * c5 * x*x*x;
        val += c4 * 12.0 * x*x;
        val += c3 * 6.0 *x;
        val += c2 * 2.0;*/
        double x2 = x * x;
        double val = 20.0 * c5 * x2 * x;
        val += 12.0 * c4 * x2;
        val += 6.0 * c3 * x;
        val += 2.0 * c2;
        return val;
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
                return dtaper(x);
            case 2:
                return secondDerivative(x);
            case 3:
                double val = 60 * c5 * x * x;
                val += 24 * c4 * x;
                val += 6 * c3;
                return val;
            case 4:
                val = 120 * c5 * x;
                val += 24 * c4;
                return val;
            case 5:
                return 120 * c5;
            default:
                return 0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return String.format("Multiplicative switch of form f(x) = %8.4g*x^5 + "
                        + "%8.4g*x^4 + %8.4g*x^3 + %8.4g*x^2 + %8.4g*x + %8.4g",
                c5, c4, c3, c2, c1, c0);
    }
}
