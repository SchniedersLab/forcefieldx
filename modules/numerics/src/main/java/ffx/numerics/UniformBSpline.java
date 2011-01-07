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
package ffx.numerics;

/**
 * Static methods to generate and differentiate uniform b-Splines.
 * <p>
 * C. de Boor, A Practical Guide to Splines. (Springer, New York, 2001)
 * 
 * @author Michael J. Schnieders
 * @see <a href="http://www.wikipedia.org/wiki/B-spline"
 *      target="_blank">b-Splines at Wikipedia</a><br>
 *      <a href="http://mathworld.wolfram.com/B-Spline.html"
 *      target="_blank">b-Splines at MathWorld</a><br>
 * @since 1.0
 */
public class UniformBSpline {

    private UniformBSpline() {
    }

    /**
     * Generate uniform b-Spline coefficients.
     *
     * @param x     A double in the range [0.0, 1.0].
     * @param order     b-Spline order (degree + 1).
     * @param coefficients  b-Spline coefficients (n coefficients for order n).
     *
     * @since 1.0
     */
    public static void bSpline(final double x, final int order,
                               final double coefficients[]) {
        // Initialization to get to a linear b-Spline (degree 1).
        coefficients[0] = 1.0 - x;
        coefficients[1] = x;
        // Apply b-Spline recursion to desired degree.
        for (int k = 2; k < order; k++) {
            bSplineRecur(x, k, coefficients, coefficients);
        }
    }

    /**
     * Uniform b-Spline recursion.
     *
     * @param x
     *            A double in the range [0.0, 1.0].
     * @param order
     *            Current b-Spline order.
     * @param coefficients
     *            Current b-Spline coefficients.
     * @param newCoefficients
     *            New b-Spline coefficients for order + 1.
     *
     * @since 1.0
     */
    private static void bSplineRecur(final double x, final int order,
                                     final double coefficients[], final double newCoefficients[]) {
        double div, k1mw;
        int i, km1, kmi;
        div = 1.0 / order;
        k1mw = order + 1 - x;
        km1 = order - 1;
        newCoefficients[order] = div * x * coefficients[km1];
        for (i = 1; i < order; i++) {
            kmi = order - i;
            newCoefficients[kmi] = div * ((x + i) * coefficients[km1 - i] + (k1mw - i) * coefficients[kmi]);
        }
        newCoefficients[0] = div * (1.0 - x) * coefficients[0];
    }

    /**
     * Generate uniform b-Spline coefficients and their derivatives.
     *
     * @param x     A double in the range [0.0, 1.0].
     * @param order     b-Spline order (degree + 1).
     * @param deriveOrder Derivative order.<br>
     *                    0 = no derivative.<br>
     *                    1 = 1rst derivative.<br>
     *                    It must not be greater than the b-Spline degree (order - 1).<br>
     *                    The method is currently limited to deriveOrder <= 5.<br>
     * @param coefficients  The b-Spline coefficient array of size [order][deriveOrder + 1].
     * @param work  A work array of size [order][order].
     * 
     * @since 1.0
     */
    public static void bSplineDerivatives(final double x, final int order,
                                          final int deriveOrder, final double coefficients[][],
                                          final double work[][]) {

        assert (deriveOrder <= order - 1 && deriveOrder <= 5);

        int j, k, o1, o2, o3, o4, o5, o6, dr_ord1;
        double tk[];
        // initialization to get to 2nd order
        work[1][0] = 1.0 - x;
        work[1][1] = x;
        // perform one pass to get to 3rd order
        work[2][0] = 0.5 * (1.0 - x) * work[1][0];
        work[2][1] = 0.5 * ((x + 1.0) * work[1][0] + (2.0 - x) * work[1][1]);
        work[2][2] = 0.5 * x * work[1][1];
        // compute standard B-spline recursion to desired order
        for (k = 3; k < order; k++) {
            bSplineRecur(x, k, work[k - 1], work[k]);
        }
        o1 = order - 1;
        // do derivatives
        if (deriveOrder > 0) {
            o2 = order - 2;
            bSplineDiff(work[o2], o1);
            if (deriveOrder > 1) {
                o3 = order - 3;
                bSplineDiff(work[o3], o2);
                bSplineDiff(work[o3], o1);
                if (deriveOrder > 2) {
                    o4 = order - 4;
                    bSplineDiff(work[o4], o3);
                    bSplineDiff(work[o4], o2);
                    bSplineDiff(work[o4], o1);
                    if (deriveOrder > 3) {
                        o5 = order - 5;
                        bSplineDiff(work[o5], o4);
                        bSplineDiff(work[o5], o3);
                        bSplineDiff(work[o5], o2);
                        bSplineDiff(work[o5], o1);
                        if (deriveOrder > 4) {
                            o6 = order - 6;
                            bSplineDiff(work[o6], o5);
                            bSplineDiff(work[o6], o4);
                            bSplineDiff(work[o6], o3);
                            bSplineDiff(work[o6], o2);
                            bSplineDiff(work[o6], o1);
                            if (deriveOrder > 5) {
                                //throw new Exception("Unsupported option: dr_ord > 5");
                            }
                        }
                    }
                }
            }
        }
        dr_ord1 = deriveOrder + 1;
        for (k = 0; k < order; k++) {
            tk = coefficients[k];
            for (j = 0; j < dr_ord1; j++) {
                tk[j] = work[o1 - j][k];
            }
        }
    }

    /**
     * Differentiate a uniform b-Spline in place.
     *
     * @param coefficients
     *            B-Spline coefficients.
     * @param order
     *            B-Spline order.
     *
     * @since 1.0
     */
    private static void bSplineDiff(final double coefficients[], final int order) {
        int km1, i;
        km1 = order - 1;
        coefficients[order] = coefficients[km1];
        for (i = km1; i > 0; i--) {
            coefficients[i] = coefficients[i - 1] - coefficients[i];
        }
        coefficients[0] = -coefficients[0];
    }
}
