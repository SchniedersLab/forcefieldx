//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.numerics.spline;

/**
 * Static methods to generate and differentiate uniform b-Splines.
 *
 * @author Michael J. Schnieders
 * @see <ul>
 * <li>
 * <a href="http://www.springer.com/mathematics/analysis/book/978-0-387-95366-3"
 * target="_blank">
 * C. de Boor, A Practical Guide to Splines. (Springer, New York, 2001)
 * </a>
 * </li>
 * <li>
 * <a href="http://www.wikipedia.org/wiki/B-spline" target="_blank">b-Splines at
 * Wikipedia</a>
 * </li>
 * <li>
 * <a href="http://mathworld.wolfram.com/B-Spline.html"
 * target="_blank">b-Splines at MathWorld</a>
 * </li>
 * </ul>
 * @since 1.0
 */
public class UniformBSpline {

    /**
     * Do not allow instantiation of UniformBSpline. All methods are static.
     */
    private UniformBSpline() {
    }

    /**
     * Generate uniform b-Spline coefficients.
     *
     * @param x            A double in the range [0.0, 1.0].
     * @param order        b-Spline order (degree + 1).
     * @param coefficients b-Spline coefficients (n coefficients for order n).
     * @since 1.0
     */
    public static void bSpline(final double x, final int order,
                               final double[] coefficients) {

        // Initialization to get to a linear b-Spline (degree 1).
        coefficients[0] = 1.0 - x;
        coefficients[1] = x;

        // Apply b-Spline recursion to desired degree.
        for (int k = 2; k < order; k++) {
            bSplineRecursion(x, k, coefficients, coefficients);
        }
    }

    /**
     * Uniform b-Spline recursion.
     *
     * @param x               A double in the range [0.0, 1.0].
     * @param order           Current b-Spline order.
     * @param coefficients    Current b-Spline coefficients.
     * @param newCoefficients New b-Spline coefficients for order + 1.
     * @since 1.0
     */
    private static void bSplineRecursion(final double x, final int order,
                                         final double[] coefficients, final double[] newCoefficients) {

        final double div = 1.0 / order;
        final double k1mw = order + 1 - x;
        final int km1 = order - 1;

        newCoefficients[order] = div * x * coefficients[km1];
        for (int i = 1; i < order; i++) {
            int kmi = order - i;
            newCoefficients[kmi] = div * ((x + i) * coefficients[km1 - i] + (k1mw - i) * coefficients[kmi]);
        }
        newCoefficients[0] = div * (1.0 - x) * coefficients[0];
    }

    /**
     * Generate uniform b-Spline coefficients and their derivatives.
     *
     * @param x            A double in the range [0.0, 1.0].
     * @param order        b-Spline order (degree + 1).
     * @param deriveOrder  Derivative order.
     *                     <br>
     *                     0 = no derivative.
     *                     <br>
     *                     1 = 1rst derivative.
     *                     <br> It must not be greater than the b-Spline degree (order - 1).
     *                     <br>
     *                     The method is currently limited to deriveOrder .LE. 5.
     *                     <br>
     * @param coefficients The b-Spline coefficient array of size
     *                     [order][deriveOrder + 1].
     * @param work         A work array of size [order][order].
     * @since 1.0
     */
    public static void bSplineDerivatives(final double x, final int order,
                                          final int deriveOrder, final double[][] coefficients,
                                          final double[][] work) {

        assert (deriveOrder <= order - 1 && deriveOrder <= 5);

        // Initialization to get to 2nd order.
        work[1][0] = 1.0 - x;
        work[1][1] = x;

        // Perform one pass to get to 3rd order.
        work[2][0] = 0.5 * (1.0 - x) * work[1][0];
        work[2][1] = 0.5 * ((x + 1.0) * work[1][0] + (2.0 - x) * work[1][1]);
        work[2][2] = 0.5 * x * work[1][1];

        // Compute standard B-spline recursion to desired order.
        for (int k = 3; k < order; k++) {
            bSplineRecursion(x, k, work[k - 1], work[k]);
        }
        int o1 = order - 1;

        // do derivatives
        try {
            if (deriveOrder > 0) {
                int o2 = order - 2;
                bSplineDiff(work[o2], o1);
                if (deriveOrder > 1) {
                    int o3 = order - 3;
                    bSplineDiff(work[o3], o2);
                    bSplineDiff(work[o3], o1);
                    if (deriveOrder > 2) {
                        int o4 = order - 4;
                        bSplineDiff(work[o4], o3);
                        bSplineDiff(work[o4], o2);
                        bSplineDiff(work[o4], o1);
                        if (deriveOrder > 3) {
                            int o5 = order - 5;
                            bSplineDiff(work[o5], o4);
                            bSplineDiff(work[o5], o3);
                            bSplineDiff(work[o5], o2);
                            bSplineDiff(work[o5], o1);
                            if (deriveOrder > 4) {
                                int o6 = order - 6;
                                bSplineDiff(work[o6], o5);
                                bSplineDiff(work[o6], o4);
                                bSplineDiff(work[o6], o3);
                                bSplineDiff(work[o6], o2);
                                bSplineDiff(work[o6], o1);
                                if (deriveOrder > 5) {
                                    throw new Exception(" Unsupported option: dr_ord > 5");
                                }
                            }
                        }
                    }
                }
            }
        } catch (Exception e) {
            // Index out of bounds if deriveOrder is too high for order.
        }

        int deriveOrder1 = deriveOrder + 1;
        for (int k = 0; k < order; k++) {
            double[] tk = coefficients[k];
            try {
                for (int j = 0; j < deriveOrder1; j++) {
                    tk[j] = work[o1 - j][k];
                }
            } catch (Exception e) {
                // Index out of bounds if deriveOrder is too high for order.
            }
        }
    }

    /**
     * Differentiate a uniform b-Spline in place.
     *
     * @param coefficients B-Spline coefficients.
     * @param order        B-Spline order.
     * @since 1.0
     */
    private static void bSplineDiff(final double[] coefficients, final int order) {
        final int order1 = order - 1;
        coefficients[order] = coefficients[order1];
        for (int i = order1; i > 0; i--) {
            coefficients[i] = coefficients[i - 1] - coefficients[i];
        }
        coefficients[0] = -coefficients[0];
    }
}
