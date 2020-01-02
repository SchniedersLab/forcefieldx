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
package ffx.numerics.special;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Static methods to evaluate erf(x) and erfc(x) for a real argument x. Rational
 * functions are used that approximate erf(x) and erfc(x) to machine precision
 * (approximately 15 decimal digits).
 * <p>
 * Adapted from an original program written by W. J. Cody, Mathematics and
 * Computer Science Division, Argonne National Laboratory, Argonne, IL 60439
 * <p>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Erf {

    private Erf() {
    }

    /**
     * Mathematical and machine-dependent constants.
     */
    private static final double sqrtPI = 1.0 / sqrt(PI);
    private static final double oneSixteenth = 1.0 / 16.0;
    private static final double thresh = 0.46875;
    /**
     * xsmall argument below which erf(x) may be represented by 2*x/sqrt(pi) and
     * above which x*x won't underflow.
     * <p>
     * A conservative value is the largest machine number X such that 1.0 + X =
     * 1.0 to machine precision
     */
    private static final double xsmall = 1.11e-16;
    /**
     * xbig largest argument acceptable for erfc.
     * <p>
     * Solution to the equation:
     * <p>
     * W(x) * (1-0.5/x**2) = XMIN, where
     * <p>
     * W(x) = exp(-x*x)/[x*sqrt(pi)]
     */
    private static final double xbig = 26.543;

    /**
     * Evaluates erf(x) for a real argument x.
     *
     * @param arg the value to evaluate erf at.
     * @return erf of the argument.
     * @since 1.0
     */
    public static double erf(double arg) {
        return erfCore(arg, false);
    }

    /**
     * Evaluate erfc(x) for a real argument x.
     *
     * @param arg the value to evaluate erfc at.
     * @return erfc of the argument.
     * @since 1.0
     */
    public static double erfc(double arg) {
        return erfCore(arg, true);
    }

    /**
     * Evaluates erf(x) or erfc(x) for a real argument x. When called with mode
     * = false, erf is returned, while with mode = true, erfc is returned.
     *
     * @param arg  the value to evaluate erf or erfc at.
     * @param mode if mode is true, evaluate erfc, otherwise evaluate erf.
     * @return if (!mode) erf(arg), else erfc(arg)
     * @since 1.0
     */
    private static double erfCore(double arg, boolean mode) {

        // Store the argument and its absolute value.
        final double x = arg;
        final double y = abs(x);
        double result = 0.0;

        // Evaluate error function for |x| less than 0.46875.
        if (y <= thresh) {
            double ysq = 0.0;
            if (y > xsmall) {
                ysq = y * y;
            }
            double xnum = 1.85777706184603153e-1 * ysq;
            double xden = ysq;
            xnum = (xnum + 3.16112374387056560e0) * ysq;
            xden = (xden + 2.36012909523441209e1) * ysq;
            xnum = (xnum + 1.13864154151050156e2) * ysq;
            xden = (xden + 2.44024637934444173e2) * ysq;
            xnum = (xnum + 3.77485237685302021e2) * ysq;
            xden = (xden + 1.28261652607737228e3) * ysq;
            result = x * (xnum + 3.20937758913846947e3)
                    / (xden + 2.84423683343917062e3);
            if (mode) {
                result = 1.0 - result;
            }
        } else if (y <= 4.0) {
            // Get complementary error function for 0.46875 <= |x| <= 4.0.
            double xnum = 2.15311535474403846e-8 * y;
            double xden = y;
            xnum = (xnum + 5.64188496988670089e-1) * y;
            xden = (xden + 1.57449261107098347e1) * y;
            xnum = (xnum + 8.88314979438837594e0) * y;
            xden = (xden + 1.17693950891312499e2) * y;
            xnum = (xnum + 6.61191906371416295e1) * y;
            xden = (xden + 5.37181101862009858e2) * y;
            xnum = (xnum + 2.98635138197400131e2) * y;
            xden = (xden + 1.62138957456669019e3) * y;
            xnum = (xnum + 8.81952221241769090e2) * y;
            xden = (xden + 3.29079923573345963e3) * y;
            xnum = (xnum + 1.71204761263407058e3) * y;
            xden = (xden + 4.36261909014324716e3) * y;
            xnum = (xnum + 2.05107837782607147e3) * y;
            xden = (xden + 3.43936767414372164e3) * y;
            result = (xnum + 1.23033935479799725e3)
                    / (xden + 1.23033935480374942e3);
            double ysq = floor(16.0 * y) * oneSixteenth;
            double del = (y - ysq) * (y + ysq);
            result = exp(-ysq * ysq - del) * result;
            if (!mode) {
                result = 1.0 - result;
                if (x < 0.0) {
                    result = -result;
                }
            } else if (x < 0.0) {
                result = 2.0 - result;
            }
        } else {
            // Get complementary error function for |x| greater than 4.0.
            if (y < xbig) {
                double iy = 1.0 / y;
                double ysq = iy * iy;
                double xnum = 1.63153871373020978e-2 * ysq;
                double xden = ysq;
                xnum = (xnum + 3.05326634961232344e-1) * ysq;
                xden = (xden + 2.56852019228982242e0) * ysq;
                xnum = (xnum + 3.60344899949804439e-1) * ysq;
                xden = (xden + 1.87295284992346047e0) * ysq;
                xnum = (xnum + 1.25781726111229246e-1) * ysq;
                xden = (xden + 5.27905102951428412e-1) * ysq;
                xnum = (xnum + 1.60837851487422766e-2) * ysq;
                xden = (xden + 6.05183413124413191e-2) * ysq;
                result = ysq * (xnum + 6.58749161529837803e-4)
                        / (xden + 2.33520497626869185e-3);
                result = (sqrtPI - result) * iy;
                ysq = floor(16.0 * y) * oneSixteenth;
                double del = (y - ysq) * (y + ysq);
                result = exp(-ysq * ysq - del) * result;
            }
            if (!mode) {
                result = 1.0 - result;
                if (x < 0.0) {
                    result = -result;
                }
            } else {
                if (x < 0.0) {
                    result = 2.0 - result;
                }
            }
        }
        return result;
    }
}
