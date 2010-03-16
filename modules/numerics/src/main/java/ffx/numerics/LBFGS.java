/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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

import static java.lang.Math.*;
import static java.lang.System.arraycopy;

import ffx.numerics.LineSearch.LineSearchResult;
import java.util.logging.Logger;

/**
 * This class implements the limited-memory Broyden-Fletcher-Goldfarb-Shanno
 * (L-BFGS) algorithm for large-scale multidimensional unconstrained
 * optimization problems.<br>
 *
 * @author Michael J. Schnieders<br>
 *         Derived from:<br>
 *         Robert Dodier's Java translation of orignal FORTRAN code by Jorge
 *         Nocedal.
 * @see <a href="http://www.jstor.org/stable/2006193">
 *      J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage",
 *      Mathematics of Computation, 35, 773-782 (1980)</a><br>
 *      <a href="http://dx.doi.org/10.1007/BF01589116" target="_blank">
 *      D. C. Liu and J. Nocedal, "On the Limited Memory BFGS Method for Large
 *      Scale Optimization", Mathematical Programming 45 (3), 503 (1989).</a><br>
 *      <a href="http://www.springer.com/math/book/978-0-387-30303-1">
 *      J. Nocedal and S. J. Wright, "Numerical Optimization",
 *      Springer-Verlag New York, 1999, Section 9.1</a><br>
 *      <a href="http://www.netlib.org/opt/lbfgs_um.shar" target="_blank">
 *      Nocedal's original FORTRAN code at Netlib</a>
 * @since 1.0
 */
public class LBFGS {

    private static final Logger logger = Logger.getLogger(LBFGS.class.getName());
    /**
     * Controls the accuracy of the line search. 
     * 
     * If the function and gradient evaluations are inexpensive with respect
     * to the cost of the iteration (which is sometimes the case when
     * solving very large problems) it may be advantageous to set 
     * <code>cappa</code> to a small value. A typical small value is 0.1.
     * Restriction: <code>cappa</code> should be greater than 1e-4.
     */
    public static final double cappa = 0.9;
    /**
     * This specifies the lower bound for the step in the line search.
     *
     * The default value is 1.0e-16. This value need not be modified unless
     * the problem is extremely badly scaled (in which case the exponent
     * should be increased).
     */
    public static final double stepMin = 1.0e-16;
    /**
     * This specifies the upper bound for the step in the line search.
     *
     * The default value is 5.0. This value need not be modified unless
     * the problem is extremely badly scaled (in which case the exponent
     * should be increased).
     */
    public static final double stepMax = 5.0;

    /**
     * This method solves the unconstrained minimization problem
     * <pre>
     *     min f(x),    x = (x1,x2,...,x_n),
     * </pre>
     * using the limited-memory BFGS method. The routine is especially
     * effective on problems involving a large number of variables. In
     * a typical iteration of this method an approximation <code>Hk</code> to 
     * the inverse of the Hessian is obtained by applying <code>m</code> BFGS
     * updates to a diagonal matrix <code>Hk0</code>, using information from the
     * previous <code>m</code> steps.
     *
     * The user specifies the number <code>m</code>, which determines the amount
     * of storage required by the routine.
     *
     * The user is required to calculate the function value <code>f</code> and 
     * its gradient <code>g</code>.
     *
     * The steplength is determined at each iteration by means of the
     * line search routine <code>lineSearch</code>, which is a slight
     * modification of the routine <code>CSRCH</code> written
     * by More' and Thuente.
     *
     * @param n The number of variables in the minimization problem.
     *		Restriction: <code>n &gt; 0</code>.
     *
     * @param m The number of corrections used in the BFGS update.
     *		Values of <code>m</code> less than 3 are not recommended;
     *		large values of <code>m</code> will result in excessive
     *		computing time. <code>3 &lt;= m &lt;= 7</code> is recommended.
     *		Restriction: <code>m &gt; 0</code>.
     *
     * @param x On initial entry this must be set by the user to the values
     *		of the initial estimate of the solution vector. On exit it
     *          contains the values of the variables at the best point found
     *          (usually a solution).
     *
     * @param f The value of the function <code>f</code> at the
     *          point <code>x</code>.
     *
     * @param g The components of the gradient <code>g</code> at the
     *          point <code>x</code>.
     *
     * @param eps Determines the accuracy with which the solution
     *		is to be found. The subroutine terminates when
     *		<code>
     *            G RMS &lt; EPS
     *		</code>
     *
     * @param optimizationSystem Implements the {@link Optimizable} interface
     *        to supply function values and gradients.
     *
     * @param listener Implements the {@link OptimizationListener} interface
     *        and will be notified after each successful step.
     *
     * @return status code (0 = success, -1 = failed)
     *
     * @since 1.0
     */
    public static int minimize(int n, int m, double[] x, double f, double[] g,
            double eps, Optimizable optimizationSystem, OptimizationListener listener) {

        assert (n > 0);
        assert (m > 0);
        assert (x != null && x.length >= n);
        assert (g != null && g.length >= n);

        boolean finish = false;
        int iterations = 0;
        int functionEvaluations = 1;
        int point = 0;
        double scaling[] = optimizationSystem.getOptimizationScaling();
        double rms = sqrt(n) / sqrt(3.0);
        double previousF = f;
        double ftol = 0.0000001;
        double previousX[] = new double[n];
        arraycopy(x, 0, previousX, 0, n);

        /**
         * Initialize inverse Hessian diagonal elements.
         */
        double diag[] = new double[n];
        for (int i = 0; i < n; i++) {
            diag[i] = 1.0;
        }
        int len = n * (2 * m + 1) + 2 * m;
        int ispt = n + 2 * m;
        int iypt = ispt + n * m;

        /**
         * Initial search direction is the steepest decent direction.
         */
        double w[] = new double[len];
        for (int i = 0; i < n; i++) {
            w[ispt + i] = -g[i];
        }
        double grms = 0.0;
        double gnorm = 0.0;
        for (int i = 0; i < n; i++) {
            double gi = g[i];
            double gis = gi * scaling[i];
            gnorm = gi * gi;
            grms += gis * gis;
        }
        gnorm = sqrt(gnorm);
        grms = sqrt(grms) / rms;
        double stp1 = 1.0 / gnorm;
        double angle = 0.0;
        if (listener != null) {
            if (!listener.optimizationUpdate(iterations, functionEvaluations, grms, 0.0, f, f - previousF, angle, null)) {
                /**
                 * Terminate the optimization.
                 */
                return 1;
            }
        } else {
            log(iterations, functionEvaluations, grms, 0.0, f, f - previousF, angle, null);
        }
        /**
         * The convergence criteria may already be satisfied.
         */
        if (grms <= eps) {
            return 0;
        }
        int npt = 0;
        LineSearchResult info[] = new LineSearchResult[1];
        int nfev[] = new int[1];
        double stp[] = new double[1];
        while (true) {
            iterations++;
            previousF = f;
            int bound = iterations - 1;
            if (iterations != 1) {
                if (iterations > m) {
                    bound = m;
                }
                double ys = XdotY(n, w, iypt + npt, 1, w, ispt + npt, 1);
                double yy = XdotY(n, w, iypt + npt, 1, w, iypt + npt, 1);
                for (int i = 0; i < n; i++) {
                    diag[i] = ys / yy;
                }
                int cp = point;
                if (point == 0) {
                    cp = m;
                }
                w[n + cp - 1] = 1.0 / ys;
                for (int i = 0; i < n; i++) {
                    w[i] = -g[i];
                }
                cp = point;
                for (int i = 0; i < bound; i++) {
                    cp = cp - 1;
                    if (cp == -1) {
                        cp = m - 1;
                    }
                    int inmc = n + m + cp;
                    int iycn = iypt + cp * n;
                    w[inmc] = w[n + cp] * XdotY(n, w, ispt + cp * n, 1, w, 0, 1);
                    aXplusY(n, -w[inmc], w, iycn, 1, w, 0, 1);
                }
                for (int i = 0; i < n; i++) {
                    w[i] = diag[i] * w[i];
                }
                for (int i = 0; i < bound; i++) {
                    double beta = w[n + cp] * XdotY(n, w, iypt + cp * n, 1, w, 0, 1);
                    int inmc = n + m + cp;
                    beta = w[inmc] - beta;
                    int iscn = ispt + cp * n;
                    aXplusY(n, beta, w, iscn, 1, w, 0, 1);
                    cp = cp + 1;
                    if (cp == m) {
                        cp = 0;
                    }
                }
                for (int i = 0; i < n; i++) {
                    w[ispt + point * n + i] = w[i];
                }
            }
            stp[0] = 1.0;
            if (iterations == 1) {
                stp[0] = stp1;
            }
            angle = 0.0;
            gnorm = 0.0;
            double snorm = 0.0;
            int searchStart = ispt + point * n;
            for (int i = 0; i < n; i++) {
                double gi = g[i];
                double si = w[searchStart + i];
                w[i] = gi;
                gnorm += gi * gi;
                snorm += si * si;
                angle += gi * si;
            }
            snorm = sqrt(snorm);
            gnorm = sqrt(gnorm);
            angle = -angle / snorm / gnorm;
            angle = min(1.0, max(-1.0, angle));
            angle = toDegrees(acos(angle));
            nfev[0] = 0;
            info[0] = null;
            f = LineSearch.search(n, x, f, g, w, searchStart, stp, ftol,
                    info, nfev, diag, optimizationSystem);
            functionEvaluations += nfev[0];
            /**
             * Compute RMS gradient and RMS coordinate change.
             */
            grms = 0.0;
            double xrms = 0.0;
            for (int i = 0; i < n; i++) {
                double gs = g[i] * scaling[i];
                double xs = (x[i] - previousX[i]) / scaling[i];
                grms += gs * gs;
                xrms += xs * xs;
            }
            grms = sqrt(grms) / rms;
            xrms = sqrt(xrms) / rms;
            if (listener != null) {
                if (!listener.optimizationUpdate(iterations, functionEvaluations, grms, xrms, f, f - previousF, angle, info[0])) {
                    /**
                     * Terminate the optimization.
                     */
                    return 1;
                }
            } else {
                log(iterations, functionEvaluations, grms, xrms, f, f - previousF, angle, info[0]);
            }

            /**
             * Terminate the optimization if the line search failed or upon
             * satisfying the convergence criteria.
             */
            if (info[0] != LineSearchResult.Success) {
                return -1;
            } else if (grms <= eps) {
                return 0;
            }

            /**
             * Store gradient and search information from this iteration.
             */
            npt = point * n;
            for (int i = 0; i < n; i++) {
                w[ispt + npt + i] = stp[0] * w[ispt + npt + i];
                w[iypt + npt + i] = g[i] - w[i];
            }
            point++;
            if (point == m) {
                point = 0;
            }
            /*
             * Cache the current solution vector.
             */
            arraycopy(x, 0, previousX, 0, n);
        }
    }

    /**
     * Print status messages for <code>LBFGS</code> if there is no listener.
     *
     * @param iter Number of iterations so far.
     * @param nfun Number of function evaluations so far.
     * @param grms Gradient RMS at current solution.
     * @param xrms Coordinate change RMS at current solution.
     * @param f Function value at current solution.
     * @param f Change in the function value compared to the previous solution.
     * @param angle Current angle between gradient and search direction.
     *
     * @since 1.0
     */
    private static void log(int iter, int nfun, double grms, double xrms, double f, double df, double angle, LineSearchResult info) {
        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" QN Iter    F Value      G RMS     F Move    X Move    Angle  FG Call  Comment\n");
        }
        if (info == null) {
            logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d",
                    iter, f, grms, df, xrms, angle, nfun));
        } else {
            logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d   %8s",
                    iter, f, grms, df, xrms, angle, nfun, info.toString()));
        }
    }

    /**
     * Compute the sum of a vector times a scalar plus another vector.
     *
     * @param n Number of points.
     * @param a scalar.
     * @param x X array.
     * @param x0 First point in the X array.
     * @param dx X increment.
     * @param y Y Array.
     * @param y0 First point in the Y array.
     * @param dy Y increment.
     *
     * @since 1.0
     */
    public static void aXplusY(int n, double a, double[] x, int x0, int dx, double[] y, int y0, int dy) {
        if (n <= 0 || a == 0) {
            return;
        }
        int stop = x0 + dx * n;
        for (int i = x0, j = y0; i != stop; i += dx, j += dy) {
            y[j] += a * x[i];
        }
        return;
    }

    /**
     * Compute the dot product of two vectors.
     *
     * @param n Number of entries to include.
     * @param x X array.
     * @param x0 First point in the X array.
     * @param dx X increment.
     * @param y Y Array.
     * @param y0 First point in the Y array.
     * @param dy Y increment.
     *
     * @return dot product
     *
     * @since 1.0
     */
    public static double XdotY(int n, double[] x, int x0, int dx, double[] y, int y0, int dy) {
        if (n <= 0) {
            return 0;
        }
        double sum = 0.0;
        int stop = x0 + dx * n;
        for (int i = x0, j = y0; i != stop; i += dx, j += dy) {
            sum += x[i] * y[j];
        }
        return sum;
    }
}
