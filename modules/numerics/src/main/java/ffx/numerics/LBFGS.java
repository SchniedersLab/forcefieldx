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
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import java.util.Arrays;
import java.util.logging.Logger;

import ffx.numerics.LineSearch.LineSearchResult;

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
    public static final double slopMax = 1.0e4;
    public static final double angMax = 180.0;
    public static final int intMax = 5;

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
     * @param mSave The number of corrections used in the BFGS update.
     *		Values of <code>mSave</code> less than 3 are not recommended;
     *		large values of <code>mSave</code> will result in excessive
     *		computing time. <code>3 &lt;= mSave &lt;= 7</code> is recommended.
     *		Restriction: <code>mSave &gt; 0</code>.
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
     * @param potential Implements the {@link Optimizable} interface
     *        to supply function values and gradients.
     *
     * @param listener Implements the {@link OptimizationListener} interface
     *        and will be notified after each successful step.
     *
     * @return status code (0 = success, -1 = failed)
     *
     * @since 1.0
     */
    public static int minimize(int n, int mSave, double[] x, double f, double[] g,
                               double eps, Potential potential,
                               OptimizationListener listener) {
        return minimize(n, mSave, x, f, g, eps, Integer.MAX_VALUE, potential, listener);

    }

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
     * @param mSave The number of corrections used in the BFGS update.
     *		Values of <code>mSave</code> less than 3 are not recommended;
     *		large values of <code>mSave</code> will result in excessive
     *		computing time. <code>3 &lt;= mSave &lt;= 7</code> is recommended.
     *		Restriction: <code>mSave &gt; 0</code>.
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
     * @param potential Implements the {@link Optimizable} interface
     *        to supply function values and gradients.
     *
     * @param listener Implements the {@link OptimizationListener} interface
     *        and will be notified after each successful step.
     *
     * @return status code (0 = success, -1 = failed)
     *
     * @since 1.0
     */
    public static int minimize(int n, int mSave, double[] x, double f, double[] g,
                               double eps, int maxIterations, Potential potential,
                               OptimizationListener listener) {

        assert (n > 0);
        assert (mSave > 0);
        assert (maxIterations > 0);
        assert (x != null && x.length >= n);
        assert (g != null && g.length >= n);

        if (mSave > n) {
            logger.warning(format("Resetting the number of saved L-BFGS vectors to %d.", n));
            mSave = n;
        }

        int iterations = 0;
        int evaluations = 1;
        int nErrors = 0;
        int maxErrors = 2;

        double rms = sqrt(n);
        double scaling[] = potential.getScaling();
        if (scaling == null) {
            scaling = new double[n];
            Arrays.fill(scaling, 1.0);
        }

        /**
         * Initial search direction is the steepest decent direction.
         */
        double s[][] = new double[mSave][n];
        double y[][] = new double[mSave][n];
        for (int i = 0; i < n; i++) {
            s[0][i] = -g[i];
        }

        double grms = 0.0;
        double gnorm = 0.0;
        for (int i = 0; i < n; i++) {
            double gi = g[i];
            if (gi == Double.NaN
                || gi == Double.NEGATIVE_INFINITY
                || gi == Double.POSITIVE_INFINITY) {
                String message = format("The gradient of variable %d is %8.3f.", i, gi);
                logger.warning(message);
                return 1;
            }
            double gis = gi * scaling[i];
            gnorm += gi * gi;
            grms += gis * gis;
        }
        gnorm = sqrt(gnorm);
        grms = sqrt(grms) / rms;

        /**
         * Notify the listeners of initial conditions.
         */
        if (listener != null) {
            if (!listener.optimizationUpdate(iterations, evaluations, grms, 0.0, f, 0.0, 0.0, null)) {
                /**
                 * Terminate the optimization.
                 */
                return 1;
            }
        } else {
            log(iterations, evaluations, grms, 0.0, f, 0.0, 0.0, null);
        }

        /**
         * The convergence criteria may already be satisfied.
         */
        if (grms <= eps) {
            return 0;
        }

        double prevF = f;
        double prevX[] = new double[n];
        double prevG[] = new double[n];

        double r[] = new double[n];
        double p[] = new double[n];
        double h0[] = new double[n];
        double q[] = new double[n];

        double alpha[] = new double[mSave];
        double rho[] = new double[mSave];
        double gamma = 1.0;

        /**
         * Line search parameters.
         */
        LineSearch lineSearch = new LineSearch(n);
        LineSearchResult info[] = {LineSearchResult.Success};
        int nFunctionEvals[] = {0};
        double angle[] = {0.0};
        double df = 0.5 * stepMax * gnorm;
        int m = -1;

        while (true) {
            iterations++;

            if (iterations > maxIterations) {
                logger.info(" Maximum number of iterations reached.");
                return 0;
            }


            int muse = min(iterations - 1, mSave);
            m++;
            if (m > mSave - 1) {
                m = 0;
            }

            /**
             * Estimate the Hessian Diagonal.
             */
            fill(h0, gamma);
            arraycopy(g, 0, q, 0, n);
            int k = m;
            for (int j = 0; j < muse; j++) {
                k--;
                if (k < 0) {
                    k = mSave - 1;
                }
                alpha[k] = XdotY(n, s[k], 0, 1, q, 0, 1);
                alpha[k] *= rho[k];
                aXplusY(n, -alpha[k], y[k], 0, 1, q, 0, 1);
            }
            for (int i = 0; i < n; i++) {
                r[i] = h0[i] * q[i];
            }
            for (int j = 0; j < muse; j++) {
                double beta = XdotY(n, r, 0, 1, y[k], 0, 1);
                beta *= rho[k];
                aXplusY(n, alpha[k] - beta, s[k], 0, 1, r, 0, 1);
                k++;
                if (k > mSave - 1) {
                    k = 0;
                }
            }

            /**
             * Set the search direction.
             */
            prevF = f;
            for (int i = 0; i < n; i++) {
                p[i] = -r[i];
            }
            arraycopy(x, 0, prevX, 0, n);
            arraycopy(g, 0, prevG, 0, n);
            /**
             * Perform the line search along the new conjugate direction
             */
            nFunctionEvals[0] = 0;
            f = lineSearch.search(n, x, f, g, p, angle, df,
                                  info, nFunctionEvals, potential);
            evaluations += nFunctionEvals[0];

            /**
             * Update variables based on the results of this iteration.
             */
            for (int i = 0; i < n; i++) {
                s[m][i] = x[i] - prevX[i];
                y[m][i] = g[i] - prevG[i];
            }
            double ys = XdotY(n, y[m], 0, 1, s[m], 0, 1);
            double yy = XdotY(n, y[m], 0, 1, y[m], 0, 1);
            gamma = abs(ys / yy);
            rho[m] = 1.0 / ys;

            /**
             * Get the sizes of the moves made during this iteration.
             */
            df = prevF - f;
            prevF = f;
            double xrms = 0.0;
            grms = 0.0;
            for (int i = 0; i < n; i++) {
                double dx = (x[i] - prevX[i]) / scaling[i];
                xrms += dx * dx;
                double gx = g[i] * scaling[i];
                grms += gx * gx;
            }
            xrms = sqrt(xrms) / rms;
            grms = sqrt(grms) / rms;

            boolean done = false;
            if (info[0] == LineSearchResult.BadIntpln
                || info[0] == LineSearchResult.IntplnErr) {
                nErrors++;
                if (nErrors >= maxErrors) {
                    done = true;
                }
            } else {
                nErrors = 0;
            }

            if (listener != null) {
                if (!listener.optimizationUpdate(iterations, evaluations,
                                                 grms, xrms, f, df, angle[0], info[0])) {
                    /**
                     * Terminate the optimization.
                     */
                    return 1;
                }
            } else {
                log(iterations, evaluations, grms, xrms, f, df, angle[0], info[0]);
            }

            /**
             * Terminate the optimization if the line search failed or upon
             * satisfying the convergence criteria.
             */
            if (done) {
                return -1;
            } else if (grms <= eps) {
                return 0;
            }
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
