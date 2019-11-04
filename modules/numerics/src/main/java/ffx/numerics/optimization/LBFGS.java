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
package ffx.numerics.optimization;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.Potential;
import ffx.numerics.optimization.LineSearch.LineSearchResult;

/**
 * This class implements the limited-memory Broyden-Fletcher-Goldfarb-Shanno
 * (L-BFGS) algorithm for large-scale multidimensional unconstrained
 * optimization problems.<br>
 *
 * @author Michael J. Schnieders<br> Derived from:
 * <br>
 * Robert Dodier's Java translation of original FORTRAN code by Jorge Nocedal.
 * <br>
 * J. Nocedal, "Updating Quasi-Newton Matrices with Limited Storage", Mathematics of Computation, 35, 773-782 (1980)
 * <br>
 * D. C. Lui and J. Nocedal, "On the Limited Memory BFGS Method for Large Scale Optimization", Mathematical Programming, 45, 503-528 (1989)
 * <br>
 * J. Nocedal and S. J. Wright, "Numerical Optimization", Springer-Verlag, New York, 1999, Section 9.1
 * @since 1.0
 */
public class LBFGS {

    private static final Logger logger = Logger.getLogger(LBFGS.class.getName());
    /**
     * Controls the accuracy of the line search.
     * <p>
     * If the function and gradient evaluations are inexpensive with respect to
     * the cost of the iteration (which is sometimes the case when solving very
     * large problems) it may be advantageous to set <code>CAPPA</code> to a
     * small value. A typical small value is 0.1. Restriction:
     * <code>CAPPA</code> should be greater than 1e-4.
     */
    static final double CAPPA = 0.9;
    /**
     * This specifies the lower bound for the step in the line search.
     * <p>
     * The default value is 1.0e-16. This value need not be modified unless the
     * problem is extremely badly scaled (in which case the exponent should be
     * increased).
     */
    static final double STEPMIN = 1.0e-16;
    /**
     * This specifies the upper bound for the step in the line search.
     * <p>
     * The default value is 5.0. This value need not be modified unless the
     * problem is extremely badly scaled (in which case the exponent should be
     * increased).
     */
    static final double STEPMAX = 5.0;
    /**
     * Constant <code>SLOPEMAX=1.0e4</code>
     */
    static final double SLOPEMAX = 1.0e4;
    /**
     * Constant <code>ANGLEMAX=180.0</code>
     */
    static final double ANGLEMAX = 180.0;
    /**
     * Constant <code>INTMAX=5</code>
     */
    static final int INTMAX = 5;

    /**
     * Make the constructor private so that the L-BFGS cannot be instantiated.
     */
    private LBFGS() {
    }

    /**
     * This method solves the unconstrained minimization problem
     * <pre>
     *     min f(x),    x = (x1,x2,...,x_n),
     * </pre> using the limited-memory BFGS method. The routine is especially
     * effective on problems involving a large number of variables. In a typical
     * iteration of this method an approximation <code>Hk</code> to the inverse
     * of the Hessian is obtained by applying <code>m</code> BFGS updates to a
     * diagonal matrix <code>Hk0</code>, using information from the previous
     * <code>m</code> steps.
     * <p>
     * The user specifies the number <code>m</code>, which determines the amount
     * of storage required by the routine.
     * <p>
     * The user is required to calculate the function value <code>f</code> and
     * its gradient <code>g</code>.
     * <p>
     * The steplength is determined at each iteration by means of the line
     * search routine <code>lineSearch</code>, which is a slight modification of
     * the routine <code>CSRCH</code> written by More' and Thuente.
     *
     * @param n             The number of variables in the minimization problem.
     *                      Restriction: <code>n &gt; 0</code>.
     * @param mSave         The number of corrections used in the BFGS update. Values of
     *                      <code>mSave</code> less than 3 are not recommended; large values of
     *                      <code>mSave</code> will result in excessive computing time.
     *                      <code>3 &lt;= mSave &lt;= 7</code> is recommended. *	Restriction:
     *                      <code>mSave &gt; 0</code>.
     * @param x             On initial entry this must be set by the user to the values of
     *                      the initial estimate of the solution vector. On exit it contains the
     *                      values of the variables at the best point found (usually a solution).
     * @param f             The value of the function <code>f</code> at the point
     *                      <code>x</code>.
     * @param g             The components of the gradient <code>g</code> at the point
     *                      <code>x</code>.
     * @param eps           Determines the accuracy with which the solution is to be
     *                      found. The subroutine terminates when <code>G RMS &lt; EPS</code>
     * @param maxIterations Maximum number of optimization steps.
     * @param potential     Implements the {@link ffx.numerics.Potential} interface to supply
     *                      function values and gradients.
     * @param listener      Implements the {@link OptimizationListener} interface and
     *                      will be notified after each successful step.
     * @return status code (0 = success, 1 = max iterations reached, -1 =
     * failed)
     * @since 1.0
     */
    public static int minimize(final int n, int mSave, final double[] x, double f, double[] g,
                               final double eps, final int maxIterations, Potential potential,
                               OptimizationListener listener) {

        assert (n > 0);
        assert (mSave > 0);
        assert (maxIterations > 0);
        assert (x != null && x.length >= n);
        assert (g != null && g.length >= n);

        if (mSave > n) {
            logger.fine(format(" Resetting the number of saved L-BFGS vectors to %d.", n));
            mSave = n;
        }

        int iterations = 0;
        int evaluations = 1;
        int nErrors = 0;
        int maxErrors = 2;

        double rms = sqrt(n);
        double[] scaling = potential.getScaling();
        if (scaling == null) {
            scaling = new double[n];
            fill(scaling, 1.0);
        }

        // Initial search direction is the steepest decent direction.
        double[][] s = new double[mSave][n];
        double[][] y = new double[mSave][n];
        for (int i = 0; i < n; i++) {
            s[0][i] = -g[i];
        }

        double grms = 0.0;
        double gnorm = 0.0;
        for (int i = 0; i < n; i++) {
            double gi = g[i];
            if (isNaN(gi) || isInfinite(gi)) {
                logger.warning(format("The gradient of variable %d is %8.3f.", i, gi));
                return 1;
            }
            double gis = gi * scaling[i];
            gnorm += gi * gi;
            grms += gis * gis;
        }
        gnorm = sqrt(gnorm);
        grms = sqrt(grms) / rms;

        // Notify the listeners of initial conditions.
        if (listener != null) {
            if (!listener.optimizationUpdate(iterations, evaluations, grms, 0.0, f, 0.0, 0.0, null)) {
                // Terminate the optimization.
                return 1;
            }
        } else {
            log(iterations, evaluations, grms, 0.0, f, 0.0, 0.0, null);
        }

        // The convergence criteria may already be satisfied.
        if (grms <= eps) {
            return 0;
        }

        final double[] prevX = new double[n];
        final double[] prevG = new double[n];
        final double[] r = new double[n];
        final double[] p = new double[n];
        final double[] h0 = new double[n];
        final double[] q = new double[n];
        final double[] alpha = new double[mSave];
        final double[] rho = new double[mSave];
        double gamma = 1.0;

        // Line search parameters.
        final LineSearch lineSearch = new LineSearch(n);
        final LineSearchResult[] info = {LineSearchResult.Success};
        final int[] nFunctionEvals = {0};
        final double[] angle = {0.0};
        double df = 0.5 * STEPMAX * gnorm;
        int m = -1;

        while (true) {
            iterations++;
            if (iterations > maxIterations) {
                logger.info(format(" Maximum number of iterations reached: %d.",
                        maxIterations));
                return 1;
            }

            int muse = min(iterations - 1, mSave);
            m++;
            if (m > mSave - 1) {
                m = 0;
            }

            // Estimate the Hessian Diagonal.
            fill(h0, gamma);
            arraycopy(g, 0, q, 0, n);
            int k = m;
            for (int j = 0; j < muse; j++) {
                k--;
                if (k < 0) {
                    k = mSave - 1;
                }
                alpha[k] = v1DotV2(n, s[k], 0, 1, q, 0, 1);
                alpha[k] *= rho[k];
                aV1PlusV2(n, -alpha[k], y[k], 0, 1, q, 0, 1);
            }
            for (int i = 0; i < n; i++) {
                r[i] = h0[i] * q[i];
            }
            for (int j = 0; j < muse; j++) {
                double beta = v1DotV2(n, r, 0, 1, y[k], 0, 1);
                beta *= rho[k];
                aV1PlusV2(n, alpha[k] - beta, s[k], 0, 1, r, 0, 1);
                k++;
                if (k > mSave - 1) {
                    k = 0;
                }
            }

            // Set the search direction.
            for (int i = 0; i < n; i++) {
                p[i] = -r[i];
            }
            arraycopy(x, 0, prevX, 0, n);
            arraycopy(g, 0, prevG, 0, n);

            // Perform the line search along the new conjugate direction.
            nFunctionEvals[0] = 0;
            double prevF = f;
            f = lineSearch.search(n, x, f, g, p, angle, df,
                    info, nFunctionEvals, potential);
            evaluations += nFunctionEvals[0];

            // Update variables based on the results of this iteration.
            for (int i = 0; i < n; i++) {
                s[m][i] = x[i] - prevX[i];
                y[m][i] = g[i] - prevG[i];
            }
            double ys = v1DotV2(n, y[m], 0, 1, s[m], 0, 1);
            double yy = v1DotV2(n, y[m], 0, 1, y[m], 0, 1);
            gamma = abs(ys / yy);
            rho[m] = 1.0 / ys;

            // Get the sizes of the moves made during this iteration.
            df = prevF - f;
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
                    logger.log(Level.OFF, " Algorithm failure: bad interpolation.");
                    done = true;
                }
            } else {
                nErrors = 0;
            }

            if (listener != null) {
                if (!listener.optimizationUpdate(iterations, evaluations,
                        grms, xrms, f, df, angle[0], info[0])) {
                    // Terminate the optimization.
                    return 1;
                }
            } else {
                log(iterations, evaluations, grms, xrms, f, df, angle[0], info[0]);
            }

            /*
              Terminate the optimization if the line search failed or upon
              satisfying the convergence criteria.
             */
            if (done) {
                return -1;
            } else if (grms <= eps) {
                return 0;
            }
        }
    }

    /**
     * This method solves the unconstrained minimization problem
     * <pre>
     *     min f(x),    x = (x1,x2,...,x_n),
     * </pre> using the limited-memory BFGS method. The routine is especially
     * effective on problems involving a large number of variables. In a typical
     * iteration of this method an approximation <code>Hk</code> to the inverse
     * of the Hessian is obtained by applying <code>m</code> BFGS updates to a
     * diagonal matrix <code>Hk0</code>, using information from the previous
     * <code>m</code> steps.
     * <p>
     * The user specifies the number <code>m</code>, which determines the amount
     * of storage required by the routine.
     * <p>
     * The user is required to calculate the function value <code>f</code> and
     * its gradient <code>g</code>.
     * <p>
     * The steplength is determined at each iteration by means of the line
     * search routine <code>lineSearch</code>, which is a slight modification of
     * the routine <code>CSRCH</code> written by More' and Thuente.
     *
     * @param n         The number of variables in the minimization problem.
     *                  Restriction: <code>n &gt; 0</code>.
     * @param mSave     The number of corrections used in the BFGS update. Values of
     *                  <code>mSave</code> less than 3 are not recommended; large values of
     *                  <code>mSave</code> will result in excessive computing time.
     *                  <code>3 &lt;= mSave &lt;= 7</code> is recommended. *	Restriction:
     *                  <code>mSave &gt; 0</code>.
     * @param x         On initial entry this must be set by the user to the values of
     *                  the initial estimate of the solution vector. On exit it contains the
     *                  values of the variables at the best point found (usually a solution).
     * @param f         The value of the function <code>f</code> at the point
     *                  <code>x</code>.
     * @param g         The components of the gradient <code>g</code> at the point
     *                  <code>x</code>.
     * @param eps       Determines the accuracy with which the solution is to be
     *                  found. The subroutine terminates when <code>G RMS &lt; EPS</code>
     * @param potential Implements the {@link ffx.numerics.Potential} interface to supply
     *                  function values and gradients.
     * @param listener  Implements the {@link OptimizationListener} interface and
     *                  will be notified after each successful step.
     * @return status code (0 = success, -1 = failed)
     * @since 1.0
     */
    public static int minimize(int n, int mSave, double[] x, double f, double[] g,
                               double eps, Potential potential,
                               OptimizationListener listener) {
        return minimize(n, mSave, x, f, g, eps, Integer.MAX_VALUE - 1, potential, listener);
    }

    /**
     * Print status messages for <code>LBFGS</code> if there is no listener.
     *
     * @param iter  Number of iterations so far.
     * @param nfun  Number of function evaluations so far.
     * @param grms  Gradient RMS at current solution.
     * @param xrms  Coordinate change RMS at current solution.
     * @param f     Function value at current solution.
     * @param df    Change in the function value compared to the previous solution.
     * @param angle Current angle between gradient and search direction.
     * @since 1.0
     */
    private static void log(int iter, int nfun, double grms, double xrms,
                            double f, double df, double angle, LineSearchResult info) {
        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n");
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
     * @param n       The number of points.
     * @param a       The scalar.
     * @param v1      The X array.
     * @param v1Start The first point in the X array.
     * @param v1Step  The X array increment.
     * @param v2      The Y array.
     * @param v2Start The first point in the Y array.
     * @param v2Step  The Y array increment.
     * @since 1.0
     */
    static void aV1PlusV2(final int n, final double a,
                          final double[] v1, final int v1Start, final int v1Step,
                          final double[] v2, final int v2Start, final int v2Step) {
        /*
          Require the number of entries (n) to be greater than zero. If the
          scalar (a) is zero, then the v2 array is unchanged.
         */
        if (n <= 0 || a == 0) {
            return;
        }

        int stop = v1Start + v1Step * n;
        for (int i = v1Start, j = v2Start; i != stop; i += v1Step, j += v2Step) {
            v2[j] += a * v1[i];
        }
    }

    /**
     * Compute the dot product of two vectors.
     *
     * @param n       Number of entries to include.
     * @param v1      The X array.
     * @param v1Start The first point in the X array.
     * @param v1Step  The X array increment.
     * @param v2      The Y array.
     * @param v2Start The first point in the Y array.
     * @param v2Step  The Y increment.
     * @return dot product
     * @since 1.0
     */
    static double v1DotV2(final int n,
                          final double[] v1, final int v1Start, final int v1Step,
                          final double[] v2, final int v2Start, final int v2Step) {

        // Require the number of entries to be greater than zero.
        if (n <= 0) {
            return 0;
        }

        double sum = 0.0;
        int stop = v1Start + v1Step * n;
        for (int i = v1Start, j = v2Start; i != stop; i += v1Step, j += v2Step) {
            sum += v1[i] * v2[j];
        }
        return sum;
    }
}
