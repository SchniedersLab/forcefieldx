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

import static ffx.numerics.LBFGS.cappa;
import static ffx.numerics.LBFGS.stepMax;
import static ffx.numerics.LBFGS.stepMin;

/**
 * This class implements an algorithm for multi-dimensional line search. This
 * file is a translation of FORTRAN code written by Jorge Nocedal.<br>
 *
 * @author Michael J. Schnieders<br>
 *         Derived from Robert Dodier's Java translation of FORTRAN code by
 *         Jorge Noceal.<br>
 *         Derived from original FORTRAN version by Jorge J. More' and David J.
 *         Thuente as part of the Minpack project, June 1983, Argonne National
 *         Laboratory.<br>
 * @see <a href="http://www.netlib.org/opt/lbfgs_um.shar" target="_blank">
 *      Nocedal's original FORTRAN code at Netlib</a>
 * @since 1.0
 */
public class LineSearch {

    public enum LineSearchResult {

        Success, GradErr, XTol, MaxEval, StepMin, StepMax, MachPrec
    };
    private static double machinePrecision = 1.0e-16;
    private static int maxFunctionEvaluations = 20;

    /**
     * Minimize a function along a search direction.
     *
     * The purpose of <code>search</code> is to find a step which satisfies
     * a sufficient decrease condition and a curvature condition.<p>
     *
     * At each stage this function updates an interval of uncertainty with
     * endpoints <code>stx</code> and <code>sty</code>. The interval of
     * uncertainty is initially chosen so that it contains a minimizer of the
     * modified function
     * <pre>
     *      f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
     * </pre>
     *
     * If a step is obtained for which the modified function has a nonpositive
     * function value and nonnegative derivative, then the interval of
     * uncertainty is chosen so that it contains a minimizer of
     * <code>f(x+stp*s)</code>.<p>
     *
     * The algorithm is designed to find a step which satisfies the sufficient
     * decrease condition
     * <pre>
     *       f(x+stp*s) &lt;= f(X) + ftol*stp*(gradf(x)'s),
     * </pre>
     * and the curvature condition
     * <pre>
     *       abs(gradf(x+stp*s)'s)) &lt;= machinePrecision*abs(gradf(x)'s).
     * </pre>
     * If <code>ftol</code> is less than <code>machinePrecision</code> and if, for example,
     * the function is bounded below, then there is always a step which
     * satisfies both conditions. If no step can be found which satisfies both
     * conditions, then the algorithm usually stops when rounding
     * errors prevent further progress. In this case <code>stp</code> only
     * satisfies the sufficient decrease condition.<p>
     *
     * @param n The number of variables.
     * @param x On entry this contains the base point for the line search.
     *		On exit it contains <code>x + stp*s</code>.
     * @param f On entry this contains the value of the objective function
     *		at <code>x</code>. On exit the method returns the value of the
     *          objective function at <code>x + stp*s</code>.
     * @param g On entry this contains the gradient of the objective function
     *		at <code>x</code>. On exit it contains the gradient at
     *		<code>x + stp*s</code>.
     * @param s The search direction.
     * @param is0 Starting index for the search direction array.
     * @param stepSize On entry this contains an initial estimate of a satifactory
     *		step length. On exit <code>stp</code> contains the final
     *          estimate.
     * @param ftol Tolerance for the sufficient decrease condition.
     * @param info This is the output status:
     *		<ul>
     *		<li><code>Success</code>  The sufficient decrease condition and
     *			the directional derivative condition hold.
     *		<li><code>Xtol</code>     Relative width of the interval of
     *                  uncertainty is at most <code>xtol</code>.
     *		<li><code>MaxEval</code>  Number of function evaluations
     *                  has reached <code>maxFunctionEvaluations</code>.
     *		<li><code>StepMin</code>  The step is at the lower bound
     *                  <code>stepMin</code>.
     *		<li><code>StepMax</code>  The step is at the upper bound
     *                  <code>stepMax</code>.
     *		<li><code>MachPrec</code> Rounding errors prevent further progress.
     *			There may not be a step which satisfies the sufficient
     *                  decrease and curvature conditions. Tolerances may be
     *                  too small.
     *		</ul>
     *	@param functionEvaluations On exit, this is set to the number of
     *           function evaluations.
     *	@param wa Temporary storage array, of length <code>n</code>.
     *  @param optimizationSystem Implements the OptimizationInterface.
     *  @return the objective function at <code>x + stp*s</code>.
     *
     *  @since 1.0
     */
    public static double search(int n, double[] x, double f, double[] g,
            double[] s, int is0, double[] stepSize, double ftol,
            LineSearchResult[] info, int[] functionEvaluations, double[] wa,
            Optimizable optimizationSystem) {

        assert (n > 0);
        assert (stepSize[0] > 0);
        assert (ftol > 0);

        /**
         * Initialization.
         */
        final double half = 0.5;
        final double twoThirds = 2.0 / 3.0;
        final double four = 4.0;
        final double initialFunctionValue = f;
        functionEvaluations[0] = 0;

        double dgx[] = new double[1];
        double dgxm[] = new double[1];
        double dgy[] = new double[1];
        double dgym[] = new double[1];
        double fx[] = new double[1];
        double fxm[] = new double[1];
        double fy[] = new double[1];
        double fym[] = new double[1];
        double stx[] = new double[1];
        double sty[] = new double[1];
        boolean bracket[] = new boolean[1];

        int stepStatus = 1;

        /**
         * Compute the initial gradient in the search direction and check
         * that s is a descent direction.
         */
        double dgInit = 0.0;
        for (int j = 0; j < n; j++) {
            dgInit += g[j] * s[is0 + j];
        }
        if (dgInit >= 0) {
            info[0] = LineSearchResult.GradErr;
            return f;
        }

        bracket[0] = false;
        boolean stage1 = true;
        double dgTest = ftol * dgInit;
        double width = stepMax - stepMin;
        double doubleWidth = width / half;
        for (int j = 0; j < n; j++) {
            wa[j] = x[j];
        }

        /**
         * The variables stx, fx, dgx contain the values of the step,
         * function, and directional derivative at the best step.
         *
         * The variables sty, fy, dgy contain the value of the step,
         * function, and derivative at the other endpoint of the interval of
         * uncertainty.
         *
         * The variables stp, f, dg contain the values of the
         * step, function, and derivative at the current step.
         */
        stx[0] = 0.0;
        fx[0] = initialFunctionValue;
        dgx[0] = dgInit;
        sty[0] = 0.0;
        fy[0] = initialFunctionValue;
        dgy[0] = dgInit;

        while (true) {
            /**
             * Set the minimum and maximum steps to correspond to the
             * present interval of uncertainty.
             */
            double stmin, stmax;
            if (bracket[0]) {
                stmin = min(stx[0], sty[0]);
                stmax = max(stx[0], sty[0]);
            } else {
                stmin = stx[0];
                stmax = stepSize[0] + four * (stepSize[0] - stx[0]);
            }

            /**
             * Force the step to be within the bounds stepMax and stepMin.
             */
            stepSize[0] = max(stepSize[0], stepMin);
            stepSize[0] = min(stepSize[0], stepMax);

            /**
             * If an unusual termination is to occur then let stepSize be the
             * lowest point obtained so far.
             */
            if ((bracket[0] && (stepSize[0] <= stmin || stepSize[0] >= stmax))
                    || functionEvaluations[0] >= maxFunctionEvaluations - 1
                    || stepStatus == 0 || (bracket[0] && stmax - stmin <= machinePrecision * stmax)) {
                stepSize[0] = stx[0];
            }

            /**
             * Evaluate f and g at stepSize and compute the directional
             * derivative.
             */
            for (int j = 0; j < n; j++) {
                x[j] = wa[j] + stepSize[0] * s[is0 + j];
            }
            f = optimizationSystem.energyAndGradient(x, g);
            functionEvaluations[0] = functionEvaluations[0] + 1;
            double dg = 0.0;
            for (int j = 0; j < n; j++) {
                dg = dg + g[j] * s[is0 + j];
            }
            double ftest1 = initialFunctionValue + stepSize[0] * dgTest;

            /**
             * Test for convergence.
             */
            if ((bracket[0] && (stepSize[0] <= stmin || stepSize[0] >= stmax)) || stepStatus == 0) {
                info[0] = LineSearchResult.MachPrec;
                return f;
            }
            if (stepSize[0] == stepMax && f <= ftest1 && dg <= dgTest) {
                info[0] = LineSearchResult.StepMax;
                return f;
            }
            if (stepSize[0] == stepMin && (f > ftest1 || dg >= dgTest)) {
                info[0] = LineSearchResult.StepMin;
                return f;
            }
            if (functionEvaluations[0] >= maxFunctionEvaluations) {
                info[0] = LineSearchResult.MaxEval;
                return f;
            }
            if (bracket[0] && stmax - stmin <= machinePrecision * stmax) {
                info[0] = LineSearchResult.XTol;
                return f;
            }
            if (f <= ftest1 && abs(dg) <= cappa * (-dgInit)) {
                info[0] = LineSearchResult.Success;
                return f;
            }

            /**
             * In the first stage we seek a step for which the modified
             * function has a nonpositive value and nonnegative derivative.
             */
            if (stage1 && f <= ftest1 && dg >= min(ftol, cappa) * dgInit) {
                stage1 = false;
            }

            /**
             * A modified function is used to predict the step only if
             * we have not obtained a step for which the modified
             * function has a nonpositive function value and nonnegative
             * derivative, and if a lower function value has been
             * obtained but the decrease is not sufficient.
             */
            if (stage1 && f <= fx[0] && f > ftest1) {
                /**
                 * Define the modified function and derivative values.
                 */
                double fm = f - stepSize[0] * dgTest;
                fxm[0] = fx[0] - stx[0] * dgTest;
                fym[0] = fy[0] - sty[0] * dgTest;
                double dgm = dg - dgTest;
                dgxm[0] = dgx[0] - dgTest;
                dgym[0] = dgy[0] - dgTest;
                /**
                 * Call step to update the interval of uncertainty
                 * and to compute the new step.
                 */
                stepStatus = step(stx, fxm, dgxm, sty, fym, dgym, stepSize, fm,
                        dgm, bracket, stmin, stmax);
                /**
                 * Reset the function and gradient values for f.
                 */
                fx[0] = fxm[0] + stx[0] * dgTest;
                fy[0] = fym[0] + sty[0] * dgTest;
                dgx[0] = dgxm[0] + dgTest;
                dgy[0] = dgym[0] + dgTest;
            } else {
                /**
                 * Call step to update the interval of uncertainty and to
                 * compute the new step.
                 */
                stepStatus = step(stx, fx, dgx, sty, fy, dgy, stepSize, f,
                        dg, bracket, stmin, stmax);
            }
            /**
             * Force a sufficient decrease in the size of the interval 
             * of uncertainty.
             */
            if (bracket[0]) {
                if (abs(sty[0] - stx[0]) >= twoThirds * doubleWidth) {
                    stepSize[0] = stx[0] + half * (sty[0] - stx[0]);
                }
                doubleWidth = width;
                width = abs(sty[0] - stx[0]);
            }
        }
    }

    /**
     * The purpose of this method is to compute a safeguarded step for
     * a search and to update an interval of uncertainty for
     * a minimizer of the function.<p>
     *
     * The parameter <code>stx</code> contains the step with the least function
     * value. The parameter <code>stp</code> contains the current step. It is
     * assumed that the derivative at <code>stx</code> is negative in the
     * direction of the step. If <code>brackt[0]</code> is <code>true</code>
     * when <code>step</code> returns then a
     * minimizer has been bracketed in an interval of uncertainty
     * with endpoints <code>stx</code> and <code>sty</code>.<p>
     *
     * Variables that must be modified by <code>mcStep</code> are
     * implemented as 1-element arrays.
     *
     * @param stx Step at the best step obtained so far.
     *   This variable is modified by <code>step</code>.
     * @param fx Function value at the best step obtained so far.
     *   This variable is modified by <code>step</code>.
     * @param dx Derivative at the best step obtained so far. The derivative
     *   must be negative in the direction of the step, that is, <code>dx</code>
     *   and <code>stp-stx</code> must have opposite signs.
     *   This variable is modified by <code>step</code>.
     *
     * @param sty Step at the other endpoint of the interval of uncertainty.
     *   This variable is modified by <code>step</code>.
     * @param fy Function value at the other endpoint of the interval of uncertainty.
     *   This variable is modified by <code>step</code>.
     * @param dy Derivative at the other endpoint of the interval of
     *   uncertainty. This variable is modified by <code>mcStep</code>.
     *
     * @param stepSize Step at the current step. If <code>brackt</code> is set
     *   then on input <code>stp</code> must be between <code>stx</code>
     *   and <code>sty</code>. On output <code>stp</code> is set to the
     *   new step.
     * @param fp Function value at the current step.
     * @param dp Derivative at the current step.
     *
     * @param bracket Tells whether a minimizer has been bracketed.
     *   If the minimizer has not been bracketed, then on input this
     *   variable must be set <code>false</code>. If the minimizer has
     *   been bracketed, then on output this variable is <code>true</code>.
     *
     * @param stepMin Lower bound for the step.
     * @param stepMax Upper bound for the step.
     *
     * @return On return from <code>step</code> a status of 1, 2, 3 or 4
     *   indicates the step has been computed successfully. Otherwise a value
     *   of 0 indicates improper input parameters.
     */
    private static int step(double[] stx, double[] fx, double[] dx,
            double[] sty, double[] fy, double[] dy, double[] stp, double fp,
            double dp, boolean[] bracket, double stepMin, double stepMax) {

        boolean bound;
        double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
        int status = 0;
        if ((bracket[0] && (stp[0] <= min(stx[0], sty[0])
                || stp[0] >= max(stx[0], sty[0])))
                || dx[0] * (stp[0] - stx[0]) >= 0.0
                || stepMax < stepMin) {
            return status;
        }
        /**
         * Determine if the derivatives have opposite sign.
         */
        sgnd = dp * (dx[0] / abs(dx[0]));

        if (fp > fx[0]) {
            /**
             * First case. A higher function value. The minimum is bracketed. 
             * If the cubic step is closer to stx than the quadratic step, 
             * the cubic step is taken, else the average of the cubic and 
             * quadratic steps is taken.
             */
            status = 1;
            bound = true;
            theta = 3.0 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
            s = max3(abs(theta), abs(dx[0]), abs(dp));
            gamma = s * sqrt(sqr(theta / s) - (dx[0] / s) * (dp / s));
            if (stp[0] < stx[0]) {
                gamma = -gamma;
            }
            p = (gamma - dx[0]) + theta;
            q = ((gamma - dx[0]) + gamma) + dp;
            r = p / q;
            stpc = stx[0] + r * (stp[0] - stx[0]);
            stpq = stx[0] + ((dx[0] / ((fx[0] - fp) / (stp[0] - stx[0]) + dx[0])) / 2.0) * (stp[0] - stx[0]);
            if (abs(stpc - stx[0]) < abs(stpq - stx[0])) {
                stpf = stpc;
            } else {
                stpf = stpc + (stpq - stpc) / 2.0;
            }
            bracket[0] = true;
        } else if (sgnd < 0.0) {
            /**
             * Second case. A lower function value and derivatives of
             * opposite sign. The minimum is bracketed. If the cubic
             * step is closer to stx than the quadratic (secant) step,
             * the cubic step is taken, else the quadratic step is taken.
             */
            status = 2;
            bound = false;
            theta = 3.0 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
            s = max3(abs(theta), abs(dx[0]), abs(dp));
            gamma = s * sqrt(sqr(theta / s) - (dx[0] / s) * (dp / s));
            if (stp[0] > stx[0]) {
                gamma = -gamma;
            }
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dx[0];
            r = p / q;
            stpc = stp[0] + r * (stx[0] - stp[0]);
            stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
            if (abs(stpc - stp[0]) > abs(stpq - stp[0])) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            bracket[0] = true;
        } else if (abs(dp) < abs(dx[0])) {
            /**
             * Third case. A lower function value, derivatives of the
             * same sign, and the magnitude of the derivative decreases.
             * The cubic step is only used if the cubic tends to infinity
             * in the direction of the step or if the minimum of the cubic
             * is beyond stp. Otherwise the cubic step is defined to be
             * either stpmin or stpmax. The quadratic (secant) step is also
             * computed and if the minimum is bracketed then the the step
             * closest to stx is taken, else the step farthest away is taken.
             */
            status = 3;
            bound = true;
            theta = 3.0 * (fx[0] - fp) / (stp[0] - stx[0]) + dx[0] + dp;
            s = max3(abs(theta), abs(dx[0]), abs(dp));
            gamma = s * sqrt(max(0, sqr(theta / s) - (dx[0] / s) * (dp / s)));
            if (stp[0] > stx[0]) {
                gamma = -gamma;
            }
            p = (gamma - dp) + theta;
            q = (gamma + (dx[0] - dp)) + gamma;
            r = p / q;
            if (r < 0.0 && gamma != 0.0) {
                stpc = stp[0] + r * (stx[0] - stp[0]);
            } else if (stp[0] > stx[0]) {
                stpc = stepMax;
            } else {
                stpc = stepMin;
            }
            stpq = stp[0] + (dp / (dp - dx[0])) * (stx[0] - stp[0]);
            if (bracket[0]) {
                if (abs(stp[0] - stpc) < abs(stp[0] - stpq)) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
            } else {
                if (abs(stp[0] - stpc) > abs(stp[0] - stpq)) {
                    stpf = stpc;
                } else {
                    stpf = stpq;
                }
            }
        } else {
            /**
             * Fourth case. A lower function value, derivatives of the
             * same sign, and the magnitude of the derivative does
             * not decrease. If the minimum is not bracketed, the step
             * is either stpmin or stpmax, else the cubic step is taken.
             */
            status = 4;
            bound = false;
            if (bracket[0]) {
                theta = 3.0 * (fp - fy[0]) / (sty[0] - stp[0]) + dy[0] + dp;
                s = max3(abs(theta), abs(dy[0]), abs(dp));
                gamma = s * sqrt(sqr(theta / s) - (dy[0] / s) * (dp / s));
                if (stp[0] > sty[0]) {
                    gamma = -gamma;
                }
                p = (gamma - dp) + theta;
                q = ((gamma - dp) + gamma) + dy[0];
                r = p / q;
                stpc = stp[0] + r * (sty[0] - stp[0]);
                stpf = stpc;
            } else if (stp[0] > stx[0]) {
                stpf = stepMax;
            } else {
                stpf = stepMin;
            }
        }
        /**
         * Update the interval of uncertainty. This update does not depend on
         * the new step or the case analysis above.
         */
        if (fp > fx[0]) {
            sty[0] = stp[0];
            fy[0] = fp;
            dy[0] = dp;
        } else {
            if (sgnd < 0.0) {
                sty[0] = stx[0];
                fy[0] = fx[0];
                dy[0] = dx[0];
            }
            stx[0] = stp[0];
            fx[0] = fp;
            dx[0] = dp;
        }
        /**
         * Compute the new step and safeguard it.
         */
        stpf = min(stepMax, stpf);
        stpf = max(stepMin, stpf);
        stp[0] = stpf;
        if (bracket[0] && bound) {
            if (sty[0] > stx[0]) {
                stp[0] = min(stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
            } else {
                stp[0] = max(stx[0] + 0.66 * (sty[0] - stx[0]), stp[0]);
            }
        }
        return status;
    }

    static double sqr(double x) {
        return x * x;
    }

    static double max3(double x, double y, double z) {
        return x < y ? (y < z ? z : y) : (x < z ? z : x);
    }
}
