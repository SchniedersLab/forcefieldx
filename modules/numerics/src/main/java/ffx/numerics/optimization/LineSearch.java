/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.numerics.optimization;

import static java.lang.System.arraycopy;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.numerics.Potential;
import static ffx.numerics.optimization.LBFGS.ANGLEMAX;
import static ffx.numerics.optimization.LBFGS.CAPPA;
import static ffx.numerics.optimization.LBFGS.INTMAX;
import static ffx.numerics.optimization.LBFGS.SLOPEMAX;
import static ffx.numerics.optimization.LBFGS.STEPMAX;
import static ffx.numerics.optimization.LBFGS.STEPMIN;
import static ffx.numerics.optimization.LBFGS.aV1PlusV2;
import static ffx.numerics.optimization.LBFGS.v1DotV2;

/**
 * This class implements an algorithm for uni-dimensional line search. This file
 * is a translation of FORTRAN code written by Jay Ponder.<br>
 *
 * @author Michael J. Schnieders
 * <br>
 * Derived from Jay Ponder's FORTRAN code (search.f).
 * @since 1.0
 */
public class LineSearch {

    /**
     * The six possible line search results (Success, WideAngle, ScaleStep, IntplnErr, ReSearch, BadIntpln).
     */
    public enum LineSearchResult {

        Success, WideAngle, ScaleStep, IntplnErr, ReSearch, BadIntpln
    }

    /**
     * Number of parameters to optimize.
     */
    private final int n;
    /**
     * Implementation of the energy and gradient for the system.
     */
    private Potential optimizationSystem;
    /**
     * Number of function evaluations (pass by reference).
     */
    private int[] functionEvaluations;
    /**
     * Line search result (pass by reference).
     */
    private LineSearchResult[] info;
    /**
     * Array of current coordinates.
     */
    private double[] x;
    /**
     * The gradient array.
     */
    private double[] g;
    /**
     * Step direction.
     */
    private final double[] s;
    /**
     * Storage for a copy of the parameters.
     */
    private final double[] x0;
    /**
     * Double step size.
     */
    private double step;
    /**
     * Interpolation.
     */
    private int interpolation;
    /**
     * Function values.
     */
    private double f0, fA, fB, fC;
    /**
     * Dot products.
     */
    private double sg0, sgA, sgB, sgC;
    /**
     * True if a restart is allowed (set to true at the beginning of the algorithm).
     */
    private boolean restart;

    /**
     * LineSearch constructor.
     *
     * @param n Number of variables.
     * @since 1.0
     */
    LineSearch(int n) {
        s = new double[n];
        x0 = new double[n];
        this.n = n;
    }

    /**
     * Begin the parabolic extrapolation procedure.
     */
    private double begin() {
        restart = true;
        interpolation = 0;
        fB = f0;
        sgB = sg0;
        return step();
    }

    /**
     * Replace last point by latest and take another step.
     */
    private double step() {
        fA = fB;
        sgA = sgB;
        aV1PlusV2(n, step, s, 0, 1, x, 0, 1);

        // Get new function and projected gradient following a step
        functionEvaluations[0]++;
        fB = optimizationSystem.energyAndGradient(x, g);
        sgB = v1DotV2(n, s, 0, 1, g, 0, 1);

        // Scale step size if initial gradient change is too large
        if (abs(sgB / sgA) >= SLOPEMAX && restart) {
            arraycopy(x0, 0, x, 0, n);
            step /= 10.0;
            info[0] = LineSearchResult.ScaleStep;
            begin();
        }
        restart = false;

        /*
          We now have an appropriate step size. Return if the gradient is small
          and function decreases.
         */
        if (abs(sgB / sg0) <= CAPPA && fB < fA) {
            if (info[0] == null) {
                info[0] = LineSearchResult.Success;
            }
            f0 = fB;
            sg0 = sgB;
            return f0;
        }

        // Interpolate if gradient changes sign or function increases.
        if (sgB * sgA < 0.0 || fB > fA) {
            return cubic();
        }

        /*
          If the finite difference curvature is negative double the step; or if
          step is less than parabolic estimate less than 4 * step use this
          estimate, otherwise truncate to step or 4 * step, respectively.
         */
        step = 2.0 * step;
        if (sgB > sgA) {
            double parab = (fA - fB) / (sgB - sgA);
            if (parab > 2.0 * step) {
                parab = 2.0 * step;
            }
            if (parab < 0.5 * step) {
                parab = 0.5 * step;
            }
            step = parab;
        }
        if (step > STEPMAX) {
            step = STEPMAX;
        }
        return step();
    }

    /**
     * Beginning of the cubic interpolation procedure.
     */
    private double cubic() {
        interpolation++;
        double sss = 3.0 * (fB - fA) / step - sgA - sgB;
        double ttt = sss * sss - sgA * sgB;
        if (ttt < 0.0) {
            info[0] = LineSearchResult.IntplnErr;
            f0 = fB;
            sg0 = sgB;
            return f0;
        }
        ttt = sqrt(ttt);
        double cube = step * (sgB + ttt + sss) / (sgB - sgA + 2.0 * ttt);
        if (cube < 0 || cube > step) {
            info[0] = LineSearchResult.IntplnErr;
            f0 = fB;
            sg0 = sgB;
            return f0;
        }
        aV1PlusV2(n, -cube, s, 0, 1, x, 0, 1);

        // Get new function and gradient, then test for termination.
        functionEvaluations[0]++;
        fC = optimizationSystem.energyAndGradient(x, g);
        sgC = v1DotV2(n, s, 0, 1, g, 0, 1);
        if (abs(sgC / sg0) <= CAPPA) {
            if (info[0] == null) {
                info[0] = LineSearchResult.Success;
            }
            f0 = fC;
            sg0 = sgC;
            return f0;
        }
        /*
          Get the next pair of bracketing points by replacing one of the
          current brackets with the interpolated point
         */
        if (fC <= fA || fC <= fB) {
            double cubstp = min(abs(cube), abs(step - cube));
            if (cubstp >= STEPMIN && interpolation < INTMAX) {
                if (sgA * sgB < 0.0) {
                    /*
                      If the current brackets have slopes of opposite sign,
                      then substitute the interpolated points for the bracket
                      point with slope of same sign as the interpolated point
                     */
                    if (sgA * sgC < 0.0) {
                        fB = fC;
                        sgB = sgC;
                        step = step - cube;
                    } else {
                        fA = fC;
                        sgA = sgC;
                        step = cube;
                        aV1PlusV2(n, cube, s, 0, 1, x, 0, 1);
                    }
                } else {
                    /*
                      If current brackets have slope of same sign, then replace
                      the far bracket if the interpolated point has a slope of
                      the opposite sign or a lower function value than the near
                      bracket, otherwise replace the far bracket point.
                     */
                    if (sgA * sgC < 0.0 || fA <= fC) {
                        fB = fC;
                        sgB = sgC;
                        step -= cube;
                    } else {
                        fA = fC;
                        sgA = sgC;
                        step = cube;
                        aV1PlusV2(n, cube, s, 0, 1, x, 0, 1);
                    }
                }
                return cubic();
            }
        }

        // Interpolation has failed; reset to best current point
        double f1 = min(fA, min(fB, fC));
        double sg1;
        if (f1 == fA) {
            sg1 = sgA;
            aV1PlusV2(n, cube - step, s, 0, 1, x, 0, 1);
        } else if (f1 == fB) {
            sg1 = sgB;
            aV1PlusV2(n, cube, s, 0, 1, x, 0, 1);
        } else {
            sg1 = sgC;
        }

        // Try to restart from best point with smaller step size.
        if (f1 > f0) {
            functionEvaluations[0]++;
            f0 = optimizationSystem.energyAndGradient(x, g);
            sg0 = v1DotV2(n, s, 0, 1, g, 0, 1);
            info[0] = LineSearchResult.IntplnErr;
            return f0;
        }
        f0 = f1;
        sg0 = sg1;
        if (sg1 > 0.0) {
            for (int i = 0; i < n; i++) {
                s[i] *= -1.0;
            }
            sg0 = -sg1;
        }
        step = max(cube, step - cube) / 10.0;
        if (step < STEPMIN) {
            step = STEPMIN;
        }

        // If already restarted once, then return with the best point.
        if (info[0] == LineSearchResult.ReSearch) {
            functionEvaluations[0]++;
            f0 = optimizationSystem.energyAndGradient(x, g);
            sg0 = v1DotV2(n, s, 0, 1, g, 0, 1);
            info[0] = LineSearchResult.BadIntpln;
            return f0;
        } else {
            // Begin again.
            info[0] = LineSearchResult.ReSearch;
            return begin();
        }
    }

    /**
     * Minimize a function along a search direction.
     * <p>
     * This is a unidimensional line search based upon parabolic extrapolation
     * and cubic interpolation using both function and gradient values; if
     * forced to search in an uphill direction, return is after the initial
     * step.
     *
     * @param n                   Number of variables.
     * @param x                   Current variable values.
     * @param f                   Current function value.
     * @param g                   Current gradient values.
     * @param p                   Search direction.
     * @param angle               Angle between the gradient and search direction.
     * @param fMove               Change in function value due to previous step.
     * @param info                Line search result.
     * @param functionEvaluations Number of function evaluations.
     * @param optimizationSystem  Instance of the {@link ffx.numerics.Potential} interface.
     * @return The final function value.
     * @since 1.0
     */
    public double search(int n, double[] x, double f, double[] g,
                         double[] p, double[] angle, double fMove, LineSearchResult[] info,
                         int[] functionEvaluations, Potential optimizationSystem) {

        assert (n > 0);

        // Initialize the line search.
        this.x = x;
        this.g = g;
        this.optimizationSystem = optimizationSystem;
        this.functionEvaluations = functionEvaluations;
        this.info = info;
        fA = 0.0;
        fB = 0.0;
        fC = 0.0;
        sgA = 0.0;
        sgB = 0.0;
        sgC = 0.0;

        // Zero out the status indicator.
        info[0] = null;

        // Copy the search direction p into a new vector s.
        arraycopy(p, 0, s, 0, n);

        // Compute the length of the gradient and search direction.
        double gNorm = sqrt(v1DotV2(n, g, 0, 1, g, 0, 1));
        double sNorm = sqrt(v1DotV2(n, s, 0, 1, s, 0, 1));

        /*
          Store the initial function, then normalize the search vector and find
          the projected gradient.
         */
        f0 = f;
        arraycopy(x, 0, x0, 0, n);
        for (int i = 0; i < n; i++) {
            s[i] /= sNorm;
        }
        sg0 = v1DotV2(n, s, 0, 1, g, 0, 1);

        /*
          Check the angle between the search direction and the negative
          gradient vector.
         */
        double cosang = -sg0 / gNorm;
        cosang = min(1.0, max(-1.0, cosang));
        angle[0] = toDegrees(acos(cosang));
        if (angle[0] > ANGLEMAX) {
            info[0] = LineSearchResult.WideAngle;
            return f;
        }

        /*
          Set the initial stepSize to the length of the passed search vector,
          or based on previous function decrease.
         */
        step = 2.0 * abs(fMove / sg0);
        step = min(step, sNorm);
        if (step > STEPMAX) {
            step = STEPMAX;
        }
        if (step < STEPMIN) {
            step = STEPMIN;
        }

        return begin();
    }
}
