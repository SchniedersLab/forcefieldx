// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
// ******************************************************************************
package ffx.numerics.optimization;

import ffx.numerics.OptimizationInterface;

import static ffx.numerics.optimization.LBFGS.DEFAULT_ANGLEMAX;
import static ffx.numerics.optimization.LBFGS.DEFAULT_CAPPA;
import static ffx.numerics.optimization.LBFGS.DEFAULT_INTMAX;
import static ffx.numerics.optimization.LBFGS.DEFAULT_SLOPEMAX;
import static ffx.numerics.optimization.LBFGS.DEFAULT_STEPMAX;
import static ffx.numerics.optimization.LBFGS.DEFAULT_STEPMIN;
import static ffx.numerics.optimization.LBFGS.aV1PlusV2;
import static ffx.numerics.optimization.LBFGS.v1DotV2;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toDegrees;

/**
 * This class implements an algorithm for uni-dimensional line search using parabolic extrapolation
 * and cubic interpolation with both function and gradient values.
 *
 * <p>The algorithm works as follows:
 * <ol>
 *   <li>Initialize with the current point, function value, gradient, and search direction</li>
 *   <li>Check if the search direction makes a reasonable angle with the negative gradient</li>
 *   <li>Set an initial step size based on previous function decrease</li>
 *   <li>Perform parabolic extrapolation to find a better step size</li>
 *   <li>If needed, perform cubic interpolation to refine the step size</li>
 *   <li>Return the best point found along the search direction</li>
 * </ol>
 *
 * <p>This implementation is a translation of FORTRAN code written by Jay Ponder (search.f).
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class LineSearch {

  /**
   * Default factor to scale step size when gradient change is too large.
   */
  private static final double STEP_SCALE_FACTOR = 10.0;

  /**
   * Default factor for parabolic step size adjustment.
   */
  private static final double PARABOLIC_UPPER_LIMIT = 2.0;

  /**
   * Default factor for parabolic step size adjustment.
   */
  private static final double PARABOLIC_LOWER_LIMIT = 0.5;

  /**
   * Default factor for step size reduction during restart.
   */
  private static final double RESTART_STEP_SCALE = 10.0;

  /**
   * Number of parameters to optimize.
   */
  private final int n;

  /**
   * Step direction.
   */
  private final double[] s;

  /**
   * Storage for a copy of the parameters.
   */
  private final double[] x0;

  /**
   * Implementation of the energy and gradient for the system.
   */
  private OptimizationInterface optimizationSystem;

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
   * Current step size.
   */
  private double step;

  /**
   * Counter for the number of interpolation attempts.
   */
  private int interpolation;

  /**
   * Function values at different points:
   * f0: Initial function value
   * fA: Function value at previous point
   * fB: Function value at current point
   * fC: Function value at interpolated point
   */
  private double f0, fA, fB, fC;

  /**
   * Dot products of search direction and gradient at different points:
   * sg0: Initial dot product
   * sgA: Dot product at previous point
   * sgB: Dot product at current point
   * sgC: Dot product at interpolated point
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
   * Minimize a function along a search direction.
   *
   * <p>This method performs a unidimensional line search along the specified search direction
   * using parabolic extrapolation and cubic interpolation with both function and gradient values.
   * The search attempts to find a step size that provides sufficient decrease in the function value
   * while maintaining a reasonable angle between the search direction and negative gradient.
   *
   * <p>If the search is forced to proceed in an uphill direction (angle > DEFAULT_ANGLEMAX),
   * the method returns after the initial step without performing the search.
   *
   * <p>The method modifies the input arrays x and g to contain the coordinates and gradient
   * at the best point found along the search direction.
   *
   * @param n                   Number of variables in the optimization problem.
   * @param x                   Current variable values (modified to contain the best point found).
   * @param f                   Current function value at point x.
   * @param g                   Current gradient values at point x (modified to contain the gradient at the best point).
   * @param p                   Search direction vector.
   * @param angle               Output parameter that will contain the angle between the gradient and search direction.
   * @param fMove               Change in function value due to the previous optimization step.
   * @param info                Output parameter that will contain the line search result status.
   * @param functionEvaluations Input/output parameter for tracking the number of function evaluations.
   * @param optimizationSystem  Implementation of the {@link ffx.numerics.OptimizationInterface} that provides
   *                            function values and gradients.
   * @return The final function value at the best point found.
   * @since 1.0
   */
  public double search(int n, double[] x, double f, double[] g, double[] p, double[] angle,
                       double fMove, LineSearchResult[] info, int[] functionEvaluations,
                       OptimizationInterface optimizationSystem) {

    // Validate input parameters
    if (n <= 0) {
      throw new IllegalArgumentException("Number of variables must be positive");
    }
    if (x == null || x.length < n) {
      throw new IllegalArgumentException("Coordinate array is null or too small");
    }
    if (g == null || g.length < n) {
      throw new IllegalArgumentException("Gradient array is null or too small");
    }
    if (p == null || p.length < n) {
      throw new IllegalArgumentException("Search direction array is null or too small");
    }
    if (angle == null || angle.length < 1) {
      throw new IllegalArgumentException("Angle array is null or too small");
    }
    if (info == null || info.length < 1) {
      throw new IllegalArgumentException("Info array is null or too small");
    }
    if (functionEvaluations == null || functionEvaluations.length < 1) {
      throw new IllegalArgumentException("Function evaluations array is null or too small");
    }
    if (optimizationSystem == null) {
      throw new IllegalArgumentException("Optimization system cannot be null");
    }

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

    // Handle the case where the search direction or gradient has zero length
    if (sNorm < Double.MIN_NORMAL) {
      info[0] = LineSearchResult.IntplnErr;
      return f;
    }

    // Store the initial function, then normalize the search vector
    f0 = f;
    arraycopy(x, 0, x0, 0, n);
    for (int i = 0; i < n; i++) {
      s[i] /= sNorm;
    }

    // Find the projected gradient
    sg0 = v1DotV2(n, s, 0, 1, g, 0, 1);

    // Check the angle between the search direction and the negative gradient vector
    double cosang = (gNorm < Double.MIN_NORMAL) ? 0.0 : -sg0 / gNorm;
    cosang = min(1.0, max(-1.0, cosang)); // Ensure cosang is in [-1, 1]
    angle[0] = toDegrees(acos(cosang));
    if (angle[0] > DEFAULT_ANGLEMAX) {
      info[0] = LineSearchResult.WideAngle;
      return f;
    }

    // Set the initial step size based on previous function decrease or search vector length
    if (sg0 != 0.0) {
      step = 2.0 * abs(fMove / sg0);
    } else {
      step = sNorm;
    }
    step = min(step, sNorm);
    step = min(max(step, DEFAULT_STEPMIN), DEFAULT_STEPMAX);

    return begin();
  }

  /**
   * Begin the parabolic extrapolation procedure.
   *
   * <p>This method initializes the line search by setting up the initial conditions
   * for the parabolic extrapolation. It marks that a restart is allowed, resets the
   * interpolation counter, and sets the initial function value and gradient projection.
   *
   * @return The final function value after completing the line search.
   */
  private double begin() {
    restart = true;
    interpolation = 0;
    fB = f0;
    sgB = sg0;
    return step();
  }

  /**
   * Take a step along the search direction and evaluate the function and gradient.
   *
   * <p>This method:
   * <ol>
   *   <li>Stores the previous function value and gradient projection</li>
   *   <li>Takes a step along the search direction</li>
   *   <li>Evaluates the function and gradient at the new point</li>
   *   <li>Checks if the step size needs to be scaled down (if gradient change is too large)</li>
   *   <li>Checks if we've found a suitable point (small gradient and decreased function)</li>
   *   <li>Decides whether to interpolate or continue extrapolation</li>
   * </ol>
   *
   * @return The final function value after completing the line search.
   */
  private double step() {
    // Store previous values and take a step
    fA = fB;
    sgA = sgB;
    aV1PlusV2(n, step, s, 0, 1, x, 0, 1);

    // Get new function and projected gradient following a step
    functionEvaluations[0]++;
    fB = optimizationSystem.energyAndGradient(x, g);
    sgB = v1DotV2(n, s, 0, 1, g, 0, 1);

    // Scale step size if initial gradient change is too large
    if (abs(sgB / sgA) >= DEFAULT_SLOPEMAX && restart) {
      arraycopy(x0, 0, x, 0, n);
      step /= STEP_SCALE_FACTOR;
      info[0] = LineSearchResult.ScaleStep;
      begin();
    }
    restart = false;

    // Return if we've found a suitable point (small gradient and decreased function)
    if (abs(sgB / sg0) <= DEFAULT_CAPPA && fB < fA) {
      if (info[0] == null) {
        info[0] = LineSearchResult.Success;
      }
      f0 = fB;
      sg0 = sgB;
      return f0;
    }

    // Interpolate if gradient changes sign or function increases
    if (sgB * sgA < 0.0 || fB > fA) {
      return cubic();
    }

    // Continue extrapolation with adjusted step size
    step = 2.0 * step;

    // If the finite difference curvature is positive, use parabolic estimate
    if (sgB > sgA) {
      double parab = (fA - fB) / (sgB - sgA);
      parab = min(PARABOLIC_UPPER_LIMIT * step, max(PARABOLIC_LOWER_LIMIT * step, parab));
      step = parab;
    }

    // Ensure step size doesn't exceed maximum
    step = min(step, DEFAULT_STEPMAX);

    return step();
  }

  /**
   * Perform cubic interpolation to refine the step size.
   *
   * <p>This method implements cubic interpolation to find a better step size when:
   * <ul>
   *   <li>The gradient changes sign between two points (indicating a minimum between them)</li>
   *   <li>The function value increases (indicating we've stepped too far)</li>
   * </ul>
   *
   * <p>The cubic interpolation uses both function values and gradients at the bracketing
   * points to estimate the location of the minimum. If the interpolation fails or doesn't
   * produce a better point, the method falls back to using the best point found so far.
   *
   * @return The final function value after completing the line search.
   */
  private double cubic() {
    interpolation++;

    // Calculate cubic interpolation coefficients
    double sss = 3.0 * (fB - fA) / step - sgA - sgB;
    double ttt = sss * sss - sgA * sgB;

    // Check if cubic interpolation is possible (discriminant must be non-negative)
    if (ttt < 0.0) {
      info[0] = LineSearchResult.IntplnErr;
      f0 = fB;
      sg0 = sgB;
      return f0;
    }

    // Calculate the cubic step size
    ttt = sqrt(ttt);
    double cube = step * (sgB + ttt + sss) / (sgB - sgA + 2.0 * ttt);

    // Check if the cubic step is valid (must be between 0 and step)
    if (cube < 0 || cube > step) {
      info[0] = LineSearchResult.IntplnErr;
      f0 = fB;
      sg0 = sgB;
      return f0;
    }

    // Move to the interpolated point
    aV1PlusV2(n, -cube, s, 0, 1, x, 0, 1);

    // Evaluate function and gradient at the interpolated point
    functionEvaluations[0]++;
    fC = optimizationSystem.energyAndGradient(x, g);
    sgC = v1DotV2(n, s, 0, 1, g, 0, 1);

    // Check if we've found a suitable point (small gradient)
    if (abs(sgC / sg0) <= DEFAULT_CAPPA) {
      if (info[0] == null) {
        info[0] = LineSearchResult.Success;
      }
      f0 = fC;
      sg0 = sgC;
      return f0;
    }

    // If the interpolated point is better than at least one of the brackets,
    // continue with further interpolation
    if (fC <= fA || fC <= fB) {
      double cubstp = min(abs(cube), abs(step - cube));

      // Continue interpolation if step size is reasonable and we haven't exceeded max iterations
      if (cubstp >= DEFAULT_STEPMIN && interpolation < DEFAULT_INTMAX) {
        if (sgA * sgB < 0.0) {
          // If the current brackets have slopes of opposite sign,
          // substitute the interpolated point for the bracket point
          // with slope of same sign as the interpolated point
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
          // If current brackets have slopes of same sign, then replace
          // the far bracket if the interpolated point has a slope of
          // the opposite sign or a lower function value than the near bracket
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

    // Interpolation has failed; reset to the best current point
    double f1 = min(fA, min(fB, fC));
    double sg1;

    // Move to the best point
    if (f1 == fA) {
      sg1 = sgA;
      aV1PlusV2(n, cube - step, s, 0, 1, x, 0, 1);
    } else if (f1 == fB) {
      sg1 = sgB;
      aV1PlusV2(n, cube, s, 0, 1, x, 0, 1);
    } else {
      sg1 = sgC;
    }

    // If the best point is worse than the initial point, return to the initial point
    if (f1 > f0) {
      functionEvaluations[0]++;
      f0 = optimizationSystem.energyAndGradient(x, g);
      sg0 = v1DotV2(n, s, 0, 1, g, 0, 1);
      info[0] = LineSearchResult.IntplnErr;
      return f0;
    }

    // Update with the best point found
    f0 = f1;
    sg0 = sg1;

    // If the gradient at the best point is positive, reverse the search direction
    if (sg1 > 0.0) {
      for (int i = 0; i < n; i++) {
        s[i] *= -1.0;
      }
      sg0 = -sg1;
    }

    // Reduce the step size
    step = max(cube, step - cube) / RESTART_STEP_SCALE;
    step = max(step, DEFAULT_STEPMIN);

    // If already restarted once, then return with the best point
    if (info[0] == LineSearchResult.ReSearch) {
      functionEvaluations[0]++;
      f0 = optimizationSystem.energyAndGradient(x, g);
      sg0 = v1DotV2(n, s, 0, 1, g, 0, 1);
      info[0] = LineSearchResult.BadIntpln;
      return f0;
    } else {
      // Begin again
      info[0] = LineSearchResult.ReSearch;
      return begin();
    }
  }

  /**
   * Enum representing the possible outcomes of a line search operation.
   * <p>
   * These status codes provide information about how the line search terminated,
   * which can be useful for diagnosing optimization issues or understanding the
   * behavior of the algorithm.
   */
  public enum LineSearchResult {
    /**
     * Successful line search.
     * <p>
     * The algorithm found a point with sufficient decrease in function value
     * and a small enough gradient projection along the search direction.
     */
    Success,

    /**
     * Angle between gradient and search direction is too wide.
     * <p>
     * The search direction makes an angle with the negative gradient that exceeds
     * DEFAULT_ANGLEMAX (180 degrees). This usually indicates a problem with the
     * search direction calculation in the optimization algorithm.
     */
    WideAngle,

    /**
     * Step size was scaled down.
     * <p>
     * The initial step resulted in a gradient change that was too large
     * (exceeding DEFAULT_SLOPEMAX), so the step size was reduced to ensure
     * more stable convergence.
     */
    ScaleStep,

    /**
     * Interpolation error occurred.
     * <p>
     * An error occurred during cubic interpolation, such as a negative discriminant
     * or an invalid step size. The algorithm falls back to using the best point found.
     */
    IntplnErr,

    /**
     * Search was restarted.
     * <p>
     * The line search was restarted with a smaller step size after interpolation
     * failed to find a better point. This is an intermediate status that may lead
     * to eventual success or failure.
     */
    ReSearch,

    /**
     * Bad interpolation result.
     * <p>
     * Interpolation failed repeatedly, even after restarting the search. This
     * typically indicates a difficult optimization landscape or numerical issues.
     */
    BadIntpln
  }
}
