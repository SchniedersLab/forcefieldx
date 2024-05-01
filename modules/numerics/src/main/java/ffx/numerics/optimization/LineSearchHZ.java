// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import static ffx.numerics.optimization.LBFGS.aV1PlusV2;
import static ffx.numerics.optimization.LBFGS.v1DotV2;
import static java.lang.System.arraycopy;
import static org.apache.commons.math3.util.FastMath.*;

import ffx.numerics.OptimizationInterface;

import java.util.Arrays;
import java.util.Collections;
import java.util.logging.Logger;

/**
 * This class implements a line search described by Hager and Zhang (2004, 2005, 2006).
 *
 * @author Michael J. Schnieders <br> Derived from Jay Ponder's FORTRAN code (search.f).
 * @since 1.0
 */
public class LineSearchHZ {
  private static final Logger logger = Logger.getLogger(LineSearchHZ.class.getName()); // TODO - remove after debugging

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
   * Array of current coordinates.
   */
  private double[] x;
  /**
   * The gradient array.
   */
  private double[] g;
  /**
   * Line search result (pass by reference).
   */
//  private LineSearchHZResult[] info;
  private LineSearch.LineSearchResult[] info;
  /**
   * Implementation of the energy and gradient for the system.
   */
  private OptimizationInterface optimizationSystem;

  private final double DELTA = 0.1; // HZ paper; Nocedal & Wright recommends 0.01?
  private final double SIGMA = 0.9; // HZ paper; Nocedal & Wright recommends 0.1 for GradientDescent?
  private final double eps = 0.000001; // 1e-6 (HZ paper used for testing)
  private final double psi0 = 0.0100;
  private final double psi1 = 0.2000;
  private final double psi2 = 2.0;
  private final double psi3 = 0.1000;
  private final double rho = 5.0; // Hager and Zhang propose rho = 5
  private final double gamma = 0.66; // gamma is the decay factor for the bracketing interval (2/3)
  private double alphaMax = 100; // not sure where this is from todo
  private final int lineSearchMax = 100; // not sure where from
  private double[] alphaArr, valueArr, slopeArr;
  private int alphaIdx, valueIdx, slopeIdx;
  private int nf[]; // number of function evaluations
  private boolean first; // flag for first iteration TODO change to just get from other values
  private boolean QuadStep;

  /**
   * LineSearch constructor.
   *
   * @param n Number of variables.
   * @since 1.0
   */
  LineSearchHZ(int n) {
    s = new double[n];
    x0 = new double[n];
    this.n = n;

    first = true;
    QuadStep = false; // user's choice? todo was false
    alphaArr = new double[10000]; // [lineSearchMax]; //todo? this size?
    valueArr = new double[10000]; // [lineSearchMax];
    slopeArr = new double[10000]; // [lineSearchMax];
    alphaIdx = 0;
    valueIdx = 0;
    slopeIdx = 0;
  }

  /**
   * Minimize a function along a search direction.
   *
   * <p>This is a line search based upon Hager and Zhang (2005).
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
   * @return The final function value.
   * @since 1.0
   */
//  search(int n, double[] x, double f, double[] g, double[] p, double[] angle,
//         double fMove, LineSearch.LineSearchResult[] info, int[] functionEvaluations,
//         OptimizationInterface optimizationSystem)
  // TODO MATLAB code has it returning double[] {alpha, phi}
  public double search(int n, double[] x, double f, double[] g, double[] p, double[] angle,
                       double fMove, LineSearch.LineSearchResult[] info, int[] functionEvaluations,
                       OptimizationInterface optimizationSystem) {
    assert(n > 0);
//    logger.info(String.format("==ENTERED LINE SEARCH: n = %d",n));

    // Initialize the line search.
    this.x = x;
    this.g = g;
    this.optimizationSystem = optimizationSystem;
    this.nf = functionEvaluations;
    this.info = info;
//    logger.info(String.format("===Other parameters: x[0] = %f, x.length=%d, g[0] = %f, g.length=%d, nf = %d",x[0],x.length,g[0],g.length,nf[0]));

    // Zero out the status indicator.
    info[0] = null;

    // Copy the search direction p into a new vector s.
    arraycopy(p, 0, s, 0, n); // s = p
    arraycopy(x, 0, x0, 0, n); // copy coordinates to x0 (start)

    // look at SATURDAY todon lines


    // L0
    // I0-I2 in HZ paper, denoted initial(k)
//    double c = initialGuess(n,x,f,g,p,angle,fMove,functionEvaluations,optimizationSystem);
    double c = initialGuess(f); // using initial f given when called from LBFGS
    // log: "initial c = %f",c

    double[] tempPhi = getPhiDPhi(0);
    double phi_0 = tempPhi[0];
    double dphi_0 = tempPhi[1];

    if (!(Double.isFinite(phi_0) && Double.isFinite(dphi_0))) {
      // warning: "Value and slope at step length = 0 must be finite."
//      info[0] = LineSearchHZResult.StartErr;
      info[0] = LineSearch.LineSearchResult.IntplnErr;
    }
    if (dphi_0 >= ulp(1)*abs(phi_0)) { // 2.2204e-16
      // warning: "Search direction is not a direction of descent; dphi_0 = %f, (eps*abs(phi_0)) = %f", dphi_0, ulp(1) * abs(phi_0)
//      info[0] = LineSearchHZResult.AngleWarn;
      info[0] = LineSearch.LineSearchResult.WideAngle;
    } else if (dphi_0 >= 0) {
//      double alphas = ulp(1); // 2.2204e-16 in MATLAB
//      double values = phi_0;
      alphaArr[0] = ulp(1);
      valueArr[0] = phi_0;
      alphaIdx++;
      valueIdx++;
//      return new double[]{alphaArr[0], valueArr[0]}; // {alphas, values} todo
//      logger.info(String.format("===Returned at line 199. Step slope >= 0. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",alphaArr[0],valueArr[0],dphi_0,x[0],g[0]));
      x = x0; // SATURDAY todo
      return valueArr[0];
    }

    // Prevent values of x_new = x+c.s that are likely to make phi(x_new) infinite
    int iterFiniteMax = (int) ceil(-log(2,2.2204*pow(10,-16)));
//    double alphas = 0; // for bisection
//    double values = phi_0;
//    double slopes = dphi_0;
    alphaArr[0] = 0; // for bisection
    valueArr[0] = phi_0;
    slopeArr[0] = dphi_0;
    alphaIdx++;
    valueIdx++;
    slopeIdx++;

    // For condition (4.3) in Hager and Zhang, 2005. ε_k = ε|f(x_k)|
    double phiLim = phi_0 + eps * abs(phi_0);

    assert( Double.isFinite(c) && c <= alphaMax );

    tempPhi = getPhiDPhi(c);
    double phi_c = tempPhi[0];
    double dphi_c = tempPhi[1];

    int iterFinite = 1;

    while ( !(Double.isFinite(phi_c) && Double.isFinite(dphi_c)) && (iterFinite < iterFiniteMax) ) {
      iterFinite++;
      c = psi3 * c; // Contracts c, since psi3 = 0.1

      tempPhi = getPhiDPhi(c);
      phi_c = tempPhi[0];
      dphi_c = tempPhi[1];
    }

    if ( !(Double.isFinite(phi_c) && Double.isFinite(dphi_c)) ) {
//      info[0] = LineSearchHZResult.FailNewEvalPt; // Failed to achieve finite new evaluation point, using alpha = 0.
      info[0] = LineSearch.LineSearchResult.IntplnErr;
      alphaArr[0] = 0; // alphas = 0
      valueArr[0] = phi_0; // values = 0
//      return new double[]{alphaArr[0], valueArr[0]}; // {alphas, values} TODO
//      logger.info(String.format("===Returned at line 241. Phi/dPhi not finite. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",alphaArr[0],valueArr[0],dphi_0,x[0],g[0]));
      return valueArr[0];
    }
//    double[] alphasArr = new double[]{alphas, c};
//    double[] valuesArr = new double[]{values, phi_c};
//    double[] slopesArr = new double[]{slopes, dphi_c};
    alphaArr[alphaIdx] = c; // use iteration counter todo check if equals one here?
    valueArr[valueIdx] = phi_c;
    slopeArr[slopeIdx] = dphi_c;
    alphaIdx++;
    valueIdx++;
    slopeIdx++;

    // If c was generated by quadratic interpolation, check whether it satisfies the Wolfe conditions
    if (satisfiesWolfe(c, phi_c, dphi_c, phi_0, dphi_0, phiLim)) {
      // Save the accepted c for future use. info_global(current_lev,column_idx) TODO
//      return new double[]{c, phi_c}; //todo return c, phi_c (step/alpha and fxn value)
      info[0] = LineSearch.LineSearchResult.Success;
//      logger.info(String.format("===Returned at line 259. Satisfies Wolfe conditions. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",c,phi_c,dphi_c,x[0],g[0]));
      return phi_c;
    }

    // [a,b] = bracket(c)
    // Initial bracketing step (HZ, B0-B3).
    // The following code is used to generate an initial interval [a,b] satisfying the opposite slope condition (2.9),
    // beginning with the initial guess [0,c].
    boolean isBracketed = false;
    int ia = 0;
    int ib = 1;
//    assert( alphasArr.length == 2 );
    assert( alphaArr[1] != 0 ); // todo not exactly true, just want to say that index 0 and 1 have been placed
    int iter = 1; // todo make sure this doesn't get used as an index
    while ( !isBracketed && iter < lineSearchMax ) {
      // B1
      if (dphi_c >= 0) {
        // Reached the upward slope, so we have b; examine previous values to find a
        ib = alphaIdx-1; //alphasArr.length-1; //TODO make sure it's minus 1
        for (int i = ib-1; i >= 0; i--) { //TODO - make sure it's minus 1 and loop works
          if (valueArr[i] <= phiLim) {
            ia = i;
            break;
          }
        }
        isBracketed = true;
      } else if (valueArr[valueIdx-1] > phiLim) { // TODO same as above
        // B2
        // The value is higher, but the slope is downward, so we must have crested over the peak. Use bisection!
        ib = alphaIdx-1; // alphasArr.length-1; // This guarantees that we enter bisect with \bar{b} = c
        ia = 0; // This guarantees that we enter bisect with \bar{a} = 0

        if ( (c != alphaArr[ib]) || (slopeArr[ib] >= 0) ) {
//          info[0] = LineSearchHZResult.IntplnErr; // error: "c = num2str(c)" ??
          info[0] = LineSearch.LineSearchResult.IntplnErr;
        }

        // This implements U3a-c
        int[] iab = bisect(ia, ib, phiLim); //, x,s);
        ia = iab[0];
        ib = iab[1];

        isBracketed = true;
      } else {
        // B3
        // Will still go downhill, expand the interval and try again. Reaching this branch means that dphi_c < 0 and
        // phi_c <= phi_0 + eps. So cold = c has a lower objective than phi_0 up to epsilon. This makes it a viable step
        // to return if bracketing fails.

        // Bracketing can fail if no cold < c <= alphamax can be found with finite phi_c and dphi_c. Going back to the
        // loop with c = cold will only result in infinite cycling. So returning (cold, phi_cold) and exiting the line
        // search is the best move
        double cold = c;
        double phi_cold = phi_c;
        if (cold >= alphaMax) {
//          return new double[]{cold,phi_cold}; TODO
          info[0] = LineSearch.LineSearchResult.ReSearch;
//          logger.info(String.format("===Returned at line 315. cold (step) >= alphaMax! alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",cold,phi_cold,-69.69,x[0],g[0]));
          return phi_cold;
        }
        // c_{j+1} = rho * c_{j}, where rho is the growth factor for the bracketing interval. Hager and Zhang suggest to
        // use rho = 5.
        c = rho * c;

        if (c > alphaMax) {
          c = alphaMax;
          // log: "B3: exceeding alphamax, using c = alphamax = %d, cold = %d", alphaMax, cold)
        }

        tempPhi = getPhiDPhi(c);
        phi_c = tempPhi[0];
        dphi_c = tempPhi[1];

        iterFinite = 1;

        while (!(Double.isFinite(phi_c) && Double.isFinite(dphi_c)) && c > (cold + ulp(1)) && iterFinite < iterFiniteMax) {
          alphaMax = c; // shrinks alphamax, assumes that steps >= c can never have finite phi_c and dphi_c
          iterFinite++;
          // log: "B3: non-finite value, bisection."
          c = (cold + c) / 2;
          tempPhi = getPhiDPhi(c);
          phi_c = tempPhi[0];
          dphi_c = tempPhi[1];
        }
        if ( !(Double.isFinite(phi_c) && Double.isFinite(dphi_c)) ) {
          // log: "Warning: failed to expand interval to bracket with finite values. If this happens frequently, check your function and gradient."
          // log: "c = %f, alphamax = %f, phi_c = %f, dphi_c = %f", c, alphamax, phi_c, dphi_c
//          alphas = cold;
//          values = phi_cold;
//          return new double[]{cold, phi_cold}; TODO
          info[0] = LineSearch.LineSearchResult.IntplnErr;
//          logger.info(String.format("===Returned at line 349. Phi/dphi not Finite. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",cold,phi_cold,dphi_c,x[0],g[0]));
          return phi_cold;
        }
        alphaArr[alphaIdx] = c;
        valueArr[valueIdx] = phi_c;
        slopeArr[slopeIdx] = dphi_c;
        alphaIdx++;
        valueIdx++;
        slopeIdx++;
      } // End of B3
      iter++;
    } // End of the while loop
    if (slopeArr[ia] * slopeArr[ib] < 0) {
      // log: "End of bracketing: The opposite slope conditions is enforced: dphi(a) = %f; dphi(b) = %f", slopeArr[ia],slopeArr[ib]
      // todo? they dont have anything and don't do anything in else
      // else {}
    }
    // End of [a,b] = bracket(c) and end of L0

    // Main loop of the linesearch algorithm (L1-L2-L3)
    while (iter < lineSearchMax) {
      double a = alphaArr[ia];
      double b = alphaArr[ib];
      // assert( b > a );
      // log: "linesearch: ia = %d, ib = %d, a = %f, b = %f, phi(a) = %f, phi(b) = %f",ia,ib,a,b,valueArr[ia],valueArr[ib]
      if (b-a <= ulp(b)) {
//        return new double[]{a, valueArr[ia]}; TODO
        info[0] = LineSearch.LineSearchResult.ReSearch;
//        logger.info(String.format("===Returned at line 378. [a,b] too close. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",a,valueArr[ia],slopeArr[ia],x[0],g[0]));
        return valueArr[ia];
      }

      // L1: Take the secant step
      // secant2(alphas,values,slopes,ia,ib,phiLim,x,s,problem)
      boolean isWolfe;
      int[] WiAB = secant2(ia,ib,phiLim);
      isWolfe = WiAB[0] == 1; // isWolfe = true if isWolfe =1; false otherwise (== 0)
      int iA = WiAB[1];
      int iB = WiAB[2];
      // log: "L1: END of secant2."
      if (isWolfe) {
//        info[0] = LineSearchHZResult.Success;d
        info[0] = LineSearch.LineSearchResult.Success;
//        return new double[]{alphaArr[iA],valueArr[iA]};TODO
//        logger.info(String.format("===Returned at line 394. Satisfies Wolfe conditions. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",alphaArr[iA],valueArr[iA],slopeArr[iA],x[0],g[0]));
        return valueArr[iA];
      }
      double A = alphaArr[iA];
      double B = alphaArr[iB];
      assert( B > A );

      if ( (B - A) < gamma * (b - a) ) {
        // log: "L3: secant succeeded."
        if ( (valueArr[ia] + ulp(1)) >= valueArr[ib] && (valueArr[iA] + ulp(1)) >= valueArr[iB] ) {
          // It's so flat, secant didn't do anything useful, time to quit
          // log: "Linesearch: secant suggests it is flat."
//          info[0] = LineSearchHZResult.Flat;
          info[0] = LineSearch.LineSearchResult.WideAngle;
//          return new double[]{A, valueArr[iA]}; // return A, values(iA) TODO
//          logger.info(String.format("===Returned at line 409. Flat function. alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",alphaArr[iA],valueArr[iA],slopeArr[iA],x[0],g[0]));
          return valueArr[iA];
        }
        ia = iA;
        ib = iB;
      } else {
        // L2: bisection step, then update
        // Secant is converging too slowly, use bisection
        // log: "L2: secant failed, using bisection."

        // Bisection:
        c = (A + B) / 2;

        tempPhi = getPhiDPhi(c);
        phi_c = tempPhi[0];
        dphi_c = tempPhi[1];

        assert( Double.isFinite(phi_c) && Double.isFinite(dphi_c) );
        alphaArr[alphaIdx] = c;
        valueArr[valueIdx] = phi_c;
        slopeArr[slopeIdx] = dphi_c;
        alphaIdx++;
        valueIdx++;
        slopeIdx++;

        int[] iab = updateHZ(iA, iB, alphaIdx-1, phiLim);
      }

      // L3: increment iter and go to L1
      iter++;
    }
    // End of main loop of linesearch algorithm
    // warning: "Linesearch failed to converge, reached maximum iterations linesearchmax = %d",lineSearchMax
//    info[0] = LineSearchHZResult.NoConvergence;
    info[0] = LineSearch.LineSearchResult.BadIntpln;

//    return new double[]{c, phi_0}; //TODO should have state too?
//    logger.info(String.format("===Returned at line 446. End of loop, no convergence (return starting values). alhpa = %f, phi = %f, dphi = %f, x[0] = %f, g[0] = %f",c,phi_0,dphi_0,x[0],g[0]));
    return phi_0;
  }

  // Implements stage U3 of the update routine (with theta=0.5).
  private int[] bisect(int ia, int ib, double phiLim) { // double[] x, double[] s) {
    double a = alphaArr[ia]; // Assuming method giving correct index (ia & ib)
    double b = alphaArr[ib];

    // Debugging (HZ, conditions shown following U3)
    // assert( slopeArr[ia] < 0 );
    // assert( valueArr[ia] <= phiLim );
    // assert( slopeArr[ib] < 0 );
    // assert( valueArr[ib] > phiLim );
    // assert( b > a );

    // This is the loop encoded in lines U3a-U3c of HZ pseudocode.
    while (b-a > ulp(b)) {
      // U3a: set d = ( \bar{a} + \bar{b} ) /2
      double d = (a + b) / 2;

      double[] tempPhi = getPhiDPhi(d);
      double phi_d = tempPhi[0];
      double dphi_d = tempPhi[1];

      assert( Double.isFinite(phi_d) && Double.isFinite(dphi_d) );

      // THEY ADD TO ARRAY HERE. CAN'T DO. USING PREDETERMINED SIZE with LINESEARCHMAX
      // do we even use the prevoius values in the array? maybe dont need to store, or only store one in another var.
      alphaArr[alphaIdx] = d;
      valueArr[valueIdx] = phi_d;
      slopeArr[slopeIdx] = dphi_d;
      alphaIdx++;
      valueIdx++;
      slopeIdx++;

      int id = alphaIdx-1; // alphaArr.length; // - 1 ?? todo - check
      if (dphi_d >= 0) { // || abs(dphi_d) < ulp(1) -> 2019 MATLAB code relaxed this condition
        ib = id; // replace b, return
        return new int[]{ia, ib};
      }
      // dphi < 0, otherwise we would not reach this line of the code
      // U3b: If phi_d <= phiLim, then \bar{a} = d, and go to U3a
      if (phi_d <= phiLim) {
        a = d; // replace a, but keep bisecting until dphi_b > 0
        ia = id;
      } else {
        // U3c: If phi_d > phiLim, then \bar{b} = d
        b = d;
        ib = id;
      }
    }

    // If bisect terminates without enforcing dphi_d >= 0, throw a warning (not implemented)
    // if dphi_d < 0
    //  warning('U3 ended without enforcing the condition dphi_d >= 0 !!!!')
    // end
    return new int[]{ia,ib};
  }

  // Implements the routine [c] = initial(k) to generate the starting guess 'c' used by the bracket routine, as
  // described in HZ paper.
//  private double initialGuess(int n, double[] x, double f, double[] g, double[] p, double[] angle,
//                              double fMove, int[] functionEvaluations,
//                              OptimizationInterface optimizationSystem) {
  private double initialGuess(double f) {

    // TODO: INITIAL GUESS -- already have initial guess when this is called?
//    nf[0] = 0; // set number of function evals. to zero
    double alpha;
    if (first) { // At first iteration, pick initial step according to H.Z. I0
//      alpha = hzStepI0(n,x,f,g,p,angle,fMove,functionEvaluations,optimizationSystem);
      alpha = hzStepI0(f);

      double[] tempPhi = getPhiDPhi(alpha);
      double dphi_0 = tempPhi[1];

      // If dphi_0 is not strictly negative, then do this to force return an alpha where dphi_0 <0.
      if (dphi_0 >= 0) {
        while (dphi_0 >= 0) {
          alpha = 0.5 * alpha;

          tempPhi = getPhiDPhi(alpha);
          dphi_0 = tempPhi[1];
          if (alpha == 0.0) {
            alpha = ulp(1); //2.2204*pow(10,-16);
            return alpha;
          }
        }
      }

      assert( dphi_0 < 0 ); // need to have this before proceeding further
      first = false; // Keep track that the first iteration has been done
    } else {
      // If here, means that this is not the first iteration of the optimization method, so use alpha from previous iteration
      alpha = 1; // TODO How to get/store previous alphas? global var?
      double[] tempPhi = getPhiDPhi(alpha);
      double phi_0 = tempPhi[0];
      double dphi_0 = tempPhi[1];

      alpha = hzStepI12(alpha,phi_0,dphi_0);
      return alpha;
    }
    return alpha; //todo why wouldn't if/else work
  }

  private double hzStepI0(double f) {
    double alpha = 1.0;
    double f_x = f; // initial function value
//    double[] gr = g;
    double gr_max = sqrt(v1DotV2(n, g, 0, 1, g, 0, 1));
    if (gr_max != 0) {
//      double x_max = max(x[0],x[1]); // Infinity norm - should be maximum absolute value (magnitude) of whole vector
      double[] sortX = x;
      Arrays.sort(sortX);
      double x_max = Math.max(abs(sortX[0]), abs(sortX[sortX.length - 1]));
      if (x_max != 0) { // I0.(a)
        alpha = psi0 * x_max / gr_max;
      } else if (f_x != 0) { // point I0.(b)
        alpha = psi0 * abs(f_x) / Math.pow(gr_max,2);
      }
    }
    alpha = min(alpha, alphaMax); // don't let alpha be above set maximum value

    return alpha;
  }

  // Pick the initial step size (Step I1-I2 in HZ paper)
  private double hzStepI12(double alpha, double phi_0, double dphi_0) {
    // Prevent values of 'x_new' that are likely to make phi(x_new) infinite
    int iterFiniteMax = (int) ceil(-log(2, ulp(1)));

    double alphaTest = psi1 * alpha; // This is psi1*alpha_{k-1} in HZ paper.
    alphaTest = min(alphaTest, alphaMax);

    double phiTest = getPhi(alphaTest);

    double[] aPhi = getFinite(alphaTest,phiTest,iterFiniteMax,phi_0);
    alphaTest = aPhi[0];
    phiTest = aPhi[1];

    boolean quadStepSuccess = false; // TODO?

    // I1
    if (QuadStep) { // If QuadStep is true... perform quadratic interpolation
      double a = ( (phiTest - phi_0) / alphaTest - dphi_0 ) / alphaTest; // quadratic fit

      // ... and if phi(psi1*alpha_{k-1}) <= phi(0), and the quadratic interpolant q() is strongly convex (i.e., a > 0), then ...
      if (Double.isFinite(a) && a > 0 && phiTest <= phi_0) {
        // ... choose minimum of quadratic interpolant q()
        double alphaTest2 = -dphi_0 / 2 / a;
        alphaTest2 = min(alphaTest2, alphaMax);

        //
        double phiTest2 = getPhi(alphaTest2);
        if (Double.isFinite(phiTest2)) {
          alphaTest = alphaTest2;
          phiTest = phiTest2;
        }
      }
    }

    // I2
    if ((QuadStep || !quadStepSuccess) && (phiTest <= phi_0)) {
      // If no QuadStep or it fails, expand the interval. While the phiTest <= phi_0 condition was not in the paper, it
      // gives a significant boost to the speed. The rationale behind it is that since the slope at alpha = 0 is
      // negative, if phiTest > phi_0 then a local minimum must be between alpha = 0 and alpha = alphaTest, so alphaTest
      // is good enough to return. todo understand
      alphaTest = psi2 * alpha; // This is c = psi2 * alpha_{k-1} in HZ paper.
      alphaTest = min(alphaTest, alphaMax);

      phiTest = getPhi(alphaTest);
      aPhi = getFinite(alphaTest, phiTest, iterFiniteMax, phi_0); //
      alphaTest = aPhi[0]; // don't use the phi this time
    }
    return alphaTest;
  }

  // Returns alpha and phi TODO: I DON'T GET WHAT THIS IS DOING
  private double[] getFinite(double alpha,double phi,int iterMax,double phi_0) {
    int iter = 1;
    while (Double.isFinite(phi)) {
      if (iter >= iterMax) {
        phi = phi_0;
        return new double[]{alpha, phi}; // return original phi
      }
      alpha = psi3 * alpha;

      phi = getPhi(alpha);

      iter++;
    }
    return null;
  }

  private double[] getPhiDPhi(double c) {
    arraycopy(x0, 0, x, 0, n); // reset X vector before adding c*s
    // TODO - do i need to store a copy of the previous coordinates, like are these just trial moves??
    // Returns phi_c = f(x + c*s) derivative ( dphi_c = f'(x+c*s).s
    // Retraction
    aV1PlusV2(n, c, s, 0, 1, x, 0, 1); // updates x; (x_new = x + c *s)

    // Compute phi_c
    double phi_c = optimizationSystem.energyAndGradient(x, g); // f(x_new)

    // Increase counter for the number of function evaluations
    nf[0]++;

    // Compute dphi_c (already calculated graient above)
//    return v1DotV2(n, s, 0, 1, g, 0, 1); // s . g
    double dphi_c = v1DotV2(n, s, 0, 1, g, 0, 1);
    return new double[]{phi_c, dphi_c};
  }

  private double getPhi(double c) {
    arraycopy(x0, 0, x, 0, n); // reset X vector before adding c*s
    aV1PlusV2(n, c, s, 0, 1, x, 0, 1); // updates x; (x_new = x + c *s)
    double phi_c = optimizationSystem.energyAndGradient(x, g); // f(x_new)
    nf[0]++;
    return phi_c;
  }


  /**
   * Check Wolfe and approximate Wolfe conditions for termination.
   * @param c location between a & b
   * @param phi_c function value at c
   * @param phi_0 initial function value
   * @param phi_lim not sure TODO
   * @return true if both Wolfe OR approximate Wolfe is satisfied.
   */
  private boolean satisfiesWolfe(double c, double phi_c, double dphi_c, double phi_0, double dphi_0, double phi_lim) {
    // Original Wolfe conditions (T1)
    boolean wolfe = (DELTA * dphi_0 >= (phi_c - phi_0) / c) && (dphi_c >= SIGMA * dphi_0);

    // Approximate Wolfe conditions (T2)
    boolean approxWolfe = ( ((2 * DELTA - 1) * dphi_0 >= dphi_c) && (dphi_c >= SIGMA * dphi_0) ) && (phi_c <= phi_lim);

    return (wolfe || approxWolfe);
  }

  // Implementation of [\bar{a},\bar{b}] = secant^{2}(a,b)
  private int[] secant2(int ia, int ib, double phiLim) { // output: [iswolfe, ia, ib] - which to do - need both?
    int isWolfe; // boolean isWolfe; 1 = true, 0 = false
    double phi_0 = valueArr[0];
    double dphi_0 = slopeArr[0];
    double a = alphaArr[ia];
    double b = alphaArr[ib];
    double dphi_a = slopeArr[ia];
    double dphi_b = slopeArr[ib];

    if ( !(dphi_a < 0 && dphi_b >= 0) ) {
      // warning: "The opposite slope condition is not satisfied: dphi(a) = %f; dphi(b) = %f",dphi_a,dphi_b
//      info[0] = LineSearchHZResult.BadSlopes;
      info[0] = LineSearch.LineSearchResult.BadIntpln;
    }

    // S1

    double c = secantFormula(a, b, dphi_a, dphi_b);
    // log: "S1: a = %f, b = %f, c = %f",a,b,c
    assert( Double.isFinite(c) );

    double[] tempPhi = getPhiDPhi(c);
    double phi_c = tempPhi[0];
    double dphi_c = tempPhi[1];

    assert( Double.isFinite(phi_c) && Double.isFinite(dphi_c) );

    alphaArr[alphaIdx] = c;
    valueArr[valueIdx] = phi_c;
    slopeArr[slopeIdx] = dphi_c;
    alphaIdx++;
    valueIdx++;
    slopeIdx++;

    int ic = alphaIdx-1;

    if (satisfiesWolfe(c,phi_c,dphi_c,phi_0,dphi_0,phiLim)) {
      // log: "S1: first c satisfied Wolfe conditions."

      // Save the accepted c for future use. TODO - ??
      // info_global(column_idx) = c;

      isWolfe = 1; // isWolfe = true;
      ia = ic;
      ib = ic;
      return new int[]{isWolfe,ia,ib}; // return true, ic, ic
    }

    // log: "S1: The opposite  slope condition is enforced: dphi(a) = %f; dphi(b) = %f", slopeArr[ia], slopeArr[ib]
    // log: "S1: Values of phi: phi(a) = %f; phi(b) = %f",valueArr[ia],valueArr[ib]
    // log: "S1: update_hz is called by secant2"

    int[] iAB = updateHZ(ia,ib,ic,phiLim);
    int iA = iAB[0];
    int iB = iAB[1];
    // log: "S1: iA = %d, iB = %d, ic = %d",iA,iB,ic
    a = alphaArr[iA];
    b = alphaArr[iB];

    // S2
    if (iB == ic) {
      // updated b, make sure also update a
      c = secant(ib,iB);
    } else if (iA == ic) {
      // S3
      // updated a, do it for b too
      c = secant(ia,iA);
    }

    // S4
    if ( (iA == ic || iB == ic) && ( (a <= c) && (c <= b) ) ) {
      // log: "S4: second c = %f",c

      tempPhi = getPhiDPhi(c);
      phi_c = tempPhi[0];
      dphi_c = tempPhi[1];

      assert( Double.isFinite(phi_c) && Double.isFinite(dphi_c) );

      alphaArr[alphaIdx] = c;
      valueArr[valueIdx] = phi_c;
      slopeArr[slopeIdx] = dphi_c;
      alphaIdx++;
      valueIdx++;
      slopeIdx++;

      ic = alphaIdx-1; // length(alphas) todo - check

      if (satisfiesWolfe(c,phi_c,dphi_c,phi_0,dphi_0,phiLim)) {
        // log: "S4: second c satisfied Wolfe conditions."

        // Save the accepted c for future use. TODO ??
        // info_global(current_lev,column_dix) = c;

        isWolfe = 1; // isWolfe = true;
        ia = ic;
        ib = ic;
        return new int[]{isWolfe,ia,ib}; // return true, ic, ic
      }

      // log: "S4: updateHZ is called by secant2"
      iAB = updateHZ(iA,iB,ic,phiLim);
      iA = iAB[0];
      iB = iAB[1];
    }
    // log: "secant2 output: a = %f, b = %f.",alphaArr[iA],alphaArr[iB]
    // log2:"                dphi(a) = %f; dphi(b) = %f.",slopeArr[iA],slopeArr[iB]
    isWolfe = 0; // isWolfe = false;
    ia = iA;
    ib = iB;
    return new int[]{isWolfe,ia,ib}; // return false, iA, iB
  }

  // call secantFormula with given indices
  private double secant(int ia, int ib) {
    return secantFormula(alphaArr[ia], alphaArr[ib], slopeArr[ia], slopeArr[ib]);
  }

  // Secant formula (p. 123 of HZ)
  private double secantFormula(double a, double b, double dphi_a, double dphi_b) {
    return (a * dphi_b - b * dphi_a) / (dphi_b - dphi_a);
  }

  // Implements [ \bar{a}, \bar{b} ] = update(a,b,c); see Hager and Zhang, stages U0-U3, p. 123
  // Given a third point, pick the best two that retain the bracket around the minimum (as defined by HZ eq. 29)
  // (opposite slope condition). b will be the upper bound, and a the lower bound.
  private int[] updateHZ(int ia, int ib, int ic, double phiLim) {
    double a = alphaArr[ia];
    double b = alphaArr[ib];
    assert( valueArr[ia] <= phiLim );
    assert( slopeArr[ib] >= 0 );
    assert( b > a );
    double c = alphaArr[ic];
    double phi_c = valueArr[ic];
    double dphi_c = slopeArr[ic];
    // log: "update: ia = %d, a = %f, ib = %d, b = %f, c = %f, phi_c = %f, dphi_c = %f",ia,a,ib,b,c,phi_c,dphi_c

    // U0
    if (c < a || c > b ) {
      // log: "U0: c is out of the bracketing interval: a = %f, b = %f, c = %f",a,b,c
      return new int[]{ia, ib}; // keep the same a and b
    }

    // U1
    if (dphi_c >= 0) {
      ib = ic; // replace b with a closer point
      return new int[]{ia, ib};
    }

    // U2
    // NB: When in U2 and U3, the opposite slope condition is naturally NOT satisfied.
    // Know dphi_c < 0. However, phi may not be monotonic between a and c, so check that the value is also smaller than
    // phi_0. (It's more dangerous to replace a than b, since we're leaving the secure environment of alhpa=0; that's
    // why didn't check this above.)
    if (phi_c <= phiLim) {
      ia = ic; // replace a
      return new int[]{ia,ib};
    } else {
      // U3
      // phi_c is bigger than phi_0, which implies that the minimum lies between a and c. Find it via bisection.

      // Debugging
      assert( dphi_c < 0 );
      assert( phi_c > phiLim );
      assert( alphaArr[ia] < alphaArr[ic] );
      assert( valueArr[ic] > valueArr[ia] );

      return bisect(ia,ic,phiLim); // ,x,s); // return new ia and ib
    }
  }

  /**
   * The six possible line search results:
   * <p>
   * Success, WideAngle, ScaleStep, IntplnErr, ReSearch, BadIntpln
   */
  public enum LineSearchHZResult {
    Success,
    NoConvergence, // failed to converge before maximum iterations allowed
    AngleWarn, // Search direction is not a direction of descent; dphi_0 >= eps*abs(phi_0)
    IntplnErr,
    StartErr, // Value and slope at step length = 0 must be finite
    FailNewEvalPt, // Failed to achieve finite new evaluation point, using alpha = 0.
    Flat, // inside bracket is flat
    BadSlopes,
  }
}