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
package ffx.algorithms;

import static java.lang.Math.sqrt;

import java.util.logging.Logger;

import java.util.logging.Level;

import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealConvergenceChecker;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.general.ConjugateGradientFormula;
import org.apache.commons.math.optimization.general.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math.analysis.solvers.BrentSolver;
import org.apache.commons.math.analysis.solvers.NewtonSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactoryImpl;


import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.potential.PotentialEnergy;
import ffx.potential.bonded.MolecularAssembly;

/**
 * Minimize the potential energy of a system to an RMS gradient per atom
 * convergence criteria.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class Minimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(Minimize.class.getName());
    private final int n;
    private final double[] x;
    private final double[] grad;
    private final double[] scaling;
    private final MolecularAssembly molecularAssembly;
    private final PotentialEnergy potentialEnergy;
    private AlgorithmListener algorithmListener;
    private NonLinearConjugateGradientOptimizer optimizer;
    private ConvergenceChecker convergenceChecker;
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    public Minimize(MolecularAssembly assembly, AlgorithmListener listener) {
        assert (assembly != null);
        molecularAssembly = assembly;
        algorithmListener = listener;
        n = molecularAssembly.getAtomList().size() * 3;
        if (molecularAssembly.getPotentialEnergy() == null) {
            potentialEnergy = new PotentialEnergy(molecularAssembly);
            molecularAssembly.setPotential(potentialEnergy);
        } else {
            potentialEnergy = molecularAssembly.getPotentialEnergy();
        }
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        for (int i = 0; i < n; i++) {
            scaling[i] = 12.0;
        }
        potentialEnergy.setOptimizationScaling(scaling);
    }

    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Exception terminating minimization.\n", e);
                }
            }
        }
    }

    public PotentialEnergy minimize() {
        return minimize(1.0);
    }

    public PotentialEnergy minimize(double eps) {
        return minimize(7, eps);
    }

    public PotentialEnergy minimize(int m, double eps) {

        time = System.nanoTime();
        potentialEnergy.getCoordinates(x);
        /**
         * Scale coordinates.
         */
        for (int i = 0; i < n; i++) {
            x[i] *= scaling[i];
        }
        done = false;
        int status = 2;
        double e = potentialEnergy.energyAndGradient(x, grad);
        status = LBFGS.minimize(n, m, x, e, grad, eps, potentialEnergy, this);
        done = true;
        

        /*
        optimizer = new NonLinearConjugateGradientOptimizer(ConjugateGradientFormula.FLETCHER_REEVES);
        convergenceChecker = new ConvergenceChecker(eps);
        optimizer.setConvergenceChecker(convergenceChecker);
        optimizer.setMaxEvaluations(Integer.MAX_VALUE);
        optimizer.setMaxIterations(Integer.MAX_VALUE);
        UnivariateRealSolver solver = UnivariateRealSolverFactoryImpl.newInstance().newBrentSolver();
        solver.setMaximalIterationCount(Integer.MAX_VALUE);
        solver.setFunctionValueAccuracy(1.0e-16);
        optimizer.setLineSearchSolver(solver);

        try {
            optimizer.optimize(potentialEnergy, GoalType.MINIMIZE, x);
            grms = convergenceChecker.getGRMS();
            done = true;
            if (grms < eps) {
                status = 0;
            } else {
                status = 1;
                nSteps = optimizer.getIterations();
            }
        } catch (Exception e) {
            String message = "Exception during optimization optimization.";
            logger.log(Level.WARNING, message, e);
        } */

        switch (status) {
            case 0:
                logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
                break;
            case 1:
                logger.info(String.format("\n Optimization terminated at step %d.\n", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.\n");
        }
        return potentialEnergy;
    }

    /**
     * Implement the OptimizationListener interface.
     *
     * @param iter Number of iterations so far.
     * @param nfun Number of function evaluations so far.
     * @param grms Gradient RMS at current solution.
     * @param xrms Coordinate change RMS at current solution.
     * @param f Function value at current solution.
     * @param df Change in the function value compared to the previous solution.
     * @param angle Current angle between gradient and search direction.
     *
     * @since 1.0
     */
    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df, double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time\n");
        }
        if (info == null) {
            logger.info(String.format("%6d%13.4f%11.4f", iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8.3f",
                                          iter, f, grms, df, xrms, angle, nfun, seconds));
            } else {
                logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8s",
                                          iter, f, grms, df, xrms, angle, nfun, info.toString()));
            }
        }
        // Update the listener and check for an termination request.
        if (algorithmListener != null) {
            algorithmListener.algorithmUpdate(molecularAssembly);
        }
        if (terminate) {
            logger.info(" The optimization recieved a termination request.");
            // Tell the L-BFGS optimizer to terminate.
            return false;
        }
        return true;
    }

    /**
     * Implement the OptimizationListener interface.
     *
     * @param iter Number of iterations so far.
     * @param nfun Number of function evaluations so far.
     * @param grms Gradient RMS at current solution.
     * @param xrms Coordinate change RMS at current solution.
     * @param f Function value at current solution.
     * @param df Change in the function value compared to the previous solution.
     *
     * @since 1.0
     */
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        if (iter == 1) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Evals     Time\n");
        }
        logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%7d %8.3f",
                                  iter, f, grms, df, xrms, nfun, seconds));
        // Update the listener and check for an termination request.
        if (algorithmListener != null) {
            algorithmListener.algorithmUpdate(molecularAssembly);
        }
        if (terminate) {
            logger.info(" The optimization recieved a termination request.");
            // Tell the optimizer to terminate.
            return false;
        }
        return true;
    }

    private class ConvergenceChecker implements RealConvergenceChecker {

        double rms = Math.sqrt(n) / Math.sqrt(3.0);
        double g[] = new double[n];
        final double eps;
        private double grms;

        public ConvergenceChecker(double eps) {
            this.eps = eps;
        }

        public double getGRMS() {
            return grms;
        }

        @Override
        public boolean converged(int k, RealPointValuePair previous, RealPointValuePair current) {
            double f = current.getValue();
            double df = f - previous.getValue();
            /**
             * Sometimes the NonLinearOptimizer calls this with a positive step?
             */
            if (df > 0) {
                return false;
            }

            int iter = k;
            int nfun = optimizer.getGradientEvaluations();

            /**
             * Compute the RMS gradient per atom.
             * Compute the RMS coordinate change per atom.
             */
            grms = 0.0;
            double xrms = 0.0;
            double prevX[] = previous.getPointRef();
            double currentX[] = current.getPointRef();
            potentialEnergy.getGradients(g);
            for (int i = 0; i < n; i++) {
                double gs = g[i];
                double xs = currentX[i] - prevX[i];
                grms += gs * gs;
                xrms += xs * xs;
            }
            grms = sqrt(grms) / rms;
            xrms = sqrt(xrms) / rms;
            /**
             * Log this step of the optimization and see if an interrupt has
             * been requested.
             */
            boolean status = optimizationUpdate(iter, nfun, grms, xrms, f, df);
            /**
             * Return true if we have converged to the RMS gradient criteria
             * or if an interrupt has been requested.
             */
            if (!status || grms < eps) {
                return true;
            }
            /**
             * Otherwise, continue the optimization. 
             */
            return false;
        }
    }
}
