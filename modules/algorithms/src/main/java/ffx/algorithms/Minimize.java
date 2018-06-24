/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.algorithms;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;

/**
 * Minimize the potential energy of a system to an RMS gradient per atom
 * convergence criteria.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class Minimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(Minimize.class.getName());
    private final int n;
    private final double[] x;
    private final double[] grad;
    private final double[] scaling;
    private final MolecularAssembly molecularAssembly;
    private final Potential potential;
    private final AlgorithmListener algorithmListener;
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double energy;
    private double grms;
    private int status;
    private int nSteps;

    /**
     * <p>
     * Constructor for Minimize.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param potential a {@link ffx.numerics.Potential} object.
     * @param algorithmListener a {@link ffx.algorithms.AlgorithmListener}
     * object.
     */
    public Minimize(MolecularAssembly molecularAssembly, Potential potential,
                    AlgorithmListener algorithmListener) {
        assert (molecularAssembly != null);
        this.molecularAssembly = molecularAssembly;
        this.algorithmListener = algorithmListener;
        this.potential = potential;
        n = potential.getNumberOfVariables();
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        Arrays.fill(scaling, 12.0);
    }

    /**
     * <p>
     * Constructor for Minimize.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param algorithmListener a {@link ffx.algorithms.AlgorithmListener}
     * object.
     */
    public Minimize(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener) {
        assert (molecularAssembly != null);
        this.molecularAssembly = molecularAssembly;
        this.algorithmListener = algorithmListener;
        if (molecularAssembly.getPotentialEnergy() == null) {
            molecularAssembly.setPotential(ForceFieldEnergy.energyFactory(molecularAssembly));
        }
        potential = molecularAssembly.getPotentialEnergy();
        n = potential.getNumberOfVariables();
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        Arrays.fill(scaling, 12.0);
    }

    public double getGRMS() {
        return grms;
    }

    public int getStatus() {
        return status;
    }

    public double getEnergy() {
        return energy;
    }

    /**
     * {@inheritDoc}
     */
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

    /**
     * <p>
     * minimize</p>
     *
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize() {
        return minimize(7,1.0, Integer.MAX_VALUE);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps The convergence criteria.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps) {
        return minimize(7, eps, Integer.MAX_VALUE);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps The convergence criteria.
     * @param maxIterations The maximum number of iterations.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps, int maxIterations) {
        return minimize(7, eps, maxIterations);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param m The number of previous steps used to estimate the Hessian.
     * @param eps The convergence criteria.
     * @param maxIterations The maximum number of iterations.
     *
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(int m, double eps, int maxIterations) {
        time = System.nanoTime();
        potential.getCoordinates(x);
        potential.setScaling(scaling);

        /**
         * Scale coordinates.
         */
        for (int i = 0; i < n; i++) {
            x[i] *= scaling[i];
        }

        done = false;
        energy = potential.energyAndGradient(x, grad);

        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Minimize initial energy: %16.8f", energy));
        }

        status = LBFGS.minimize(n, m, x, energy, grad, eps, maxIterations, potential, this);
        done = true;

        switch (status) {
            case 0:
                logger.info(format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
                break;
            case 1:
                logger.info(format("\n Optimization terminated at step %d.\n", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.\n");
        }

        potential.setScaling(null);
        return potential;
    }

    /**
     * {@inheritDoc}
     *
     * Implement the OptimizationListener interface.
     *
     * @since 1.0
     */
    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df,
                                      double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;
        this.energy = f;

        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time\n");
        }
        if (info == null) {
            logger.info(format("%6d%13.4f%11.4f", iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8.3f",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
            } else {
                logger.info(format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8s",
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
     * @since 1.0
     * @return a boolean.
     */
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;
        this.energy = f;

        if (iter == 1) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Evals     Time\n");
        }
        logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%7d %8.3f",
                iter, f, grms, df, xrms, nfun, seconds));
        // Update the listener and check for a termination request.
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
}
