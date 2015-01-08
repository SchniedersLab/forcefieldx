/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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

import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.XtalEnergy;

/**
 * Minimize the potential energy of a system to an RMS gradient per atom
 * convergence criteria.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class CrystalMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(CrystalMinimize.class.getName());
    private final int n;
    private final double[] x;
    private final double[] grad;
    private final double[] scaling;
    private final MolecularAssembly molecularAssembly;
    private final XtalEnergy potential;
    private AlgorithmListener algorithmListener;
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;
    private Crystal crystal;

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
    public CrystalMinimize(MolecularAssembly molecularAssembly, XtalEnergy potential,
            AlgorithmListener algorithmListener) {
        assert (molecularAssembly != null);
        this.molecularAssembly = molecularAssembly;
        this.algorithmListener = algorithmListener;
        this.potential = potential;
        n = potential.getNumberOfVariables();
        x = new double[n];
        grad = new double[n];
        crystal = molecularAssembly.getCrystal();
        scaling = new double[n];
        for (int i = 0; i < n - 6; i += 3) {
            scaling[i] = 12.0 * crystal.a;
            scaling[i + 1] = 12.0 * crystal.b;
            scaling[i + 2] = 12.0 * crystal.c;
        }
        scaling[n - 6] = 4.0 * sqrt(crystal.a);
        scaling[n - 5] = 4.0 * sqrt(crystal.b);
        scaling[n - 4] = 4.0 * sqrt(crystal.c);
        scaling[n - 3] = 0.02 * sqrt(crystal.alpha);
        scaling[n - 2] = 0.02 * sqrt(crystal.beta);
        scaling[n - 1] = 0.02 * sqrt(crystal.gamma);

        potential.setScaling(scaling);
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
    public CrystalMinimize(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener) {
        assert (molecularAssembly != null);
        this.molecularAssembly = molecularAssembly;
        this.algorithmListener = algorithmListener;
        if (molecularAssembly.getPotentialEnergy() == null) {
            molecularAssembly.setPotential(new ForceFieldEnergy(molecularAssembly));
        }
        potential = new XtalEnergy(molecularAssembly.getPotentialEnergy(), molecularAssembly);
        n = potential.getNumberOfVariables();
        x = new double[n];
        grad = new double[n];
        crystal = molecularAssembly.getCrystal();
        scaling = new double[n];
        for (int i = 0; i < n - 6; i += 3) {
            scaling[i] = 12.0 * crystal.a;
            scaling[i + 1] = 12.0 * crystal.b;
            scaling[i + 2] = 12.0 * crystal.c;
        }
        scaling[n - 6] = 4.0 * sqrt(crystal.a);
        scaling[n - 5] = 4.0 * sqrt(crystal.b);
        scaling[n - 4] = 4.0 * sqrt(crystal.c);
        scaling[n - 3] = 0.02 * sqrt(crystal.alpha);
        scaling[n - 2] = 0.02 * sqrt(crystal.beta);
        scaling[n - 1] = 0.02 * sqrt(crystal.gamma);

        potential.setScaling(scaling);
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
        return minimize(1.0);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps a double.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps) {
        return minimize(7, eps);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param m a int.
     * @param eps a double.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(int m, double eps) {
        time = System.nanoTime();
        potential.getCoordinates(x);
        /**
         * Scale coordinates.
         */
        for (int i = 0; i < n; i++) {
            x[i] *= scaling[i];
        }

        done = false;
        int status = 2;
        double e = potential.energyAndGradient(x, grad);
        status = LBFGS.minimize(n, m, x, e, grad, eps, potential, this);
        done = true;

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
        crystal = molecularAssembly.getCrystal();
        logger.info(String.format("\n Final lattice parameters" + crystal));

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
     * @since 1.0
     * @return a boolean.
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
