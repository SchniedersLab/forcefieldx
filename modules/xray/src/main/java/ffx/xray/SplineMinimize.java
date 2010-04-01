/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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
package ffx.xray;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.xray.SplineEnergy.Type;

/**
 *
 * @author fennt
 */
public class SplineMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final Crystal crystal;
    private final SplineEnergy splineenergy;
    private final int n;
    private final double x[];
    private final double grad[];
    private final double scaling[];
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    public SplineMinimize(ReflectionList reflectionlist,
            RefinementData refinementdata, double x[], int type) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.x = x;

        n = x.length;
        splineenergy = new SplineEnergy(reflectionlist, refinementdata,
                n, type);
        grad = new double[n];
        scaling = new double[n];
        for (int i = 0; i < n; i++) {
            if (type == Type.FOTOESQ
                    || type == Type.FCTOESQ) {
                x[i] = 0.1;
            } else {
                x[i] = 1.0;
            }
            scaling[i] = 1.0;
        }
        splineenergy.setOptimizationScaling(scaling);
    }

    public SplineEnergy minimize() {
        return minimize(0.5);
    }

    public SplineEnergy minimize(double eps) {
        return minimize(5, eps);
    }

    public SplineEnergy minimize(int m, double eps) {

        double e = splineenergy.energyAndGradient(x, grad);

        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, splineenergy, this);
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
        return splineenergy;
    }

    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df, double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        /*
        if (iter == 0) {
        logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
        logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time\n");
        }
        if (info == null) {
        logger.info(String.format("%6d %13.4g %11.4g\n",
        iter, f, grms));
        } else {
        if (info == LineSearchResult.Success) {
        logger.info(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8.3g\n",
        iter, f, grms, df, xrms, angle, nfun, seconds));
        } else {
        logger.info(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8s\n",
        iter, f, grms, df, xrms, angle, nfun, info.toString()));
        }
        }
         */

        if (terminate) {
            logger.info(" The optimization recieved a termination request.");
            // Tell the L-BFGS optimizer to terminate.
            return false;
        }
        return true;
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
}
