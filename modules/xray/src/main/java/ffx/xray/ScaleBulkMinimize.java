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

import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author fennt
 */
public class ScaleBulkMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(SplineOptimizer.class.getName());
    private static final double eightpi2 = 8.0 * Math.PI * Math.PI;
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final Crystal crystal;
    private final ScaleBulkOptimizer bulksolventoptimizer;
    private final int n;
    private final double x[];
    private final double grad[];
    private final double scaling[];
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    public ScaleBulkMinimize(ReflectionList reflectionlist,
            RefinementData refinementdata) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;

        n = 9;
        bulksolventoptimizer = new ScaleBulkOptimizer(reflectionlist, refinementdata);
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        x[0] = refinementdata.model_k;
        x[7] = refinementdata.solvent_k;
        x[8] = refinementdata.solvent_ueq;
        for (int i = 0; i < n; i++) {
            scaling[i] = 1.0;
        }
        bulksolventoptimizer.setOptimizationScaling(scaling);


        int scale_flag = crystal.scale_flag;

        if ((scale_flag & Crystal.SCALE_B11) == Crystal.SCALE_B11) {
            System.out.println("scale B11");
        }
        if ((scale_flag & Crystal.SCALE_B22) == Crystal.SCALE_B22) {
            System.out.println("scale B22");
        }
        if ((scale_flag & Crystal.SCALE_B33) == Crystal.SCALE_B33) {
            System.out.println("scale B33");
        }
        if ((scale_flag & Crystal.SCALE_B12) == Crystal.SCALE_B12) {
            System.out.println("scale B12");
        }
        if ((scale_flag & Crystal.SCALE_B13) == Crystal.SCALE_B13) {
            System.out.println("scale B13");
        }
        if ((scale_flag & Crystal.SCALE_B23) == Crystal.SCALE_B23) {
            System.out.println("scale B23");
        }
    }

    public ScaleBulkOptimizer minimize() {
        return minimize(0.5);
    }

    public ScaleBulkOptimizer minimize(double eps) {
        return minimize(5, eps);
    }

    public ScaleBulkOptimizer minimize(int m, double eps) {

        double e = bulksolventoptimizer.energyAndGradient(x, grad);

        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, bulksolventoptimizer, this);
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
        refinementdata.model_k = x[0];
        refinementdata.aniso_b[0] = x[1];
        refinementdata.aniso_b[1] = x[2];
        refinementdata.aniso_b[2] = x[3];
        refinementdata.aniso_b[3] = x[4];
        refinementdata.aniso_b[4] = x[5];
        refinementdata.aniso_b[5] = x[6];
        refinementdata.solvent_k = x[7];
        refinementdata.solvent_ueq = x[8];

        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            sb.append(String.format("\n final scale:\n"));
            sb.append(String.format("  overall scale: %g\n", x[0]));
            sb.append(String.format("  aniso B tensor:\n"));
            sb.append(String.format("    %g %g %g\n", x[1], x[4], x[5]));
            sb.append(String.format("    %g %g %g\n", x[4], x[2], x[6]));
            sb.append(String.format("    %g %g %g\n", x[5], x[6], x[3]));
            sb.append(String.format("  bulk solvent scale: %g  B: %g\n\n",
                    x[7], x[8] * eightpi2));
            logger.info(sb.toString());
        }

        return bulksolventoptimizer;
    }

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
