/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
 */
package ffx.xray;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

/**
 * <p>ScaleBulkMinimize class.</p>
 *
 * @author Timothy Fenn
 *
 */
public class ScaleBulkMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(ScaleBulkMinimize.class.getName());
    private static double toSeconds = 1.0e-9;
    private static final double eightpi2 = 8.0 * Math.PI * Math.PI;
    private final ReflectionList reflectionlist;
    private final DiffractionRefinementData refinementData;
    private final Crystal crystal;
    private final CrystalReciprocalSpace crs;
    private final ScaleBulkEnergy bulkSolventEnergy;
    private final int solvent_n;
    private final int n;
    private final double x[];
    private final double grad[];
    private final double scaling[];
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    /**
     * <p>Constructor for ScaleBulkMinimize.</p>
     *
     * @param reflectionlist a {@link ffx.crystal.ReflectionList} object.
     * @param refinementdata a {@link ffx.xray.DiffractionRefinementData}
     * object.
     * @param crs a {@link ffx.xray.CrystalReciprocalSpace} object.
     */
    public ScaleBulkMinimize(ReflectionList reflectionlist,
            DiffractionRefinementData refinementdata, CrystalReciprocalSpace crs) {
        this.reflectionlist = reflectionlist;
        this.refinementData = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.crs = crs;

        if (crs.solventmodel == SolventModel.NONE) {
            solvent_n = 1;
        } else {
            solvent_n = 3;
        }
        n = solvent_n + refinementdata.scale_n;
        bulkSolventEnergy = new ScaleBulkEnergy(reflectionlist, refinementdata, n);

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        x[0] = refinementdata.model_k;
        if (solvent_n > 1) {
            x[1] = refinementdata.solvent_k;
            x[2] = refinementdata.solvent_ueq;
        }
        for (int i = 0; i < 6; i++) {
            if (crystal.scale_b[i] >= 0) {
                x[solvent_n + crystal.scale_b[i]] =
                        refinementdata.model_b[i];
            }
        }

        for (int i = 0; i < n; i++) {
            scaling[i] = 1.0;
        }

        bulkSolventEnergy.setScaling(scaling);

        setInitialScale();
    }

    private void setInitialScale() {
        double fc[][] = refinementData.fc;
        double fo[][] = refinementData.fsigf;

        double e = bulkSolventEnergy.energyAndGradient(x, grad);
        double fct, sumfofc, sumfc;
        sumfofc = sumfc = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])
                    || fo[i][1] <= 0.0) {
                continue;
            }

            fct = refinementData.fctot_f(i);
            sumfofc += fo[i][0] * fct;
            sumfc += fct * fct;
        }

        x[0] = Math.log(4.0 * sumfofc / sumfc);
    }

    /**
     * <p>ksbsGridOptimize</p>
     */
    public void ksbsGridOptimize() {
        if (solvent_n < 3) {
            return;
        }

        double min = Double.POSITIVE_INFINITY;
        double k = refinementData.solvent_k;
        double kmin = 0.05;
        double kmax = 0.9;
        double kstep = 0.05;
        double b = refinementData.solvent_ueq;
        double bmin = 10.0;
        double bmax = 200.0;
        double bstep = 10.0;
        for (double i = kmin; i <= kmax; i += kstep) {
            for (double j = bmin; j <= bmax; j += bstep) {

                x[1] = i;
                x[2] = j;
                double sum = bulkSolventEnergy.energyAndGradient(x, grad);

                System.out.println("ks: " + i + " bs: " + j + " sum: " + sum);
                if (sum < min) {
                    min = sum;
                    k = i;
                    b = j;
                }
            }
        }
        System.out.println("minks: " + k + " minbs: " + b + " min: " + min);
        refinementData.solvent_k = k;
        refinementData.solvent_ueq = b;
    }

    /**
     * <p>GridOptimize</p>
     */
    public void GridOptimize() {
        if (crs == null) {
            return;
        }

        double min = Double.POSITIVE_INFINITY;
        double a = crs.solvent_a;
        double b = crs.solvent_b;
        double amin = a - 1.0;
        double amax = (a + 1.0) / 0.9999;
        double astep = 0.25;
        double bmin = b - 0.2;
        double bmax = (b + 0.2) / 0.9999;
        double bstep = 0.05;
        if (crs.solventmodel == SolventModel.BINARY) {
            amin = a - 0.2;
            amax = (a + 0.2) / 0.9999;
            astep = 0.05;
        }
        for (double i = amin; i <= amax; i += astep) {
            for (double j = bmin; j <= bmax; j += bstep) {
                crs.setSolventAB(i, j);

                crs.computeDensity(refinementData.fs);
                double sum = bulkSolventEnergy.energyAndGradient(x, grad);

                System.out.println("a: " + i + " b: " + j + " sum: " + sum);
                if (sum < min) {
                    min = sum;
                    a = i;
                    b = j;
                }
            }
        }
        System.out.println("mina: " + a + " minb: " + b + " min: " + min);
        crs.setSolventAB(a, b);
        refinementData.solvent_a = a;
        refinementData.solvent_b = b;
        crs.computeDensity(refinementData.fs);
    }

    /**
     * <p>minimize</p>
     *
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    public ScaleBulkEnergy minimize() {
        return minimize(0.5);
    }

    /**
     * <p>minimize</p>
     *
     * @param eps a double.
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    public ScaleBulkEnergy minimize(double eps) {
        return minimize(5, eps);
    }

    /**
     * <p>minimize</p>
     *
     * @param m a int.
     * @param eps a double.
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    public ScaleBulkEnergy minimize(int m, double eps) {

        double e = bulkSolventEnergy.energyAndGradient(x, grad);

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, bulkSolventEnergy, this);
        done = true;
        switch (status) {
            case 0:
                logger.info(String.format("\n Optimization achieved convergence criteria: %10.5f\n", grms));
                break;
            case 1:
                logger.info(String.format("\n Optimization terminated at step %d.\n", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.\n");
        }
        refinementData.model_k = x[0] / scaling[0];
        if (solvent_n > 1) {
            refinementData.solvent_k = x[1] / scaling[1];
            refinementData.solvent_ueq = x[2] / scaling[2];

            if (crs != null) {
                refinementData.solvent_a = crs.solvent_a;
                refinementData.solvent_b = crs.solvent_b;
            }
        }
        for (int i = 0; i < 6; i++) {
            if (crystal.scale_b[i] >= 0) {
                refinementData.model_b[i] =
                        x[solvent_n + crystal.scale_b[i]]
                        / scaling[solvent_n + crystal.scale_b[i]];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            mtime += System.nanoTime();
            logger.info(String.format(" Optimization time: %8.3f (sec)\n", mtime * toSeconds));
        }

        return bulkSolventEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df, double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization of Fc to Fo Scaling Parameters\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time");
        }
        if (info == null) {
            logger.info(String.format("%6d %12.5f %10.6f",
                    iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(String.format("%6d %12.5f %10.6f %10.6f %9.5f %8.2f %6d %8.3f",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
            } else {
                logger.info(String.format("%6d %12.5f %10.6f %10.6f %9.5f %8.2f %6d %8s",
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
}
