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

/**
 *
 * @author fennt
 */
public class ScaleBulkMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());
    private static double toSeconds = 0.000000001;
    private static final double eightpi2 = 8.0 * Math.PI * Math.PI;
    private final ReflectionList reflectionlist;
    private final RefinementData refinementdata;
    private final Crystal crystal;
    private final CrystalReciprocalSpace crs;
    private final ScaleBulkEnergy bulksolventenergy;
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
        this(reflectionlist, refinementdata, null);
    }

    public ScaleBulkMinimize(ReflectionList reflectionlist,
            RefinementData refinementdata, CrystalReciprocalSpace crs) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;
        this.crs = crs;

        n = refinementdata.solvent_n + refinementdata.scale_n;
        bulksolventenergy = new ScaleBulkEnergy(reflectionlist, refinementdata, n);

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        x[0] = refinementdata.model_k;
        if (refinementdata.solvent_n > 1) {
            x[1] = refinementdata.solvent_k;
            x[2] = refinementdata.solvent_ueq;
        }
        for (int i = 0; i < 6; i++) {
            if (crystal.scale_b[i] >= 0) {
                x[refinementdata.solvent_n + crystal.scale_b[i]] =
                        refinementdata.model_b[i];
            }
        }

        for (int i = 0; i < n; i++) {
            scaling[i] = 1.0;
        }

        bulksolventenergy.setOptimizationScaling(scaling);
    }

    public void ksbsGridOptimize() {
        if (refinementdata.solvent_n < 3) {
            return;
        }

        double min = Double.MAX_VALUE;
        double k = refinementdata.solvent_k;
        double kmin = 0.1;
        double kmax = 0.8;
        double kstep = 0.1;
        double b = refinementdata.solvent_ueq;
        double bmin = 10.0;
        double bmax = 200.0;
        double bstep = 10.0;
        for (double i = kmin; i <= kmax; i += kstep) {
            for (double j = bmin; j <= bmax; j += bstep) {

                x[1] = i;
                x[2] = j;
                double sum = bulksolventenergy.energyAndGradient(x, grad);

                System.out.println("ks: " + i + " bs: " + j + " sum: " + sum);
                if (sum < min) {
                    min = sum;
                    k = i;
                    b = j;
                }
            }
        }
        System.out.println("minks: " + k + " minbs: " + b + " min: " + min);
        refinementdata.solvent_k = k;
        refinementdata.solvent_ueq = b;
    }

    public void GridOptimize() {
        if (crs == null) {
            return;
        }

        if (refinementdata.binarysolvent){
            binaryradGridOptimize();
        } else {
            asdGridOptimize();
        }
    }

    public void asdGridOptimize() {
        if (crs == null) {
            return;
        }

        double min = Double.MAX_VALUE;
        double a = refinementdata.solvent_a;
        double amin = a - 1.0;
        double amax = (a + 1.0) / 0.9999;
        double astep = 0.5;
        double sd = refinementdata.solvent_sd;
        double sdmin = sd - 0.15;
        double sdmax = (sd + 0.15) / 0.9999;
        double sdstep = 0.05;
        for (double i = amin; i <= amax; i += astep) {
            crs.setSolventA(i);
            for (double j = sdmin; j <= sdmax; j += sdstep) {
                crs.setSolventsd(j);

                crs.computeDensity(refinementdata.fs);
                double sum = bulksolventenergy.energyAndGradient(x, grad);

                System.out.println("a: " + i + " sd: " + j + " sum: " + sum);
                if (sum < min) {
                    min = sum;
                    a = i;
                    sd = j;
                }
            }
        }
        System.out.println("mina: " + a + " minsd: " + sd + " min: " + min);
        crs.setSolventA(a);
        crs.setSolventsd(sd);
        refinementdata.solvent_a = a;
        refinementdata.solvent_sd = sd;
        crs.computeDensity(refinementdata.fs);
    }

    public void binaryradGridOptimize() {
        if (crs == null) {
            return;
        }

        double min = Double.MAX_VALUE;
        double r = 2.0;
        double rmin = r - 1.5;
        double rmax = (r + 1.5) / 0.9999;
        double rstep = 0.5;
        for (double i = rmin; i <= rmax; i += rstep) {
            crs.setSolventBinaryrad(i);
            crs.computeDensity(refinementdata.fs);
            double sum = bulksolventenergy.energyAndGradient(x, grad);

            System.out.println("rad: " + i + " sum: " + sum);
            if (sum < min) {
                min = sum;
                r = i;
            }
        }
        System.out.println("minrad: " + r + " min: " + min);
        crs.setSolventBinaryrad(r);
        crs.computeDensity(refinementdata.fs);
    }

    public ScaleBulkEnergy minimize() {
        return minimize(0.5);
    }

    public ScaleBulkEnergy minimize(double eps) {
        return minimize(5, eps);
    }

    public ScaleBulkEnergy minimize(int m, double eps) {

        double e = bulksolventenergy.energyAndGradient(x, grad);

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, bulksolventenergy, this);
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
        refinementdata.model_k = x[0] / scaling[0];
        if (refinementdata.solvent_n > 1) {
            refinementdata.solvent_k = x[1] / scaling[1];
            refinementdata.solvent_ueq = x[2] / scaling[2];
        }
        for (int i = 0; i < 6; i++) {
            if (crystal.scale_b[i] >= 0) {
                refinementdata.model_b[i] =
                        x[refinementdata.solvent_n + crystal.scale_b[i]]
                        / scaling[refinementdata.solvent_n + crystal.scale_b[i]];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuffer sb = new StringBuffer();
            mtime += System.nanoTime();
            sb.append(String.format("minimizer time: %g\n", mtime * toSeconds));
            logger.info(sb.toString());
        }

        return bulksolventenergy;
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
