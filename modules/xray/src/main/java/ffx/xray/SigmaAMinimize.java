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

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.ComplexNumber;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;

/**
 *
 * @author fennt
 */
public class SigmaAMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(SplineEnergy.class.getName());
    private static double toSeconds = 0.000000001;
    private final ReflectionList reflectionlist;
    protected final RefinementData refinementdata;
    private final Crystal crystal;
    private final SigmaAEnergy sigmaaenergy;
    private final int n;
    private final double x[];
    private final double grad[];
    private final double scaling[];
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    public SigmaAMinimize(ReflectionList reflectionlist,
            RefinementData refinementdata) {
        this.reflectionlist = reflectionlist;
        this.refinementdata = refinementdata;
        this.crystal = reflectionlist.crystal;

        n = refinementdata.nparams * 2;
        sigmaaenergy = new SigmaAEnergy(reflectionlist, refinementdata);
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];

        for (int i = 0; i < refinementdata.nparams; i++) {
            // for optimizationscaling, best to move to 0.0
            x[i] = refinementdata.sigmaa[i] - 1.0;
            scaling[i] = 1.0;
            x[i + refinementdata.nparams] = refinementdata.sigmaw[i];
            scaling[i + refinementdata.nparams] = 2.0;
        }

        // generate Es
        int type = SplineEnergy.Type.FCTOESQ;
        SplineMinimize splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.fcesq, type);
        splineminimize.minimize(7, 1.0);

        type = SplineEnergy.Type.FOTOESQ;
        splineminimize = new SplineMinimize(reflectionlist,
                refinementdata, refinementdata.foesq, type);
        splineminimize.minimize(7, 1.0);

        // generate initial w estimate
        ReflectionSpline spline = new ReflectionSpline(reflectionlist,
                refinementdata.nparams);
        int nmean[] = new int[refinementdata.nparams];
        for (int i = 0; i < refinementdata.nparams; i++) {
            nmean[i] = 0;
        }
        double mean = 0.0, tot = 0.0;
        double fc[][] = refinementdata.fctot;
        double fo[][] = refinementdata.fsigf;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (ih.allowed() == 0.0
                    || Double.isNaN(fc[i][0])
                    || Double.isNaN(fo[i][0])) {
                continue;
            }

            double s2 = Crystal.invressq(crystal, ih);
            double epsc = ih.epsilonc();
            ComplexNumber fct = new ComplexNumber(fc[i][0], fc[i][1]);
            double ecscale = spline.f(s2, refinementdata.fcesq);
            double eoscale = spline.f(s2, refinementdata.foesq);
            double ec = fct.times(sqrt(ecscale)).abs();
            double eo = fo[i][0] * sqrt(eoscale);
            double wi = pow(eo - ec, 2.0) / epsc;

            nmean[spline.i1()]++;
            tot++;

            x[spline.i1() + refinementdata.nparams] += (wi
                    - x[spline.i1() + refinementdata.nparams])
                    / nmean[spline.i1()];
            mean += (wi - mean) / tot;
        }

        logger.info("starting mean w: " + mean + " w scaling: " + 1.0 / mean);
        for (int i = 0; i < refinementdata.nparams; i++) {
            x[i] -= x[i + refinementdata.nparams];
            x[i] *= scaling[i];
            scaling[i + refinementdata.nparams] = 1.0 / mean;
            x[i + refinementdata.nparams] *= scaling[i + refinementdata.nparams];
        }

        sigmaaenergy.setScaling(scaling);
    }

    public double calculateLikelihood() {
        sigmaaenergy.energyAndGradient(x, grad);
        return refinementdata.llkr;
    }

    public double calculateLikelihoodFree() {
        return sigmaaenergy.energyAndGradient(x, grad);
    }

    public SigmaAEnergy minimize() {
        return minimize(0.5);
    }

    public SigmaAEnergy minimize(double eps) {
        return minimize(7, eps);
    }

    public SigmaAEnergy minimize(int m, double eps) {

        double e = sigmaaenergy.energyAndGradient(x, grad);

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, sigmaaenergy, this);
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

        for (int i = 0; i < refinementdata.nparams; i++) {
            refinementdata.sigmaa[i] = 1.0 + x[i] / scaling[i];
            refinementdata.sigmaw[i] = x[i + refinementdata.nparams]
                    / scaling[i + refinementdata.nparams];
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            mtime += System.nanoTime();
            sb.append(String.format("minimizer time: %g\n", mtime * toSeconds));
            logger.info(sb.toString());
        }

        return sigmaaenergy;
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
