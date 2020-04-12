//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
//******************************************************************************
package ffx.xray;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.ParallelTeam;

import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.crystal.ReflectionSpline;
import ffx.numerics.math.ComplexNumber;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch.LineSearchResult;
import ffx.numerics.optimization.OptimizationListener;

/**
 * <p>
 * SigmaAMinimize class.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class SigmaAMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(SigmaAMinimize.class.getName());
    private static final double toSeconds = 1.0e-9;
    protected final DiffractionRefinementData refinementData;
    private final ReflectionList reflectionList;
    private final Crystal crystal;
    private final SigmaAEnergy sigmaAEnergy;
    private final int n;
    private final double[] x;
    private final double[] grad;
    private final double[] scaling;
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;

    /**
     * <p>
     * Constructor for SigmaAMinimize.</p>
     *
     * @param reflectionList a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData a {@link ffx.xray.DiffractionRefinementData} object.
     * @param parallelTeam   the ParallelTeam to execute the SigmaAMinimize.
     */
    SigmaAMinimize(ReflectionList reflectionList, DiffractionRefinementData refinementData, ParallelTeam parallelTeam) {
        this.reflectionList = reflectionList;
        this.refinementData = refinementData;
        this.crystal = reflectionList.crystal;

        n = refinementData.nBins * 2;
        sigmaAEnergy = new SigmaAEnergy(reflectionList, refinementData, parallelTeam);
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];

        for (int i = 0; i < refinementData.nBins; i++) {
            // For optimization scaling, best to move to 0.0.
            x[i] = refinementData.sigmaA[i] - 1.0;
            scaling[i] = 1.0;
            x[i + refinementData.nBins] = refinementData.sigmaW[i];
            scaling[i + refinementData.nBins] = 2.0;
        }

        // Generate Es
        int type = SplineEnergy.Type.FCTOESQ;
        SplineMinimize splineMinimize = new SplineMinimize(reflectionList, refinementData, refinementData.esqFc, type);
        splineMinimize.minimize(7, 1.0);

        type = SplineEnergy.Type.FOTOESQ;
        splineMinimize = new SplineMinimize(reflectionList, refinementData, refinementData.esqFo, type);
        splineMinimize.minimize(7, 1.0);

        setWEstimate();
    }

    /**
     * <p>
     * calculateLikelihoodFree</p>
     *
     * @return a double.
     */
    public double calculateLikelihoodFree() {
        sigmaAEnergy.setScaling(scaling);
        double energy = sigmaAEnergy.energyAndGradient(x, grad);
        sigmaAEnergy.setScaling(null);
        return energy;
    }

    /**
     * <p>getCoordinates.</p>
     *
     * @param x an array of {@link double} objects.
     * @return an array of {@link double} objects.
     */
    public double[] getCoordinates(double[] x) {
        if (x == null) {
            x = new double[this.x.length];
        }
        arraycopy(this.x, 0, x, 0, this.x.length);
        return x;
    }

    /**
     * <p>getNumberOfVariables.</p>
     *
     * @return a int.
     */
    public int getNumberOfVariables() {
        return x.length;
    }

    /**
     * <p>
     * minimize</p>
     *
     * @return a {@link ffx.xray.SigmaAEnergy} object.
     */
    public SigmaAEnergy minimize() {
        return minimize(0.5);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps a double.
     * @return a {@link ffx.xray.SigmaAEnergy} object.
     */
    public SigmaAEnergy minimize(double eps) {
        return minimize(7, eps);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param m   a int.
     * @param eps a double.
     * @return a {@link ffx.xray.SigmaAEnergy} object.
     */
    public SigmaAEnergy minimize(int m, double eps) {

        sigmaAEnergy.setScaling(scaling);

        double e = sigmaAEnergy.energyAndGradient(x, grad);

        long mTime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, sigmaAEnergy, this);
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

        for (int i = 0; i < refinementData.nBins; i++) {
            refinementData.sigmaA[i] = 1.0 + x[i] / scaling[i];
            int index = i + refinementData.nBins;
            refinementData.sigmaW[i] = x[index] / scaling[index];
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            mTime += System.nanoTime();
            sb.append(format(" Optimization time: %8.3f (sec)\n", mTime * toSeconds));
            logger.info(sb.toString());
        }

        sigmaAEnergy.setScaling(null);
        return sigmaAEnergy;
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
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization of SigmaA Parameters\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time");
        }
        if (info == null) {
            logger.info(format("%6d %12.2f %10.2f",
                    iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(format("%6d %12.2f %10.2f %10.5f %9.5f %8.2f %6d %8.3f",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
            } else {
                logger.info(format("%6d %12.2f %10.2f %10.5f %9.5f %8.2f %6d %8s",
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

    private void setWEstimate() {
        // generate initial w estimate
        ReflectionSpline spline = new ReflectionSpline(reflectionList, refinementData.nBins);
        int[] nMean = new int[refinementData.nBins];
        for (int i = 0; i < refinementData.nBins; i++) {
            nMean[i] = 0;
        }
        double mean = 0.0;
        double tot = 0.0;
        double[][] fcTot = refinementData.fcTot;
        double[][] fSigF = refinementData.fSigF;
        for (HKL ih : reflectionList.hkllist) {
            int i = ih.index();
            if (ih.allowed() == 0.0 || isNaN(fcTot[i][0]) || isNaN(fSigF[i][0])) {
                continue;
            }

            double s2 = Crystal.invressq(crystal, ih);
            double epsc = ih.epsilonc();
            ComplexNumber fct = new ComplexNumber(fcTot[i][0], fcTot[i][1]);
            double ecscale = spline.f(s2, refinementData.esqFc);
            double eoscale = spline.f(s2, refinementData.esqFo);
            double ec = fct.times(sqrt(ecscale)).abs();
            double eo = fSigF[i][0] * sqrt(eoscale);
            double wi = pow(eo - ec, 2.0) / epsc;

            nMean[spline.i1()]++;
            tot++;

            x[spline.i1() + refinementData.nBins] += (wi - x[spline.i1() + refinementData.nBins]) / nMean[spline.i1()];
            mean += (wi - mean) / tot;
        }
        logger.info(format(" Starting mean w:    %8.3f", mean));
        logger.info(format(" Starting w scaling: %8.3f", 1.0 / mean));
        for (int i = 0; i < refinementData.nBins; i++) {
            x[i] -= x[i + refinementData.nBins];
            x[i] *= scaling[i];
            scaling[i + refinementData.nBins] = 1.0 / mean;
            x[i + refinementData.nBins] *= scaling[i + refinementData.nBins];
        }
    }

    /**
     * <p>Getter for the field <code>sigmaAEnergy</code>.</p>
     *
     * @return a {@link ffx.xray.SigmaAEnergy} object.
     */
    SigmaAEnergy getSigmaAEnergy() {
        return sigmaAEnergy;
    }

    /**
     * <p>
     * calculateLikelihood</p>
     *
     * @return a double.
     */
    double calculateLikelihood() {
        sigmaAEnergy.setScaling(scaling);
        sigmaAEnergy.energyAndGradient(x, grad);
        sigmaAEnergy.setScaling(null);

        return refinementData.llkR;
    }
}
