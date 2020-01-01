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
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.log;

import edu.rit.pj.ParallelTeam;

import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch.LineSearchResult;
import ffx.numerics.optimization.OptimizationListener;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

/**
 * <p>
 * ScaleBulkMinimize class.</p>
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class ScaleBulkMinimize implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(ScaleBulkMinimize.class.getName());
    private static double toSeconds = 1.0e-9;
    private final ReflectionList reflectionlist;
    private final DiffractionRefinementData refinementData;
    private final Crystal crystal;
    private final CrystalReciprocalSpace crystalReciprocalSpace;
    private final ScaleBulkEnergy bulkSolventEnergy;
    private final int solventN;
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
     * Constructor for ScaleBulkMinimize.</p>
     *
     * @param reflectionList         a {@link ffx.crystal.ReflectionList} object.
     * @param refinementData         a {@link ffx.xray.DiffractionRefinementData} object.
     * @param crystalReciprocalSpace a {@link ffx.xray.CrystalReciprocalSpace} object.
     * @param parallelTeam           the ParallelTeam to execute the ScaleBulkMinimize.
     */
    public ScaleBulkMinimize(ReflectionList reflectionList, DiffractionRefinementData refinementData,
                             CrystalReciprocalSpace crystalReciprocalSpace, ParallelTeam parallelTeam) {
        this.reflectionlist = reflectionList;
        this.refinementData = refinementData;
        this.crystal = reflectionList.crystal;
        this.crystalReciprocalSpace = crystalReciprocalSpace;

        if (crystalReciprocalSpace.solventModel == SolventModel.NONE) {
            solventN = 1;
        } else {
            solventN = 3;
        }
        n = solventN + refinementData.nScale;
        bulkSolventEnergy = new ScaleBulkEnergy(reflectionList, refinementData, n, parallelTeam);

        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        fill(scaling, 1.0);

        x[0] = refinementData.modelScaleK;
        if (solventN > 1) {
            x[1] = refinementData.bulkSolventK;
            x[2] = refinementData.bulkSolventUeq;
        }
        for (int i = 0; i < 6; i++) {
            if (crystal.scaleB[i] >= 0) {
                x[solventN + crystal.scaleB[i]] = refinementData.modelAnisoB[i];
            }
        }

        setInitialScale();
    }

    /**
     * <p>getScaleBulkEnergy.</p>
     *
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    ScaleBulkEnergy getScaleBulkEnergy() {
        return bulkSolventEnergy;
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

    private void setInitialScale() {
        double[][] fc = refinementData.fc;
        double[][] fSigf = refinementData.fSigF;

        bulkSolventEnergy.setScaling(scaling);
        double e = bulkSolventEnergy.energyAndGradient(x, grad);
        bulkSolventEnergy.setScaling(null);

        double sumfofc = 0.0;
        double sumfc = 0.0;
        for (HKL ih : reflectionlist.hkllist) {
            int i = ih.index();
            if (isNaN(fc[i][0]) || isNaN(fSigf[i][0]) || fSigf[i][1] <= 0.0) {
                continue;
            }

            double fcTotF = refinementData.fcTotF(i);
            sumfofc += fSigf[i][0] * fcTotF;
            sumfc += fcTotF * fcTotF;
        }

        x[0] = log(4.0 * sumfofc / sumfc);
    }

    /**
     * <p>
     * ksbsGridOptimize</p>
     */
    public void ksbsGridOptimize() {
        if (solventN < 3) {
            return;
        }

        bulkSolventEnergy.setScaling(scaling);

        double min = Double.POSITIVE_INFINITY;
        double bulkSolventK = refinementData.bulkSolventK;
        double bulkSolventUeq = refinementData.bulkSolventUeq;
        double kmin = 0.05;
        double kmax = 0.9;
        double kstep = 0.05;
        double bmin = 10.0;
        double bmax = 200.0;
        double bstep = 10.0;
        for (double i = kmin; i <= kmax; i += kstep) {
            for (double j = bmin; j <= bmax; j += bstep) {
                x[1] = i;
                x[2] = j;
                double sum = bulkSolventEnergy.energyAndGradient(x, grad);
                logger.info(" ks: " + i + " bs: " + j + " sum: " + sum);
                if (sum < min) {
                    min = sum;
                    bulkSolventK = i;
                    bulkSolventUeq = j;
                }
            }
        }
        logger.info(" minks: " + bulkSolventK + " minbs: " + bulkSolventUeq + " min: " + min);
        refinementData.bulkSolventK = bulkSolventK;
        refinementData.bulkSolventUeq = bulkSolventUeq;
        bulkSolventEnergy.setScaling(null);
    }

    /**
     * <p>
     * GridOptimize</p>
     */
    void gridOptimize() {
        if (crystalReciprocalSpace == null) {
            return;
        }

        bulkSolventEnergy.setScaling(scaling);

        double min = Double.POSITIVE_INFINITY;
        double solventA = crystalReciprocalSpace.solventA;
        double solventB = crystalReciprocalSpace.solventB;
        double amin = solventA - 1.0;
        double amax = (solventA + 1.0) / 0.9999;
        double astep = 0.25;
        double bmin = solventB - 0.2;
        double bmax = (solventB + 0.2) / 0.9999;
        double bstep = 0.05;
        if (crystalReciprocalSpace.solventModel == SolventModel.BINARY) {
            amin = solventA - 0.2;
            amax = (solventA + 0.2) / 0.9999;
            astep = 0.05;
        }

        logger.info(" Bulk Solvent Grid Search");
        for (double i = amin; i <= amax; i += astep) {
            for (double j = bmin; j <= bmax; j += bstep) {
                crystalReciprocalSpace.setSolventAB(i, j);
                crystalReciprocalSpace.computeDensity(refinementData.fs);
                double sum = bulkSolventEnergy.energy(x);
                logger.info(format(" A: %6.3f B: %6.3f Sum: %12.8f", i, j, sum));
                if (sum < min) {
                    min = sum;
                    solventA = i;
                    solventB = j;
                }
            }
        }

        logger.info(format("\n Minimum at\n A: %6.3f B: %6.3f Sum: %12.8f", solventA, solventB, min));
        crystalReciprocalSpace.setSolventAB(solventA, solventB);
        refinementData.solventA = solventA;
        refinementData.solventB = solventB;
        crystalReciprocalSpace.computeDensity(refinementData.fs);
        bulkSolventEnergy.setScaling(null);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    public ScaleBulkEnergy minimize() {
        return minimize(0.5);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps a double.
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    public ScaleBulkEnergy minimize(double eps) {
        return minimize(5, eps);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param m   a int.
     * @param eps a double.
     * @return a {@link ffx.xray.ScaleBulkEnergy} object.
     */
    public ScaleBulkEnergy minimize(int m, double eps) {
        bulkSolventEnergy.setScaling(scaling);
        double e = bulkSolventEnergy.energyAndGradient(x, grad);

        long mtime = -System.nanoTime();
        time = -System.nanoTime();
        done = false;
        int status = LBFGS.minimize(n, m, x, e, grad, eps, bulkSolventEnergy, this);
        done = true;
        switch (status) {
            case 0:
                logger.info(format("\n Optimization achieved convergence criteria: %10.5f\n", grms));
                break;
            case 1:
                logger.info(format("\n Optimization terminated at step %d.\n", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.\n");
        }
        refinementData.modelScaleK = x[0] / scaling[0];
        if (solventN > 1) {
            refinementData.bulkSolventK = x[1] / scaling[1];
            refinementData.bulkSolventUeq = x[2] / scaling[2];
            if (crystalReciprocalSpace != null) {
                refinementData.solventA = crystalReciprocalSpace.solventA;
                refinementData.solventB = crystalReciprocalSpace.solventB;
            }
        }
        for (int i = 0; i < 6; i++) {
            int offset = crystal.scaleB[i];
            if (offset >= 0) {
                int index = solventN + offset;
                refinementData.modelAnisoB[i] = x[index] / scaling[index];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            mtime += System.nanoTime();
            logger.info(format(" Optimization time: %8.3f (sec)\n", mtime * toSeconds));
        }

        bulkSolventEnergy.setScaling(null);

        return bulkSolventEnergy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms,
                                      double f, double df, double angle, LineSearchResult info) {
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
            logger.info(format("%6d %12.5f %10.6f",
                    iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(format("%6d %12.5f %10.6f %10.6f %9.5f %8.2f %6d %8.3f",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
            } else {
                logger.info(format("%6d %12.5f %10.6f %10.6f %9.5f %8.2f %6d %8s",
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
