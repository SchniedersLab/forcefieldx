// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
// ******************************************************************************
package ffx.xray;

import edu.rit.pj.ParallelTeam;
import ffx.algorithms.Terminatable;
import ffx.crystal.Crystal;
import ffx.crystal.HKL;
import ffx.crystal.ReflectionList;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch.LineSearchResult;
import ffx.numerics.optimization.OptimizationListener;
import ffx.xray.solvent.SolventModel;

import javax.annotation.Nullable;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.ScalarMath.b2u;
import static ffx.numerics.math.ScalarMath.u2b;
import static ffx.utilities.Constants.NS2SEC;
import static java.lang.Double.isNaN;
import static java.lang.Math.max;
import static java.lang.String.format;
import static java.lang.System.arraycopy;
import static java.util.Arrays.fill;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.log;

/**
 * ScaleBulkMinimize class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class ScaleBulkMinimize implements OptimizationListener, Terminatable {

  private static final Logger logger = Logger.getLogger(ScaleBulkMinimize.class.getName());

  private final ReflectionList reflectionList;
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
   * Constructor for ScaleBulkMinimize.
   *
   * @param reflectionList         a {@link ffx.crystal.ReflectionList} object.
   * @param refinementData         a {@link ffx.xray.DiffractionRefinementData} object.
   * @param crystalReciprocalSpace a {@link ffx.xray.CrystalReciprocalSpace} object.
   * @param parallelTeam           the ParallelTeam to execute the ScaleBulkMinimize.
   */
  public ScaleBulkMinimize(ReflectionList reflectionList,
                           DiffractionRefinementData refinementData,
                           CrystalReciprocalSpace crystalReciprocalSpace,
                           ParallelTeam parallelTeam) {
    this.reflectionList = reflectionList;
    this.refinementData = refinementData;
    this.crystal = reflectionList.crystal;
    this.crystalReciprocalSpace = crystalReciprocalSpace;

    if (crystalReciprocalSpace.getSolventModel() == SolventModel.NONE ||
        refinementData.bulkSolventFixed) {
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

    // Initialize K overall if its zero.
    if (refinementData.modelScaleK == 0.0) {
      x[0] = getInitialKOverall(x);
      logger.info(format("\n Initial K_Overall (Fc to Fo Scale): %4.2f", exp(0.25 * x[0])));
    } else {
      x[0] = refinementData.modelScaleK;
    }

    // Refinement of Bulk Solvent Parameters.
    if (solventN > 1) {
      x[1] = refinementData.bulkSolventK;
      x[2] = refinementData.bulkSolventUeq;
    }
    for (int i = 0; i < 6; i++) {
      if (crystal.scaleB[i] >= 0) {
        x[solventN + crystal.scaleB[i]] = refinementData.modelAnisoB[i];
      }
    }
  }

  /**
   * getCoordinates.
   *
   * @param x the array to populate with parameters or null to create a new array.
   * @return an array containing the parameters.
   */
  public double[] getCoordinates(@Nullable double[] x) {
    if (x == null) {
      x = new double[this.x.length];
    }
    arraycopy(this.x, 0, x, 0, this.x.length);
    return x;
  }

  /**
   * Retrieves the number of variables used in the optimization process.
   *
   * @return the number of variables.
   */
  public int getNumberOfVariables() {
    return x.length;
  }

  /**
   * Performs a minimization procedure with the specified convergence criterion (epsilon).
   *
   * @param eps the convergence threshold for the minimization process; smaller values indicate tighter convergence criteria.
   * @return a {@link ScaleBulkEnergy} object resulting from the minimization process.
   */
  public ScaleBulkEnergy minimize(double eps) {
    return minimize(7, eps);
  }

  /**
   * Performs a minimization process to optimize scaling and bulk solvent energy parameters.
   *
   * @param m   the number of prior gradient evaluations used by the L-BFGS algorithm;
   *            larger values can improve convergence for challenging optimizations,
   *            but may increase memory usage.
   * @param eps the convergence threshold for the minimization process; smaller values indicate
   *            tighter convergence criteria.
   * @return the {@link ScaleBulkEnergy} object representing the refined bulk solvent energy after
   * the minimization process.
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
    }
    for (int i = 0; i < 6; i++) {
      int offset = crystal.scaleB[i];
      if (offset >= 0) {
        int index = solventN + offset;
        refinementData.modelAnisoB[i] = x[index] / scaling[index];
      }
    }

    mtime += System.nanoTime();
    logger.info(format(" Optimization time: %8.3f (sec)\n", mtime * NS2SEC));

    bulkSolventEnergy.setScaling(null);

    return bulkSolventEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean optimizationUpdate(int iter, int nBFGS, int nfun, double grms, double xrms,
                                    double f, double df, double angle, LineSearchResult info) {
    long currentTime = System.nanoTime();
    Double seconds = (currentTime - time) * NS2SEC;
    time = currentTime;
    this.grms = grms;
    this.nSteps = iter;

    if (iter == 0) {
      if (nBFGS > 0) {
        if (solventN > 1) {
          logger.info("\n Limited Memory BFGS Quasi-Newton Optimization of Scaling and Solvent Parameters\n");
        } else {
          logger.info("\n Limited Memory BFGS Quasi-Newton Optimization of Overall Scaling Parameters\n");
        }
      } else {
        if (solventN > 1) {
          logger.info("\n Steepest Decent Optimization of Scaling and Solvent Parameters\n");
        } else {
          logger.info("\n Steepest Decent Optimization of Overall Scaling Parameters\n");
        }
      }
      logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time      R  Rfree");
    }
    if (info == null) {
      logger.info(format("%6d %12.5f %10.6f", iter, f, grms));
    } else {
      if (info == LineSearchResult.Success) {
        double R = bulkSolventEnergy.getR();
        double Rfree = bulkSolventEnergy.getRfree();
        logger.info(format("%6d %12.5f %10.6f %10.6f %9.5f %8.2f %6d %8.3f %5.3f %5.3f",
            iter, f, grms, df, xrms, angle, nfun, seconds, R, Rfree));
      } else {
        logger.info(format("%6d %12.5f %10.6f %10.6f %9.5f %8.2f %6d %8s",
            iter, f, grms, df, xrms, angle, nfun, info));
      }
    }
    if (terminate) {
      logger.info(" The optimization received a termination request.");
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

  /**
   * getScaleBulkEnergy.
   *
   * @return a {@link ffx.xray.ScaleBulkEnergy} object.
   */
  ScaleBulkEnergy getScaleBulkEnergy() {
    return bulkSolventEnergy;
  }

  /**
   * Calculates the initial overall scaling factor (KOverall) for the refinement process
   * based on the provided parameters.
   *
   * @param x an array of double values representing current refinement parameters.
   * @return a double value representing the initial overall scaling factor.
   */
  private double getInitialKOverall(double[] x) {
    double[][] fc = refinementData.fc;
    double[][] fSigf = refinementData.fSigF;
    double[] grad = new double[x.length];

    bulkSolventEnergy.setScaling(scaling);
    double e = bulkSolventEnergy.energyAndGradient(x, grad);
    bulkSolventEnergy.setScaling(null);

    double sumfofc = 0.0;
    double sumfc = 0.0;
    for (HKL ih : reflectionList.hklList) {
      int i = ih.getIndex();
      if (isNaN(fc[i][0]) || isNaN(fSigf[i][0]) || fSigf[i][1] <= 0.0) {
        continue;
      }

      double fcTotF = refinementData.fcTotF(i);
      sumfofc += fSigf[i][0] * fcTotF;
      sumfc += fcTotF * fcTotF;
    }

    // x[0] = log(4.0 * sumfofc / sumfc);

    // logger.info(" Setting initial scale factor.");
    // logger.info(format(" Sum FoFc: %16.8f", sumfofc));
    // logger.info(format(" Sum FcFc: %16.8f", sumfc));
    // logger.info(format(" 4.0 * Log(Sum FoFc / Sum FcFc): %16.8f", x[0]));

    return 4.0 * log(sumfofc / sumfc);
  }

  /**
   * Perform a grid optimization of the bulk solvent model to minimize the residual energy and improve
   * refinement parameters. This method adjusts the bulk solvent parameters (A and B) within a specified
   * range using a grid search algorithm, evaluates the associated energy, and identifies optimal
   * parameters based on the minimum residual.
   */
  public void gridOptimizeBulkSolventModel() {
    if (crystalReciprocalSpace == null) {
      return;
    }

    bulkSolventEnergy.setScaling(scaling);

    // Reset solvent A & B parameters to their default values.
    crystalReciprocalSpace.setDefaultSolventAB();

    // Initial target values.
    crystalReciprocalSpace.computeDensity(refinementData.fs);
    double initialTarget = bulkSolventEnergy.energy(x);
    double min = initialTarget;
    double R = bulkSolventEnergy.getR();
    double Rfree = bulkSolventEnergy.getRfree();
    double minR = R;
    double minRfree = Rfree;

    double solventA = crystalReciprocalSpace.getSolventA();
    double solventB = crystalReciprocalSpace.getSolventB();
    double amin = solventA - 2.0;
    double amax = solventA + 2.0;
    double astep = 0.25;
    double bmin = solventB - 0.2;
    double bmax = solventB + 0.2;
    double bstep = 0.05;
    if (crystalReciprocalSpace.getSolventModel() == SolventModel.BINARY) {
      amin = solventA - 0.4;
      amax = solventA + 0.4;
      astep = 0.05;
    }

    logger.info(" Grid Search for Bulk Solvent Model Parameters");
    logger.info("               A      B    Target      R  Rfree");
    logger.info(format(" Initial: %6.3f %6.3f %9.6f %6.3f %6.3f ",
        solventA, solventB, min, R, Rfree));

    int index = 0;
    for (double a = amin; a <= amax; a += astep) {
      for (double b = bmin; b <= bmax; b += bstep) {
        crystalReciprocalSpace.setSolventAB(a, b);
        crystalReciprocalSpace.computeDensity(refinementData.fs);
        double sum = bulkSolventEnergy.energy(x);
        R = bulkSolventEnergy.getR();
        Rfree = bulkSolventEnergy.getRfree();

        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" %8d %6.3f %6.3f %9.6f %6.3f %6.3f ", index, a, b, sum, R, Rfree));
        } else if (sum < min) {
          logger.info(format(" %8d %6.3f %6.3f %9.6f %6.3f %6.3f ", index, a, b, sum, R, Rfree));
        }

        index++;
        if (sum < min) {
          min = sum;
          minR = R;
          minRfree = Rfree;
          solventA = a;
          solventB = b;
        }
      }
    }

    logger.info(format(" Minimum: %6.3f %6.3f %9.6f %6.3f %6.3f ", solventA, solventB, min, minR, minRfree));
    if (min < initialTarget) {
      crystalReciprocalSpace.setSolventAB(solventA, solventB);
      crystalReciprocalSpace.computeDensity(refinementData.fs);
    }

    bulkSolventEnergy.setScaling(null);
  }

  /**
   * Perform a grid search optimization for bulk solvent and scaling parameters.
   * The method adjusts the scaling factor and bulk solvent parameters (Ks and Bs)
   * to minimize the energy function and improve the refinement results.
   */
  public void gridOptimizeKsBs() {
    if (solventN < 3) {
      return;
    }

    bulkSolventEnergy.setScaling(scaling);

    double bulkSolventK = refinementData.bulkSolventK;
    double bulkSolventUeq = refinementData.bulkSolventUeq;

    // Starting point.
    x[0] = refinementData.modelScaleK;
    x[1] = bulkSolventK;
    x[2] = bulkSolventUeq;
    double initialSum = bulkSolventEnergy.energy(x);
    double initialR = bulkSolventEnergy.getR();
    double initialRfree = bulkSolventEnergy.getRfree();
    double min = initialSum;

    // Prepare for the grid search.
    // Make sure Ks is at least 0.3.
    if (bulkSolventK < 0.3) {
      bulkSolventK = 0.3;
    }
    // Make sure Bs is at least 40.0.
    if (u2b(bulkSolventUeq) < 40.0) {
      bulkSolventUeq = b2u(40.0);
    }

    double kmin = bulkSolventK - 0.3;
    double kmax = bulkSolventK + 0.3;
    double kstep = 0.02;
    double initialB = u2b(bulkSolventUeq);
    double umin = b2u(max(0.0, initialB - 40.0));
    double umax = b2u(initialB + 40.0);
    double ustep = b2u(2.0);

    logger.info("\n Grid Search for Bulk Solvent Ks and Bs");
    logger.info("             Ks       Bs   Target      R  Rfree");
    logger.info(format(" Initial: %5.3f %8.3f %8.5f %5.3f %5.3f",
        bulkSolventK, u2b(bulkSolventUeq), initialSum, initialR, initialRfree));

    int index = 0;
    for (double ks = kmin; ks <= kmax; ks += kstep) {
      for (double ku = umin; ku <= umax; ku += ustep) {
        x[1] = ks;
        x[2] = ku;
        double sum = bulkSolventEnergy.energy(x);
        double R = bulkSolventEnergy.getR();
        double Rfree = bulkSolventEnergy.getRfree();
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" %8d %5.3f %8.3f %8.5f %5.3f %5.3f", index, ks, u2b(ku), sum, R, Rfree));
        } else if (sum < min) {
          logger.info(format(" %8d %5.3f %8.3f %8.5f %5.3f %5.3f", index, ks, u2b(ku), sum, R, Rfree));
        }
        index++;
        if (sum < min) {
          min = sum;
          bulkSolventK = ks;
          bulkSolventUeq = ku;
        }
      }
    }
    x[1] = bulkSolventK;
    x[2] = bulkSolventUeq;
    double sum = bulkSolventEnergy.energy(x);
    double R = bulkSolventEnergy.getR();
    double Rfree = bulkSolventEnergy.getRfree();
    logger.info(format(" Minimum: %5.3f %8.3f %8.5f %5.3f %5.3f",
        bulkSolventK, u2b(bulkSolventUeq), sum, R, Rfree));
    if (sum < initialSum) {
      refinementData.bulkSolventK = bulkSolventK;
      refinementData.bulkSolventUeq = bulkSolventUeq;
    }

    bulkSolventEnergy.setScaling(null);
  }
}
