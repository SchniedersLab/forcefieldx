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

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch.LineSearchResult;
import ffx.numerics.optimization.OptimizationListener;
import ffx.potential.MolecularAssembly;
import ffx.realspace.RealSpaceData;
import ffx.xray.refine.RefinementMode;
import ffx.xray.refine.RefinementModel;

import javax.annotation.Nullable;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.utilities.Constants.NS2SEC;
import static java.lang.String.format;

/**
 * The RefinementMinimize class is responsible for performing energy minimization for refinement
 * tasks. It provides multiple methods for minimization, allowing customization of gradient
 * root mean square (RMS) values, maximum iterations, and matrix conditioning settings if necessary.
 * The class enables iterative optimization with monitoring capabilities through an optional listener,
 * and supports graceful termination during the process.
 * <p>
 * This class implements both the OptimizationListener and Terminatable interfaces and
 * acts as an intermediary for managing refinement energy computation and convergence behavior.
 *
 * @author Timothy D. Fenn
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RefinementMinimize implements OptimizationListener, Terminatable {

  private static final Logger logger = Logger.getLogger(RefinementMinimize.class.getName());

  /**
   * A {@link DataContainer} instance that serves as the primary data model for the
   * refinement process in {@link RefinementMinimize}.
   * <p>
   * This data container must contain a valid {@link RefinementModel} and
   * at least one type of data representation, such as {@link DiffractionData}
   * or {@link RealSpaceData}. It is utilized to manage parameters, weights,
   * and additional metadata necessary for optimization and energy calculations during
   * refinement.
   * <p>
   * The data within this container supports methods for managing weights, printing
   * optimization updates, and handling data-specific configurations.
   */
  private final DataContainer dataContainer;
  /**
   * The {@code listener} field is an instance of {@link AlgorithmListener}.
   * It serves as an observer to handle updates during the progress of the refinement algorithm.
   * The listener provides mechanisms for monitoring the algorithm's state or terminating
   * the process based on user interaction or external conditions.
   */
  private final AlgorithmListener listener;
  /**
   * The {@code refinementModel} field represents the {@link RefinementModel}
   * used in the optimization process during refinement. This model encapsulates
   * structural and energy-related parameters crucial for refining the input data
   * during the minimization workflows.
   * <p>
   * It is a required component of the {@link RefinementMinimize} class and is
   * initialized during object construction to ensure accurate and efficient
   * refinement operations.
   * <p>
   * The field is immutable and must be present in the associated
   * {@link DataContainer} provided during the instantiation of
   * {@link RefinementMinimize}.
   */
  private final RefinementModel refinementModel;
  /**
   * Represents the energy of the refinement process, utilized during the refinement and optimization
   * steps in the {@link RefinementMinimize} class. This variable encapsulates the current state of
   * the refinement energy, which may include contributions from diffraction or real-space data.
   * <p>
   * The {@link RefinementEnergy} instance is central to the refinement workflow, enabling evaluation,
   * gradients, and convergence determination based on the energy values during minimization stages.
   */
  private final RefinementEnergy refinementEnergy;

  /**
   * The variable n represents the number of data points or observations
   * involved in the refinement process. It is a constant value, initialized
   * during the creation of the RefinementMinimize object, and does not change
   * during the lifetime of the object.
   */
  private final int n;
  /**
   * Stores the current set of variables for the minimization process.
   * Represents an array of doubles used to describe the system's state
   * during the refinement or optimization procedure.
   * <p>
   * This is immutable and finalized to maintain the integrity of the optimization
   * algorithm and prevent unintended modifications during runtime.
   */
  private final double[] x;
  /**
   * The `grad` array stores the current gradient vector for the refinement process.
   * It represents the gradient of the objective function with respect to the parameters
   * being optimized during the minimization procedure.
   * <p>
   * This variable is crucial in determining the direction and magnitude of parameter adjustments
   * to refine the model toward a minimized energy state. It is updated iteratively during the optimization process.
   */
  private final double[] grad;
  /**
   * This array defines scaling factors applied to specific problem variables during optimization.
   * The scaling factors are used to normalize variable gradients, ensuring convergence and
   * stability during numerical minimization. Scaling may depend on the problem's dimensionality
   * or the relative magnitude of specific variable contributions to the objective function.
   * <p>
   * Initialized and used internally by the optimization algorithm to adjust step sizes or
   * condition the optimization landscape.
   */
  private final double[] scaling;
  /**
   * A flag indicating whether the refinement process has completed.
   * <p>
   * The `done` variable is used within the `RefinementMinimize` class to track the
   * completion status of the refinement procedure. It is set to `false` initially
   * and updated to `true` when the process is finished.
   */
  private boolean done = false;
  /**
   * Flag indicating whether the optimization or refinement process should be terminated.
   * <p>
   * This variable is used to signal that the current operation should stop, typically in response
   * to an external or internal condition. It is managed by the class and leveraged by relevant methods
   * to ensure termination of iterative processes when required.
   */
  private boolean terminate = false;
  /**
   * A variable that represents the elapsed time or duration associated
   * with the refinement process in milliseconds.
   * <p>
   * This variable may be used internally to measure or limit the
   * computation time of the refinement operations, such as minimization
   * or optimization processes.
   */
  private long time;
  /**
   * The grms variable represents the gradient root mean square (RMS) of the
   * current optimization step. It is used as a measure of convergence in
   * optimization algorithms, indicating how close the current system is to a
   * locally optimal state. A lower value of grms typically corresponds to a
   * better convergence during the refinement process.
   */
  private double grms;
  /**
   * The variable nSteps represents the number of steps or iterations completed during
   * the refinement minimization process. This value is used to track the progress of
   * the optimization algorithm and ensure it adheres to a defined maximum iteration limit.
   */
  private int nSteps;

  /**
   * Constructor for the RefinementMinimize class, which initializes the refinement process
   * with a specified data container and no algorithm listener.
   *
   * @param dataContainer the data container object containing refinement data and methods
   */
  public RefinementMinimize(DataContainer dataContainer) {
    this(dataContainer, null);
  }

  /**
   * Constructor for the RefinementMinimize class, initializing the refinement process
   * with a specified data container and an optional algorithm listener for updates.
   *
   * @param dataContainer     the data container object containing refinement data and methods
   * @param algorithmListener an optional algorithm listener that provides updates during processing (can be null)
   */
  public RefinementMinimize(DataContainer dataContainer, @Nullable AlgorithmListener algorithmListener) {
    this.dataContainer = dataContainer;
    this.listener = algorithmListener;
    this.refinementModel = dataContainer.getRefinementModel();

    // Create the target potential to optimize.
    refinementEnergy = new RefinementEnergy(dataContainer);

    // Define the number of parameters to optimize.
    n = refinementModel.getNumParameters();
    x = new double[n];
    grad = new double[n];
    scaling = new double[n];
    refinementModel.loadOptimizationScaling(scaling);

    refinementEnergy.getCoordinates(x);
    refinementEnergy.setScaling(scaling);
  }

  /**
   * Minimizes the refinement energy using a default gradient root mean square (RMS) value of 1.0.
   *
   * @return a {@link RefinementEnergy} object representing the result of the minimization process.
   */
  public RefinementEnergy minimize() {
    return minimize(1.0);
  }

  /**
   * Minimizes the refinement energy using a specified gradient root mean square (RMS) value.
   *
   * @param eps the desired gradient RMS value for the refinement process
   * @return a {@link RefinementEnergy} object representing the result of the minimization process
   */
  public RefinementEnergy minimize(double eps) {
    return minimize(7, eps, Integer.MAX_VALUE);
  }

  /**
   * Minimizes the refinement energy with a specified maximum number of iterations,
   * using a gradient RMS value of 1.0 and default matrix conditioning cycles.
   *
   * @param maxIterations the maximum number of iterations allowed for the minimization process
   * @return a {@link RefinementEnergy} object representing the result of the minimization process
   */
  public RefinementEnergy minimize(int maxIterations) {
    return minimize(7, 1.0, maxIterations);
  }

  /**
   * Minimizes the refinement energy using a specified gradient root mean square (RMS) value
   * and maximum number of iterations.
   *
   * @param eps           the desired gradient RMS value for the refinement process
   * @param maxIterations the maximum number of iterations allowed for the minimization process
   * @return a {@link RefinementEnergy} object representing the result of the minimization process
   */
  public RefinementEnergy minimize(double eps, int maxIterations) {
    return minimize(7, eps, maxIterations);
  }

  /**
   * Minimizes the refinement energy using the provided matrix conditioning cycles, gradient RMS value,
   * and a maximum number of iterations. The default number of matrix conditioning steps is 7, while
   * 0 is for steepest decent.
   *
   * @param m       the number of matrix conditioning cycles for the optimization process
   * @param eps     the desired gradient root mean square (RMS) value for convergence
   * @param maxiter the maximum number of iterations allowed for the minimization process
   * @return a {@link RefinementEnergy} object representing the outcome of the minimization process
   */
  public RefinementEnergy minimize(int m, double eps, int maxiter) {
    if (dataContainer instanceof DiffractionData) {
      logger.info(" Beginning X-ray Refinement");
    } else if (dataContainer instanceof RealSpaceData) {
      logger.info(" Beginning Real Space Refinement");
    } else {
      logger.info(" Beginning Refinement");
    }

    RefinementMode refinementMode = refinementModel.getRefinementMode();
    logger.info(refinementMode.toString());
    logger.info(refinementModel.toString());

    refinementEnergy.getCoordinates(x);

    // Scale coordinates.
    for (int i = 0; i < n; i++) {
      x[i] *= scaling[i];
    }

    long mtime = -System.nanoTime();
    time = -System.nanoTime();
    done = false;
    int status;
    double e = refinementEnergy.energyAndGradient(x, grad);
    status = LBFGS.minimize(n, m, x, e, grad, eps, maxiter, refinementEnergy, this);
    done = true;
    switch (status) {
      case 0:
        logger.info(format("\n Optimization achieved convergence criteria: %8.5f", grms));
        break;
      case 1:
        logger.info(format("\n Optimization terminated at step %d.", nSteps));
        break;
      default:
        logger.warning("\n Optimization failed.");
    }

    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder();
      mtime += System.nanoTime();
      sb.append(format(" Optimization time: %g (sec)", mtime * NS2SEC));
      logger.info(sb.toString());
    }

    return refinementEnergy;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean optimizationUpdate(int iter, int nBFGS, int nfun, double grms, double xrms, double f,
                                    double df, double angle, @Nullable LineSearchResult info) {

    long currentTime = System.nanoTime();
    double seconds = (currentTime - time) * NS2SEC;
    time = currentTime;
    this.grms = grms;
    this.nSteps = iter;

    // Update display.
    if (listener != null) {
      RefinementModel refinementModel = dataContainer.getRefinementModel();
      MolecularAssembly[] molecularAssembly = refinementModel.getMolecularAssemblies();
      for (MolecularAssembly ma : molecularAssembly) {
        listener.algorithmUpdate(ma);
      }
    }

    if (iter == 0) {
      if (nBFGS > 0) {
        logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n");
      } else {
        logger.info("\n Steepest Decent Optimization: \n");
      }
      logger.info(
          " Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time      "
              + dataContainer.printOptimizationHeader());
    }
    if (info == null) {
      logger.info(format("%6d %12.3f %10.3f", iter, f, grms));
    } else if (info == LineSearchResult.Success) {
      StringBuilder sb = new StringBuilder();
      sb.append(format("%6d %12.3f %10.3f %10.3f %9.4f %8.2f %6d %8.3f ",
          iter, f, grms, df, xrms, angle, nfun, seconds));
      sb.append(dataContainer.printOptimizationUpdate());
      logger.info(sb.toString());
    } else {
      logger.info(
          format("%6d %12.3f %10.3f %10.3f %9.4f %8.2f %6d %8s",
              iter, f, grms, df, xrms, angle, nfun, info));
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
          logger.log(Level.WARNING, " Exception terminating minimization.\n", e);
        }
      }
    }
  }

}
