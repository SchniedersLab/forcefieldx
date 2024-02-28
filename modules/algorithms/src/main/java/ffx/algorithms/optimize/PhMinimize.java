// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.algorithms.optimize;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.numerics.Potential;
import ffx.numerics.optimization.LBFGS;
import ffx.numerics.optimization.LineSearch;
import ffx.numerics.optimization.OptimizationListener;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.extended.ExtendedSystem;

import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.fill;

/**
 * Minimize the potential energy of a system to an RMS gradient per atom convergence criteria.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PhMinimize implements OptimizationListener, Terminatable {

  private static final Logger logger = Logger.getLogger(PhMinimize.class.getName());

  /** The MolecularAssembly being operated on. */
  protected final MolecularAssembly molecularAssembly;
  /** The potential energy to optimize. */
  protected final Potential potential;
  /** The extended system that contains the fictitious particle */
  protected final ExtendedSystem esvSystem;
  /** The AlgorithmListener to update the UI. */
  protected final AlgorithmListener algorithmListener;
  /** Number of variables. */
  protected final int n;
  /** Number of Extended System variables. */
  protected final int nESV;
  /** Current value of each variable. */
  protected final double[] x;
  /** Current value of each Extended System variable. */
  protected double[] theta;
  /** The gradient. */
  protected final double[] grad;
  /** The gradient. */
  protected double[] gradESV;
  /** Scaling applied to each variable. */
  protected final double[] scaling;
  /** A flag to indicate the algorithm is done. */
  protected boolean done = false;
  /** A flag to indicate the algorithm should be terminated. */
  protected boolean terminate = false;
  /** Minimization time in nanoseconds. */
  protected long time;
  /** The final potential energy. */
  protected double energy;
  /** The return status of the optimization. */
  protected int status;
  /** The number of optimization steps taken. */
  protected int nSteps;
  /** The final RMS gradient. */
  double rmsGradient;

  /**
   * Constructor for Minimize.
   *
   * @param molecularAssembly a {@link MolecularAssembly} object.
   * @param potential a {@link Potential} object.
   * @param algorithmListener a {@link AlgorithmListener} object.
   */
  public PhMinimize(MolecularAssembly molecularAssembly, Potential potential,
      AlgorithmListener algorithmListener, ExtendedSystem esvSystem) {
    this.molecularAssembly = molecularAssembly;
    this.algorithmListener = algorithmListener;
    this.potential = potential;
    this.esvSystem = esvSystem;
    n = potential.getNumberOfVariables();
    nESV = esvSystem.getNumberOfVariables();
    x = new double[n];
    theta = new double[nESV];
    grad = new double[n];
    gradESV = new double[nESV];
    scaling = new double[n];
    fill(scaling, 12.0);
  }

  /**
   * Constructor for Minimize.
   *
   * @param molecularAssembly a {@link MolecularAssembly} object.
   * @param algorithmListener a {@link AlgorithmListener} object.
   */
  public PhMinimize(MolecularAssembly molecularAssembly, AlgorithmListener algorithmListener,
      ExtendedSystem esvSystem) {
    this.molecularAssembly = molecularAssembly;
    this.algorithmListener = algorithmListener;
    this.esvSystem = esvSystem;
    if (molecularAssembly.getPotentialEnergy() == null) {
      molecularAssembly.setPotential(ForceFieldEnergy.energyFactory(molecularAssembly));
    }
    potential = molecularAssembly.getPotentialEnergy();
    n = potential.getNumberOfVariables();
    nESV = esvSystem.getNumberOfVariables();
    x = new double[n];
    theta = new double[nESV];
    grad = new double[n];
    gradESV = new double[nESV];
    scaling = new double[n];
    fill(scaling, 12.0);
  }

  /**
   * Getter for the field <code>energy</code>.
   *
   * @return a double.
   */
  public double getEnergy() {
    return energy;
  }

  /**
   * getRMSGradient.
   *
   * @return a double.
   */
  public double getRMSGradient() {
    return rmsGradient;
  }

  /**
   * Getter for the field <code>status</code>.
   *
   * @return The status of the optimization.
   */
  public int getStatus() {
    return status;
  }

  /**
   * minimize
   *
   * @return a {@link Potential} object.
   */
  public Potential minimizeCoordinates() {
    return minimizeCoordinates(7, 1.0, Integer.MAX_VALUE);
  }

  /**
   * minimize
   *
   * @param eps The convergence criteria.
   * @return a {@link Potential} object.
   */
  public Potential minimizeCoordinates(double eps) {
    return minimizeCoordinates(7, eps, Integer.MAX_VALUE);
  }

  /**
   * minimize
   *
   * @param eps The convergence criteria.
   * @param maxIterations The maximum number of iterations.
   * @return a {@link Potential} object.
   */
  public Potential minimizeCoordinates(double eps, int maxIterations) {
    return minimizeCoordinates(7, eps, maxIterations);
  }

  /**
   * minimize
   *
   * @param m The number of previous steps used to estimate the Hessian.
   * @param eps The convergence criteria.
   * @param maxIterations The maximum number of iterations.
   * @return a {@link Potential} object.
   */
  public Potential minimizeCoordinates(int m, double eps, int maxIterations) {
    time = System.nanoTime();
    potential.getCoordinates(x);
    potential.setScaling(scaling);

    // Scale coordinates.
    for (int i = 0; i < n; i++) {
      x[i] *= scaling[i];
    }

    done = false;
    energy = potential.energyAndGradient(x, grad);

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format(" Minimize initial energy: %16.8f", energy));
    }

    status = LBFGS.minimize(n, m, x, energy, grad, eps, maxIterations, potential, this);
    done = true;

    switch (status) {
      case 0 -> logger.info(format("\n Optimization achieved convergence criteria: %8.5f", rmsGradient));
      case 1 -> logger.info(format("\n Optimization terminated at step %d.", nSteps));
      default -> logger.warning("\n Optimization failed.");
    }

    potential.setScaling(null);
    return potential;
  }

  /**
   * minimize
   *
   * @param eps The convergence criteria.
   * @param maxIterations The maximum number of iterations.
   * @return a {@link Potential} object.
   */
  public Potential minimizeTitration(double eps, int maxIterations) {
    return minimizeTitration(7, eps, maxIterations);
  }

  /**
   * minimize
   *
   * @param m The number of previous steps used to estimate the Hessian.
   * @param eps The convergence criteria.
   * @param maxIterations The maximum number of iterations.
   * @return a {@link Potential} object.
   */
  public Potential minimizeTitration(int m, double eps, int maxIterations) {
    time = System.nanoTime();
    theta = esvSystem.getThetaPosition();

    done = false;
    potential.getCoordinates(x);
    theta = esvSystem.getThetaPosition();
    energy = esvSystem.energyAndGradient(theta, gradESV);

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format(" Minimize initial energy: %16.8f", energy));
    }

    status = LBFGS.minimize(nESV, m, theta, energy, gradESV, eps, maxIterations, esvSystem, this);
    done = true;

    switch (status) {
      case 0 -> {
        logger.info(format("\n Optimization achieved convergence criteria: %8.5f", rmsGradient));
        for (Atom atom : molecularAssembly.getAtomList()) {
          int atomIndex = atom.getIndex() - 1;
          atom.setOccupancy(esvSystem.getTitrationLambda(atomIndex));
          atom.setTempFactor(esvSystem.getTautomerLambda(atomIndex));
        }
      }
      case 1 -> logger.info(format("\n Optimization terminated at step %d.", nSteps));
      default -> logger.warning("\n Optimization failed.");
    }

    potential.setScaling(null);
    return potential;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Implement the OptimizationListener interface.
   *
   * @since 1.0
   */
  @Override
  public boolean optimizationUpdate(int iteration, int nBFGS, int functionEvaluations,
      double rmsGradient, double rmsCoordinateChange, double energy, double energyChange,
      double angle, LineSearch.LineSearchResult lineSearchResult) {
    long currentTime = System.nanoTime();
    Double seconds = (currentTime - time) * 1.0e-9;
    time = currentTime;
    this.rmsGradient = rmsGradient;
    this.nSteps = iteration;
    this.energy = energy;

    if (iteration == 0) {
      if (nBFGS > 0) {
        logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n");
      } else {
        logger.info("\n Steepest Decent Optimization: \n");
      }
      logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time\n");
    }
    if (lineSearchResult == null) {
      logger.info(format("%6d%13.4f%11.4f", iteration, energy, rmsGradient));
    } else {
      if (lineSearchResult == LineSearch.LineSearchResult.Success) {
        logger.info(
            format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8.3f", iteration, energy, rmsGradient,
                energyChange, rmsCoordinateChange, angle, functionEvaluations, seconds));
      } else {
        logger.info(format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8s", iteration, energy, rmsGradient,
            energyChange, rmsCoordinateChange, angle, functionEvaluations, lineSearchResult));
      }
    }
    // Update the listener and check for a termination request.
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

  /** {@inheritDoc} */
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
