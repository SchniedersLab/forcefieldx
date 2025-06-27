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
package ffx.algorithms.dynamics;

import ffx.numerics.math.RunningStatistics;

import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * NonEquilbriumDynamics is a class that contains methods to control
 * non-equilibrium molecular dynamics simulations.
 */
public class NonEquilbriumDynamics {

  private static final Logger logger = Logger.getLogger(NonEquilbriumDynamics.class.getName());

  /**
   * Begin at L=1 and decrease to L=0.
   */
  private final boolean reverseNEQ;
  /**
   * The number of non-equilibrium lambda steps.
   */
  private int nonEquilibriumLambdaSteps;
  /**
   * The non-equilibrium work values.
   */
  private final RunningStatistics nonEquilibriumWorkValues;
  /**
   * The total number of MD steps.
   */
  private long totalMDSteps = 1;
  /**
   * The non-equilibrium lambda update frequency.
   * <p>
   * The number of MD steps between updates of the non-equilibrium lambda value, which is a function of
   * the total number of MD times steps and the number of non-equilibrium lambda steps.
   */
  private long nonEquilibiumLambdaUpdateFrequency = Long.MAX_VALUE;
  /**
   * The non-equilibrium lambda step size.
   */
  private double lambdaStepSize;
  /**
   * The lambda value to restart non-equilibrium dynamics.
   * Set to 0.0 or 1.0 if doing a new forward or reverse simulation, respectively.
   */
  private double restartLambda;
  /**
   * Restarted from a DYN file. Will cause a change in the number of MD steps, NEQ steps, and what lambda to start at.
   */
  private boolean restartFromDYN = false;

  /**
   * Constructor for NonEquilbriumDynamics.
   *
   * @param nonEquilibriumLambdaSteps The number of non-equilibrium lambda steps.
   * @param reverseNEQ                If true, lambda values should decrease from 1 to 0.
   */
  public NonEquilbriumDynamics(int nonEquilibriumLambdaSteps, boolean reverseNEQ) {
    if (nonEquilibriumLambdaSteps < 1) {
      this.nonEquilibriumLambdaSteps = 100;
    } else {
      this.nonEquilibriumLambdaSteps = nonEquilibriumLambdaSteps;
    }
    this.reverseNEQ = reverseNEQ;
    nonEquilibriumWorkValues = new RunningStatistics();
  }

  /**
   * Get the number of non-equilibrium lambda steps.
   *
   * @return The number of non-equilibrium lambda steps.
   */
  public int getNonEquilibriumLambdaSteps() {
    return nonEquilibriumLambdaSteps;
  }

  /**
   * Get the initial lambda value.
   * @return The initial lambda value.
   */
  public double getInitialLambda() {
    if (!restartFromDYN) {
      restartLambda = reverseNEQ ? 1.0 : 0.0;
    }
    return restartLambda;
  }

  /**
   * Set the lambda value read from a DYN file. Will be read as the initial lambda value.
   *
   * @param dynLambda The lambda value read from DYN file.
   */
  public void setRestartLambda(double dynLambda) {
    restartLambda = dynLambda;
    restartFromDYN = true;
  }

  /**
   * Configure increments of the non-equilibrium lambda values based on the total number of MD steps.
   *
   * @param nSteps The total number of MD steps.
   * @return The total number of MD steps may be adjusted to be a multiple of the non-equilibrium lambda steps.
   */
  public long setMDSteps(long nSteps) {
    if (nSteps < 1) {
      long defaultSteps = 100L * nonEquilibriumLambdaSteps;
      logger.info(format(" Invalid number of MD steps %d. Setting the number of steps to %d.", nSteps, defaultSteps));
      nSteps = defaultSteps;
    }

    if (nSteps % nonEquilibriumLambdaSteps != 0) {
      logger.info(format(" Non-equilibrium lambda steps (%d) is not a multiple of total steps (%d).",
          nonEquilibriumLambdaSteps, nSteps));
      nSteps = nSteps - (nSteps % nonEquilibriumLambdaSteps);
      logger.info(format(" Number of steps adjusted to %d.", nSteps));
    }

    nonEquilibiumLambdaUpdateFrequency = nSteps / nonEquilibriumLambdaSteps;
    lambdaStepSize = 1.0 / nonEquilibriumLambdaSteps;
    if (restartFromDYN) {
      int binReached; // last bin from the previous run
      if (reverseNEQ) {
        binReached = (int) Math.round((1.0-restartLambda) / lambdaStepSize);
      } else {
        binReached = (int) Math.round(restartLambda / lambdaStepSize);
      }
      nonEquilibriumLambdaSteps = nonEquilibriumLambdaSteps - binReached;
      nSteps = nonEquilibiumLambdaUpdateFrequency * nonEquilibriumLambdaSteps - (nonEquilibiumLambdaUpdateFrequency-1); // remove 9 steps
    }
    totalMDSteps = nSteps;
    return nSteps;
  }

  /**
   * Check if the non-equilibrium lambda value should be updated at a given MD step.
   *
   * @param step The MD step number.
   * @return True if the non-equilibrium lambda value should be updated.
   */
  public boolean isUpdateStep(long step) {
    if (step < 1 || step > totalMDSteps) {
      logger.severe(format(" Invalid MD step number %d. Must be between 1 and %d.", step, totalMDSteps));
      return false;
    }
    // The last step is a special case.
    if (step == totalMDSteps) {
      return true;
    }
    return (step - 1) % nonEquilibiumLambdaUpdateFrequency == 0;
  }

  /**
   * Add a work contribution.
   *
   * @param work The work value.
   */
  public void addWork(double work) {
    nonEquilibriumWorkValues.addValue(work);
  }

  /**
   * Get the total work for a given range of lambda bins.
   *
   * @return The total work.
   */
  public double getWork() {
    return nonEquilibriumWorkValues.getSum();
  }

  /**
   * Get the non-equilibrium lambda value for a given MD step.
   *
   * @param step          The MD step number.
   * @param currentLambda The current lambda value.
   * @return The lambda value.
   */
  public double getNextLambda(long step, double currentLambda) {
    if (isUpdateStep(step)) {
      int lambdaBin = getCurrentLambdaBin(step);
      if (reverseNEQ) {
        return restartLambda - lambdaBin * lambdaStepSize;
      } else {
        return restartLambda + lambdaBin * lambdaStepSize;
      }
    } else {
      logger.warning(format(" Non-equilibrium lambda update frequency is %d, but step %d is not a multiple of this frequency.",
          nonEquilibiumLambdaUpdateFrequency, step - 1));
      logger.warning(format(" Returning the current lambda value %6.4f.", currentLambda));
      return currentLambda;
    }
  }

  /**
   * Get the current lambda bin for a given MD step.
   *
   * @param step The MD step number.
   * @return The lambda bin.
   */
  public int getCurrentLambdaBin(long step) {
    if (step == totalMDSteps) {
      return nonEquilibriumLambdaSteps;
    } else if (isUpdateStep(step)) {
      int bin = (int) ((step - 1) / nonEquilibiumLambdaUpdateFrequency);
      if (restartFromDYN) {
        return bin + 1;
      } else {
        return bin;
      }
    } else {
      logger.warning(format(" Non-equilibrium lambda update frequency is %d, but step %d is not a multiple of this frequency.",
          nonEquilibiumLambdaUpdateFrequency, step - 1));
      return 0;
    }
  }

}
