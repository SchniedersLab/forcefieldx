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
package ffx.openmm;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_getFriction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_getTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_setFriction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_setTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BrownianIntegrator_step;

/**
 * This is an Integrator which simulates a System using Brownian dynamics.
 */
public class BrownianIntegrator extends Integrator {

  /**
   * Create a BrownianIntegrator.
   *
   * @param temperature   the temperature of the heat bath (in Kelvin)
   * @param frictionCoeff the friction coefficient which couples the system to the heat bath, measured in 1/ps
   * @param stepSize      the step size with which to integrate the system (in picoseconds)
   */
  public BrownianIntegrator(double temperature, double frictionCoeff, double stepSize) {
    super(OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, stepSize));
  }

  /**
   * Destroy the integrator.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_BrownianIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the friction coefficient which determines how strongly the system is coupled to
   * the heat bath (in inverse ps).
   *
   * @return the friction coefficient, measured in 1/ps
   */
  public double getFriction() {
    return OpenMM_BrownianIntegrator_getFriction(pointer);
  }

  /**
   * Get the random number seed.  See setRandomNumberSeed() for details.
   */
  public int getRandomNumberSeed() {
    return OpenMM_BrownianIntegrator_getRandomNumberSeed(pointer);
  }

  /**
   * Get the temperature of the heat bath (in Kelvin).
   *
   * @return the temperature of the heat bath (in Kelvin).
   */
  public double getTemperature() {
    return OpenMM_BrownianIntegrator_getTemperature(pointer);
  }

  /**
   * Set the friction coefficient which determines how strongly the system is coupled to
   * the heat bath (in inverse ps).
   *
   * @param coeff the friction coefficient, measured in 1/ps
   */
  public void setFriction(double coeff) {
    OpenMM_BrownianIntegrator_setFriction(pointer, coeff);
  }

  /**
   * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
   * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
   * are run with different random number seeds, the sequence of random forces will be different.  On
   * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
   * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
   * results on successive runs, even if those runs were initialized identically.
   * <p>
   * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
   * is created from this Force. This is done to ensure that each Context receives unique random seeds
   * without you needing to set them explicitly.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_BrownianIntegrator_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Set the temperature of the heat bath (in Kelvin).
   *
   * @param temp the temperature of the heat bath, measured in Kelvin.
   */
  public void setTemperature(double temp) {
    OpenMM_BrownianIntegrator_setTemperature(pointer, temp);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps the number of time steps to take
   */
  public void step(int steps) {
    OpenMM_BrownianIntegrator_step(pointer, steps);
  }
}