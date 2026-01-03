// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_getFriction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_getTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_setFriction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_setTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_step;

/**
 * This is an Integrator which simulates a System using Langevin dynamics, with
 * the LFMiddle discretization (J. Phys. Chem. A 2019, 123, 28, 6056-6079).
 * This method tend to produce more accurate configurational sampling than other
 * discretizations, such as the one used in LangevinIntegrator.
 * <p>
 * The algorithm is closely related to the BAOAB discretization
 * (Proc. R. Soc. A. 472: 20160138).  Both methods produce identical trajectories,
 * but LFMiddle returns half step (leapfrog) velocities, while BAOAB returns
 * on-step velocities.  The former provide a much more accurate sampling of the
 * thermal ensemble.
 */
public class LangevinMiddleIntegrator extends Integrator {

  /**
   * Create a LangevinMiddleIntegrator.
   *
   * @param dt    the step size with which to integrate the system (in picoseconds)
   * @param temp  the temperature of the heat bath (in Kelvin)
   * @param gamma the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
   */
  public LangevinMiddleIntegrator(double dt, double temp, double gamma) {
    super(OpenMM_LangevinMiddleIntegrator_create(temp, gamma, dt));
  }

  /**
   * Destroy the integrator.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_LangevinMiddleIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the friction coefficient.
   *
   * @return The friction coefficient in inverse picoseconds.
   */
  public double getFriction() {
    return OpenMM_LangevinMiddleIntegrator_getFriction(pointer);
  }

  /**
   * Get the random number seed.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_LangevinMiddleIntegrator_getRandomNumberSeed(pointer);
  }

  /**
   * Get the temperature.
   *
   * @return The temperature in Kelvin.
   */
  public double getTemperature() {
    return OpenMM_LangevinMiddleIntegrator_getTemperature(pointer);
  }

  /**
   * Set the friction coefficient.
   *
   * @param gamma The friction coefficient in inverse picoseconds.
   */
  public void setFriction(double gamma) {
    OpenMM_LangevinMiddleIntegrator_setFriction(pointer, gamma);
  }

  /**
   * Set the random number seed.
   *
   * @param seed The random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_LangevinMiddleIntegrator_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Set the temperature.
   *
   * @param temp The temperature in Kelvin.
   */
  public void setTemperature(double temp) {
    OpenMM_LangevinMiddleIntegrator_setTemperature(pointer, temp);
  }

  /**
   * Step the integrator.
   *
   * @param steps The number of steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMM_LangevinMiddleIntegrator_step(pointer, steps);
  }
}