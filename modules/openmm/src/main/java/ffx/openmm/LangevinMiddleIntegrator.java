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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinMiddleIntegrator_step;

/**
 * Langevin Integrator.
 */
public class LangevinMiddleIntegrator extends Integrator {

  /**
   * Constructor.
   *
   * @param dt    The time step.
   * @param temp  The temperature.
   * @param gamma The friction coefficient.
   */
  public LangevinMiddleIntegrator(double dt, double temp, double gamma) {
    pointer = OpenMM_LangevinMiddleIntegrator_create(temp, gamma, dt);
  }

  /**
   * Step the integrator.
   *
   * @param steps The number of steps to take.
   */
  public void step(int steps) {
    OpenMM_LangevinMiddleIntegrator_step(pointer, steps);
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
   * Destroy the integrator.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_LangevinMiddleIntegrator_destroy(pointer);
      pointer = null;
    }
  }
}
