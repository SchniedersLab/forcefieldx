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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_getErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_getMaximumStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_setErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_setMaximumStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VariableVerletIntegrator_stepTo;

/**
 * This class implements a Verlet integrator with variable time stepping.
 * The integrator automatically adjusts the step size to maintain a specified
 * error tolerance, making it suitable for systems with widely varying time scales.
 * <p>
 * Unlike the standard Verlet integrator which uses a fixed step size, this
 * variable step size algorithm monitors the local truncation error and adjusts
 * the step size accordingly. This can lead to more efficient integration for
 * systems where different parts evolve on different time scales.
 */
public class VariableVerletIntegrator extends Integrator {

  /**
   * Create a VariableVerletIntegrator.
   *
   * @param errorTol The error tolerance for adaptive step sizing.
   */
  public VariableVerletIntegrator(double errorTol) {
    super(OpenMM_VariableVerletIntegrator_create(errorTol));
  }

  /**
   * Destroy the integrator.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_VariableVerletIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the error tolerance for adaptive step sizing.
   *
   * @return The error tolerance.
   */
  public double getErrorTolerance() {
    return OpenMM_VariableVerletIntegrator_getErrorTolerance(pointer);
  }

  /**
   * Get the maximum step size the integrator is allowed to use (in ps).
   *
   * @return The maximum step size.
   */
  public double getMaximumStepSize() {
    return OpenMM_VariableVerletIntegrator_getMaximumStepSize(pointer);
  }

  /**
   * Set the error tolerance for adaptive step sizing.
   *
   * @param tol The error tolerance.
   */
  public void setErrorTolerance(double tol) {
    OpenMM_VariableVerletIntegrator_setErrorTolerance(pointer, tol);
  }

  /**
   * Set the maximum step size the integrator is allowed to use (in ps).
   *
   * @param size The maximum step size.
   */
  public void setMaximumStepSize(double size) {
    OpenMM_VariableVerletIntegrator_setMaximumStepSize(pointer, size);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps The number of time steps to take.
   */
  public void step(int steps) {
    OpenMM_VariableVerletIntegrator_step(pointer, steps);
  }

  /**
   * Advance the simulation by integrating until a specified time is reached.
   *
   * @param time The time to which the simulation should be advanced (in ps).
   */
  public void stepTo(double time) {
    OpenMM_VariableVerletIntegrator_stepTo(pointer, time);
  }
}