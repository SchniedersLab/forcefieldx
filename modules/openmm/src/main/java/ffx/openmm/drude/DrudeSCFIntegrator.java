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
package ffx.openmm.drude;

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_step;

/**
 * This is a leap-frog Verlet Integrator that simulates systems with Drude particles.  It uses the
 * self-consistent field (SCF) method: at every time step, the positions of Drude particles are
 * adjusted to minimize the potential energy.
 * <p>
 * This Integrator requires the System to include a DrudeForce, which it uses to identify the Drude
 * particles.
 */
public class DrudeSCFIntegrator extends DrudeIntegrator {

  /**
   * Create a DrudeSCFIntegrator.
   *
   * @param stepSize the step size with which to integrator the system (in picoseconds)
   */
  public DrudeSCFIntegrator(double stepSize) {
    super(OpenMM_DrudeSCFIntegrator_create(stepSize));
  }

  /**
   * Destroy the integrator.
   * <p>
   * This method releases the memory associated with the DrudeSCFIntegrator object.
   * After calling this method, the integrator should not be used.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_DrudeSCFIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the error tolerance to use when minimizing the potential energy.  This roughly corresponds
   * to the maximum allowed force magnitude on the Drude particles after minimization.
   *
   * @return the error tolerance to use, measured in kJ/mol/nm
   */
  public double getMinimizationErrorTolerance() {
    return OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(pointer);
  }

  /**
   * Set the error tolerance to use when minimizing the potential energy.  This roughly corresponds
   * to the maximum allowed force magnitude on the Drude particles after minimization.
   *
   * @param tolerance the error tolerance to use, measured in kJ/mol/nm
   */
  public void setMinimizationErrorTolerance(double tolerance) {
    OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(pointer, tolerance);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps the number of time steps to take
   */
  @Override
  public void step(int steps) {
    OpenMM_DrudeSCFIntegrator_step(pointer, steps);
  }
}