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
package ffx.openmm.drude;

import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMDrudeLibrary.OpenMM_DrudeSCFIntegrator_step;

/**
 * DrudeSCFIntegrator implements a self-consistent field (SCF) integrator for systems with Drude oscillators.
 * <p>
 * This integrator extends the basic DrudeIntegrator by implementing a self-consistent field approach
 * for handling Drude polarizable systems. The SCF method iteratively solves for the equilibrium
 * positions of Drude particles at each time step, ensuring that the induced dipoles are consistent
 * with the local electric field environment.
 * <p>
 * The SCF approach is particularly useful for:
 * <ul>
 * <li>Systems requiring highly accurate polarization responses</li>
 * <li>Cases where Drude particles need to be at their instantaneous equilibrium positions</li>
 * <li>Simulations where the polarization energy must be minimized at each step</li>
 * </ul>
 * <p>
 * The integrator performs iterative minimization of the Drude particle positions until
 * convergence is achieved within a specified error tolerance. This ensures that the
 * polarization is always in equilibrium with the instantaneous configuration of charges,
 * providing the most accurate representation of induced dipole interactions.
 */
public class DrudeSCFIntegrator extends DrudeIntegrator {

  /**
   * Create a new DrudeSCFIntegrator.
   * <p>
   * This constructor initializes a Drude SCF integrator with the specified step size.
   * The integrator will use self-consistent field iterations to ensure that Drude
   * particles are at their equilibrium positions at each time step.
   *
   * @param stepSize The integration step size in picoseconds.
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
   * Get the error tolerance for the SCF minimization.
   * <p>
   * This method returns the convergence criterion used in the self-consistent field
   * iterations. The SCF procedure continues until the change in Drude particle
   * positions between iterations falls below this tolerance value.
   *
   * @return The minimization error tolerance.
   */
  public double getMinimizationErrorTolerance() {
    return OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(pointer);
  }

  /**
   * Set the error tolerance for the SCF minimization.
   * <p>
   * This method sets the convergence criterion for the self-consistent field
   * iterations. Smaller tolerance values lead to more accurate polarization
   * but require more iterations. Typical values range from 1e-4 to 1e-6.
   *
   * @param tolerance The minimization error tolerance.
   */
  public void setMinimizationErrorTolerance(double tolerance) {
    OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(pointer, tolerance);
  }

  /**
   * Integrate the system forward in time by the specified number of time steps.
   * <p>
   * This method advances the simulation using the SCF approach. At each time step,
   * the real atoms are propagated using standard dynamics, while the Drude particles
   * are iteratively minimized to their equilibrium positions consistent with the
   * instantaneous electric field environment.
   *
   * @param steps The number of steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMM_DrudeSCFIntegrator_step(pointer, steps);
  }
}