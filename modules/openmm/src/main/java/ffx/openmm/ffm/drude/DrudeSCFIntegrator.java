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
package ffx.openmm.ffm.drude;

import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * Leap-frog Verlet integrator that uses the self-consistent-field method for Drude particles.
 *
 * <p>At each time step, the positions of Drude particles are adjusted to minimize potential
 * energy. The System must contain a DrudeForce so the integrator can identify Drude particles.</p>
 */
public class DrudeSCFIntegrator extends DrudeIntegrator {

  /**
   * Create a Drude SCF integrator.
   *
   * @param stepSize Integration time step in picoseconds.
   */
  public DrudeSCFIntegrator(double stepSize) {
    super(create(stepSize));
  }

  /**
   * Release the native integrator handle.
   *
   * <p>The integrator must not be used after it is destroyed.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_DrudeSCFIntegrator_destroy);
  }

  /**
   * Get the force tolerance used when minimizing Drude coordinates.
   *
   * <p>The tolerance roughly corresponds to the maximum allowed force magnitude on Drude
   * particles after minimization.</p>
   *
   * @return Minimization force tolerance in kJ/mol/nm.
   */
  public double getMinimizationErrorTolerance() {
    return OpenMMNative.OpenMM_DrudeSCFIntegrator_getMinimizationErrorTolerance(getPointer());
  }

  /**
   * Set the force tolerance used when minimizing Drude coordinates.
   *
   * <p>The tolerance roughly corresponds to the maximum allowed force magnitude on Drude
   * particles after minimization.</p>
   *
   * @param tolerance Minimization force tolerance in kJ/mol/nm.
   */
  public void setMinimizationErrorTolerance(double tolerance) {
    OpenMMNative.OpenMM_DrudeSCFIntegrator_setMinimizationErrorTolerance(
        getPointer(), tolerance);
  }

  /**
   * Advance the simulation by the requested number of time steps.
   *
   * @param steps Number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_DrudeSCFIntegrator_step(getPointer(), steps);
  }

  private static MemorySegment create(double stepSize) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_DrudeSCFIntegrator_create(stepSize);
  }
}
