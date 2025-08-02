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

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_addIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getCurrentIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getIntegrationForceGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getNumIntegrators;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_getStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setCurrentIntegrator;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setIntegrationForceGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_setStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CompoundIntegrator_step;

/**
 * This class implements a compound integrator that allows you to use multiple
 * integration algorithms within a single simulation, and to switch between them
 * during the simulation. This is useful for implementing multiple time step
 * algorithms, where different forces are evaluated with different frequencies.
 * <p>
 * To use this class, create a CompoundIntegrator, then call addIntegrator() to
 * add child integrators to it. You can then call setCurrentIntegrator() to
 * specify which child integrator should be used for the next integration step.
 * <p>
 * The compound integrator maintains its own step size, which may be different
 * from the step sizes of the child integrators. When you call step(), it will
 * invoke the current child integrator with the appropriate number of steps to
 * advance the simulation by the compound integrator's step size.
 */
public class CompoundIntegrator extends Integrator {

  /**
   * Create a CompoundIntegrator.
   */
  public CompoundIntegrator() {
    pointer = OpenMM_CompoundIntegrator_create();
  }

  /**
   * Add a child integrator to this CompoundIntegrator.
   *
   * @param integrator The integrator to add.
   * @return The index of the integrator that was added.
   */
  public int addIntegrator(Integrator integrator) {
    return OpenMM_CompoundIntegrator_addIntegrator(pointer, integrator.getPointer());
  }

  /**
   * Destroy the integrator.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CompoundIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the tolerance within which constraints must be satisfied during the simulation.
   *
   * @return The constraint tolerance in nm.
   */
  @Override
  public double getConstraintTolerance() {
    return OpenMM_CompoundIntegrator_getConstraintTolerance(pointer);
  }

  /**
   * Get the index of the current child integrator.
   *
   * @return The index of the current child integrator.
   */
  public int getCurrentIntegrator() {
    return OpenMM_CompoundIntegrator_getCurrentIntegrator(pointer);
  }

  /**
   * Get the set of force groups this integrator acts on.
   *
   * @return The bit flags indicating which force groups this integrator acts on.
   */
  @Override
  public int getIntegrationForceGroups() {
    return OpenMM_CompoundIntegrator_getIntegrationForceGroups(pointer);
  }

  /**
   * Get a child integrator by index.
   *
   * @param index The index of the integrator to get.
   * @return The integrator at the specified index.
   */
  public Integrator getIntegrator(int index) {
    PointerByReference integratorPointer = OpenMM_CompoundIntegrator_getIntegrator(pointer, index);
    Integrator integrator = new Integrator() {
    };
    integrator.setPointer(integratorPointer);
    return integrator;
  }

  /**
   * Get the number of child integrators that have been added.
   *
   * @return The number of child integrators.
   */
  public int getNumIntegrators() {
    return OpenMM_CompoundIntegrator_getNumIntegrators(pointer);
  }

  /**
   * Get the size of each time step, in picoseconds.
   *
   * @return The step size in ps.
   */
  @Override
  public double getStepSize() {
    return OpenMM_CompoundIntegrator_getStepSize(pointer);
  }

  /**
   * Set the tolerance within which constraints must be satisfied during the
   * simulation. The default value is 1e-5 nm.
   *
   * @param tolerance The tolerance within which constraints must be satisfied.
   */
  @Override
  public void setConstraintTolerance(double tolerance) {
    OpenMM_CompoundIntegrator_setConstraintTolerance(pointer, tolerance);
  }

  /**
   * Set the current child integrator.
   *
   * @param index The index of the child integrator to use.
   */
  public void setCurrentIntegrator(int index) {
    OpenMM_CompoundIntegrator_setCurrentIntegrator(pointer, index);
  }

  /**
   * Set the force groups this integrator acts on.
   *
   * @param groups The bit flags indicating which force groups this integrator acts on.
   */
  @Override
  public void setIntegrationForceGroups(int groups) {
    OpenMM_CompoundIntegrator_setIntegrationForceGroups(pointer, groups);
  }

  /**
   * Set the size of each time step, in picoseconds.
   *
   * @param stepSize The step size in ps.
   */
  @Override
  public void setStepSize(double stepSize) {
    OpenMM_CompoundIntegrator_setStepSize(pointer, stepSize);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps The number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMM_CompoundIntegrator_step(pointer, steps);
  }
}