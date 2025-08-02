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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_getConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_getIntegrationForceGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_getStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setIntegrationForceGroups;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setStepSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;

/**
 * An Integrator defines a method for simulating a System by integrating the equations of motion.
 * This is an abstract class. Subclasses define particular integration methods.
 * <p>
 * Each Integrator object is bound to a particular Context which it integrates. This connection
 * is specified by passing the Integrator as an argument to the constructor of the Context.
 */
public abstract class Integrator {

  /**
   * OpenMM Integrator pointer.
   */
  protected PointerByReference pointer;

  /**
   * Constructor.
   */
  public Integrator() {
    pointer = null;
  }

  /**
   * This method will be called by subclasses when the integrator is destroyed.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_Integrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the tolerance within which constraints must be satisfied during the simulation.
   *
   * @return The constraint tolerance in nm.
   */
  public double getConstraintTolerance() {
    return OpenMM_Integrator_getConstraintTolerance(pointer);
  }

  /**
   * Get the set of force groups this integrator acts on.
   *
   * @return The bit flags indicating which force groups this integrator acts on.
   */
  public int getIntegrationForceGroups() {
    return OpenMM_Integrator_getIntegrationForceGroups(pointer);
  }

  /**
   * Get the OpenMM Integrator pointer.
   *
   * @return The OpenMM Integrator pointer.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get the size of each time step, in picoseconds.
   *
   * @return The step size in ps.
   */
  public double getStepSize() {
    return OpenMM_Integrator_getStepSize(pointer);
  }

  /**
   * Set the tolerance within which constraints must be satisfied during the
   * simulation. The default value is 1e-5 nm.
   *
   * @param tolerance The tolerance within which constraints must be satisfied.
   */
  public void setConstraintTolerance(double tolerance) {
    OpenMM_Integrator_setConstraintTolerance(pointer, tolerance);
  }

  /**
   * Set the force groups this integrator acts on.
   *
   * @param groups The bit flags indicating which force groups this integrator acts on.
   */
  public void setIntegrationForceGroups(int groups) {
    OpenMM_Integrator_setIntegrationForceGroups(pointer, groups);
  }

  /**
   * Set the OpenMM Integrator pointer.
   *
   * @param pointer The OpenMM Integrator pointer.
   */
  public void setPointer(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Set the size of each time step, in picoseconds.
   *
   * @param stepSize The step size in ps.
   */
  public void setStepSize(double stepSize) {
    OpenMM_Integrator_setStepSize(pointer, stepSize);
  }

  /**
   * Integrate the system forward in time by the specified number of time steps.
   *
   * @param steps The number of steps to take.
   */
  public void step(int steps) {
    OpenMM_Integrator_step(pointer, steps);
  }
}