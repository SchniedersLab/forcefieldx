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
import edu.uiowa.jopenmm.OpenMM_Vec3;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_applyConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_reinitialize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocities;
import static ffx.openmm.Vec3Array.toVec3Array;

/**
 * A Context stores the complete state of a simulation.  More specifically, it includes:
 *
 * <ul>
 * <li>The current time</li>
 * <li>The position of each particle</li>
 * <li>The velocity of each particle</li>
 * <li>The values of configurable parameters defined by Force objects in the System</li>
 * </ul>
 *
 * You can retrieve a snapshot of the current state at any time by calling getState().  This
 * allows you to record the state of the simulation at various points, either for analysis
 * or for checkpointing.  getState() can also be used to retrieve the current forces on each
 * particle and the current energy of the System.
 */
public class Context {

  /**
   * Context pointer.
   */
  private PointerByReference pointer;

  /**
   * Constructor.
   */
  public Context() {
    pointer = null;
  }

  /**
   * Constructor.
   *
   * @param system     The system to simulate.
   * @param integrator The integrator to use for simulating the system.
   * @param platform   The platform to use for performing computations.
   */
  public Context(System system, Integrator integrator, Platform platform) {
    pointer = OpenMM_Context_create_2(system.getPointer(), integrator.getPointer(), platform.getPointer());
  }

  /**
   * Update the context.
   *
   * @param system     The system to simulate.
   * @param integrator The integrator to use for simulating the system.
   * @param platform   The platform to use for performing computations.
   */
  public void updateContext(System system, Integrator integrator, Platform platform) {
    // Destroy the old context.
    destroy();
    pointer = OpenMM_Context_create_2(system.getPointer(), integrator.getPointer(), platform.getPointer());
  }

  /**
   * Get the pointer to the context.
   *
   * @return The pointer to the context.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Does the context have a pointer?
   *
   * @return True if the context pointer is not null.
   */
  public boolean hasContextPointer() {
    return pointer != null;
  }

  /**
   * Get the state of the context.
   *
   * @param types              The bit-flag of properties to include in the state.
   * @param enforcePeriodicBox If true, the positions will be adjusted so atoms are inside the main periodic box.
   * @return The state of the context.
   */
  public State getState(int types, int enforcePeriodicBox) {
    return new State(OpenMM_Context_getState(pointer, types, enforcePeriodicBox));
  }

  /**
   * Reinitialize the context.
   *
   * <p>When a Context is created, it may cache information about the System being simulated and
   * the Force objects contained in it. This means that, if the System or Forces are then modified,
   * the Context might not see all changes. Call reinitialize() to force the Context to rebuild its
   * internal representation of the System and pick up any changes that have been made.
   *
   * <p>This is an expensive operation, so you should try to avoid calling it too frequently.
   *
   * @param preserveState If true, the state will be restored to the same state it had before the call.
   */
  public void reinitialize(int preserveState) {
    if (pointer != null) {
      OpenMM_Context_reinitialize(pointer, preserveState);
    }
  }

  /**
   * Set a parameter value.
   *
   * @param name  The name of the parameter.
   * @param value The value of the parameter.
   */
  public void setParameter(String name, double value) {
    OpenMM_Context_setParameter(pointer, name, value);
  }

  /**
   * Set the periodic box vectors.
   *
   * @param a The first vector.
   * @param b The second vector.
   * @param c The third vector.
   */
  public void setPeriodicBoxVectors(OpenMM_Vec3 a, OpenMM_Vec3 b, OpenMM_Vec3 c) {
    OpenMM_Context_setPeriodicBoxVectors(pointer, a, b, c);
  }

  /**
   * Set the atomic positions.
   *
   * @param positions The atomic positions.
   */
  public void setPositions(double[] positions) {
    Vec3Array vec3Array = toVec3Array(positions);
    OpenMM_Context_setPositions(pointer, vec3Array.getPointer());
    vec3Array.destroy();
  }

  /**
   * Set the atomic velocities.
   *
   * @param velocities The atomic velocities.
   */
  public void setVelocities(double[] velocities) {
    Vec3Array velArray = toVec3Array(velocities);
    OpenMM_Context_setVelocities(pointer, velArray.getPointer());
    velArray.destroy();
  }

  /**
   * Apply constraints to the current positions.
   *
   * @param tol The constraint tolerance.
   */
  public void applyConstraints(double tol) {
    OpenMM_Context_applyConstraints(pointer, tol);
  }

  /**
   * Destroy the context.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_Context_destroy(pointer);
      pointer = null;
    }
  }

}
