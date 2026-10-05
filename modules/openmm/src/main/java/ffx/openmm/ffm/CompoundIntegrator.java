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
package ffx.openmm.ffm;

import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * Allows several integration algorithms to be used within one simulation, switching between them.
 *
 * <p>Create the other integrators, add each with {@link #addIntegrator(Integrator)}, create a {@link Context}
 * with this compound integrator, and select the active integrator with {@link #setCurrentIntegrator(int)}. The active
 * integrator handles all {@link #step(int)} calls until it is changed. All integrators must be added before the context is
 * created, per the OpenMM header.</p>
 *
 * <p>Switching integrators requires that they interpret positions and velocities identically. Leapfrog-style
 * integrators assume velocities are offset from positions by half a time step, so velocities must be adjusted when
 * switching between a leapfrog and a non-leapfrog integrator, or between leapfrog integrators with different step sizes,
 * to avoid introducing error.</p>
 *
 * <p>The inherited step-size, constraint-tolerance and integration-force-group accessors act on the current
 * integrator, except that {@link #setIntegrationForceGroups(int)} sets the force groups of all contained integrators.</p>
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>Contained integrators are owned by this compound integrator once added: OpenMM deletes them when this object is
 * deleted, and this wrapper invalidates the wrapper of each added integrator.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.CompoundIntegrator} is not. This
 * class invalidates the added integrator's wrapper to reflect the ownership transfer; the JNA {@code addIntegrator}
 * does not. {@link #getIntegrator(int)} returns a borrowed {@link MemorySegment}, whereas the JNA method returns a
 * {@code PointerByReference}.</p>
 */
public class CompoundIntegrator extends Integrator {

  /**
   * Create an empty compound integrator. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   */
  public CompoundIntegrator() {
    super(create());
  }

  /**
   * Add an integrator to this compound integrator.
   *
   * <p>Native ownership passes to the compound integrator, which deletes it when it is itself deleted. After the call
   * the passed wrapper is invalidated (its handle is cleared without calling the native destructor), so it must not be used,
   * destroyed or added elsewhere. All integrators must be added before the {@link Context} is created.</p>
   *
   * @param integrator live integrator to add; must not be null and must not already be owned elsewhere.
   * @return index of the integrator that was added.
   * @throws NullPointerException if {@code integrator} is null.
   * @throws IllegalStateException if {@code integrator} or this object has been destroyed.
   */
  public int addIntegrator(Integrator integrator) {
    int index = OpenMMNative.OpenMM_CompoundIntegrator_addIntegrator(
        getPointer(), integrator.getPointer());
    integrator.invalidate();
    return index;
  }

  /**
   * Destroy the native integrator.
   *
   * <p>The native handle is released once and this wrapper is invalidated, along with the contained integrators owned by it; repeated calls have no effect. If the
   * integrator was passed to a {@link Context}, that context's destruction also destroys the integrator, so do not call
   * this for an integrator that a live context owns.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CompoundIntegrator_destroy);
  }

  /**
   * Get the distance tolerance within which constraints are maintained, as a fraction of the constrained
   * distance. This calls the native method for the current integrator.
   *
   * @return constraint tolerance (dimensionless fraction of the constrained distance, per the OpenMM header).
   */
  @Override
  public double getConstraintTolerance() {
    return OpenMMNative.OpenMM_CompoundIntegrator_getConstraintTolerance(getPointer());
  }

  /**
   * Get the index of the current integrator.
   *
   * @return index of the integrator used by {@link #step(int)}.
   */
  public int getCurrentIntegrator() {
    return OpenMMNative.OpenMM_CompoundIntegrator_getCurrentIntegrator(getPointer());
  }

  /**
   * Get the force groups used for integration, as bit flags (group i is included if {@code (groups & (1 << i)) != 0};
   * all groups by default). This returns the value for the current integrator.
   *
   * @return force-group bit flags.
   */
  @Override
  public int getIntegrationForceGroups() {
    return OpenMMNative.OpenMM_CompoundIntegrator_getIntegrationForceGroups(getPointer());
  }

  /**
   * Get a contained integrator.
   *
   * <p>The returned segment is a borrowed, non-owning native handle that is valid only while this compound integrator is
   * alive; it is not wrapped in a Java class and must not be destroyed by the caller.</p>
   *
   * @param index index of the integrator, as returned by {@link #addIntegrator(Integrator)}.
   * @return borrowed native integrator handle.
   */
  public MemorySegment getIntegrator(int index) {
    return OpenMMNative.OpenMM_CompoundIntegrator_getIntegrator(getPointer(), index);
  }

  /**
   * Get the number of integrators that have been added.
   *
   * @return number of contained integrators.
   */
  public int getNumIntegrators() {
    return OpenMMNative.OpenMM_CompoundIntegrator_getNumIntegrators(getPointer());
  }

  /**
   * Get the size of each time step, from the current integrator.
   *
   * @return step size, in ps.
   */
  @Override
  public double getStepSize() {
    return OpenMMNative.OpenMM_CompoundIntegrator_getStepSize(getPointer());
  }

  /**
   * Set the constraint tolerance on the current integrator.
   *
   * @param tolerance distance tolerance within which constraints are maintained, as a fraction of the constrained
   *     distance (per the OpenMM header).
   */
  @Override
  public void setConstraintTolerance(double tolerance) {
    OpenMMNative.OpenMM_CompoundIntegrator_setConstraintTolerance(getPointer(), tolerance);
  }

  /**
   * Select the integrator used by subsequent {@link #step(int)} calls. Make sure the integrators are compatible; see
   * the class description.
   *
   * @param index index of the integrator to use.
   */
  public void setCurrentIntegrator(int index) {
    OpenMMNative.OpenMM_CompoundIntegrator_setCurrentIntegrator(getPointer(), index);
  }

  /**
   * Set the force groups used for integration for all contained integrators, as bit flags (group i is included if
   * {@code (groups & (1 << i)) != 0}).
   *
   * @param groups force-group bit flags.
   */
  @Override
  public void setIntegrationForceGroups(int groups) {
    OpenMMNative.OpenMM_CompoundIntegrator_setIntegrationForceGroups(getPointer(), groups);
  }

  /**
   * Set the size of each time step on the current integrator.
   *
   * @param stepSize step size, in ps.
   */
  @Override
  public void setStepSize(double stepSize) {
    OpenMMNative.OpenMM_CompoundIntegrator_setStepSize(getPointer(), stepSize);
  }

  /**
   * Advance the simulation by calling step on the current integrator.
   *
   * @param steps number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_CompoundIntegrator_step(getPointer(), steps);
  }

  /**
   * Create the native compound integrator after loading the FFM runtime.
   *
   * @return native integrator handle owned by the new wrapper.
   */
  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_CompoundIntegrator_create.makeInvoker().apply();
  }
}
