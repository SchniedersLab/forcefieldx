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
 * An error-controlled, variable-step integrator using the leap-frog Verlet algorithm. OpenMM compares the Verlet
 * result with an explicit Euler result, takes the difference as the integration error of each step, and continuously
 * adjusts the step size to keep it below the specified tolerance.
 *
 * <p>The error tolerance has no absolute meaning; it is an adjustable parameter that affects step size and accuracy,
 * and OpenMM suggests 0.001 as a common starting point. Unlike fixed-step Verlet, this method is not symplectic, so it
 * is most appropriate for constant-temperature simulations rather than long constant-energy simulations. An optional
 * maximum step size can be set.</p>
 *
 * <p>Because the step size varies, the inherited {@link Integrator#getStepSize()} returns the size of the most
 * recent step, and OpenMM states that the effect of {@link Integrator#setStepSize(double)} is undefined (it may be
 * ignored).</p>
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.VariableVerletIntegrator} is not.</p>
 */
public class VariableVerletIntegrator extends Integrator {

  /**
   * Create a variable-step Verlet integrator. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   *
   * @param errorTol error tolerance (dimensionless adjustable accuracy parameter).
   */
  public VariableVerletIntegrator(double errorTol) {
    super(create(errorTol));
  }

  /**
   * Destroy the native integrator.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. If the
   * integrator was passed to a {@link Context}, that context's destruction also destroys the integrator, so do not call
   * this for an integrator that a live context owns.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_VariableVerletIntegrator_destroy);
  }

  /**
   * Get the error tolerance.
   *
   * @return error tolerance (dimensionless adjustable accuracy parameter).
   */
  public double getErrorTolerance() {
    return OpenMMNative.OpenMM_VariableVerletIntegrator_getErrorTolerance(getPointer());
  }

  /**
   * Get the maximum step size the integrator will ever use.
   *
   * @return maximum step size, in ps; 0 (the default) means no limit.
   */
  public double getMaximumStepSize() {
    return OpenMMNative.OpenMM_VariableVerletIntegrator_getMaximumStepSize(getPointer());
  }

  /**
   * Set the error tolerance.
   *
   * @param tol error tolerance (dimensionless adjustable accuracy parameter).
   */
  public void setErrorTolerance(double tol) {
    OpenMMNative.OpenMM_VariableVerletIntegrator_setErrorTolerance(getPointer(), tol);
  }

  /**
   * Set the maximum step size the integrator will ever use, which prevents excessively large steps such as at a
   * local energy minimum.
   *
   * @param size maximum step size, in ps; 0 means no limit.
   */
  public void setMaximumStepSize(double size) {
    OpenMMNative.OpenMM_VariableVerletIntegrator_setMaximumStepSize(getPointer(), size);
  }

  /**
   * Advance the simulation by taking the specified number of adaptive steps.
   *
   * @param steps number of time steps to take.
   *
   * <p>This overrides {@link Integrator#step(int)} with the variable-Verlet native step call.</p>
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_VariableVerletIntegrator_step(getPointer(), steps);
  }

  /**
   * Advance the simulation by adaptive steps until the specified time is reached. On return the simulation time
   * exactly equals the requested time. A time earlier than the current time returns without integrating.
   *
   * @param time target simulation time, in ps.
   */
  public void stepTo(double time) {
    OpenMMNative.OpenMM_VariableVerletIntegrator_stepTo(getPointer(), time);
  }

  /**
   * Create the native integrator after loading the FFM runtime.
   *
   * @param errorTol error tolerance.
   * @return native integrator handle owned by the new wrapper.
   */
  private static MemorySegment create(double errorTol) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_VariableVerletIntegrator_create(errorTol);
  }
}
