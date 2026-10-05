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

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.util.Objects;

/**
 * Executes an OpenMM simulation for a system on a selected platform.
 *
 * <p>A context refers to its system, integrator, and platform. The caller must keep the native
 * system alive for the context's lifetime. This façade retains the integrator wrapper and destroys
 * it when the context is destroyed, matching the JNA façade lifecycle; the system and platform
 * wrappers are not destroyed by the context.</p>
 *
 * <p>This façade does not expose OpenMM's constructor overload that accepts a map of
 * platform-specific properties.</p>
 */
public class Context extends OpenMMHandle {

  private Integrator integrator;
  private Platform platform;

  /**
   * Create an uninitialized context wrapper without a native context.
   *
   * <p>Call {@link #updateContext(System, Integrator, Platform)} before using native context
   * operations.</p>
   */
  public Context() {
    super(MemorySegment.NULL, true);
    integrator = null;
    platform = null;
  }

  /**
   * Create a context using OpenMM's selected platform.
   *
   * @param system     live system definition; its native handle must remain alive while this
   *                   context exists.
   * @param integrator live integrator; this façade destroys its wrapper's native integrator when
   *                   the context is destroyed.
   */
  public Context(System system, Integrator integrator) {
    super(create(system, integrator));
    this.integrator = integrator;
    platform = null;
  }

  /**
   * Create a context on a specified platform. The platform wrapper is borrowed and is not
   * destroyed by the context.
   *
   * @param system     live system definition; its native handle must remain alive while this
   *                   context exists.
   * @param integrator live integrator; this façade destroys its wrapper's native integrator when
   *                   the context is destroyed.
   * @param platform   live registered platform to use; its native handle is borrowed.
   */
  public Context(System system, Integrator integrator, Platform platform) {
    super(create(system, integrator, platform));
    this.integrator = integrator;
    this.platform = platform;
  }

  /**
   * Adjust particle positions to satisfy distance constraints, and recompute virtual-site
   * positions.
   *
   * @param tolerance distance tolerance for satisfying the constraints, in nanometers.
   */
  public void applyConstraints(double tolerance) {
    OpenMMNative.OpenMM_Context_applyConstraints(getPointer(), tolerance);
  }

  /**
   * Adjust particle velocities so that the velocity along every constrained distance is zero.
   *
   * @param tolerance velocity tolerance in nanometers per picosecond.
   */
  public void applyVelocityConstraints(double tolerance) {
    OpenMMNative.OpenMM_Context_applyVelocityConstraints(getPointer(), tolerance);
  }

  /**
   * Recompute positions of all virtual sites without applying distance constraints.
   */
  public void computeVirtualSites() {
    OpenMMNative.OpenMM_Context_computeVirtualSites(getPointer());
  }

  /**
   * Destroy the native context and this façade's associated integrator wrapper.
   *
   * <p>The system and platform handles are not destroyed. Repeated calls have no effect.</p>
   */
  @Override
  public void destroy() {
    Integrator currentIntegrator = integrator;
    integrator = null;
    if (currentIntegrator != null) {
      currentIntegrator.destroy();
    }
    destroy(OpenMMNative::OpenMM_Context_destroy);
  }

  /**
   * Get the integrator wrapper associated with this context.
   *
   * @return associated integrator wrapper, or {@code null} if none is associated or the context
   *         has been destroyed.
   */
  public Integrator getIntegrator() {
    return integrator;
  }

  /**
   * Get the current value of a global context parameter.
   *
   * <p>The name is encoded as temporary NUL-terminated UTF-8 for the native call.</p>
   *
   * @param name non-null parameter name.
   * @return parameter value.
   */
  public double getParameter(String name) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_Context_getParameter(getPointer(), value));
  }

  /**
   * Get the current value of a context parameter using a caller-provided native name.
   *
   * <p>The caller retains ownership of {@code name}, which must be a readable NUL-terminated
   * UTF-8 string for the duration of the call.</p>
   *
   * @param name caller-owned NUL-terminated UTF-8 parameter name.
   * @return parameter value.
   */
  public double getParameter(MemorySegment name) {
    return OpenMMNative.OpenMM_Context_getParameter(getPointer(), name);
  }

  /**
   * Get the borrowed native parameter-array view for this context.
   *
   * @return borrowed opaque {@code OpenMM_ParameterArray} handle, valid only while this context is
   *         alive; this façade does not convert it to a Java map.
   */
  public MemorySegment getParameters() {
    return OpenMMNative.OpenMM_Context_getParameters(getPointer());
  }

  /**
   * Get the explicitly supplied platform wrapper, if any.
   *
   * @return platform wrapper, or {@code null} when OpenMM selected the platform automatically.
   */
  public Platform getPlatform() {
    return platform;
  }

  /**
   * Create an owned snapshot of selected context data.
   *
   * <p>{@code types} is a bitwise combination of OpenMM state-data flags. When periodic-box
   * enforcement is requested, OpenMM translates molecule positions so each molecule's center lies
   * in one periodic box.</p>
   *
   * @param types              bitwise combination of OpenMM state-data flags.
   * @param enforcePeriodicBox whether to translate molecule positions into one periodic box.
   * @return owned state snapshot; caller must close it.
   */
  public State getState(int types, boolean enforcePeriodicBox) {
    return new State(OpenMMNative.OpenMM_Context_getState(
        getPointer(), types, OpenMMBooleans.toNative(enforcePeriodicBox)));
  }

  /**
   * Create an owned snapshot using the native integer representation of the periodic-box flag.
   *
   * @param types              bitwise combination of OpenMM state-data flags.
   * @param enforcePeriodicBox nonzero to translate molecule positions into one periodic box.
   * @return owned state snapshot; caller must close it.
   */
  public State getState(int types, int enforcePeriodicBox) {
    return new State(OpenMMNative.OpenMM_Context_getState(getPointer(), types, enforcePeriodicBox));
  }

  /**
   * Create an owned snapshot of selected data and force groups.
   *
   * <p>{@code groups} is a bit mask: force group {@code i} is included when bit {@code i} is set.
   * The group selection applies when forces or energies are requested.</p>
   *
   * @param types              bitwise combination of OpenMM state-data flags.
   * @param enforcePeriodicBox whether to translate molecule positions into one periodic box.
   * @param groups             bit mask selecting force groups for force and energy calculations.
   * @return owned state snapshot; caller must close it.
   */
  public State getState(int types, boolean enforcePeriodicBox, int groups) {
    return new State(OpenMMNative.OpenMM_Context_getState_2(
        getPointer(), types, OpenMMBooleans.toNative(enforcePeriodicBox), groups));
  }

  /**
   * Create an owned snapshot of selected data and force groups using the native integer flag.
   *
   * @param types              bitwise combination of OpenMM state-data flags.
   * @param enforcePeriodicBox nonzero to translate molecule positions into one periodic box.
   * @param groups             bit mask selecting force groups for force and energy calculations.
   * @return owned state snapshot; caller must close it.
   */
  public State getState(int types, int enforcePeriodicBox, int groups) {
    return new State(OpenMMNative.OpenMM_Context_getState_2(
        getPointer(), types, enforcePeriodicBox, groups));
  }

  /**
   * Get a non-owning wrapper for the native system simulated by this context.
   *
   * @return non-owning system wrapper; closing it only invalidates that wrapper. It is usable only
   *         while the original owning system remains alive.
   */
  public System getSystem() {
    return new System(OpenMMNative.OpenMM_Context_getSystem(getPointer()), false);
  }

  /**
   * @return {@code true} if this wrapper currently holds a non-null native context handle.
   */
  public boolean hasContextPointer() {
    return !isDestroyed();
  }

  /**
   * @return number of completed integration steps.
   */
  public long getStepCount() {
    return OpenMMNative.OpenMM_Context_getStepCount(getPointer());
  }

  /**
   * @return simulation time in picoseconds.
   */
  public double getTime() {
    return OpenMMNative.OpenMM_Context_getTime(getPointer());
  }

  /**
   * Reinitialize this context after changing its system or forces.
   *
   * <p>This is expensive. By default, reinitialization discards the current state; setting
   * {@code preserveState} asks OpenMM to restore state through a checkpoint and can fail if the
   * modified system is incompatible with that checkpoint.</p>
   *
   * @param preserveState whether to attempt to preserve positions, velocities, parameters, and time.
   */
  public void reinitialize(boolean preserveState) {
    OpenMMNative.OpenMM_Context_reinitialize(
        getPointer(), OpenMMBooleans.toNative(preserveState));
  }

  /**
   * Reinitialize this context using the native integer representation of the state flag.
   *
   * @param preserveState nonzero to attempt to preserve positions, velocities, parameters, and
   *                     time.
   */
  public void reinitialize(int preserveState) {
    OpenMMNative.OpenMM_Context_reinitialize(getPointer(), preserveState);
  }

  /**
   * Set the value of a global context parameter.
   *
   * <p>The name is encoded as temporary NUL-terminated UTF-8 for the native call.</p>
   *
   * @param name  non-null parameter name.
   * @param value parameter value.
   */
  public void setParameter(String name, double value) {
    OpenMMStrings.withUtf8String(name,
        pointer -> OpenMMNative.OpenMM_Context_setParameter(getPointer(), pointer, value));
  }

  /**
   * Set a context parameter using a caller-provided native name.
   *
   * <p>The caller retains ownership of {@code name}, which must be a readable NUL-terminated
   * UTF-8 string for the duration of the call.</p>
   *
   * @param name  caller-owned NUL-terminated UTF-8 parameter name.
   * @param value new parameter value.
   */
  public void setParameter(MemorySegment name, double value) {
    OpenMMNative.OpenMM_Context_setParameter(getPointer(), name, value);
  }

  /**
   * Set the periodic box vectors in nanometers. Vectors must satisfy OpenMM's periodic-box
   * requirements.
   *
   * @param a first box vector in nm.
   * @param b second box vector in nm.
   * @param c third box vector in nm.
   */
  public void setPeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c) {
    try (Arena arena = Arena.ofConfined()) {
      OpenMMNative.OpenMM_Context_setPeriodicBoxVectors(
          getPointer(), a.toNative(arena), b.toNative(arena), c.toNative(arena));
    }
  }

  /**
   * Set particle positions from packed {@code x,y,z} coordinates, with one vector per system
   * particle.
   *
   * @param positions non-null packed coordinates of length three times the system particle count,
   *                  in nanometers.
   */
  public void setPositions(double[] positions) {
    try (Vec3Array vectors = Vec3Array.toVec3Array(positions)) {
      setPositions(vectors);
    }
  }

  /**
   * Set particle positions from a native FFM vector array with one vector per system particle.
   *
   * @param positions live position vectors in nanometers; ownership remains with the caller.
   */
  public void setPositions(Vec3Array positions) {
    Objects.requireNonNull(positions, "Positions cannot be null.");
    OpenMMNative.OpenMM_Context_setPositions(getPointer(), positions.getPointer());
  }

  /**
   * Copy the data included in a state snapshot into this context.
   *
   * <p>Information not included in the state is left unchanged. A state snapshot is not a full
   * checkpoint and does not restore all internal simulation state, such as random-number generator
   * state.</p>
   *
   * @param state live state snapshot to apply; ownership remains with the caller.
   */
  public void setState(State state) {
    Objects.requireNonNull(state, "State cannot be null.");
    OpenMMNative.OpenMM_Context_setState(getPointer(), state.getPointer());
  }

  /**
   * Set the integration step count.
   *
   * @param count new integration step count.
   */
  public void setStepCount(long count) {
    OpenMMNative.OpenMM_Context_setStepCount(getPointer(), count);
  }

  /**
   * Set simulation time.
   *
   * @param time simulation time in picoseconds.
   */
  public void setTime(double time) {
    OpenMMNative.OpenMM_Context_setTime(getPointer(), time);
  }

  /**
   * Set particle velocities from packed {@code x,y,z} components, with one vector per system
   * particle.
   *
   * @param velocities non-null packed velocities of length three times the system particle count,
   *                   in nanometers per picosecond.
   */
  public void setVelocities(double[] velocities) {
    try (Vec3Array vectors = Vec3Array.toVec3Array(velocities)) {
      setVelocities(vectors);
    }
  }

  /**
   * Set particle velocities from a native FFM vector array with one vector per system particle.
   *
   * @param velocities live velocity vectors in nanometers per picosecond; ownership remains with
   *                   the caller.
   */
  public void setVelocities(Vec3Array velocities) {
    Objects.requireNonNull(velocities, "Velocities cannot be null.");
    OpenMMNative.OpenMM_Context_setVelocities(getPointer(), velocities.getPointer());
  }

  /**
   * Assign velocities sampled for the specified temperature.
   *
   * @param temperature temperature in kelvin.
   * @param randomSeed  random-number seed passed to OpenMM.
   */
  public void setVelocitiesToTemperature(double temperature, int randomSeed) {
    OpenMMNative.OpenMM_Context_setVelocitiesToTemperature(
        getPointer(), temperature, randomSeed);
  }

  /**
   * Replace this wrapper's context, destroying its current context and associated integrator first.
   *
   * @param system     live system to simulate; its native handle must remain alive while the
   *                   context exists.
   * @param integrator live integrator to bind to the new context.
   * @param platform   live registered platform on which to execute.
   */
  public void updateContext(System system, Integrator integrator, Platform platform) {
    destroy();
    replacePointer(create(system, integrator, platform));
    this.integrator = integrator;
    this.platform = platform;
  }

  /**
   * Create a context with OpenMM-selected platform.
   *
   * @param system     system definition.
   * @param integrator integrator.
   * @return owned native context handle.
   */
  private static MemorySegment create(System system, Integrator integrator) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_Context_create(
        Objects.requireNonNull(system, "System cannot be null.").getPointer(),
        Objects.requireNonNull(integrator, "Integrator cannot be null.").getPointer());
  }

  /**
   * Create a context on an explicit platform.
   *
   * @param system     system definition.
   * @param integrator integrator.
   * @param platform   platform.
   * @return owned native context handle.
   */
  private static MemorySegment create(System system, Integrator integrator, Platform platform) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_Context_create_2(
        Objects.requireNonNull(system, "System cannot be null.").getPointer(),
        Objects.requireNonNull(integrator, "Integrator cannot be null.").getPointer(),
        Objects.requireNonNull(platform, "Platform cannot be null.").getPointer());
  }
}
