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
import java.lang.foreign.ValueLayout;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * Owns or borrows an OpenMM molecular system containing particles, constraints, forces, virtual
 * sites, and default periodic box vectors.
 *
 * <p>A system defines particles, constraints, forces, virtual sites, and default periodic box
 * vectors. Adding a force or virtual site transfers its native ownership to the system; this
 * wrapper tracks and invalidates the corresponding Java wrappers when the native system releases
 * those objects. A non-owning system wrapper, such as one returned by {@link Context#getSystem()},
 * only invalidates itself when closed.</p>
 */
public class System extends OpenMMHandle {

  private final Map<Long, Force> forces = new HashMap<>();
  private final Map<Long, VirtualSite> virtualSites = new HashMap<>();
  private boolean ownsHandle;

  /**
   * Create an empty OpenMM molecular system.
   */
  public System() {
    super(create());
    ownsHandle = true;
  }

  /**
   * Wrap a native system handle and assume ownership of it.
   *
   * @param pointer native system handle.
   */
  public System(MemorySegment pointer) {
    super(pointer);
    ownsHandle = true;
  }

  /**
   * Wrap a native system handle with explicit ownership semantics.
   *
   * @param pointer    native system handle.
   * @param ownsHandle whether this wrapper should destroy the native system.
   */
  System(MemorySegment pointer, boolean ownsHandle) {
    super(pointer);
    this.ownsHandle = ownsHandle;
  }

  /**
   * Add a distance constraint between two particles. Both particles must be in this system and
   * have nonzero mass; the distance is specified in nanometers.
   *
   * @param particle1 zero-based index of the first particle.
   * @param particle2 zero-based index of the second particle.
   * @param distance  constrained distance in nm.
   * @return constraint index.
   */
  public int addConstraint(int particle1, int particle2, double distance) {
    return OpenMMNative.OpenMM_System_addConstraint(
        getPointer(), particle1, particle2, distance);
  }

  /**
   * Add a force to this system and transfer native ownership to it.
   *
   * <p>The force wrapper remains usable while the system owns the native force, but must not be
   * destroyed separately. Removing the force or destroying the system invalidates the wrapper.</p>
   *
   * @param force live, non-null FFM force to add.
   * @return force index.
   */
  public int addForce(Force force) {
    if (force != null) {
      int index = OpenMMNative.OpenMM_System_addForce(getPointer(), force.getPointer());
      force.setForceIndex(index);
      forces.put(force.getPointer().address(), force);
      return index;
    }
    return -1;
  }

  /**
   * Add a particle to the system.
   *
   * @param mass particle mass in daltons (atomic mass units); zero denotes a particle whose
   *             position and velocity are not modified by integrators.
   * @return particle index.
   */
  public int addParticle(double mass) {
    return OpenMMNative.OpenMM_System_addParticle(getPointer(), mass);
  }

  /**
   * Destroy the native system if this wrapper owns it, and invalidate this wrapper's tracked force
   * and virtual-site wrappers.
   *
   * <p>Closing a non-owning wrapper does not destroy the native system or its contents.</p>
   */
  @Override
  public void destroy() {
    forces.values().forEach(Force::invalidate);
    virtualSites.values().forEach(VirtualSite::invalidate);
    forces.clear();
    virtualSites.clear();
    if (ownsHandle) {
      destroy(OpenMMNative::OpenMM_System_destroy);
    } else {
      invalidate();
    }
  }

  /**
   * Get the two particles and distance defining a constraint.
   *
   * @param index zero-based constraint index.
   * @return independent Java value with the particle indices and distance in nanometers.
   */
  public ConstraintParameters getConstraintParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment particle1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment particle2 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment distance = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_System_getConstraintParameters(
          getPointer(), index, particle1, particle2, distance);
      return new ConstraintParameters(
          particle1.get(ValueLayout.JAVA_INT, 0),
          particle2.get(ValueLayout.JAVA_INT, 0),
          distance.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the default periodic box vectors, copied into Java values.
   *
   * @return copied box vectors in nanometers.
   */
  public PeriodicBoxVectors getDefaultPeriodicBoxVectors() {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment a = new Vec3(0.0, 0.0, 0.0).toNative(arena);
      MemorySegment b = new Vec3(0.0, 0.0, 0.0).toNative(arena);
      MemorySegment c = new Vec3(0.0, 0.0, 0.0).toNative(arena);
      OpenMMNative.OpenMM_System_getDefaultPeriodicBoxVectors(getPointer(), a, b, c);
      return new PeriodicBoxVectors(Vec3.fromNative(a), Vec3.fromNative(b), Vec3.fromNative(c));
    }
  }

  /**
   * Get a force previously added through this FFM system wrapper.
   *
   * @param index zero-based force index.
   * @return tracked force wrapper, or {@code null} if the force was not added through this
   *         wrapper.
   */
  public Force getForce(int index) {
    return forces.get(OpenMMNative.OpenMM_System_getForce(getPointer(), index).address());
  }

  /**
   * @return number of constraints.
   */
  public int getNumConstraints() {
    return OpenMMNative.OpenMM_System_getNumConstraints(getPointer());
  }

  /**
   * @return number of forces.
   */
  public int getNumForces() {
    return OpenMMNative.OpenMM_System_getNumForces(getPointer());
  }

  /**
   * @return number of particles.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_System_getNumParticles(getPointer());
  }

  /**
   * Get a particle's mass in daltons (atomic mass units).
   *
   * @param index zero-based particle index.
   * @return particle mass in daltons.
   */
  public double getParticleMass(int index) {
    return OpenMMNative.OpenMM_System_getParticleMass(getPointer(), index);
  }

  /**
   * Rebind this wrapper to a native system handle without destroying the previously referenced
   * system.
   *
   * <p>The previously referenced handle is not destroyed, matching the pointer-setter behavior
   * of the JNA façade. The caller remains responsible for that displaced handle.</p>
   *
   * @param pointer replacement native system handle.
   */
  public void setPointer(MemorySegment pointer) {
    rebindPointer(pointer);
    ownsHandle = true;
  }

  /**
   * Get the virtual-site wrapper previously installed at a particle through this FFM system.
   *
   * @param index zero-based particle index, which must identify a virtual site.
   * @return tracked virtual-site wrapper, or {@code null} if it was not installed through this
   *         wrapper.
   */
  public VirtualSite getVirtualSite(int index) {
    return virtualSites.get(OpenMMNative.OpenMM_System_getVirtualSite(getPointer(), index).address());
  }

  /**
   * Determine whether a particle is a virtual site.
   *
   * @param index zero-based particle index.
   * @return {@code true} if the particle is a virtual site.
   */
  public boolean isVirtualSite(int index) {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_System_isVirtualSite(getPointer(), index));
  }

  /**
   * Remove a distance constraint.
   *
   * @param index zero-based constraint index.
   */
  public void removeConstraint(int index) {
    OpenMMNative.OpenMM_System_removeConstraint(getPointer(), index);
  }

  /**
   * Remove and destroy a force.
   *
   * <p>If this wrapper tracks the removed force, its Java wrapper is invalidated.</p>
   *
   * @param index zero-based force index.
   */
  public void removeForce(int index) {
    MemorySegment pointer = OpenMMNative.OpenMM_System_getForce(getPointer(), index);
    OpenMMNative.OpenMM_System_removeForce(getPointer(), index);
    Force force = forces.remove(pointer.address());
    if (force != null) {
      force.invalidate();
    }
  }

  /**
   * Set the particles and distance defining a constraint. Particles must be in this system, have
   * nonzero mass, and the distance is in nanometers.
   *
   * @param index     zero-based constraint index.
   * @param particle1 zero-based index of the first particle.
   * @param particle2 zero-based index of the second particle.
   * @param distance  constrained distance in nanometers.
   */
  public void setConstraintParameters(int index, int particle1, int particle2, double distance) {
    OpenMMNative.OpenMM_System_setConstraintParameters(
        getPointer(), index, particle1, particle2, distance);
  }

  /**
   * Set the default periodic box vectors in nanometers. Contexts created from this system use these
   * vectors initially.
   *
   * @param a first box vector in nm.
   * @param b second box vector in nm.
   * @param c third box vector in nm.
   */
  public void setDefaultPeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c) {
    try (Arena arena = Arena.ofConfined()) {
      OpenMMNative.OpenMM_System_setDefaultPeriodicBoxVectors(
          getPointer(), a.toNative(arena), b.toNative(arena), c.toNative(arena));
    }
  }

  /**
   * Set a particle's mass in daltons (atomic mass units).
   *
   * @param index zero-based particle index.
   * @param mass  mass in daltons.
   */
  public void setParticleMass(int index, double mass) {
    OpenMMNative.OpenMM_System_setParticleMass(getPointer(), index, mass);
  }

  /**
   * Set a particle's virtual-site definition and transfer native ownership of the site to this
   * system.
   *
   * @param index       zero-based particle index to designate as a virtual site.
   * @param virtualSite live, non-null FFM virtual site.
   */
  public void setVirtualSite(int index, VirtualSite virtualSite) {
    Objects.requireNonNull(virtualSite, "Virtual site cannot be null.");
    OpenMMNative.OpenMM_System_setVirtualSite(getPointer(), index, virtualSite.getPointer());
    virtualSites.put(virtualSite.getPointer().address(), virtualSite);
  }

  /**
   * Determine whether any force in this system uses periodic boundaries.
   *
   * <p>OpenMM reports an error if a force does not implement the periodic-boundary query.</p>
   *
   * @return {@code true} if at least one force uses periodic boundary conditions.
   */
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_System_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Immutable Java value containing the parameters of one distance constraint.
   *
   * @param particle1 zero-based first particle index.
   * @param particle2 zero-based second particle index.
   * @param distance  constrained distance in nanometers.
   */
  public record ConstraintParameters(int particle1, int particle2, double distance) {
  }

  /**
   * Immutable Java value containing the three periodic box vectors.
   *
   * @param a first box vector in nanometers.
   * @param b second box vector in nanometers.
   * @param c third box vector in nanometers.
   */
  public record PeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c) {
  }

  /**
   * Create the native system after initializing the FFM runtime.
   *
   * @return owned native system handle.
   */
  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_System_create.makeInvoker().apply();
  }
}
