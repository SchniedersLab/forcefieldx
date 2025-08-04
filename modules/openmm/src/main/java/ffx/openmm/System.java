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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.HashMap;
import java.util.Map;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addConstraint;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getConstraintParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getDefaultPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getNumConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getNumForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getParticleMass;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getVirtualSite;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_isVirtualSite;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_removeConstraint;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_removeForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setConstraintParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setDefaultPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setParticleMass;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setVirtualSite;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_usesPeriodicBoundaryConditions;

/**
 * This class represents a molecular system.  The definition of a System involves
 * four elements:
 *
 * <ul>
 * <li>The set of particles in the system</li>
 * <li>The forces acting on them</li>
 * <li>Pairs of particles whose separation should be constrained to a fixed value</li>
 * <li>For periodic systems, the dimensions of the periodic box</li>
 * </ul>
 * <p>
 * The particles and constraints are defined directly by the System object, while
 * forces are defined by objects that extend the Force class.  After creating a
 * System, call addParticle() once for each particle, addConstraint() for each constraint,
 * and addForce() for each Force.
 * <p>
 * In addition, particles may be designated as "virtual sites".  These are particles
 * whose positions are computed automatically based on the positions of other particles.
 * To define a virtual site, call setVirtualSite(), passing in a VirtualSite object
 * that defines the rules for computing its position.
 */
public class System {

  /**
   * System pointer.
   */
  private PointerByReference pointer;

  /**
   * Map to store Force instances using their PointerByReference as the key.
   */
  private final Map<PointerByReference, Force> forceMap = new HashMap<>();

  /**
   * Map to store VirtualSite instances using their PointerByReference as the key.
   */
  private final Map<PointerByReference, VirtualSite> virtualSiteMap = new HashMap<>();

  /**
   * Constructor.
   */
  public System() {
    pointer = OpenMM_System_create();
  }

  /**
   * Constructor.
   *
   * @param pointer The OpenMM System pointer.
   */
  public System(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Add a constraint to the system.
   *
   * @param particle1 The first particle.
   * @param particle2 The second particle.
   * @param distance  The distance between the particles.
   * @return The index of the constraint that was added.
   */
  public int addConstraint(int particle1, int particle2, double distance) {
    return OpenMM_System_addConstraint(pointer, particle1, particle2, distance);
  }

  /**
   * Add a force to the system.
   *
   * @param force The force to add.
   * @return The index of the force that was added.
   */
  public int addForce(Force force) {
    if (force != null) {
      int forceIndex = OpenMM_System_addForce(pointer, force.getPointer());
      force.setForceIndex(forceIndex);
      forceMap.put(force.getPointer(), force);
      return forceIndex;
    }
    return -1;
  }

  /**
   * Add a particle to the system.
   *
   * @param mass The mass of the particle.
   * @return The index of the particle that was added.
   */
  public int addParticle(double mass) {
    return OpenMM_System_addParticle(pointer, mass);
  }

  /**
   * Destroy the system.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_System_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters defining a constraint.
   *
   * @param index     The index of the constraint.
   * @param particle1 The index of the first particle involved in the constraint (output).
   * @param particle2 The index of the second particle involved in the constraint (output).
   * @param distance  The required distance between the two particles (output).
   */
  public void getConstraintParameters(int index, IntByReference particle1, IntByReference particle2, DoubleByReference distance) {
    OpenMM_System_getConstraintParameters(pointer, index, particle1, particle2, distance);
  }

  /**
   * Get the parameters defining a constraint.
   *
   * @param index     The index of the constraint.
   * @param particle1 The index of the first particle involved in the constraint (output).
   * @param particle2 The index of the second particle involved in the constraint (output).
   * @param distance  The required distance between the two particles (output).
   */
  public void getConstraintParameters(int index, IntBuffer particle1, IntBuffer particle2, DoubleBuffer distance) {
    OpenMM_System_getConstraintParameters(pointer, index, particle1, particle2, distance);
  }

  /**
   * Get the default periodic box vectors.
   *
   * @param a The first vector (output).
   * @param b The second vector (output).
   * @param c The third vector (output).
   */
  public void getDefaultPeriodicBoxVectors(OpenMM_Vec3 a, OpenMM_Vec3 b, OpenMM_Vec3 c) {
    OpenMM_System_getDefaultPeriodicBoxVectors(pointer, a, b, c);
  }

  /**
   * Get a force in the system.
   *
   * @param index The index of the force to get.
   * @return The force object.
   */
  public Force getForce(int index) {
    PointerByReference forcePointer = OpenMM_System_getForce(pointer, index);
    return forceMap.get(forcePointer);
  }

  /**
   * Get the number of constraints in the system.
   *
   * @return The number of constraints in the system.
   */
  public int getNumConstraints() {
    return OpenMM_System_getNumConstraints(pointer);
  }

  /**
   * Get the number of forces in the system.
   *
   * @return The number of forces in the system.
   */
  public int getNumForces() {
    return OpenMM_System_getNumForces(pointer);
  }

  /**
   * Get the number of particles in the system.
   *
   * @return The number of particles in the system.
   */
  public int getNumParticles() {
    return OpenMM_System_getNumParticles(pointer);
  }

  /**
   * Get the mass of a particle.
   *
   * @param index The index of the particle.
   * @return The mass of the particle.
   */
  public double getParticleMass(int index) {
    return OpenMM_System_getParticleMass(pointer, index);
  }

  /**
   * Get the pointer to the system.
   *
   * @return The pointer to the system.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Get the virtual site for a particle.
   *
   * @param index The index of the particle.
   * @return The virtual site object.
   */
  public VirtualSite getVirtualSite(int index) {
    PointerByReference virtualSitePointer = OpenMM_System_getVirtualSite(pointer, index);
    return virtualSiteMap.get(virtualSitePointer);
  }

  /**
   * Check if a particle is a virtual site.
   *
   * @param index The index of the particle.
   * @return True if the particle is a virtual site.
   */
  public boolean isVirtualSite(int index) {
    int result = OpenMM_System_isVirtualSite(pointer, index);
    return result != 0;
  }

  /**
   * Remove a constraint from the system.
   *
   * @param index The index of the constraint to remove.
   */
  public void removeConstraint(int index) {
    OpenMM_System_removeConstraint(pointer, index);
  }

  /**
   * Remove a force from the system.
   *
   * @param index The index of the force to remove.
   */
  public void removeForce(int index) {
    Force force = getForce(index);
    forceMap.remove(force.getPointer());
    OpenMM_System_removeForce(pointer, index);
  }

  /**
   * Set the parameters defining a constraint.
   *
   * @param index     The index of the constraint.
   * @param particle1 The index of the first particle involved in the constraint.
   * @param particle2 The index of the second particle involved in the constraint.
   * @param distance  The required distance between the two particles.
   */
  public void setConstraintParameters(int index, int particle1, int particle2, double distance) {
    OpenMM_System_setConstraintParameters(pointer, index, particle1, particle2, distance);
  }

  /**
   * Set the default periodic box vectors.
   *
   * @param a The first vector.
   * @param b The second vector.
   * @param c The third vector.
   */
  public void setDefaultPeriodicBoxVectors(OpenMM_Vec3 a, OpenMM_Vec3 b, OpenMM_Vec3 c) {
    OpenMM_System_setDefaultPeriodicBoxVectors(pointer, a, b, c);
  }

  /**
   * Set the mass of a particle.
   *
   * @param index The index of the particle.
   * @param mass  The mass of the particle.
   */
  public void setParticleMass(int index, double mass) {
    OpenMM_System_setParticleMass(pointer, index, mass);
  }

  /**
   * Set the pointer to the system.
   *
   * @param pointer The pointer to the system.
   */
  public void setPointer(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Set a particle to be a virtual site.
   *
   * @param index       The index of the particle.
   * @param virtualSite The virtual site object.
   */
  public void setVirtualSite(int index, VirtualSite virtualSite) {
    if (virtualSite != null) {
      virtualSiteMap.put(virtualSite.getPointer(), virtualSite);
    }
    OpenMM_System_setVirtualSite(pointer, index, virtualSite.getPointer());
  }

  /**
   * Check if the system uses periodic boundary conditions.
   *
   * @return True if the system uses periodic boundary conditions.
   */
  public boolean usesPeriodicBoundaryConditions() {
    int result = OpenMM_System_usesPeriodicBoundaryConditions(pointer);
    return result == OpenMM_True;
  }
}