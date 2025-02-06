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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addConstraint;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getNumConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_removeForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setDefaultPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setParticleMass;

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
 *
 * The particles and constraints are defined directly by the System object, while
 * forces are defined by objects that extend the Force class.  After creating a
 * System, call addParticle() once for each particle, addConstraint() for each constraint,
 * and addForce() for each Force.
 *
 * In addition, particles may be designated as "virtual sites".  These are particles
 * whose positions are computed automatically based on the positions of other particles.
 * To define a virtual site, call setVirtualSite(), passing in a VirtualSite object
 * that defines the rules for computing its position.
 */
public class System {

  /**
   * Context pointer.
   */
  PointerByReference pointer;

  /**
   * Constructor.
   */
  public System() {
    pointer = OpenMM_System_create();
  }

  /**
   * Add a constraint to the system.
   *
   * @param particle1 The first particle.
   * @param particle2 The second particle.
   * @param distance  The distance between the particles.
   */
  public void addConstraint(int particle1, int particle2, double distance) {
    OpenMM_System_addConstraint(pointer, particle1, particle2, distance);
  }

  /**
   * Add a force to the system.
   *
   * @param force The force to add.
   */
  public void addForce(Force force) {
    if (force != null) {
      int forceIndex = OpenMM_System_addForce(pointer, force.getPointer());
      force.setForceIndex(forceIndex);
    }
  }

  /**
   * Add a particle to the system.
   *
   * @param mass The mass of the particle.
   */
  public void addParticle(double mass) {
    OpenMM_System_addParticle(pointer, mass);
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
   * Remove a force from the system.
   *
   * @param index The index of the force to remove.
   */
  public void removeForce(int index) {
    OpenMM_System_removeForce(pointer, index);
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
   * Get the pointer to the system.
   *
   * @return The pointer to the system.
   */
  public PointerByReference getPointer() {
    return pointer;
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

}
