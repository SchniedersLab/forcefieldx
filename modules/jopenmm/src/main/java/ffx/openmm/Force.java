// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_getForceGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_getName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setForceGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_usesPeriodicBoundaryConditions;

/**
 * OpenMM Force.
 */
public abstract class Force {

  /**
   * The forcePointer is allocated and deallocated by classes that extend OpenMMForce.
   */
  protected PointerByReference pointer = null;

  /**
   * The forceIndex is returned by OpenMMSystem.addForce and is used to remove the force.
   */
  private int forceIndex = -1;

  /**
   * Get the pointer to the OpenMM Force.
   *
   * @return The pointer to the OpenMM Force.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Set the force group.
   *
   * @param forceGroup The force group.
   */
  public void setForceGroup(int forceGroup) {
    OpenMM_Force_setForceGroup(pointer, forceGroup);
  }

  /**
   * Get the force group.
   *
   * @return The force group.
   */
  public int getForceGroup() {
    return OpenMM_Force_getForceGroup(pointer);
  }

  /**
   * Set the name of the force.
   *
   * @param name The name of the force.
   */
  public void setName(String name) {
    OpenMM_Force_setName(pointer, name);
  }

  /**
   * Get the name of the force.
   *
   * @return The name of the force.
   */
  public String getName() {
    Pointer pointer = OpenMM_Force_getName(this.pointer);
    if (pointer == null) {
      return null;
    }
    return pointer.getString(0);
  }

  /**
   * Set the force index.
   *
   * @param forceIndex The force index.
   */
  public void setForceIndex(int forceIndex) {
    this.forceIndex = forceIndex;
  }

  /**
   * Get the force index.
   *
   * @return The force index.
   */
  public int getForceIndex() {
    return forceIndex;
  }

  /**
   * Check if the force use periodic boundary conditions. This is a virtual method
   * that must be implemented by classes that extend OpenMMForce.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_Force_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }

}
