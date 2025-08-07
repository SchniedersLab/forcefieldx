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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VirtualSite_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VirtualSite_getParticle;

/**
 * A VirtualSite describes the rules for computing a particle's position based on
 * other particles. This is an abstract class. Subclasses define particular rules.
 * To define a virtual site, create an instance of a VirtualSite subclass and then
 * call setVirtualSite() on the System.
 */
public abstract class VirtualSite {

  /**
   * The pointer is allocated and deallocated by classes that extend VirtualSite.
   */
  protected PointerByReference pointer;

  /**
   * Create a VirtualSite from an existing pointer.
   *
   * @param pointer The pointer to the OpenMM VirtualSite.
   */
  public VirtualSite(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Destroy the virtual site.
   */
  public abstract void destroy();

  /**
   * Get the number of particles this virtual site depends on.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_VirtualSite_getNumParticles(pointer);
  }

  /**
   * Get the index of a particle this virtual site depends on.
   *
   * @param particle the particle to get (between 0 and getNumParticles())
   * @return the index of the particle in the System
   */
  public int getParticle(int particle) {
    return OpenMM_VirtualSite_getParticle(pointer, particle);
  }

  /**
   * Get the pointer to the OpenMM VirtualSite.
   *
   * @return The pointer to the OpenMM VirtualSite.
   */
  public PointerByReference getPointer() {
    return pointer;
  }
}