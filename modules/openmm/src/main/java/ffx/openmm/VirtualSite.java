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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VirtualSite_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VirtualSite_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VirtualSite_getParticle;

/**
 * A VirtualSite describes a method for computing a particle's position based on other particles.
 * This is an abstract class. Subclasses define particular algorithms.
 * <p>
 * Some Force objects make use of virtual sites. They are used for two purposes:
 * <ul>
 * <li>They can serve as the location of a force center, such as a point charge.</li>
 * <li>They can serve as extra particles with well defined positions that can be used in
 * calculating forces applied to other particles.</li>
 * </ul>
 * <p>
 * In either case, the position of a virtual site is calculated based on the positions of
 * other particles. The calculation is done outside the integration loop, so the virtual
 * site locations are always consistent with the current particle positions.
 * <p>
 * When a Context is created, virtual sites are automatically added as extra particles.
 * Their masses and charges are set to zero. Forces that act on virtual sites can be
 * applied by specifying the virtual site's particle index, which equals the number of
 * ordinary particles plus the virtual site's index.
 */
public abstract class VirtualSite {

  /**
   * The pointer is allocated and deallocated by classes that extend VirtualSite.
   */
  protected PointerByReference pointer = null;

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
  public void destroy() {
    if (pointer != null) {
      OpenMM_VirtualSite_destroy(pointer);
      pointer = null;
    }
  }

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
   * @param particle The particle whose index to get.
   * @return The index of the particle.
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