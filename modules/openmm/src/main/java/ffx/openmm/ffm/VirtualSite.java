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
 * Base class for FFM-backed OpenMM virtual sites.
 *
 * <p>A virtual site's position is computed from positions of its defining particles rather than
 * integrated directly. Before installation in a {@link System}, the wrapper can release its
 * native site. Installing it in a system transfers native ownership to that system.</p>
 */
public abstract class VirtualSite extends OpenMMHandle {

  /**
   * Wrap a native virtual-site handle.
   *
   * @param pointer native virtual-site handle.
   */
  public VirtualSite(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Get the number of particles on which this virtual site depends.
   *
   * @return number of particles defining this site.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_VirtualSite_getNumParticles(getPointer());
  }

  /**
   * Get the System particle index of one particle defining this site.
   *
   * @param particle zero-based index among this site's defining particles, in the range
   *                 {@code [0, getNumParticles())}.
   * @return zero-based System particle index.
   */
  public int getParticle(int particle) {
    return OpenMMNative.OpenMM_VirtualSite_getParticle(getPointer(), particle);
  }
}
