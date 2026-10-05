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
 * Base class for FFM-backed OpenMM forces.
 *
 * <p>A force contributes energy and forces to an OpenMM {@link System}. Before it is added to a
 * system, the wrapper is responsible for releasing the native force. Adding it transfers native
 * ownership to the system; the system tracks this wrapper and invalidates it when the system
 * removes or releases the force. Do not destroy a force while it is attached to a system.</p>
 */
public abstract class Force extends OpenMMHandle {

  private int forceIndex = -1;

  /**
   * Wrap a native force handle.
   *
   * <p>The wrapper may release the force until ownership is transferred to a {@link System}.</p>
   *
   * @param pointer native force handle.
   */
  public Force(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Get the force group to which this force belongs. OpenMM supports group indices from 0 through
   * 31.
   *
   * @return force-group index.
   */
  public int getForceGroup() {
    return OpenMMNative.OpenMM_Force_getForceGroup(getPointer());
  }

  /**
   * Get the index most recently recorded by this wrapper when the force was added to an FFM
   * {@link System}.
   *
   * <p>This is Java-side bookkeeping, not a live native query. It is not adjusted if forces are
   * subsequently removed or reordered through other native code.</p>
   *
   * @return last recorded system force index, or {@code -1} before it is added through this
   * wrapper's system API.
   */
  public int getForceIndex() {
    return forceIndex;
  }

  /**
   * Get the user-assigned force name.
   *
   * @return force name.
   */
  public String getName() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_Force_getName(getPointer()));
  }

  /**
   * Assign this force to a force group. OpenMM supports group indices from 0 through 31.
   *
   * @param forceGroup group index in the range 0 through 31.
   */
  public void setForceGroup(int forceGroup) {
    OpenMMNative.OpenMM_Force_setForceGroup(getPointer(), forceGroup);
  }

  /**
   * Record the index assigned when this force is added to a system. This changes only Java-side
   * bookkeeping and does not modify the native system.
   *
   * @param forceIndex system force index recorded by the wrapper.
   */
  public final void setForceIndex(int forceIndex) {
    this.forceIndex = forceIndex;
  }

  /**
   * Set the user-assigned force name. The Java string is encoded as temporary UTF-8 for the native
   * call; OpenMM stores its own copy.
   *
   * @param name non-null force name.
   */
  public void setName(String name) {
    OpenMMStrings.withUtf8String(name, string -> {
      OpenMMNative.OpenMM_Force_setName(getPointer(), string);
    });
  }

  /**
   * Determine whether this force's implementation uses periodic boundary conditions.
   *
   * @return {@code true} when the force depends on periodic boundary conditions.
   */
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_Force_usesPeriodicBoundaryConditions(getPointer()));
  }
}
