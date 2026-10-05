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

import static java.lang.foreign.ValueLayout.JAVA_INT;

/**
 * Owns a native OpenMM array of pairs of particle indices.
 *
 * <p>Indices passed to or returned by this array are zero-based System particle indices. The
 * wrapper owns the native array; {@link #close()} and {@link #destroy()} release it.</p>
 */
public class BondArray extends OpenMMHandle {

  /**
   * Create a bond array with the specified initial size.
   *
   * @param size Initial number of bonds; must be nonnegative.
   */
  public BondArray(int size) {
    super(create(size));
  }

  /**
   * An immutable bond value containing two zero-based particle indices.
   *
   * @param particle1 zero-based index of the first particle.
   * @param particle2 zero-based index of the second particle.
   */
  public record Bond(int particle1, int particle2) {
  }

  /**
   * Append a bond to the end of this array.
   *
   * @param particle1 zero-based index of the first particle.
   * @param particle2 zero-based index of the second particle.
   */
  public void append(int particle1, int particle2) {
    OpenMMNative.OpenMM_BondArray_append(requirePointer(), particle1, particle2);
  }

  /**
   * Release this native array. Repeated calls have no effect.
   */
  public void destroy() {
    destroy(OpenMMNative::OpenMM_BondArray_destroy);
  }

  /**
   * Copy one bond from this array into a Java value.
   *
   * @param index zero-based array index in the range {@code [0, getSize())}.
   * @return copied pair of particle indices.
   */
  public Bond get(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment particle1 = arena.allocate(JAVA_INT);
      MemorySegment particle2 = arena.allocate(JAVA_INT);
      OpenMMNative.OpenMM_BondArray_get(getPointer(), index, particle1, particle2);
      return new Bond(particle1.get(JAVA_INT, 0), particle2.get(JAVA_INT, 0));
    }
  }

  /**
   * Get the number of bonds currently in this array.
   *
   * @return number of bonds.
   */
  public int getSize() {
    return OpenMMNative.OpenMM_BondArray_getSize(getPointer());
  }

  /**
   * Set the number of bonds in this array.
   *
   * @param size new number of bonds; must be nonnegative.
   */
  public void resize(int size) {
    OpenMMNative.OpenMM_BondArray_resize(getPointer(), size);
  }

  /**
   * Replace one bond in this array.
   *
   * @param index     zero-based array index in the range {@code [0, getSize())}.
   * @param particle1 zero-based index of the first particle.
   * @param particle2 zero-based index of the second particle.
   */
  public void set(int index, int particle1, int particle2) {
    OpenMMNative.OpenMM_BondArray_set(getPointer(), index, particle1, particle2);
  }

  private static MemorySegment create(int size) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_BondArray_create(size);
  }
}
