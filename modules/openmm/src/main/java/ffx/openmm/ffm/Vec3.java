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

import ffx.openmm.ffm.bindings.OpenMM_Vec3;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;

/**
 * Immutable three-component OpenMM vector value.
 *
 * <p>This Java record owns no native memory. Its component units depend on the OpenMM operation
 * using the vector.</p>
 *
 * @param x x component.
 * @param y y component.
 * @param z z component.
 */
public record Vec3(double x, double y, double z) {

  /**
   * Copy an OpenMM {@code Vec3} struct into a Java value.
   *
   * <p>The struct may be passed directly or through a pointer to it. The returned record does not
   * retain the native segment.</p>
   *
   * @param vector live native vector struct or pointer to one, readable for the duration of this
   *               call.
   * @return independent Java value with the copied components.
   */
  public static Vec3 fromNative(MemorySegment vector) {
    MemorySegment struct = vector.reinterpret(OpenMM_Vec3.layout().byteSize());
    return new Vec3(OpenMM_Vec3.x(struct), OpenMM_Vec3.y(struct), OpenMM_Vec3.z(struct));
  }

  /**
   * Allocate and populate an OpenMM {@code Vec3} struct in the supplied arena.
   *
   * <p>The returned segment is valid only while {@code arena} remains open.</p>
   *
   * @param arena arena that owns the native struct.
   * @return native vector struct allocated in {@code arena}.
   */
  public MemorySegment toNative(Arena arena) {
    MemorySegment struct = arena.allocate(OpenMM_Vec3.layout());
    OpenMM_Vec3.x(struct, x);
    OpenMM_Vec3.y(struct, y);
    OpenMM_Vec3.z(struct, z);
    return struct;
  }
}
