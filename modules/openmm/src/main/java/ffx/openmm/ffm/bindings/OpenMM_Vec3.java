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
package ffx.openmm.ffm.bindings;

import java.lang.foreign.Arena;
import java.lang.foreign.GroupLayout;
import java.lang.foreign.MemoryLayout;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.SegmentAllocator;
import java.util.function.Consumer;

import static java.lang.foreign.MemoryLayout.PathElement.groupElement;
import static java.lang.foreign.ValueLayout.OfDouble;

/**
 * {@snippet lang = c:
 * struct {
 *     double x;
 *     double y;
 *     double z;
 * }
 *}
 */
public class OpenMM_Vec3 {

  OpenMM_Vec3() {
    // Should not be called directly
  }

  private static final GroupLayout $LAYOUT = MemoryLayout.structLayout(
      OpenMMNative.C_DOUBLE.withName("x"),
      OpenMMNative.C_DOUBLE.withName("y"),
      OpenMMNative.C_DOUBLE.withName("z")
  ).withName("$anon$93:9");

  /**
   * The layout of this struct
   */
  public static final GroupLayout layout() {
    return $LAYOUT;
  }

  private static final OfDouble x$LAYOUT = (OfDouble) $LAYOUT.select(groupElement("x"));

  /**
   * Layout for field:
   * {@snippet lang = c:
   * double x
   *}
   */
  public static final OfDouble x$layout() {
    return x$LAYOUT;
  }

  private static final long x$OFFSET = $LAYOUT.byteOffset(groupElement("x"));

  /**
   * Offset for field:
   * {@snippet lang = c:
   * double x
   *}
   */
  public static final long x$offset() {
    return x$OFFSET;
  }

  /**
   * Getter for field:
   * {@snippet lang = c:
   * double x
   *}
   */
  public static double x(MemorySegment struct) {
    return struct.get(x$LAYOUT, x$OFFSET);
  }

  /**
   * Setter for field:
   * {@snippet lang = c:
   * double x
   *}
   */
  public static void x(MemorySegment struct, double fieldValue) {
    struct.set(x$LAYOUT, x$OFFSET, fieldValue);
  }

  private static final OfDouble y$LAYOUT = (OfDouble) $LAYOUT.select(groupElement("y"));

  /**
   * Layout for field:
   * {@snippet lang = c:
   * double y
   *}
   */
  public static final OfDouble y$layout() {
    return y$LAYOUT;
  }

  private static final long y$OFFSET = $LAYOUT.byteOffset(groupElement("y"));

  /**
   * Offset for field:
   * {@snippet lang = c:
   * double y
   *}
   */
  public static final long y$offset() {
    return y$OFFSET;
  }

  /**
   * Getter for field:
   * {@snippet lang = c:
   * double y
   *}
   */
  public static double y(MemorySegment struct) {
    return struct.get(y$LAYOUT, y$OFFSET);
  }

  /**
   * Setter for field:
   * {@snippet lang = c:
   * double y
   *}
   */
  public static void y(MemorySegment struct, double fieldValue) {
    struct.set(y$LAYOUT, y$OFFSET, fieldValue);
  }

  private static final OfDouble z$LAYOUT = (OfDouble) $LAYOUT.select(groupElement("z"));

  /**
   * Layout for field:
   * {@snippet lang = c:
   * double z
   *}
   */
  public static final OfDouble z$layout() {
    return z$LAYOUT;
  }

  private static final long z$OFFSET = $LAYOUT.byteOffset(groupElement("z"));

  /**
   * Offset for field:
   * {@snippet lang = c:
   * double z
   *}
   */
  public static final long z$offset() {
    return z$OFFSET;
  }

  /**
   * Getter for field:
   * {@snippet lang = c:
   * double z
   *}
   */
  public static double z(MemorySegment struct) {
    return struct.get(z$LAYOUT, z$OFFSET);
  }

  /**
   * Setter for field:
   * {@snippet lang = c:
   * double z
   *}
   */
  public static void z(MemorySegment struct, double fieldValue) {
    struct.set(z$LAYOUT, z$OFFSET, fieldValue);
  }

  /**
   * Obtains a slice of {@code arrayParam} which selects the array element at {@code index}.
   * The returned segment has address {@code arrayParam.address() + index * layout().byteSize()}
   */
  public static MemorySegment asSlice(MemorySegment array, long index) {
    return array.asSlice(layout().byteSize() * index);
  }

  /**
   * The size (in bytes) of this struct
   */
  public static long sizeof() {
    return layout().byteSize();
  }

  /**
   * Allocate a segment of size {@code layout().byteSize()} using {@code allocator}
   */
  public static MemorySegment allocate(SegmentAllocator allocator) {
    return allocator.allocate(layout());
  }

  /**
   * Allocate an array of size {@code elementCount} using {@code allocator}.
   * The returned segment has size {@code elementCount * layout().byteSize()}.
   */
  public static MemorySegment allocateArray(long elementCount, SegmentAllocator allocator) {
    return allocator.allocate(MemoryLayout.sequenceLayout(elementCount, layout()));
  }

  /**
   * Reinterprets {@code addr} using target {@code arena} and {@code cleanupAction} (if any).
   * The returned segment has size {@code layout().byteSize()}
   */
  public static MemorySegment reinterpret(MemorySegment addr, Arena arena, Consumer<MemorySegment> cleanup) {
    return reinterpret(addr, 1, arena, cleanup);
  }

  /**
   * Reinterprets {@code addr} using target {@code arena} and {@code cleanupAction} (if any).
   * The returned segment has size {@code elementCount * layout().byteSize()}
   */
  public static MemorySegment reinterpret(MemorySegment addr, long elementCount, Arena arena, Consumer<MemorySegment> cleanup) {
    return addr.reinterpret(layout().byteSize() * elementCount, arena, cleanup);
  }
}

