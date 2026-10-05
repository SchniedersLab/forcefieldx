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

/**
 * Owns a native OpenMM array of three-component vectors.
 *
 * <p>OpenMM does not associate units with these generic vectors; units are determined by the
 * operation or data being represented.</p>
 */
public class Vec3Array extends OpenMMHandle {

  /**
   * Create a native vector array.
   *
   * @param size initial number of vectors; must be nonnegative.
   */
  public Vec3Array(int size) {
    super(create(size));
  }

  /**
   * Wrap an existing native vector array and assume responsibility for releasing it.
   *
   * @param pointer native vector-array handle whose lifetime is transferred to this wrapper.
   */
  public Vec3Array(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Append a copy of a vector to the end of this array.
   *
   * @param value vector to append.
   */
  public void append(Vec3 value) {
    try (Arena arena = Arena.ofConfined()) {
      OpenMMNative.OpenMM_Vec3Array_append(getPointer(), value.toNative(arena));
    }
  }

  /**
   * Release this native array. Repeated calls have no effect.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_Vec3Array_destroy);
  }

  /**
   * Copy one vector from this array into a Java record.
   *
   * @param index zero-based array index in the range {@code [0, getSize())}.
   * @return copied vector.
   */
  public Vec3 get(int index) {
    return Vec3.fromNative(OpenMMNative.OpenMM_Vec3Array_get(getPointer(), index));
  }

  /**
   * Copy this vector array into packed {@code x,y,z} components, with three consecutive values per
   * vector.
   *
   * @return independent array of length {@code 3 * getSize()}.
   */
  public double[] getArray() {
    int size = getSize();
    double[] values = new double[size * 3];
    for (int index = 0; index < size; index++) {
      Vec3 vector = get(index);
      int offset = index * 3;
      values[offset] = vector.x();
      values[offset + 1] = vector.y();
      values[offset + 2] = vector.z();
    }
    return values;
  }

  /**
   * Get the number of vectors currently in this array.
   *
   * @return number of vectors.
   */
  public int getSize() {
    return OpenMMNative.OpenMM_Vec3Array_getSize(getPointer());
  }

  /**
   * Set the number of vectors in this array.
   *
   * @param size new number of vectors; must be nonnegative.
   */
  public void resize(int size) {
    OpenMMNative.OpenMM_Vec3Array_resize(getPointer(), size);
  }

  /**
   * Replace one vector in this array with a copy of the supplied value.
   *
   * @param index zero-based array index in the range {@code [0, getSize())}.
   * @param value vector to store.
   */
  public void set(int index, Vec3 value) {
    try (Arena arena = Arena.ofConfined()) {
      OpenMMNative.OpenMM_Vec3Array_set(getPointer(), index, value.toNative(arena));
    }
  }

  /**
   * Convert packed {@code x,y,z} components to a newly allocated OpenMM vector array.
   *
   * @param values non-null packed vector components, with three values per vector.
   * @return new owned vector array; the caller must close or destroy it.
   * @throws IllegalArgumentException if {@code values.length} is not divisible by three.
   * @throws NullPointerException if {@code values} is null.
   */
  public static Vec3Array toVec3Array(double[] values) {
    if (values.length % 3 != 0) {
      throw new IllegalArgumentException("Vec3 array length must be divisible by three.");
    }

    Vec3Array vectors = new Vec3Array(0);
    for (int index = 0; index < values.length; index += 3) {
      vectors.append(new Vec3(values[index], values[index + 1], values[index + 2]));
    }
    return vectors;
  }

  private static MemorySegment create(int size) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_Vec3Array_create(size);
  }
}
