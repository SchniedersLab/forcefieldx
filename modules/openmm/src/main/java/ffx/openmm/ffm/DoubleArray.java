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
 * Owns a native OpenMM array of {@code double} values.
 */
public class DoubleArray extends OpenMMHandle {

  /**
   * Create a native double array.
   *
   * @param size initial number of values; must be nonnegative.
   */
  public DoubleArray(int size) {
    super(create(size));
  }

  /**
   * Append a value to the end of this array.
   *
   * @param value value to append.
   */
  public void append(double value) {
    OpenMMNative.OpenMM_DoubleArray_append(getPointer(), value);
  }

  /**
   * Release this native array. Repeated calls have no effect.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_DoubleArray_destroy);
  }

  /**
   * Get one value from this array.
   *
   * @param index zero-based array index in the range {@code [0, getSize())}.
   * @return value at the index.
   */
  public double get(int index) {
    return OpenMMNative.OpenMM_DoubleArray_get(getPointer(), index);
  }

  /**
   * Get the number of values currently in this array.
   *
   * @return number of values.
   */
  public int getSize() {
    return OpenMMNative.OpenMM_DoubleArray_getSize(getPointer());
  }

  /**
   * Set the number of values in this array.
   *
   * @param size new number of values; must be nonnegative.
   */
  public void resize(int size) {
    OpenMMNative.OpenMM_DoubleArray_resize(getPointer(), size);
  }

  /**
   * Replace one value in this array.
   *
   * @param index zero-based array index in the range {@code [0, getSize())}.
   * @param value value to store.
   */
  public void set(int index, double value) {
    OpenMMNative.OpenMM_DoubleArray_set(getPointer(), index, value);
  }

  private static MemorySegment create(int size) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_DoubleArray_create(size);
  }
}
