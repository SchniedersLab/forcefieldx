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
 * Owns a native OpenMM array of UTF-8 strings.
 *
 * <p>Strings passed to native calls are encoded as temporary NUL-terminated UTF-8 data. Strings
 * returned from the array are copied into Java strings before the native call's borrowed result is
 * released. Embedded NUL characters in input strings terminate the value seen by OpenMM.</p>
 */
public class StringArray extends OpenMMHandle {

  /**
   * Create a native string array.
   *
   * @param size initial number of strings; must be nonnegative.
   */
  public StringArray(int size) {
    super(create(size));
  }

  /**
   * Wrap an existing native string array and assume responsibility for releasing it.
   *
   * @param pointer native string-array handle whose lifetime is transferred to this wrapper.
   */
  public StringArray(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Append a UTF-8 string to the end of this array.
   *
   * @param value non-null string to append.
   */
  public void append(String value) {
    OpenMMStrings.withUtf8String(value, string -> {
      OpenMMNative.OpenMM_StringArray_append(getPointer(), string);
    });
  }

  /**
   * Release this native array. Repeated calls have no effect.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_StringArray_destroy);
  }

  /**
   * Copy one string from this array into Java.
   *
   * @param index zero-based array index.
   * @return copied string, or {@code null} when {@code index} is outside the current array.
   */
  public String get(int index) {
    if (index < 0 || index >= getSize()) {
      return null;
    }
    return OpenMMStrings.copy(OpenMMNative.OpenMM_StringArray_get(getPointer(), index));
  }

  /**
   * Get the number of strings currently in this array.
   *
   * @return number of strings.
   */
  public int getSize() {
    return OpenMMNative.OpenMM_StringArray_getSize(getPointer());
  }

  /**
   * Set the number of strings in this array.
   *
   * @param size new number of strings; must be nonnegative.
   */
  public void resize(int size) {
    OpenMMNative.OpenMM_StringArray_resize(getPointer(), size);
  }

  /**
   * Replace one string in this array with a UTF-8 value.
   *
   * @param index zero-based array index in the range {@code [0, getSize())}.
   * @param value non-null string to store.
   */
  public void set(int index, String value) {
    OpenMMStrings.withUtf8String(value, string -> {
      OpenMMNative.OpenMM_StringArray_set(getPointer(), index, string);
    });
  }

  private static MemorySegment create(int size) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_StringArray_create(size);
  }
}
