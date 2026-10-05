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

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Converts between Java strings and NUL-terminated UTF-8 strings used by the OpenMM C API.
 *
 * <p>Native pointers accepted by {@link #copy(MemorySegment)} are borrowed and must remain valid
 * until copying completes. Strings passed through the {@code withUtf8String} methods are allocated
 * in a confined arena and are valid only for the duration of the callback. Since OpenMM consumes
 * NUL-terminated C strings, an embedded NUL in a Java input string terminates the native value at
 * that character.</p>
 */
public final class OpenMMStrings {

  private OpenMMStrings() {
  }

  /**
   * Copy a NUL-terminated UTF-8 string into Java.
   *
   * <p>This method does not release or retain the native string. The caller must ensure that its
   * pointer remains valid and that a terminating NUL byte is accessible.</p>
   *
   * @param string borrowed native character pointer; {@link MemorySegment#NULL} is allowed.
   * @return copied Java string, or {@code null} for a null native pointer.
   */
  public static String copy(MemorySegment string) {
    if (string.address() == 0) {
      return null;
    }
    return string.reinterpret(Long.MAX_VALUE).getString(0);
  }

  /**
   * Allocate a temporary NUL-terminated UTF-8 string and pass it to a callback.
   *
   * <p>The segment is confined to this call and must not be retained by the callback.</p>
   *
   * @param value  non-null Java string to encode.
   * @param action non-null callback that uses the temporary string before returning.
   * @throws NullPointerException if {@code value} or {@code action} is null.
   */
  public static void withUtf8String(String value, Consumer<MemorySegment> action) {
    withUtf8StringResult(value, string -> {
      action.accept(string);
      return null;
    });
  }

  /**
   * Allocate a temporary NUL-terminated UTF-8 string, pass it to a callback, and return the
   * callback's result.
   *
   * <p>The segment is confined to this call and must not escape through the callback result.</p>
   *
   * @param value  non-null Java string to encode.
   * @param action non-null callback that uses the temporary string before returning.
   * @param <T>    callback result type.
   * @return callback result.
   * @throws NullPointerException if {@code value} or {@code action} is null.
   */
  public static <T> T withUtf8StringResult(String value, Function<MemorySegment, T> action) {
    Objects.requireNonNull(value, "String cannot be null.");
    Objects.requireNonNull(action, "String action cannot be null.");
    try (Arena arena = Arena.ofConfined()) {
      return action.apply(arena.allocateFrom(value));
    }
  }
}
