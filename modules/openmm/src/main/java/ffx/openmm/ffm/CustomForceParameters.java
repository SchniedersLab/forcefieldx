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

/**
 * Shared implementation helpers for custom-force parameter vectors and context updates.
 *
 * <p>These package-private helpers have no direct JNA façade counterpart: they describe the FFM
 * adaptation only. Array conversion creates independently owned native arrays; copying produces
 * independent Java values; context forwarding is intentionally conditional on a live context
 * pointer.
 */
final class CustomForceParameters {

  /** Prevent construction of this static helper class. */
  private CustomForceParameters() {}

  /**
   * Allocate a native double array and copy every Java value into it.
   *
   * @param parameters source Java values; must not be {@code null}
   * @return newly allocated native array owned by the caller, who must close it
   */
  static DoubleArray toNative(double[] parameters) {
    DoubleArray result = new DoubleArray(parameters.length);
    for (int i = 0; i < parameters.length; i++) {
      result.set(i, parameters[i]);
    }
    return result;
  }

  /**
   * Copy every value in a native array to a newly allocated Java array.
   *
   * @param parameters native array to read; it remains owned by its caller
   * @return independent Java copy
   */
  static double[] copy(DoubleArray parameters) {
    double[] result = new double[parameters.getSize()];
    for (int i = 0; i < result.length; i++) {
      result[i] = parameters.get(i);
    }
    return result;
  }

  /**
   * Invoke an OpenMM context update only when the context has a native handle.
   *
   * <p>A context without a handle is silently ignored, matching the legacy façade's guard.
   *
   * @param context context to inspect
   * @param update operation receiving the context's borrowed native handle
   */
  static void updateContext(Context context, java.util.function.Consumer<java.lang.foreign.MemorySegment> update) {
    if (context.hasContextPointer()) {
      update.accept(context.getPointer());
    }
  }
}
