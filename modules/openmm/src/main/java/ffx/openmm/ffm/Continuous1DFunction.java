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
import java.lang.foreign.ValueLayout;

/**
 * A one-dimensional tabulated function evaluated by natural cubic spline interpolation.
 *
 * <p>The table values are uniformly spaced between {@code min} and {@code max}; OpenMM documents
 * the function as zero outside this range. The constructor's periodicity flag is forwarded to
 * native OpenMM. Values and bounds use the units of the tabulated expression; this wrapper
 * performs no unit conversion.
 */
public class Continuous1DFunction extends TabulatedFunction {

  /**
   * Create a continuous function from values at uniformly spaced points over {@code [min,max]}.
   * A natural cubic spline interpolates between adjacent values. Unless {@code periodic} is true,
   * the function is zero outside this interval. The input array is copied to native storage.
   *
   * @param values uniformly spaced function values, with the first at {@code min} and last at
   *     {@code max}
   * @param min x coordinate of the first table element
   * @param max x coordinate of the last table element
   * @param periodic whether the native function treats the interpolated function as periodic
   */
  public Continuous1DFunction(double[] values, double min, double max, boolean periodic) {
    super(create(values, min, max, periodic));
  }

  /** Destroy the native continuous function. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_Continuous1DFunction_destroy);
  }

  /**
   * Get the current table values and bounds as Java-owned copies.
   *
   * @return copied values and the coordinates of the first and last table elements
   */
  public Parameters getFunctionParameters() {
    try (Arena arena = Arena.ofConfined(); DoubleArray values = new DoubleArray(0)) {
      MemorySegment min = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment max = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_Continuous1DFunction_getFunctionParameters(
          getPointer(), values.getPointer(), min, max);
      return new Parameters(copy(values), min.get(ValueLayout.JAVA_DOUBLE, 0),
          max.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace the table values and bounds. Values are copied to native storage; the existing
   * periodicity setting is not changed by this method.
   *
   * @param values uniformly spaced function values, with the first at {@code min} and last at
   *     {@code max}
   * @param min x coordinate of the first table element
   * @param max x coordinate of the last table element
   */
  public void setFunctionParameters(double[] values, double min, double max) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMNative.OpenMM_Continuous1DFunction_setFunctionParameters(
          getPointer(), nativeValues.getPointer(), min, max);
    }
  }

  /**
   * Copied parameters defining a one-dimensional continuous function.
   *
   * @param values Java-owned copy of the uniformly spaced function values
   * @param min x coordinate of the first table element
   * @param max x coordinate of the last table element
   */
  public record Parameters(double[] values, double min, double max) {
  }

  private static MemorySegment create(double[] values, double min, double max, boolean periodic) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_Continuous1DFunction_create(
          nativeValues.getPointer(), min, max, OpenMMBooleans.toNative(periodic));
    }
  }

  private static DoubleArray toNative(double[] values) {
    DoubleArray array = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) {
      array.set(i, values[i]);
    }
    return array;
  }

  private static double[] copy(DoubleArray values) {
    double[] result = new double[values.getSize()];
    for (int i = 0; i < result.length; i++) result[i] = values.get(i);
    return result;
  }
}
