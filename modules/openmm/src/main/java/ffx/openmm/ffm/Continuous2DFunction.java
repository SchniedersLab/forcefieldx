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
 * A two-dimensional tabulated function evaluated by natural cubic spline interpolation.
 *
 * <p>The table has {@code xsize} uniformly spaced x coordinates over {@code [xmin,xmax]} and
 * {@code ysize} uniformly spaced y coordinates over {@code [ymin,ymax]}. Values are stored with
 * x varying fastest: {@code values[i + xsize*j] = f(x_i,y_j)}. OpenMM documents the function as
 * zero when either coordinate is outside its specified range; the constructor's periodicity flag
 * is also forwarded to native OpenMM. This wrapper does not convert units.
 */
public class Continuous2DFunction extends TabulatedFunction {

  /**
   * Create a continuous function from values on uniformly spaced x and y grids.
   * The values array is copied to native storage and must contain one value per grid point
   * ({@code xsize*ysize} values).
   *
   * @param values table values in x-fastest order, {@code values[i + xsize*j] = f(x_i,y_j)}
   * @param xsize number of table points along x
   * @param ysize number of table points along y
   * @param xmin x coordinate of the first x table point
   * @param xmax x coordinate of the last x table point
   * @param ymin y coordinate of the first y table point
   * @param ymax y coordinate of the last y table point
   * @param periodic whether the native function treats the interpolated function as periodic
   */
  public Continuous2DFunction(double[] values, int xsize, int ysize, double xmin, double xmax,
                              double ymin, double ymax, boolean periodic) {
    super(create(values, xsize, ysize, xmin, xmax, ymin, ymax, periodic));
  }

  /** Destroy the native continuous function. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_Continuous2DFunction_destroy);
  }

  /**
   * Get the current table values, dimensions, and bounds as Java-owned copies.
   *
   * @return copied parameters; values use x-fastest order,
   *     {@code values[i + xsize*j] = f(x_i,y_j)}
   */
  public Parameters getFunctionParameters() {
    try (Arena arena = Arena.ofConfined(); DoubleArray values = new DoubleArray(0)) {
      MemorySegment xs = arena.allocate(ValueLayout.JAVA_INT), ys = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment xmin = arena.allocate(ValueLayout.JAVA_DOUBLE), xmax = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment ymin = arena.allocate(ValueLayout.JAVA_DOUBLE), ymax = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_Continuous2DFunction_getFunctionParameters(
          getPointer(), xs, ys, values.getPointer(), xmin, xmax, ymin, ymax);
      return new Parameters(copy(values), xs.get(ValueLayout.JAVA_INT, 0), ys.get(ValueLayout.JAVA_INT, 0),
          xmin.get(ValueLayout.JAVA_DOUBLE, 0), xmax.get(ValueLayout.JAVA_DOUBLE, 0),
          ymin.get(ValueLayout.JAVA_DOUBLE, 0), ymax.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace the table values, dimensions, and bounds. Values are copied to native storage; the
   * existing periodicity setting is not changed by this method. The array must contain
   * {@code xsize*ysize} values in x-fastest order.
   *
   * @param values table values with x varying fastest
   * @param xsize number of table points along x
   * @param ysize number of table points along y
   * @param xmin x coordinate of the first x table point
   * @param xmax x coordinate of the last x table point
   * @param ymin y coordinate of the first y table point
   * @param ymax y coordinate of the last y table point
   */
  public void setFunctionParameters(double[] values, int xsize, int ysize, double xmin,
                                    double xmax, double ymin, double ymax) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMNative.OpenMM_Continuous2DFunction_setFunctionParameters(
          getPointer(), xsize, ysize, nativeValues.getPointer(), xmin, xmax, ymin, ymax);
    }
  }

  /**
   * Copied parameters defining a two-dimensional continuous function.
   *
   * @param values Java-owned copy of the table values, with x varying fastest
   * @param xsize number of table points along x
   * @param ysize number of table points along y
   * @param xmin x coordinate of the first x table point
   * @param xmax x coordinate of the last x table point
   * @param ymin y coordinate of the first y table point
   * @param ymax y coordinate of the last y table point
   */
  public record Parameters(double[] values, int xsize, int ysize, double xmin, double xmax,
                           double ymin, double ymax) {
  }

  private static MemorySegment create(double[] values, int xsize, int ysize, double xmin, double xmax,
                                      double ymin, double ymax, boolean periodic) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_Continuous2DFunction_create(xsize, ysize, nativeValues.getPointer(),
          xmin, xmax, ymin, ymax, OpenMMBooleans.toNative(periodic));
    }
  }

  private static DoubleArray toNative(double[] values) {
    DoubleArray array = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) array.set(i, values[i]);
    return array;
  }

  private static double[] copy(DoubleArray values) {
    double[] result = new double[values.getSize()];
    for (int i = 0; i < result.length; i++) result[i] = values.get(i);
    return result;
  }
}
