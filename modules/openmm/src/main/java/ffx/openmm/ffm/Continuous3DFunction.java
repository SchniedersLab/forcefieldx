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
 * A three-dimensional tabulated function evaluated by natural cubic spline interpolation.
 *
 * <p>The table has uniformly spaced coordinates over each inclusive bound pair. Values are
 * stored with x varying fastest, then y, then z:
 * {@code values[i + xsize*j + xsize*ysize*k] = f(x_i,y_j,z_k)}. OpenMM documents the function as
 * zero when any coordinate is outside its specified range; the constructor's periodicity flag
 * is also forwarded to native OpenMM. This wrapper does not convert units.
 */
public class Continuous3DFunction extends TabulatedFunction {

  /**
   * Create a continuous function from values on uniformly spaced x, y, and z grids.
   * The values array is copied to native storage and must contain
   * {@code xsize*ysize*zsize} values in x-fastest order.
   *
   * @param values table values in x-fastest order,
   *     {@code values[i + xsize*j + xsize*ysize*k] = f(x_i,y_j,z_k)}
   * @param xsize number of table points along x
   * @param ysize number of table points along y
   * @param zsize number of table points along z
   * @param xmin x coordinate of the first x table point
   * @param xmax x coordinate of the last x table point
   * @param ymin y coordinate of the first y table point
   * @param ymax y coordinate of the last y table point
   * @param zmin z coordinate of the first z table point
   * @param zmax z coordinate of the last z table point
   * @param periodic whether the native function treats the interpolated function as periodic
   */
  public Continuous3DFunction(double[] values, int xsize, int ysize, int zsize, double xmin,
                              double xmax, double ymin, double ymax, double zmin, double zmax, boolean periodic) {
    super(create(values, xsize, ysize, zsize, xmin, xmax, ymin, ymax, zmin, zmax, periodic));
  }

  /** Destroy the native continuous function. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_Continuous3DFunction_destroy);
  }

  /**
   * Get the current table values, dimensions, and bounds as Java-owned copies.
   *
   * @return copied parameters; values use x-fastest order,
   *     {@code values[i + xsize*j + xsize*ysize*k] = f(x_i,y_j,z_k)}
   */
  public Parameters getFunctionParameters() {
    try (Arena arena = Arena.ofConfined(); DoubleArray values = new DoubleArray(0)) {
      MemorySegment xs = arena.allocate(ValueLayout.JAVA_INT), ys = arena.allocate(ValueLayout.JAVA_INT),
          zs = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment xmin = arena.allocate(ValueLayout.JAVA_DOUBLE), xmax = arena.allocate(ValueLayout.JAVA_DOUBLE),
          ymin = arena.allocate(ValueLayout.JAVA_DOUBLE), ymax = arena.allocate(ValueLayout.JAVA_DOUBLE),
          zmin = arena.allocate(ValueLayout.JAVA_DOUBLE), zmax = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_Continuous3DFunction_getFunctionParameters(
          getPointer(), xs, ys, zs, values.getPointer(), xmin, xmax, ymin, ymax, zmin, zmax);
      return new Parameters(copy(values), xs.get(ValueLayout.JAVA_INT, 0), ys.get(ValueLayout.JAVA_INT, 0),
          zs.get(ValueLayout.JAVA_INT, 0), xmin.get(ValueLayout.JAVA_DOUBLE, 0), xmax.get(ValueLayout.JAVA_DOUBLE, 0),
          ymin.get(ValueLayout.JAVA_DOUBLE, 0), ymax.get(ValueLayout.JAVA_DOUBLE, 0),
          zmin.get(ValueLayout.JAVA_DOUBLE, 0), zmax.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace the table values, dimensions, and bounds. Values are copied to native storage; the
   * existing periodicity setting is not changed by this method. The array must contain
   * {@code xsize*ysize*zsize} values in x-fastest order.
   *
   * @param values table values with x varying fastest, then y, then z
   * @param xsize number of table points along x
   * @param ysize number of table points along y
   * @param zsize number of table points along z
   * @param xmin x coordinate of the first x table point
   * @param xmax x coordinate of the last x table point
   * @param ymin y coordinate of the first y table point
   * @param ymax y coordinate of the last y table point
   * @param zmin z coordinate of the first z table point
   * @param zmax z coordinate of the last z table point
   */
  public void setFunctionParameters(double[] values, int xsize, int ysize, int zsize, double xmin,
                                    double xmax, double ymin, double ymax, double zmin, double zmax) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMNative.OpenMM_Continuous3DFunction_setFunctionParameters(
          getPointer(), xsize, ysize, zsize, nativeValues.getPointer(), xmin, xmax, ymin, ymax, zmin, zmax);
    }
  }

  /**
   * Copied parameters defining a three-dimensional continuous function.
   *
   * @param values Java-owned copy of the table values, with x varying fastest
   * @param xsize number of table points along x
   * @param ysize number of table points along y
   * @param zsize number of table points along z
   * @param xmin x coordinate of the first x table point
   * @param xmax x coordinate of the last x table point
   * @param ymin y coordinate of the first y table point
   * @param ymax y coordinate of the last y table point
   * @param zmin z coordinate of the first z table point
   * @param zmax z coordinate of the last z table point
   */
  public record Parameters(double[] values, int xsize, int ysize, int zsize, double xmin, double xmax,
                           double ymin, double ymax, double zmin, double zmax) {
  }

  private static MemorySegment create(double[] values, int xsize, int ysize, int zsize, double xmin,
                                      double xmax, double ymin, double ymax, double zmin, double zmax, boolean periodic) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_Continuous3DFunction_create(xsize, ysize, zsize,
          nativeValues.getPointer(), xmin, xmax, ymin, ymax, zmin, zmax, OpenMMBooleans.toNative(periodic));
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
