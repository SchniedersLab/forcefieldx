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
 * A three-dimensional discrete tabulated function. Each input coordinate is rounded to the
 * nearest integer; the table element at those indices is returned. If any index is outside its
 * dimension's valid range, the result is undefined. Values are stored with x varying fastest,
 * then y, then z: {@code values[i + xsize*j + xsize*ysize*k] = f(i,j,k)}.
 */
public class Discrete3DFunction extends TabulatedFunction {

  /**
   * Create a discrete function from dimensions and x-fastest values. The array is copied to
   * native storage and must contain {@code xsize*ysize*zsize} elements.
   *
   * @param xsize number of table elements along x
   * @param ysize number of table elements along y
   * @param zsize number of table elements along z
   * @param values table values in x-fastest order,
   *     {@code values[i + xsize*j + xsize*ysize*k] = f(i,j,k)}
   */
  public Discrete3DFunction(int xsize, int ysize, int zsize, double[] values) {
    super(create(xsize, ysize, zsize, values));
  }

  /** Destroy the native discrete function. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_Discrete3DFunction_destroy);
  }

  /**
   * Get the table values and dimensions as Java-owned copies.
   *
   * @return copied table values in x-fastest order, plus the x, y, and z dimensions
   */
  public Parameters getFunctionParameters() {
    try (Arena arena = Arena.ofConfined(); DoubleArray values = new DoubleArray(0)) {
      MemorySegment xs = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment ys = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment zs = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_Discrete3DFunction_getFunctionParameters(
          getPointer(), xs, ys, zs, values.getPointer());
      return new Parameters(copy(values), xs.get(ValueLayout.JAVA_INT, 0),
          ys.get(ValueLayout.JAVA_INT, 0), zs.get(ValueLayout.JAVA_INT, 0));
    }
  }

  /**
   * Replace the dimensions and x-fastest table values. The array is copied to native storage and
   * must contain {@code xsize*ysize*zsize} elements.
   *
   * @param xsize number of table elements along x
   * @param ysize number of table elements along y
   * @param zsize number of table elements along z
   * @param values table values in x-fastest order,
   *     {@code values[i + xsize*j + xsize*ysize*k] = f(i,j,k)}
   */
  public void setFunctionParameters(int xsize, int ysize, int zsize, double[] values) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMNative.OpenMM_Discrete3DFunction_setFunctionParameters(
          getPointer(), xsize, ysize, zsize, nativeValues.getPointer());
    }
  }

  /**
   * Copied table values and dimensions.
   *
   * @param values Java-owned copy of table values in x-fastest order
   * @param xsize number of table elements along x
   * @param ysize number of table elements along y
   * @param zsize number of table elements along z
   */
  public record Parameters(double[] values, int xsize, int ysize, int zsize) {
  }

  private static MemorySegment create(int xsize, int ysize, int zsize, double[] values) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_Discrete3DFunction_create(
          xsize, ysize, zsize, nativeValues.getPointer());
    }
  }

  private static DoubleArray toNative(double[] values) {
    DoubleArray result = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) result.set(i, values[i]);
    return result;
  }

  private static double[] copy(DoubleArray array) {
    double[] result = new double[array.getSize()];
    for (int i = 0; i < result.length; i++) result[i] = array.get(i);
    return result;
  }
}
