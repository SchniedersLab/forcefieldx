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
 * A one-dimensional discrete tabulated function. To evaluate it, x is rounded to the nearest
 * integer and the table element with that index is returned. If the resulting index is outside
 * {@code [0,size)}, the result is undefined. The meaning and units of table values are determined
 * by the expression using this function; this wrapper performs no unit conversion.
 */
public class Discrete1DFunction extends TabulatedFunction {

  /**
   * Create a function from discrete values. The input array is copied to native storage.
   *
   * @param values table values indexed by integer x from 0 through {@code values.length-1}
   */
  public Discrete1DFunction(double[] values) {
    super(create(values));
  }

  /** Destroy the native discrete function. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_Discrete1DFunction_destroy);
  }

  /**
   * Get a Java-owned copy of the discrete table values.
   *
   * @return copied values indexed by integer x
   */
  public double[] getFunctionParameters() {
    try (DoubleArray values = new DoubleArray(0)) {
      OpenMMNative.OpenMM_Discrete1DFunction_getFunctionParameters(getPointer(), values.getPointer());
      return copy(values);
    }
  }

  /**
   * Replace the discrete table values. The input array is copied to native storage.
   *
   * @param values replacement values indexed by integer x
   */
  public void setFunctionParameters(double[] values) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMNative.OpenMM_Discrete1DFunction_setFunctionParameters(
          getPointer(), nativeValues.getPointer());
    }
  }

  private static MemorySegment create(double[] values) {
    try (DoubleArray nativeValues = toNative(values)) {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_Discrete1DFunction_create(nativeValues.getPointer());
    }
  }

  private static DoubleArray toNative(double[] values) {
    DoubleArray result = new DoubleArray(values.length);
    for (int i = 0; i < values.length; i++) result.set(i, values[i]);
    return result;
  }

  private static double[] copy(DoubleArray values) {
    double[] result = new double[values.getSize()];
    for (int i = 0; i < result.length; i++) result[i] = values.get(i);
    return result;
  }
}
