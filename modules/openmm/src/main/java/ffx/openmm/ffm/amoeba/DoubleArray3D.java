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
package ffx.openmm.ffm.amoeba;

import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.OpenMMHandle;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;
import java.util.Objects;

/**
 * An owned native AMOEBA three-dimensional array of double arrays used to define torsion-torsion
 * energy grids.
 *
 * <p>The native array has three dimensions. {@link #set(int, int, DoubleArray)} sets the
 * innermost values at the specified first- and second-dimension indices.</p>
 */
public class DoubleArray3D extends OpenMMHandle {

  /**
   * Create a three-dimensional array with the specified dimensions.
   *
   * @param d1 size of the first dimension.
   * @param d2 size of the second dimension.
   * @param d3 size of the third dimension.
   */
  public DoubleArray3D(int d1, int d2, int d3) {
    super(create(d1, d2, d3));
  }

  /**
   * Destroy this array and release its native storage.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_3D_DoubleArray_destroy);
  }

  /**
   * Set the innermost array at the specified indices.
   *
   * @param d1     first-dimension index.
   * @param d2     second-dimension index.
   * @param values caller-owned values for the third dimension; its size must equal {@code d3}.
   *               The native array copies the supplied elements.
   * @throws NullPointerException if {@code values} is null.
   */
  public void set(int d1, int d2, DoubleArray values) {
    OpenMMNative.OpenMM_3D_DoubleArray_set(getPointer(), d1, d2,
        Objects.requireNonNull(values, "Values cannot be null.").getPointer());
  }

  private static MemorySegment create(int d1, int d2, int d3) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_3D_DoubleArray_create(d1, d2, d3);
  }
}
