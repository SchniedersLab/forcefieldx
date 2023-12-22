// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.openmm;

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_set;

/**
 * DoubleArray3D.
 */
public class DoubleArray3D {

  private PointerByReference pointer;

  /**
   * Constructor.
   *
   * @param d1 The size of the first dimension.
   * @param d2 The size of the second dimension.
   * @param d3 The size of the third dimension.
   */
  public DoubleArray3D(int d1, int d2, int d3) {
    pointer = OpenMM_3D_DoubleArray_create(d1, d2, d3);
  }

  /**
   * Set the value of the array at the given index.
   *
   * @param d1    The first dimension index.
   * @param d2    The second dimension index.
   * @param value The value to set.
   */
  public void set(int d1, int d2, DoubleArray value) {
    OpenMM_3D_DoubleArray_set(pointer, d1, d2, value.getPointer());
  }

  /**
   * Destroy the array.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_3D_DoubleArray_destroy(pointer);
      pointer = null;
    }
  }


  /**
   * Get the pointer to the array.
   *
   * @return The pointer.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

}
