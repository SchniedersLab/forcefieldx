// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_get;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_getSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_resize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_set;

/**
 * Int Array.
 */
public class IntArray {

  private PointerByReference pointer;

  /**
   * Constructor.
   *
   * @param size The size of the array.
   */
  public IntArray(int size) {
    pointer = OpenMM_IntArray_create(size);
  }

  /**
   * Append a int value to the array.
   *
   * @param value The value to append.
   */
  public void append(int value) {
    OpenMM_IntArray_append(pointer, value);
  }

  /**
   * Set a value in the array.
   *
   * @param index The index.
   * @param value The value.
   */
  public void set(int index, int value) {
    OpenMM_IntArray_set(pointer, index, value);
  }

  /**
   * Destroy the array.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_IntArray_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Resize the array.
   *
   * @param size The new size.
   */
  public void resize(int size) {
    OpenMM_IntArray_resize(pointer, size);
  }


  /**
   * Get the size of the array.
   *
   * @return The size.
   */
  public int getSize() {
    return OpenMM_IntArray_getSize(pointer);
  }

  /**
   * Get a value from the array.
   *
   * @param index The index.
   * @return The value.
   */
  public int get(int index) {
    return OpenMM_IntArray_get(pointer, index);
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
