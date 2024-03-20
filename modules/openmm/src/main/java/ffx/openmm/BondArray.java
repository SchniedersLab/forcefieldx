// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_get;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_getSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_resize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_set;

/**
 * Bond Array.
 */
public class BondArray {

  private PointerByReference pointer;

  /**
   * Create a new bond array.
   *
   * @param size The size of the bond array.
   */
  public BondArray(int size) {
    OpenMM_BondArray_create(size);
  }

  /**
   * Append a bond to the bond array.
   *
   * @param i1 The first atom index.
   * @param i2 The second atom index.
   */
  public void append(int i1, int i2) {
    OpenMM_BondArray_append(pointer, i1, i2);
  }

  /**
   * Set the bond at index to i1 and i2.
   *
   * @param index The index of the bond to set.
   * @param i1    The first atom index.
   * @param i2    The second atom index.
   */
  public void set(int index, int i1, int i2) {
    OpenMM_BondArray_set(pointer, index, i1, i2);
  }

  /**
   * Get the size of the bond array.
   *
   * @return The size of the bond array.
   */
  public int getSize() {
    return OpenMM_BondArray_getSize(pointer);
  }

  /**
   * Resize the bond array.
   *
   * @param size The new size of the bond array.
   */
  public void resize(int size) {
    OpenMM_BondArray_resize(pointer, size);
  }

  /**
   * Get the bond at index.
   *
   * @param index The index of the bond to get.
   * @param i1    The first atom index.
   * @param i2    The second atom index.
   */
  public void get(int index, IntBuffer i1, IntBuffer i2) {
    OpenMM_BondArray_get(pointer, index, i1, i2);
  }

  /**
   * Get the bond at index.
   *
   * @param index The index of the bond to get.
   * @param i1    The first atom index.
   * @param i2    The second atom index.
   */
  public void get(int index, IntByReference i1, IntByReference i2) {
    OpenMM_BondArray_get(pointer, index, i1, i2);
  }

  /**
   * Get the pointer to the bond array.
   *
   * @return The pointer to the bond array.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Destroy the bond array.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_BondArray_destroy(pointer);
      pointer = null;
    }
  }

}
