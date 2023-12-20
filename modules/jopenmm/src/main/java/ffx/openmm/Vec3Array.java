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
import edu.uiowa.jopenmm.OpenMM_Vec3;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_get;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_getSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_resize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_set;

/**
 * Vec3 Array.
 */
public class Vec3Array {

  /**
   * String vec3 array pointer.
   */
  private PointerByReference pointer;

  /**
   * OpenMM Vec3 Array constructor.
   *
   * @param size The size of the String Array.
   */
  public Vec3Array(int size) {
    pointer = OpenMM_Vec3Array_create(size);
  }

  /**
   * OpenMM Vec3 Array constructor.
   *
   * @param pointer The Vec3 Array pointer.
   */
  public Vec3Array(PointerByReference pointer) {
    this.pointer = pointer;
  }

  /**
   * Get the pointer to the vec3 array pointer.
   *
   * @return The pointer to the vec3 array.
   */
  public PointerByReference getPointer() {
    return pointer;
  }

  /**
   * Append a Vec3 to the Vec3Array.
   *
   * @param vec3 The Vec3 to append.
   */
  public void append(OpenMM_Vec3.ByValue vec3) {
    OpenMM_Vec3Array_append(pointer, vec3);
  }

  /**
   * Get a Vec3 from the Vec3Array.
   *
   * @return The Vec3 at index i.
   */
  public OpenMM_Vec3 get(int i) {
    return OpenMM_Vec3Array_get(pointer, i);
  }

  /**
   * Get the size of the Vec3Array.
   *
   * @return The size of the Vec3Array.
   */
  public int getSize() {
    return OpenMM_Vec3Array_getSize(pointer);
  }

  /**
   * Resize the Vec3Array.
   *
   * @param size The new size of the Vec3Array.
   */
  public void resize(int size) {
    OpenMM_Vec3Array_resize(pointer, size);
  }

  /**
   * Set a Vec3 in the Vec3Array.
   *
   * @param i    The index of the Vec3 to set.
   * @param vec3 The Vec3 to set.
   */
  public void set(int i, OpenMM_Vec3.ByValue vec3) {
    OpenMM_Vec3Array_set(pointer, i, vec3);
  }

  /**
   * Destroy the Vec3Array.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_Vec3Array_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Convert a double array to a Vec3Array.
   *
   * @param array The double array.
   * @return The Vec3Array.
   */
  public static Vec3Array toVec3Array(double[] array) {
    Vec3Array vec3Array = new Vec3Array(0);
    OpenMM_Vec3.ByValue vec3 = new OpenMM_Vec3.ByValue();
    for (int i = 0; i < array.length; i += 3) {
      vec3.x = array[i];
      vec3.y = array[i + 1];
      vec3.z = array[i + 2];
      vec3Array.append(vec3);
    }
    return vec3Array;
  }

  /**
   * Convert the Vec3Array to a double array.
   *
   * @return The double array.
   */
  public double[] getArray() {
    int size = getSize();
    double[] array = new double[size * 3];
    for (int i = 0; i < size; i++) {
      OpenMM_Vec3 vec3 = get(i);
      int index = i * 3;
      array[index] = vec3.x;
      array[index + 1] = vec3.y;
      array[index + 2] = vec3.z;
    }
    return array;
  }

}
