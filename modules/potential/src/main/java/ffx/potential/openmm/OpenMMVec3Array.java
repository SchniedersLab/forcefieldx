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
package ffx.potential.openmm;

import com.sun.jna.ptr.PointerByReference;
import edu.uiowa.jopenmm.OpenMM_Vec3;

import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_get;

/**
 * OpenMM Vec3 Array.
 */
public class OpenMMVec3Array {

  private static final Logger logger = Logger.getLogger(OpenMMStringArray.class.getName());

  /**
   * String vec3 array pointer.
   */
  private PointerByReference vec3ArrayPointer;

  /**
   * OpenMM Vec3 Array constructor.
   *
   * @param size The size of the String Array.
   */
  public OpenMMVec3Array(int size) {
    vec3ArrayPointer = OpenMM_Vec3Array_create(size);
  }

  /**
   * OpenMM Vec3 Array constructor.
   *
   * @param vec3ArrayPointer The Vec3 Array pointer.
   */
  public OpenMMVec3Array(PointerByReference vec3ArrayPointer) {
    this.vec3ArrayPointer = vec3ArrayPointer;
  }

  /**
   * Get the pointer to the vec3 array pointer.
   *
   * @return The pointer to the vec3 array.
   */
  public PointerByReference getPointer() {
    return vec3ArrayPointer;
  }

  /**
   * Append a Vec3 to the Vec3Array.
   *
   * @param vec3 The Vec3 to append.
   */
  public void append(OpenMM_Vec3.ByValue vec3) {
    OpenMM_Vec3Array_append(vec3ArrayPointer, vec3);
  }

  /**
   * Get a Vec3 from the Vec3Array.
   *
   * @return The Vec3 at index i.
   */
  public OpenMM_Vec3 get(int i) {
    return OpenMM_Vec3Array_get(vec3ArrayPointer, i);
  }

  /**
   * Destroy the Vec3Array.
   */
  public void destroy() {
    if (vec3ArrayPointer != null) {
      OpenMM_Vec3Array_destroy(vec3ArrayPointer);
      vec3ArrayPointer = null;
    }
  }

}
