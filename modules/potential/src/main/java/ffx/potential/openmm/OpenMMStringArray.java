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

import com.sun.jna.Pointer;
import com.sun.jna.ptr.PointerByReference;

import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_get;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_getSize;

public class OpenMMStringArray {

  private static final Logger logger = Logger.getLogger(OpenMMStringArray.class.getName());

  /**
   * String Array Platform pointer.
   */
  private PointerByReference stringArrayPointer;

  /**
   * OpenMM String Array constructor.
   *
   * @param stringArrayPointer The String Array pointer.
   */
  public OpenMMStringArray(PointerByReference stringArrayPointer) {
    this.stringArrayPointer = stringArrayPointer;
  }

  /**
   * Set the String Array pointer.
   *
   * @param stringArrayPointer The String Array pointer.
   */
  public void setStringArrayPointer(PointerByReference stringArrayPointer) {
    this.stringArrayPointer = stringArrayPointer;
  }

  /**
   * Get the number of strings in the String Array.
   *
   * @return The number of strings in the String Array.
   */
  public int getSize() {
    return OpenMM_StringArray_getSize(stringArrayPointer);
  }

  /**
   * Return the String at index i.
   *
   * @param i The index of the String to return.
   * @return String The requested String.
   */
  public String get(int i) {
    int size = getSize();
    if (i < 0 || i >= size) {
      return null;
    }
    Pointer string = OpenMM_StringArray_get(stringArrayPointer, i);
    if (string == null) {
      return null;
    }
    return string.getString(0);
  }

  /**
   * Destroy the String Array.
   */
  public void destroy() {
    if (stringArrayPointer != null) {
      OpenMM_StringArray_destroy(stringArrayPointer);
      stringArrayPointer = null;
    }
  }

}
