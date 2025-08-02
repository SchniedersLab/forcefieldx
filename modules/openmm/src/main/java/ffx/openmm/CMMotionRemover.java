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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_getFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_setFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions;

/**
 * Center of Mass Motion Remover.
 */
public class CMMotionRemover extends Force {

  /**
   * OpenMM CMMotionRemover constructor.
   *
   * @param frequency The frequency to apply the CMMotionRemover.
   */
  public CMMotionRemover(int frequency) {
    pointer = OpenMM_CMMotionRemover_create(frequency);
  }

  /**
   * Get the frequency at which center of mass motion should be removed.
   *
   * @return The frequency at which center of mass motion should be removed.
   */
  public int getFrequency() {
    return OpenMM_CMMotionRemover_getFrequency(pointer);
  }

  /**
   * Destroy the OpenMM CMMotionRemover.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CMMotionRemover_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Set the frequency at which center of mass motion should be removed.
   *
   * @param frequency The frequency at which center of mass motion should be removed.
   */
  public void setFrequency(int frequency) {
    OpenMM_CMMotionRemover_setFrequency(pointer, frequency);
  }

  /**
   * Returns whether this force makes use of periodic boundary conditions.
   *
   * @return The CMMotionRemover always returns false.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(pointer);
    return pbc != OpenMM_False;
  }
}