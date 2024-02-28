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
package ffx.openmm.amoeba;

import com.sun.jna.ptr.PointerByReference;
import ffx.openmm.Force;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid;

/**
 * Torsion-Torsion Force.
 */
public class TorsionTorsionForce extends Force {

  /**
   * Create an OpenMM TorsionTorsion Force.
   */
  public TorsionTorsionForce() {
    pointer = OpenMM_AmoebaTorsionTorsionForce_create();
  }

  /**
   * Add a torsion to the TorsionTorsionForce.
   *
   * @param atom1           The index of the first atom.
   * @param atom2           The index of the second atom.
   * @param atom3           The index of the third atom.
   * @param atom4           The index of the fourth atom.
   * @param atom5           The index of the fifth atom.
   * @param chiralCheckAtom The index of the chiral check atom.
   * @param gridIndex       The index of the grid.
   */
  public void addTorsionTorsion(int atom1, int atom2, int atom3, int atom4, int atom5, int chiralCheckAtom, int gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(pointer, atom1, atom2, atom3, atom4, atom5, chiralCheckAtom, gridIndex);
  }

  /**
   * Set the grid for a torsion-torsion.
   *
   * @param gridIndex The index of the grid.
   * @param grid      The grid.
   */
  public void setTorsionTorsionGrid(int gridIndex, PointerByReference grid) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(pointer, gridIndex, grid);
  }

  /**
   * Destroy the Amoeba Torsion-Torsion Force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaTorsionTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

}
