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
package ffx.openmm.amoeba;

import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import ffx.openmm.Force;

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;

/**
 * This class implements the Amoeba torsion-torsion interaction.
 * <p>
 * To use it, create an AmoebaTorsionTorsionForce object then call addTorsionTorsion() once for each torsion-torsion.  After
 * a torsion-torsion has been added, you can modify its force field parameters by calling setTorsionTorsionParameters().
 */
public class TorsionTorsionForce extends Force {

  /**
   * Create an AmoebaTorsionTorsionForce.
   */
  public TorsionTorsionForce() {
    super(OpenMM_AmoebaTorsionTorsionForce_create());
  }

  /**
   * Add a torsion-torsion term to the force field.
   *
   * @param particle1            the index of the first particle connected by the torsion-torsion
   * @param particle2            the index of the second particle connected by the torsion-torsion
   * @param particle3            the index of the third particle connected by the torsion-torsion
   * @param particle4            the index of the fourth particle connected by the torsion-torsion
   * @param particle5            the index of the fifth particle connected by the torsion-torsion
   * @param chiralCheckAtomIndex the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
   * @param gridIndex            the index to the grid to be used
   * @return the index of the torsion-torsion that was added
   */
  public int addTorsionTorsion(int particle1, int particle2, int particle3, int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) {
    return OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(pointer, particle1, particle2, particle3, particle4, particle5, chiralCheckAtomIndex, gridIndex);
  }

  /**
   * Destroy the Amoeba Torsion-Torsion Force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_AmoebaTorsionTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the number of torsion-torsion terms in the potential function
   *
   * @return the number of torsion-torsion terms
   */
  public int getNumTorsionTorsions() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsions(pointer);
  }

  /**
   * Get the number of torsion-torsion grids
   *
   * @return the number of torsion-torsion grids
   */
  public int getNumTorsionTorsionGrids() {
    return OpenMM_AmoebaTorsionTorsionForce_getNumTorsionTorsionGrids(pointer);
  }

  /**
   * Get the force field parameters for a torsion-torsion term.
   *
   * @param index the index of the torsion-torsion for which to get parameters
   * @param particle1                 the index of the first particle connected by the torsion-torsion
   * @param particle2                 the index of the second particle connected by the torsion-torsion
   * @param particle3                 the index of the third particle connected by the torsion-torsion
   * @param particle4                 the index of the fourth particle connected by the torsion-torsion
   * @param particle5                 the index of the fifth particle connected by the torsion-torsion
   * @param chiralCheckAtomIndex      the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
   * @param gridIndex                 the grid index
   */
  public void getTorsionTorsionParameters(int index, IntByReference particle1, IntByReference particle2,
                                          IntByReference particle3, IntByReference particle4, IntByReference particle5,
                                          IntByReference chiralCheckAtomIndex, IntByReference gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionParameters(pointer, index, particle1, particle2, particle3,
        particle4, particle5, chiralCheckAtomIndex, gridIndex);
  }

  /**
   * Get the torsion-torsion grid at the specified index
   *
   * @param index the grid index
   * @return grid         return grid reference
   */
  public PointerByReference getTorsionTorsionGrid(int index) {
    return OpenMM_AmoebaTorsionTorsionForce_getTorsionTorsionGrid(pointer, index);
  }

  /**
   * Set the force field parameters for a torsion-torsion term.
   *
   * @param index                the index of the torsion-torsion for which to set parameters
   * @param particle1            the index of the first particle connected by the torsion-torsion
   * @param particle2            the index of the second particle connected by the torsion-torsion
   * @param particle3            the index of the third particle connected by the torsion-torsion
   * @param particle4            the index of the fourth particle connected by the torsion-torsion
   * @param particle5            the index of the fifth particle connected by the torsion-torsion
   * @param chiralCheckAtomIndex the index of the particle connected to particle3, but not particle2 or particle4 to be used in chirality check
   * @param gridIndex            the grid index
   */
  public void setTorsionTorsionParameters(int index, int particle1, int particle2, int particle3,
                                          int particle4, int particle5, int chiralCheckAtomIndex, int gridIndex) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionParameters(pointer, index, particle1, particle2, particle3,
        particle4, particle5, chiralCheckAtomIndex, gridIndex);
  }

  /**
   * Set the torsion-torsion grid at the specified index
   *
   * @param gridIndex the index of the torsion-torsion for which to get parameters
   * @param grid      either 3 or 6 values may be specified per grid point.  If the derivatives
   *                  are omitted, they are calculated automatically by fitting a 2D spline to
   *                  the energies.
   *                  grid[x][y][0] = x value
   *                  grid[x][y][1] = y value
   *                  grid[x][y][2] = energy
   *                  grid[x][y][3] = dEdx value
   *                  grid[x][y][4] = dEdy value
   *                  grid[x][y][5] = dEd(xy) value
   */
  public void setTorsionTorsionGrid(int gridIndex, PointerByReference grid) {
    OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(pointer, gridIndex, grid);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
   *
   * @param periodic if true, periodic boundary conditions will be used
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_AmoebaTorsionTorsionForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Returns whether or not this force makes use of periodic boundary
   * conditions.
   *
   * @return true if force uses PBC and false otherwise
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_AmoebaTorsionTorsionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}