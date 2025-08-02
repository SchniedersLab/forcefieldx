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

import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_addMap;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_addTorsion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_getMapParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_getNumMaps;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_getNumTorsions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_getTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_setMapParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_setTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions;

/**
 * This class implements a force field term based on the CMAP torsion interaction.
 * It is used to apply corrections to the potential energy based on pairs of dihedral angles.
 * This is especially important for modeling protein backbone conformations accurately.
 */
public class CMAPTorsionForce extends Force {

  /**
   * Create a new CMAPTorsionForce.
   */
  public CMAPTorsionForce() {
    pointer = OpenMM_CMAPTorsionForce_create();
  }

  /**
   * Add a map to the force.
   *
   * @param size   The size of the map (number of grid points along each axis).
   * @param energy The energy values for the map.
   * @return The index of the map that was added.
   */
  public int addMap(int size, PointerByReference energy) {
    return OpenMM_CMAPTorsionForce_addMap(pointer, size, energy);
  }

  /**
   * Add a torsion to the force.
   *
   * @param map The index of the map to use for this torsion.
   * @param a1  The index of the first atom in the first torsion.
   * @param a2  The index of the second atom in the first torsion.
   * @param a3  The index of the third atom in the first torsion.
   * @param a4  The index of the fourth atom in the first torsion.
   * @param b1  The index of the first atom in the second torsion.
   * @param b2  The index of the second atom in the second torsion.
   * @param b3  The index of the third atom in the second torsion.
   * @param b4  The index of the fourth atom in the second torsion.
   * @return The index of the torsion that was added.
   */
  public int addTorsion(int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    return OpenMM_CMAPTorsionForce_addTorsion(pointer, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CMAPTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for a map.
   *
   * @param index  The index of the map.
   * @param size   The size of the map (output).
   * @param energy The energy values for the map (output).
   */
  public void getMapParameters(int index, IntByReference size, PointerByReference energy) {
    OpenMM_CMAPTorsionForce_getMapParameters(pointer, index, size, energy);
  }

  /**
   * Get the parameters for a map.
   *
   * @param index  The index of the map.
   * @param size   The size of the map (output).
   * @param energy The energy values for the map (output).
   */
  public void getMapParameters(int index, IntBuffer size, PointerByReference energy) {
    OpenMM_CMAPTorsionForce_getMapParameters(pointer, index, size, energy);
  }

  /**
   * Get the number of maps.
   *
   * @return The number of maps.
   */
  public int getNumMaps() {
    return OpenMM_CMAPTorsionForce_getNumMaps(pointer);
  }

  /**
   * Get the number of torsions.
   *
   * @return The number of torsions.
   */
  public int getNumTorsions() {
    return OpenMM_CMAPTorsionForce_getNumTorsions(pointer);
  }

  /**
   * Get the parameters for a torsion.
   *
   * @param index The index of the torsion.
   * @param map   The index of the map used for this torsion (output).
   * @param a1    The index of the first atom in the first torsion (output).
   * @param a2    The index of the second atom in the first torsion (output).
   * @param a3    The index of the third atom in the first torsion (output).
   * @param a4    The index of the fourth atom in the first torsion (output).
   * @param b1    The index of the first atom in the second torsion (output).
   * @param b2    The index of the second atom in the second torsion (output).
   * @param b3    The index of the third atom in the second torsion (output).
   * @param b4    The index of the fourth atom in the second torsion (output).
   */
  public void getTorsionParameters(int index, IntByReference map, IntByReference a1, IntByReference a2,
                                   IntByReference a3, IntByReference a4, IntByReference b1,
                                   IntByReference b2, IntByReference b3, IntByReference b4) {
    OpenMM_CMAPTorsionForce_getTorsionParameters(pointer, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Get the parameters for a torsion.
   *
   * @param index The index of the torsion.
   * @param map   The index of the map used for this torsion (output).
   * @param a1    The index of the first atom in the first torsion (output).
   * @param a2    The index of the second atom in the first torsion (output).
   * @param a3    The index of the third atom in the first torsion (output).
   * @param a4    The index of the fourth atom in the first torsion (output).
   * @param b1    The index of the first atom in the second torsion (output).
   * @param b2    The index of the second atom in the second torsion (output).
   * @param b3    The index of the third atom in the second torsion (output).
   * @param b4    The index of the fourth atom in the second torsion (output).
   */
  public void getTorsionParameters(int index, IntBuffer map, IntBuffer a1, IntBuffer a2,
                                   IntBuffer a3, IntBuffer a4, IntBuffer b1,
                                   IntBuffer b2, IntBuffer b3, IntBuffer b4) {
    OpenMM_CMAPTorsionForce_getTorsionParameters(pointer, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Set the parameters for a map.
   *
   * @param index  The index of the map.
   * @param size   The size of the map.
   * @param energy The energy values for the map.
   */
  public void setMapParameters(int index, int size, PointerByReference energy) {
    OpenMM_CMAPTorsionForce_setMapParameters(pointer, index, size, energy);
  }

  /**
   * Set the parameters for a torsion.
   *
   * @param index The index of the torsion.
   * @param map   The index of the map to use for this torsion.
   * @param a1    The index of the first atom in the first torsion.
   * @param a2    The index of the second atom in the first torsion.
   * @param a3    The index of the third atom in the first torsion.
   * @param a4    The index of the fourth atom in the first torsion.
   * @param b1    The index of the first atom in the second torsion.
   * @param b2    The index of the second atom in the second torsion.
   * @param b3    The index of the third atom in the second torsion.
   * @param b4    The index of the fourth atom in the second torsion.
   */
  public void setTorsionParameters(int index, int map, int a1, int a2, int a3, int a4,
                                   int b1, int b2, int b3, int b4) {
    OpenMM_CMAPTorsionForce_setTorsionParameters(pointer, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CMAPTorsionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}