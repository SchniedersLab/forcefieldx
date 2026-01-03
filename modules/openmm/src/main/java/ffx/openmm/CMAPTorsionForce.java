// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
 * This class implements an interaction between pairs of dihedral angles.  The interaction energy is
 * defined by an "energy correction map" (CMAP), which is simply a set of tabulated energy values
 * on a regular grid of (phi, psi) angles.  Natural cubic spline interpolation is used to compute
 * forces and energies at arbitrary values of the two angles.
 * <p>
 * To use this class, first create one or more energy correction maps by calling addMap().  For each
 * one, you provide an array of energies at uniformly spaced values of the two angles.  Next,
 * add interactions by calling addTorsion().  For each one, you specify the sequence of particles used
 * to calculate each of the two dihedral angles, and the index of the map used to calculate their
 * interaction energy.
 */
public class CMAPTorsionForce extends Force {

  /**
   * Create a CMAPTorsionForce.
   */
  public CMAPTorsionForce() {
    super(OpenMM_CMAPTorsionForce_create());
  }

  /**
   * Create a new map that can be used for torsion pairs.
   *
   * @param size   the size of the map along each dimension
   * @param energy the energy values for the map.  This must be of length size*size.
   *               The element energy[i+size*j] contains the energy when the first
   *               torsion angle equals i*2*PI/size and the second torsion angle
   *               equals j*2*PI/size.
   * @return the index of the map that was added
   */
  public int addMap(int size, PointerByReference energy) {
    return OpenMM_CMAPTorsionForce_addMap(pointer, size, energy);
  }

  /**
   * Add a CMAP torsion term to the force field.
   *
   * @param map the index of the map to use for this term
   * @param a1  the index of the first particle forming the first torsion
   * @param a2  the index of the second particle forming the first torsion
   * @param a3  the index of the third particle forming the first torsion
   * @param a4  the index of the fourth particle forming the first torsion
   * @param b1  the index of the first particle forming the second torsion
   * @param b2  the index of the second particle forming the second torsion
   * @param b3  the index of the third particle forming the second torsion
   * @param b4  the index of the fourth particle forming the second torsion
   * @return the index of the torsion that was added
   */
  public int addTorsion(int map, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4) {
    return OpenMM_CMAPTorsionForce_addTorsion(pointer, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CMAPTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the energy values of a map.
   *
   * @param index  the index of the map for which to get energy values
   * @param size   the size of the map along each dimension
   * @param energy the energy values for the map.  This must be of length size*size.
   *               The element energy[i+size*j] contains the energy when the first
   *               torsion angle equals i*2*PI/size and the second torsion angle
   *               equals j*2*PI/size.
   */
  public void getMapParameters(int index, IntByReference size, PointerByReference energy) {
    OpenMM_CMAPTorsionForce_getMapParameters(pointer, index, size, energy);
  }

  /**
   * Get the energy values of a map.
   *
   * @param index  the index of the map for which to get energy values
   * @param size   the size of the map along each dimension
   * @param energy the energy values for the map.  This must be of length size*size.
   *               The element energy[i+size*j] contains the energy when the first
   *               torsion angle equals i*2*PI/size and the second torsion angle
   *               equals j*2*PI/size.
   */
  public void getMapParameters(int index, IntBuffer size, PointerByReference energy) {
    OpenMM_CMAPTorsionForce_getMapParameters(pointer, index, size, energy);
  }

  /**
   * Get the number of maps that have been defined.
   */
  public int getNumMaps() {
    return OpenMM_CMAPTorsionForce_getNumMaps(pointer);
  }

  /**
   * Get the number of CMAP torsion terms in the potential function
   */
  public int getNumTorsions() {
    return OpenMM_CMAPTorsionForce_getNumTorsions(pointer);
  }

  /**
   * Get the force field parameters for a CMAP torsion term.
   *
   * @param index the index of the torsion for which to get parameters
   * @param map   the index of the map to use for this term
   * @param a1    the index of the first particle forming the first torsion
   * @param a2    the index of the second particle forming the first torsion
   * @param a3    the index of the third particle forming the first torsion
   * @param a4    the index of the fourth particle forming the first torsion
   * @param b1    the index of the first particle forming the second torsion
   * @param b2    the index of the second particle forming the second torsion
   * @param b3    the index of the third particle forming the second torsion
   * @param b4    the index of the fourth particle forming the second torsion
   */
  public void getTorsionParameters(int index, IntByReference map, IntByReference a1, IntByReference a2,
                                   IntByReference a3, IntByReference a4, IntByReference b1,
                                   IntByReference b2, IntByReference b3, IntByReference b4) {
    OpenMM_CMAPTorsionForce_getTorsionParameters(pointer, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Get the force field parameters for a CMAP torsion term.
   *
   * @param index the index of the torsion for which to get parameters
   * @param map   the index of the map to use for this term
   * @param a1    the index of the first particle forming the first torsion
   * @param a2    the index of the second particle forming the first torsion
   * @param a3    the index of the third particle forming the first torsion
   * @param a4    the index of the fourth particle forming the first torsion
   * @param b1    the index of the first particle forming the second torsion
   * @param b2    the index of the second particle forming the second torsion
   * @param b3    the index of the third particle forming the second torsion
   * @param b4    the index of the fourth particle forming the second torsion
   */
  public void getTorsionParameters(int index, IntBuffer map, IntBuffer a1, IntBuffer a2,
                                   IntBuffer a3, IntBuffer a4, IntBuffer b1,
                                   IntBuffer b2, IntBuffer b3, IntBuffer b4) {
    OpenMM_CMAPTorsionForce_getTorsionParameters(pointer, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Set the energy values of a map.
   *
   * @param index  the index of the map for which to set energy values
   * @param size   the size of the map along each dimension
   * @param energy the energy values for the map.  This must be of length size*size.
   *               The element energy[i+size*j] contains the energy when the first
   *               torsion angle equals i*2*PI/size and the second torsion angle
   *               equals j*2*PI/size.
   */
  public void setMapParameters(int index, int size, PointerByReference energy) {
    OpenMM_CMAPTorsionForce_setMapParameters(pointer, index, size, energy);
  }

  /**
   * Set the force field parameters for a CMAP torsion term.
   *
   * @param index the index of the torsion for which to set parameters
   * @param map   the index of the map to use for this term
   * @param a1    the index of the first particle forming the first torsion
   * @param a2    the index of the second particle forming the first torsion
   * @param a3    the index of the third particle forming the first torsion
   * @param a4    the index of the fourth particle forming the first torsion
   * @param b1    the index of the first particle forming the second torsion
   * @param b2    the index of the second particle forming the second torsion
   * @param b3    the index of the third particle forming the second torsion
   * @param b4    the index of the fourth particle forming the second torsion
   */
  public void setTorsionParameters(int index, int map, int a1, int a2, int a3, int a4,
                                   int b1, int b2, int b3, int b4) {
    OpenMM_CMAPTorsionForce_setTorsionParameters(pointer, index, map, a1, a2, a3, a4, b1, b2, b3, b4);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   * Usually this is not appropriate for bonded forces, but there are situations when it can be useful.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CMAPTorsionForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the map and torsion parameters in a Context to match those stored in this Force object.  This method provides
   * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
   * Simply call setMapParameters() and setTorsionParameters() to modify this object's parameters, then call updateParametersInContext()
   * to copy them over to the Context.
   * <p>
   * The only information that can be updated with this method is the energy values for a map, and the map index
   * for a torsion.  The size of a map and the set of particles involved in a torsion cannot be changed.  Also,
   * new bonds and torsions cannot be added.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CMAPTorsionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Returns whether or not this force makes use of periodic boundary
   * conditions.
   *
   * @return true if force uses PBC and false otherwise
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CMAPTorsionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}