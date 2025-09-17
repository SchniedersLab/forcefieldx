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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_addTorsion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_getNumTorsions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_getTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_setTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RBTorsionForce_usesPeriodicBoundaryConditions;

/**
 * This class implements an interaction between groups of four particles that varies with the torsion angle between them
 * according to the Ryckaert-Bellemans potential. To use it, create an RBTorsionForce object then call addTorsion() once
 * for each torsion. After a torsion has been added, you can modify its force field parameters by calling setTorsionParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 */
public class RBTorsionForce extends Force {

  /**
   * Create a new RBTorsionForce.
   */
  public RBTorsionForce() {
    super(OpenMM_RBTorsionForce_create());
  }

  /**
   * Add a torsion to the force.
   *
   * @param particle1 The index of the first particle forming the torsion.
   * @param particle2 The index of the second particle forming the torsion.
   * @param particle3 The index of the third particle forming the torsion.
   * @param particle4 The index of the fourth particle forming the torsion.
   * @param c0        The C0 RB parameter.
   * @param c1        The C1 RB parameter.
   * @param c2        The C2 RB parameter.
   * @param c3        The C3 RB parameter.
   * @param c4        The C4 RB parameter.
   * @param c5        The C5 RB parameter.
   * @return The index of the torsion that was added.
   */
  public int addTorsion(int particle1, int particle2, int particle3, int particle4,
                        double c0, double c1, double c2, double c3, double c4, double c5) {
    return OpenMM_RBTorsionForce_addTorsion(pointer, particle1, particle2, particle3, particle4,
        c0, c1, c2, c3, c4, c5);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_RBTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the number of torsions in the force.
   *
   * @return The number of torsions.
   */
  public int getNumTorsions() {
    return OpenMM_RBTorsionForce_getNumTorsions(pointer);
  }

  /**
   * Get the parameters for a torsion.
   *
   * @param index     The index of the torsion for which to get parameters.
   * @param particle1 The index of the first particle forming the torsion (output).
   * @param particle2 The index of the second particle forming the torsion (output).
   * @param particle3 The index of the third particle forming the torsion (output).
   * @param particle4 The index of the fourth particle forming the torsion (output).
   * @param c0        The C0 RB parameter (output).
   * @param c1        The C1 RB parameter (output).
   * @param c2        The C2 RB parameter (output).
   * @param c3        The C3 RB parameter (output).
   * @param c4        The C4 RB parameter (output).
   * @param c5        The C5 RB parameter (output).
   */
  public void getTorsionParameters(int index, IntByReference particle1, IntByReference particle2,
                                   IntByReference particle3, IntByReference particle4,
                                   DoubleByReference c0, DoubleByReference c1, DoubleByReference c2,
                                   DoubleByReference c3, DoubleByReference c4, DoubleByReference c5) {
    OpenMM_RBTorsionForce_getTorsionParameters(pointer, index, particle1, particle2, particle3, particle4,
        c0, c1, c2, c3, c4, c5);
  }

  /**
   * Get the parameters for a torsion.
   *
   * @param index     The index of the torsion for which to get parameters.
   * @param particle1 The index of the first particle forming the torsion (output).
   * @param particle2 The index of the second particle forming the torsion (output).
   * @param particle3 The index of the third particle forming the torsion (output).
   * @param particle4 The index of the fourth particle forming the torsion (output).
   * @param c0        The C0 RB parameter (output).
   * @param c1        The C1 RB parameter (output).
   * @param c2        The C2 RB parameter (output).
   * @param c3        The C3 RB parameter (output).
   * @param c4        The C4 RB parameter (output).
   * @param c5        The C5 RB parameter (output).
   */
  public void getTorsionParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                   IntBuffer particle3, IntBuffer particle4,
                                   DoubleBuffer c0, DoubleBuffer c1, DoubleBuffer c2,
                                   DoubleBuffer c3, DoubleBuffer c4, DoubleBuffer c5) {
    OpenMM_RBTorsionForce_getTorsionParameters(pointer, index, particle1, particle2, particle3, particle4,
        c0, c1, c2, c3, c4, c5);
  }

  /**
   * Set the parameters for a torsion.
   *
   * @param index     The index of the torsion for which to set parameters.
   * @param particle1 The index of the first particle forming the torsion.
   * @param particle2 The index of the second particle forming the torsion.
   * @param particle3 The index of the third particle forming the torsion.
   * @param particle4 The index of the fourth particle forming the torsion.
   * @param c0        The C0 RB parameter.
   * @param c1        The C1 RB parameter.
   * @param c2        The C2 RB parameter.
   * @param c3        The C3 RB parameter.
   * @param c4        The C4 RB parameter.
   * @param c5        The C5 RB parameter.
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3, int particle4,
                                   double c0, double c1, double c2, double c3, double c4, double c5) {
    OpenMM_RBTorsionForce_setTorsionParameters(pointer, index, particle1, particle2, particle3, particle4,
        c0, c1, c2, c3, c4, c5);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_RBTorsionForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_RBTorsionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_RBTorsionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}