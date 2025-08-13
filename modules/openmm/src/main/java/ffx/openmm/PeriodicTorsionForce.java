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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_addTorsion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_getNumTorsions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_getTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_setTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_usesPeriodicBoundaryConditions;

/**
 * This class implements an interaction between groups of four particles that varies periodically with the torsion angle
 * between them. To use it, create a PeriodicTorsionForce object then call addTorsion() once for each torsion. After
 * a torsion has been added, you can modify its force field parameters by calling setTorsionParameters(). This will
 * have no effect on Contexts that already exist unless you call updateParametersInContext().
 */
public class PeriodicTorsionForce extends Force {

  /**
   * Create a new PeriodicTorsionForce.
   */
  public PeriodicTorsionForce() {
    super(OpenMM_PeriodicTorsionForce_create());
  }

  /**
   * Add a torsion to the PeriodicTorsionForce.
   *
   * @param particle1   Index of the first atom.
   * @param particle2   Index of the second atom.
   * @param particle3   Index of the third atom.
   * @param particle4   Index of the fourth atom.
   * @param periodicity The periodicity of the torsion.
   * @param phase       The phase of the torsion.
   * @param k           The force constant for the torsion.
   * @return The index of the torsion that was added.
   */
  public int addTorsion(int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) {
    return OpenMM_PeriodicTorsionForce_addTorsion(pointer, particle1, particle2, particle3, particle4, periodicity, phase, k);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_PeriodicTorsionForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the number of torsions in the force.
   *
   * @return The number of torsions.
   */
  public int getNumTorsions() {
    return OpenMM_PeriodicTorsionForce_getNumTorsions(pointer);
  }

  /**
   * Get the parameters for a torsion.
   *
   * @param index       The index of the torsion.
   * @param particle1   The index of the first atom (output).
   * @param particle2   The index of the second atom (output).
   * @param particle3   The index of the third atom (output).
   * @param particle4   The index of the fourth atom (output).
   * @param periodicity The periodicity of the torsion (output).
   * @param phase       The phase of the torsion (output).
   * @param k           The force constant for the torsion (output).
   */
  public void getTorsionParameters(int index, IntByReference particle1, IntByReference particle2,
                                   IntByReference particle3, IntByReference particle4,
                                   IntByReference periodicity, DoubleByReference phase,
                                   DoubleByReference k) {
    OpenMM_PeriodicTorsionForce_getTorsionParameters(pointer, index, particle1, particle2,
        particle3, particle4, periodicity, phase, k);
  }

  /**
   * Get the parameters for a torsion.
   *
   * @param index       The index of the torsion.
   * @param particle1   The index of the first atom (output).
   * @param particle2   The index of the second atom (output).
   * @param particle3   The index of the third atom (output).
   * @param particle4   The index of the fourth atom (output).
   * @param periodicity The periodicity of the torsion (output).
   * @param phase       The phase of the torsion (output).
   * @param k           The force constant for the torsion (output).
   */
  public void getTorsionParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                   IntBuffer particle3, IntBuffer particle4,
                                   IntBuffer periodicity, DoubleBuffer phase,
                                   DoubleBuffer k) {
    OpenMM_PeriodicTorsionForce_getTorsionParameters(pointer, index, particle1, particle2,
        particle3, particle4, periodicity, phase, k);
  }

  /**
   * Set the parameters for a torsion.
   *
   * @param index       The index of the torsion for which to set parameters.
   * @param particle1   The index of the first atom in the torsion.
   * @param particle2   The index of the second atom in the torsion.
   * @param particle3   The index of the third atom in the torsion.
   * @param particle4   The index of the fourth atom in the torsion.
   * @param periodicity The periodicity of the torsion.
   * @param phase       The phase of the torsion.
   * @param k           The force constant for the torsion.
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3,
                                   int particle4, int periodicity, double phase, double k) {
    OpenMM_PeriodicTorsionForce_setTorsionParameters(pointer, index, particle1, particle2,
        particle3, particle4, periodicity, phase, k);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_PeriodicTorsionForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_PeriodicTorsionForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_PeriodicTorsionForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}