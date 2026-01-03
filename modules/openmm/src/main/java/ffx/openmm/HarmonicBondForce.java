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

import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_getBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_getNumBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_usesPeriodicBoundaryConditions;

/**
 * This class implements an interaction between pairs of particles that varies harmonically with the distance
 * between them.  To use it, create a HarmonicBondForce object then call addBond() once for each bond.  After
 * a bond has been added, you can modify its force field parameters by calling setBondParameters().  This will
 * have no effect on Contexts that already exist unless you call updateParametersInContext().
 */
public class HarmonicBondForce extends Force {

  /**
   * Create a new HarmonicBondForce.
   */
  public HarmonicBondForce() {
    super(OpenMM_HarmonicBondForce_create());
  }

  /**
   * Add a bond term to the force field.
   *
   * @param i1     The index of the first particle connected by the bond.
   * @param i2     The index of the second particle connected by the bond.
   * @param length The equilibrium length of the bond, measured in nm.
   * @param k      The harmonic force constant for the bond, measured in kJ/mol/nm&circ;2.
   * @return The index of the bond that was added.
   */
  public int addBond(int i1, int i2, double length, double k) {
    return OpenMM_HarmonicBondForce_addBond(pointer, i1, i2, length, k);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_HarmonicBondForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the force field parameters for a bond term.
   *
   * @param index  The index of the bond for which to get parameters.
   * @param i1     The index of the first particle connected by the bond (output).
   * @param i2     The index of the second particle connected by the bond (output).
   * @param length The equilibrium length of the bond, measured in nm (output).
   * @param k      The harmonic force constant for the bond, measured in kJ/mol/nm&circ;2 (output).
   */
  public void getBondParameters(int index, IntByReference i1, IntByReference i2,
                                DoubleByReference length, DoubleByReference k) {
    OpenMM_HarmonicBondForce_getBondParameters(pointer, index, i1, i2, length, k);
  }

  /**
   * Get the force field parameters for a bond term.
   *
   * @param index  The index of the bond for which to get parameters.
   * @param i1     The index of the first particle connected by the bond (output).
   * @param i2     The index of the second particle connected by the bond (output).
   * @param length The equilibrium length of the bond, measured in nm (output).
   * @param k      The harmonic force constant for the bond, measured in kJ/mol/nm&circ;2 (output).
   */
  public void getBondParameters(int index, IntBuffer i1, IntBuffer i2,
                                DoubleBuffer length, DoubleBuffer k) {
    OpenMM_HarmonicBondForce_getBondParameters(pointer, index, i1, i2, length, k);
  }

  /**
   * Get the number of bonds.
   *
   * @return The number of bonds.
   */
  public int getNumBonds() {
    return OpenMM_HarmonicBondForce_getNumBonds(pointer);
  }

  /**
   * Set the force field parameters for a bond term.
   *
   * @param i      The index of the bond for which to set parameters.
   * @param i1     The index of the first particle connected by the bond.
   * @param i2     The index of the second particle connected by the bond.
   * @param length The equilibrium length of the bond, measured in nm.
   * @param k      The harmonic force constant for the bond, measured in kJ/mol/nm&circ;2.
   */
  public void setBondParameters(int i, int i1, int i2, double length, double k) {
    OpenMM_HarmonicBondForce_setBondParameters(pointer, i, i1, i2, length, k);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_HarmonicBondForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_HarmonicBondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_HarmonicBondForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}
