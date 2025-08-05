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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_addAngle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_getAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_getNumAngles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_setAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicAngleForce_usesPeriodicBoundaryConditions;

/**
 * This class implements a harmonic angle force. The force is applied between three particles
 * to maintain a specific angle. The energy is computed as:
 * <p>
 * E = k * (theta - theta0)^2
 * <p>
 * where k is the force constant, theta is the current angle, and theta0 is the equilibrium angle.
 */
public class HarmonicAngleForce extends Force {

  /**
   * Create a new HarmonicAngleForce.
   */
  public HarmonicAngleForce() {
    super(OpenMM_HarmonicAngleForce_create());
  }

  /**
   * Add an angle to the force.
   *
   * @param particle1 The index of the first particle forming the angle.
   * @param particle2 The index of the second particle forming the angle (vertex).
   * @param particle3 The index of the third particle forming the angle.
   * @param angle     The equilibrium angle in radians.
   * @param k         The force constant for the angle.
   * @return The index of the angle that was added.
   */
  public int addAngle(int particle1, int particle2, int particle3, double angle, double k) {
    return OpenMM_HarmonicAngleForce_addAngle(pointer, particle1, particle2, particle3, angle, k);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_HarmonicAngleForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for an angle.
   *
   * @param index     The index of the angle for which to get parameters.
   * @param particle1 The index of the first particle forming the angle (output).
   * @param particle2 The index of the second particle forming the angle (output).
   * @param particle3 The index of the third particle forming the angle (output).
   * @param angle     The equilibrium angle in radians (output).
   * @param k         The force constant for the angle (output).
   */
  public void getAngleParameters(int index, IntByReference particle1, IntByReference particle2,
                                 IntByReference particle3, DoubleByReference angle, DoubleByReference k) {
    OpenMM_HarmonicAngleForce_getAngleParameters(pointer, index, particle1, particle2, particle3, angle, k);
  }

  /**
   * Get the parameters for an angle.
   *
   * @param index     The index of the angle for which to get parameters.
   * @param particle1 The index of the first particle forming the angle (output).
   * @param particle2 The index of the second particle forming the angle (output).
   * @param particle3 The index of the third particle forming the angle (output).
   * @param angle     The equilibrium angle in radians (output).
   * @param k         The force constant for the angle (output).
   */
  public void getAngleParameters(int index, IntBuffer particle1, IntBuffer particle2,
                                 IntBuffer particle3, DoubleBuffer angle, DoubleBuffer k) {
    OpenMM_HarmonicAngleForce_getAngleParameters(pointer, index, particle1, particle2, particle3, angle, k);
  }

  /**
   * Get the number of angles in the force.
   *
   * @return The number of angles.
   */
  public int getNumAngles() {
    return OpenMM_HarmonicAngleForce_getNumAngles(pointer);
  }

  /**
   * Set the parameters for an angle.
   *
   * @param index     The index of the angle for which to set parameters.
   * @param particle1 The index of the first particle forming the angle.
   * @param particle2 The index of the second particle forming the angle (vertex).
   * @param particle3 The index of the third particle forming the angle.
   * @param angle     The equilibrium angle in radians.
   * @param k         The force constant for the angle.
   */
  public void setAngleParameters(int index, int particle1, int particle2, int particle3, double angle, double k) {
    OpenMM_HarmonicAngleForce_setAngleParameters(pointer, index, particle1, particle2, particle3, angle, k);
  }

  /**
   * Set whether this force should apply periodic boundary conditions when calculating displacements.
   *
   * @param periodic If true, periodic boundary conditions will be applied.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_HarmonicAngleForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_HarmonicAngleForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_HarmonicAngleForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}