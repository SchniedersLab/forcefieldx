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

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_getParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_getReferencePositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_setParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_setReferencePositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_RMSDForce_usesPeriodicBoundaryConditions;

/**
 * This class implements a force that computes the root mean square deviation (RMSD)
 * between the current particle positions and a set of reference positions. The RMSD
 * force can be used for structural restraints or biasing simulations toward specific
 * conformations.
 * <p>
 * The RMSD is calculated as the minimum RMSD over all possible rotations and translations
 * that align the current positions with the reference positions. This force is commonly
 * used in enhanced sampling methods and structure-based drug design.
 */
public class RMSDForce extends Force {

  /**
   * Create a new RMSDForce.
   *
   * @param particles          The particles to include in the RMSD calculation.
   * @param referencePositions The reference positions for the particles.
   */
  public RMSDForce(PointerByReference particles, PointerByReference referencePositions) {
    pointer = OpenMM_RMSDForce_create(particles, referencePositions);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_RMSDForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the particles included in the RMSD calculation.
   *
   * @return The particles included in the RMSD calculation.
   */
  public PointerByReference getParticles() {
    return OpenMM_RMSDForce_getParticles(pointer);
  }

  /**
   * Get the reference positions for the particles.
   *
   * @return The reference positions for the particles.
   */
  public PointerByReference getReferencePositions() {
    return OpenMM_RMSDForce_getReferencePositions(pointer);
  }

  /**
   * Set the particles to include in the RMSD calculation.
   *
   * @param particles The particles to include in the RMSD calculation.
   */
  public void setParticles(PointerByReference particles) {
    OpenMM_RMSDForce_setParticles(pointer, particles);
  }

  /**
   * Set the reference positions for the particles.
   *
   * @param referencePositions The reference positions for the particles.
   */
  public void setReferencePositions(PointerByReference referencePositions) {
    OpenMM_RMSDForce_setReferencePositions(pointer, referencePositions);
  }

  /**
   * Update the parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_RMSDForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_RMSDForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}