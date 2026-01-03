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
 * This is a force whose energy equals the root mean squared deviation (RMSD)
 * between the current coordinates and a reference structure. It is intended for
 * use with CustomCVForce. You will not normally want a force that exactly equals
 * the RMSD, but there are many situations where it is useful to have a restraining
 * or biasing force that depends on the RMSD in some way.
 * <p>
 * The force is computed by first aligning the particle positions to the reference
 * structure, then computing the RMSD between the aligned positions and the reference.
 * The computation can optionally be done based on only a subset of the particles
 * in the system.
 */
public class RMSDForce extends Force {

  /**
   * Create an RMSDForce.
   *
   * @param referencePositions the reference positions to compute the deviation
   *                           from. The length of this vector must equal the
   *                           number of particles in the system, even if not
   *                           all particles are used in computing the RMSD.
   * @param particles          the indices of the particles to use when computing
   *                           the RMSD. If this is empty (the default), all
   *                           particles in the system will be used.
   */
  public RMSDForce(PointerByReference particles, PointerByReference referencePositions) {
    super(OpenMM_RMSDForce_create(particles, referencePositions));
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_RMSDForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the indices of the particles to use when computing the RMSD. If this
   * is empty, all particles in the system will be used.
   *
   * @return the indices of the particles to use when computing the RMSD.
   */
  public PointerByReference getParticles() {
    return OpenMM_RMSDForce_getParticles(pointer);
  }

  /**
   * Get the reference positions to compute the deviation from.
   *
   * @return the reference positions to compute the deviation from.
   */
  public PointerByReference getReferencePositions() {
    return OpenMM_RMSDForce_getReferencePositions(pointer);
  }

  /**
   * Set the indices of the particles to use when computing the RMSD. If this
   * is empty, all particles in the system will be used.
   *
   * @param particles the indices of the particles to use when computing the RMSD.
   */
  public void setParticles(PointerByReference particles) {
    OpenMM_RMSDForce_setParticles(pointer, particles);
  }

  /**
   * Set the reference positions to compute the deviation from.
   *
   * @param referencePositions the reference positions to compute the deviation from.
   */
  public void setReferencePositions(PointerByReference referencePositions) {
    OpenMM_RMSDForce_setReferencePositions(pointer, referencePositions);
  }

  /**
   * Update the reference positions and particle indices in a Context to match those stored
   * in this Force object. This method provides an efficient method to update certain parameters
   * in an existing Context without needing to reinitialize it. Simply call setReferencePositions()
   * and setParticles() to modify this object's parameters, then call updateParametersInContext()
   * to copy them over to the Context.
   *
   * @param context the Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_RMSDForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Returns whether or not this force makes use of periodic boundary
   * conditions.
   *
   * @return true if force uses PBC and false otherwise.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_RMSDForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}