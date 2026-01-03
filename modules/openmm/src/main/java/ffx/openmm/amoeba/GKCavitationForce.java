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
package ffx.openmm.amoeba;

import ffx.openmm.Context;
import ffx.openmm.Force;

/*
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_updateParametersInContext;
 */

/**
 * Implicit Solvent Cavitation Force.
 */
public class GKCavitationForce extends Force {

  /**
   * Constructor.
   */
  public GKCavitationForce() {
    super(null);
    System.out.println(" The GKCavitationForce is not currently supported.");
    System.exit(-1);
    // pointer = OpenMM_AmoebaGKCavitationForce_create();
  }

  /**
   * Add an atom to the Cavitation force.
   *
   * @param radius         Atomic radius.
   * @param surfaceTension Surface tension.
   * @param isHydrogen     Is this a hydrogen atom?
   */
  public void addParticle(double radius, double surfaceTension, int isHydrogen) {
    // OpenMM_AmoebaGKCavitationForce_addParticle(pointer, radius, surfaceTension, isHydrogen);
  }

  /**
   * Set the parameters for an atom in the Cavitation force.
   *
   * @param index          Atom index.
   * @param radius         Atomic radius.
   * @param surfaceTension Surface tension.
   * @param isHydrogen     Is this a hydrogen atom?
   */
  public void setParticleParameters(int index, double radius, double surfaceTension, int isHydrogen) {
    // OpenMM_AmoebaGKCavitationForce_setParticleParameters(pointer, index, radius, surfaceTension, isHydrogen);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method Nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    // OpenMM_AmoebaGKCavitationForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      // OpenMM_AmoebaGKCavitationForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Update the parameters in the context.
   *
   * @param context OpenMM context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      // OpenMM_AmoebaGKCavitationForce_updateParametersInContext(pointer, context.getPointer());
    }
  }
}
