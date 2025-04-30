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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addComputePerDof;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addGlobalVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addPerDofVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addUpdateContextState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_create;

/**
 * Custom Integrator.
 */
public class CustomIntegrator extends Integrator {

  /**
   * Constructor.
   *
   * @param dt The time step.
   */
  public CustomIntegrator(double dt) {
    pointer = OpenMM_CustomIntegrator_create(dt);
  }

  /**
   * Add a per-DOF computation to this Integrator.
   *
   * @param name       The name of the variable to create.
   * @param expression The expression to evaluate.
   */
  public void addComputePerDof(String name, String expression) {
    OpenMM_CustomIntegrator_addComputePerDof(pointer, name, expression);
  }

  /**
   * Add a position constraint to this Integrator.
   */
  public void addConstrainPositions() {
    OpenMM_CustomIntegrator_addConstrainPositions(pointer);
  }

  /**
   * Add a velocity constraint to this Integrator.
   */
  public void addConstrainVelocities() {
    OpenMM_CustomIntegrator_addConstrainVelocities(pointer);
  }

  /**
   * Add a global variable to this Integrator.
   *
   * @param name         The name of the variable to create.
   * @param initialValue The initial value of the variable.
   */
  public void addGlobalVariable(String name, double initialValue) {
    OpenMM_CustomIntegrator_addGlobalVariable(pointer, name, initialValue);
  }

  /**
   * Add a per-DOF variable to this Integrator.
   *
   * @param name         The name of the variable to create.
   * @param initialValue The initial value of the variable.
   */
  public void addPerDofVariable(String name, double initialValue) {
    OpenMM_CustomIntegrator_addPerDofVariable(pointer, name, initialValue);
  }

  /**
   * Add an update context state to this Integrator.
   */
  public void addUpdateContextState() {
    OpenMM_CustomIntegrator_addUpdateContextState(pointer);
  }

}
