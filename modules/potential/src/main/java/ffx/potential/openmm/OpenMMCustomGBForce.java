// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.openmm;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addComputedValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyTerm;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_updateParametersInContext;

/**
 * OpenMM CustomGBForce.
 */
public class OpenMMCustomGBForce extends OpenMMForce {

  /**
   * OpenMM CustomGBForce constructor.
   */
  public OpenMMCustomGBForce() {
    forcePointer = OpenMM_CustomGBForce_create();
  }

  /**
   * Add per particle parameter.
   *
   * @param name The name of the parameter.
   */
  public void addPerParticleParameter(String name) {
    OpenMM_CustomGBForce_addPerParticleParameter(forcePointer, name);
  }

  /**
   * Add a global parameter.
   *
   * @param name  The parameter name.
   * @param value The parameter value.
   */
  public void addGlobalParameter(String name, double value) {
    OpenMM_CustomGBForce_addGlobalParameter(forcePointer, name, value);
  }

  /**
   * Add a computed value.
   *
   * @param name       The computed value name.
   * @param expression The computed value expression.
   * @param type       The computed value type.
   */
  public void addComputedValue(String name, String expression, int type) {
    OpenMM_CustomGBForce_addComputedValue(forcePointer, name, expression, type);
  }

  /**
   * Add an energy term.
   *
   * @param expression The energy term expression.
   * @param type       The energy term type.
   */
  public void addEnergyTerm(String expression, int type) {
    OpenMM_CustomGBForce_addEnergyTerm(forcePointer, expression, type);
  }

  /**
   * Add a particle to the force.
   *
   * @param particleParameters The particle parameters.
   */
  public void addParticle(OpenMMDoubleArray particleParameters) {
    OpenMM_CustomGBForce_addParticle(forcePointer, particleParameters.getPointer());
  }

  /**
   * Set the particle parameters.
   *
   * @param index              The particle index.
   * @param particleParameters The particle parameters.
   */
  public void setParticleParameters(int index, OpenMMDoubleArray particleParameters) {
    OpenMM_CustomGBForce_setParticleParameters(forcePointer, index, particleParameters.getPointer());
  }

  /**
   * Set the cutoff distance.
   *
   * @param off The cutoff distance.
   */
  public void setCutoffDistance(double off) {
    OpenMM_CustomGBForce_setCutoffDistance(forcePointer, off);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The context.
   */
  public void updateParametersInContext(OpenMMContext context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomGBForce_updateParametersInContext(forcePointer, context.getContextPointer());
    }
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (forcePointer != null) {
      OpenMM_CustomGBForce_destroy(forcePointer);
      forcePointer = null;
    }
  }

}
