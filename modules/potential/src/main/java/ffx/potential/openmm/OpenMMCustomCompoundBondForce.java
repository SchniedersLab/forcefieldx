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

import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_updateParametersInContext;

/**
 * OpenMM Custom Compound Bond Force.
 */
public class OpenMMCustomCompoundBondForce extends OpenMMForce {

  public OpenMMCustomCompoundBondForce(int i, String energy) {
    forcePointer = OpenMM_CustomCompoundBondForce_create(i, energy);
  }

  /**
   * Add a per-bond parameter to the OpenMM System.
   *
   * @param name The name of the parameter.
   */
  public void addPerBondParameter(String name) {
    OpenMM_CustomCompoundBondForce_addPerBondParameter(forcePointer, name);
  }

  /**
   * Add a global parameter.
   *
   * @param name  The parameter name.
   * @param value The parameter value.
   */
  public void addGlobalParameter(String name, double value) {
    OpenMM_CustomCompoundBondForce_addGlobalParameter(forcePointer, name, value);
  }

  /**
   * Add a Custom Compound Bond to the OpenMM System.
   *
   * @param particles  The indices of the particles.
   * @param parameters The bond parameters.
   */
  public void addBond(OpenMMIntArray particles, OpenMMDoubleArray parameters) {
    OpenMM_CustomCompoundBondForce_addBond(forcePointer, particles.getPointer(), parameters.getPointer());
  }

  /**
   * Set the parameters for a Custom Compound Bond.
   *
   * @param index      The index of the bond.
   * @param particles  The indices of the particles.
   * @param parameters The bond parameters.
   */
  public void setBondParameters(int index, OpenMMIntArray particles, OpenMMDoubleArray parameters) {
    OpenMM_CustomCompoundBondForce_setBondParameters(forcePointer, index, particles.getPointer(), parameters.getPointer());
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param openMMContext The OpenMM Context.
   */
  public void updateParametersInContext(OpenMMContext openMMContext) {
    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomCompoundBondForce_updateParametersInContext(forcePointer, openMMContext.getContextPointer());
    }
  }

  /**
   * Clean up.
   */
  public void destroy() {
    OpenMM_CustomCompoundBondForce_destroy(forcePointer);
  }

}
