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

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_updateParametersInContext;

/**
 * OpenMM CustomBondForce wrapper.
 */
public class OpenMMCustomBondForce extends OpenMMForce {

  /**
   * Custom Bond Force Constructor.
   *
   * @param energy The definition of the Energy.
   */
  public OpenMMCustomBondForce(String energy) {
    pointer = OpenMM_CustomBondForce_create(energy);
  }

  /**
   * Add a bond to the OpenMM System.
   *
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param parameters The bond parameters.
   */
  public void addBond(int i1, int i2, OpenMMDoubleArray parameters) {
    OpenMM_CustomBondForce_addBond(pointer, i1, i2, parameters.getPointer());
  }

  /**
   * Add a per-bond parameter to the CustomBondForce.
   *
   * @param name The name of the parameter.
   */
  public void addPerBondParameter(String name) {
    OpenMM_CustomBondForce_addPerBondParameter(pointer, name);
  }

  /**
   * Add a global parameter to the CustomBondForce.
   *
   * @param name  The name of the parameter.
   * @param value The value of the parameter.
   */
  public void addGlobalParameter(String name, double value) {
    OpenMM_CustomBondForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Set the parameters for one bond in the OpenMM System.
   *
   * @param index      The index of the bond.
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param parameters The bond parameters.
   */
  public void setBondParameters(int index, int i1, int i2, OpenMMDoubleArray parameters) {
    OpenMM_CustomBondForce_setBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param openMMContext The OpenMM Context.
   */
  public void updateParametersInContext(OpenMMContext openMMContext) {
    if (openMMContext.hasContextPointer()) {
      OpenMM_CustomBondForce_updateParametersInContext(pointer, openMMContext.getPointer());
    }
  }

  /**
   * Destroy the OpenMM CustomBondForce.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomBondForce_destroy(pointer);
      pointer = null;
    }
  }

}
