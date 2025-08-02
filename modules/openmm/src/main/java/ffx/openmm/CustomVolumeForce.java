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

import com.sun.jna.Pointer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomVolumeForce_usesPeriodicBoundaryConditions;

/**
 * This class implements forces that depend on the volume of the periodic box.
 * The energy is computed as a user-supplied algebraic expression that may depend
 * on the volume of the periodic box and on global parameters.
 * <p>
 * To use this class, create a CustomVolumeForce object, passing an algebraic expression to the
 * constructor that defines the energy as a function of the volume. Then call addGlobalParameter()
 * to define global parameters that the energy may depend on.
 * <p>
 * The energy expression may use the variable "volume" to refer to the volume of the periodic box.
 * It may also use any global parameters that have been defined.
 */
public class CustomVolumeForce extends Force {

  /**
   * Create a CustomVolumeForce.
   *
   * @param energy The algebraic expression that gives the energy as a function of the volume.
   */
  public CustomVolumeForce(String energy) {
    pointer = OpenMM_CustomVolumeForce_create(energy);
  }

  /**
   * Add a new global parameter that the energy may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomVolumeForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a new global parameter that the energy may depend on.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_CustomVolumeForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Destroy the force.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomVolumeForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the algebraic expression that gives the energy as a function of the volume.
   *
   * @return The energy expression.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomVolumeForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter for which to get the default value.
   * @return The parameter default value.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomVolumeForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter for which to get the name.
   * @return The parameter name.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomVolumeForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of global parameters that the energy depends on.
   *
   * @return The number of parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomVolumeForce_getNumGlobalParameters(pointer);
  }

  /**
   * Set the algebraic expression that gives the energy as a function of the volume.
   *
   * @param energy The energy expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomVolumeForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the algebraic expression that gives the energy as a function of the volume.
   *
   * @param energy The energy expression.
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_CustomVolumeForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter for which to set the default value.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomVolumeForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomVolumeForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter for which to set the name.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_CustomVolumeForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomVolumeForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}