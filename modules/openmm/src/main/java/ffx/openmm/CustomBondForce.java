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
import com.sun.jna.ptr.IntByReference;

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getNumPerBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_getPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setPerBondParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_usesPeriodicBoundaryConditions;

/**
 * Custom Bond Force.
 */
public class CustomBondForce extends Force {

  /**
   * Create a CustomBondForce.
   *
   * @param energy The energy expression for the force.
   */
  public CustomBondForce(String energy) {
    pointer = OpenMM_CustomBondForce_create(energy);
  }

  /**
   * Add a bond to the OpenMM System.
   *
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param parameters The bond parameters.
   * @return The index of the bond that was added.
   */
  public int addBond(int i1, int i2, DoubleArray parameters) {
    return OpenMM_CustomBondForce_addBond(pointer, i1, i2, parameters.getPointer());
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative of the energy with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomBondForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a global parameter to the CustomBondForce.
   *
   * @param name  The name of the parameter.
   * @param value The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double value) {
    return OpenMM_CustomBondForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Add a per-bond parameter to the CustomBondForce.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerBondParameter(String name) {
    return OpenMM_CustomBondForce_addPerBondParameter(pointer, name);
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

  /**
   * Get the parameters for a specific bond.
   *
   * @param index      The index of the bond.
   * @param i1         The index of the first atom (output).
   * @param i2         The index of the second atom (output).
   * @param parameters The parameters for the bond (output).
   */
  public void getBondParameters(int index, IntBuffer i1, IntBuffer i2, DoubleArray parameters) {
    OpenMM_CustomBondForce_getBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Get the parameters for a specific bond.
   *
   * @param index      The index of the bond.
   * @param i1         The index of the first atom (output).
   * @param i2         The index of the second atom (output).
   * @param parameters The parameters for the bond (output).
   */
  public void getBondParameters(int index, IntByReference i1, IntByReference i2, DoubleArray parameters) {
    OpenMM_CustomBondForce_getBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Get the energy function expression.
   *
   * @return The energy function expression.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomBondForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of a parameter to compute the derivative of the energy with respect to.
   *
   * @param index The index of the parameter derivative.
   * @return The name of the parameter.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomBondForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomBondForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomBondForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of bonds.
   *
   * @return The number of bonds.
   */
  public int getNumBonds() {
    return OpenMM_CustomBondForce_getNumBonds(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomBondForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomBondForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of per-bond parameters.
   *
   * @return The number of per-bond parameters.
   */
  public int getNumPerBondParameters() {
    return OpenMM_CustomBondForce_getNumPerBondParameters(pointer);
  }

  /**
   * Get the name of a per-bond parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerBondParameterName(int index) {
    Pointer p = OpenMM_CustomBondForce_getPerBondParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for one bond in the OpenMM System.
   *
   * @param index      The index of the bond.
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param parameters The bond parameters.
   */
  public void setBondParameters(int index, int i1, int i2, DoubleArray parameters) {
    OpenMM_CustomBondForce_setBondParameters(pointer, index, i1, i2, parameters.getPointer());
  }

  /**
   * Set the energy function expression.
   *
   * @param energy The energy function expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomBondForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @param value The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMM_CustomBondForce_setGlobalParameterDefaultValue(pointer, index, value);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomBondForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-bond parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerBondParameterName(int index, String name) {
    OpenMM_CustomBondForce_setPerBondParameterName(pointer, index, name);
  }

  /**
   * Set whether this force uses periodic boundary conditions.
   *
   * @param periodic If true, periodic boundary conditions will be used.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomBondForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomBondForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}