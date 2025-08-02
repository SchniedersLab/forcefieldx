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
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addAngle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addPerAngleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumAngles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getNumPerAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_getPerAngleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setPerAngleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions;

/**
 * Custom Angle Force.
 */
public class CustomAngleForce extends Force {

  /**
   * Create a CustomAngleForce.
   *
   * @param energy The energy expression for the force.
   */
  public CustomAngleForce(String energy) {
    super(OpenMM_CustomAngleForce_create(energy));
  }

  /**
   * Add an angle force to the OpenMM System.
   *
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param i3         The index of the third atom.
   * @param parameters The parameters for the angle.
   * @return The index of the angle that was added.
   */
  public int addAngle(int i1, int i2, int i3, DoubleArray parameters) {
    return OpenMM_CustomAngleForce_addAngle(pointer, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative of the energy with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomAngleForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name  The name of the parameter.
   * @param value The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double value) {
    return OpenMM_CustomAngleForce_addGlobalParameter(pointer, name, value);
  }

  /**
   * Add a per-angle parameter to the OpenMM System.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerAngleParameter(String name) {
    return OpenMM_CustomAngleForce_addPerAngleParameter(pointer, name);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomAngleForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for a specific angle.
   *
   * @param index      The index of the angle.
   * @param i1         The index of the first atom (output).
   * @param i2         The index of the second atom (output).
   * @param i3         The index of the third atom (output).
   * @param parameters The parameters for the angle (output).
   */
  public void getAngleParameters(int index, IntBuffer i1, IntBuffer i2, IntBuffer i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_getAngleParameters(pointer, index, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Get the parameters for a specific angle.
   *
   * @param index      The index of the angle.
   * @param i1         The index of the first atom (output).
   * @param i2         The index of the second atom (output).
   * @param i3         The index of the third atom (output).
   * @param parameters The parameters for the angle (output).
   */
  public void getAngleParameters(int index, IntByReference i1, IntByReference i2, IntByReference i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_getAngleParameters(pointer, index, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Get the energy function expression.
   *
   * @return The energy function expression.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomAngleForce_getEnergyFunction(pointer);
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
    Pointer p = OpenMM_CustomAngleForce_getEnergyParameterDerivativeName(pointer, index);
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
    return OpenMM_CustomAngleForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomAngleForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of angles.
   *
   * @return The number of angles.
   */
  public int getNumAngles() {
    return OpenMM_CustomAngleForce_getNumAngles(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomAngleForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of per-angle parameters.
   *
   * @return The number of per-angle parameters.
   */
  public int getNumPerAngleParameters() {
    return OpenMM_CustomAngleForce_getNumPerAngleParameters(pointer);
  }

  /**
   * Get the name of a per-angle parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerAngleParameterName(int index) {
    Pointer p = OpenMM_CustomAngleForce_getPerAngleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for one angle in the OpenMM System.
   *
   * @param index      The index of the angle.
   * @param i1         The index of the first atom.
   * @param i2         The index of the second atom.
   * @param i3         The index of the third atom.
   * @param parameters The angle parameters.
   */
  public void setAngleParameters(int index, int i1, int i2, int i3, DoubleArray parameters) {
    OpenMM_CustomAngleForce_setAngleParameters(pointer, index, i1, i2, i3, parameters.getPointer());
  }

  /**
   * Set the energy function expression.
   *
   * @param energy The energy function expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomAngleForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @param value The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMM_CustomAngleForce_setGlobalParameterDefaultValue(pointer, index, value);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomAngleForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a per-angle parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerAngleParameterName(int index, String name) {
    OpenMM_CustomAngleForce_setPerAngleParameterName(pointer, index, name);
  }

  /**
   * Set whether this force uses periodic boundary conditions.
   *
   * @param periodic If true, periodic boundary conditions will be used.
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions(pointer, periodic ? 1 : 0);
  }

  /**
   * Update the parameters in the OpenMM Context.
   *
   * @param context The OpenMM Context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomAngleForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}