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
import com.sun.jna.ptr.PointerByReference;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_addCollectiveVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getCollectiveVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getCollectiveVariableName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getCollectiveVariableValues;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getInnerContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getNumCollectiveVariables;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_setEnergyFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCVForce_usesPeriodicBoundaryConditions;

/**
 * This class supports energy functions that depend on collective variables. To use it,
 * you define a set of collective variables (scalar valued functions that depend on the
 * particle positions), and an algebraic expression for the energy as a function of the
 * collective variables. The expression also may involve tabulated functions, and may
 * depend on arbitrary global parameters.
 * <p>
 * Each collective variable is defined by a Force object. The Force's potential energy
 * is computed, and that becomes the value of the variable. This provides enormous
 * flexibility in defining collective variables, especially by using custom forces.
 * Anything that can be computed as a potential function can also be used as a collective
 * variable.
 * <p>
 * To use this class, create a CustomCVForce object, passing an algebraic expression to the
 * constructor that defines the potential energy. Then call addCollectiveVariable() to define
 * collective variables and addGlobalParameter() to define global parameters. The values
 * of global parameters may be modified during a simulation by calling Context::setParameter().
 * <p>
 * This class also has the ability to compute derivatives of the potential energy with respect to global parameters.
 * Call addEnergyParameterDerivative() to request that the derivative with respect to a particular parameter be
 * computed. You can then query its value in a Context by calling getState() on it.
 * <p>
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta, select. All trigonometric functions
 * are defined in radians, and log is the natural logarithm. step(x) = 0 if x is less than 0, 1 otherwise. delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 * <p>
 * In addition, you can call addTabulatedFunction() to define a new function based on tabulated values. You specify the function by
 * creating a TabulatedFunction object. That function can then appear in the expression.
 */
public class CustomCVForce extends Force {

  /**
   * Create a new CustomCVForce.
   *
   * @param energy The energy function as an algebraic expression.
   */
  public CustomCVForce(String energy) {
    super(OpenMM_CustomCVForce_create(energy));
  }

  /**
   * Add a collective variable to the force.
   *
   * @param name  The name of the collective variable.
   * @param force The Force object that defines the collective variable.
   * @return The index of the collective variable that was added.
   */
  public int addCollectiveVariable(String name, PointerByReference force) {
    return OpenMM_CustomCVForce_addCollectiveVariable(pointer, name, force);
  }

  /**
   * Add a collective variable to the force.
   *
   * @param name  The name of the collective variable.
   * @param force The Force object that defines the collective variable.
   * @return The index of the collective variable that was added.
   */
  public int addCollectiveVariable(Pointer name, PointerByReference force) {
    return OpenMM_CustomCVForce_addCollectiveVariable(pointer, name, force);
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomCVForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative with respect to.
   */
  public void addEnergyParameterDerivative(Pointer name) {
    OpenMM_CustomCVForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomCVForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(Pointer name, double defaultValue) {
    return OpenMM_CustomCVForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a tabulated function to the force.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    return OpenMM_CustomCVForce_addTabulatedFunction(pointer, name, function.getPointer());
  }

  /**
   * Add a tabulated function to the force.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(Pointer name, TabulatedFunction function) {
    return OpenMM_CustomCVForce_addTabulatedFunction(pointer, name, function.getPointer());
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomCVForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get a collective variable by index.
   *
   * @param index The index of the collective variable.
   * @return The Force object that defines the collective variable.
   */
  public PointerByReference getCollectiveVariable(int index) {
    return OpenMM_CustomCVForce_getCollectiveVariable(pointer, index);
  }

  /**
   * Get the name of a collective variable.
   *
   * @param index The index of the collective variable.
   * @return The name of the collective variable.
   */
  public String getCollectiveVariableName(int index) {
    Pointer p = OpenMM_CustomCVForce_getCollectiveVariableName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the values of all collective variables in a given context.
   *
   * @param context The context for which to get the values.
   * @param values  The values of the collective variables (output).
   */
  public void getCollectiveVariableValues(Context context, PointerByReference values) {
    OpenMM_CustomCVForce_getCollectiveVariableValues(pointer, context.getPointer(), values);
  }

  /**
   * Get the energy function.
   *
   * @return The energy function as an algebraic expression.
   */
  public String getEnergyFunction() {
    Pointer p = OpenMM_CustomCVForce_getEnergyFunction(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the name of an energy parameter derivative.
   *
   * @param index The index of the parameter derivative.
   * @return The name of the parameter derivative.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomCVForce_getEnergyParameterDerivativeName(pointer, index);
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
    return OpenMM_CustomCVForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomCVForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the inner context used for evaluating collective variables.
   *
   * @param context The main context.
   * @return The inner context.
   */
  public PointerByReference getInnerContext(Context context) {
    return OpenMM_CustomCVForce_getInnerContext(pointer, context.getPointer());
  }

  /**
   * Get the number of collective variables.
   *
   * @return The number of collective variables.
   */
  public int getNumCollectiveVariables() {
    return OpenMM_CustomCVForce_getNumCollectiveVariables(pointer);
  }

  /**
   * Get the number of energy parameter derivatives.
   *
   * @return The number of energy parameter derivatives.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomCVForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomCVForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomCVForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get a tabulated function by index.
   *
   * @param index The index of the function.
   * @return The TabulatedFunction object.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomCVForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function.
   *
   * @param index The index of the function.
   * @return The name of the function as it appears in expressions.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomCVForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the energy function.
   *
   * @param energy The energy function as an algebraic expression.
   */
  public void setEnergyFunction(String energy) {
    OpenMM_CustomCVForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the energy function.
   *
   * @param energy The energy function as an algebraic expression.
   */
  public void setEnergyFunction(Pointer energy) {
    OpenMM_CustomCVForce_setEnergyFunction(pointer, energy);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomCVForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomCVForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, Pointer name) {
    OpenMM_CustomCVForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Update the parameters in a Context to match those stored in this Force object.
   *
   * @param context The Context in which to update the parameters.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomCVForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomCVForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}