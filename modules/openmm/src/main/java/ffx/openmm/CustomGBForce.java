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
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addComputedValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyParameterDerivative;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyTerm;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getComputedValueParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getEnergyParameterDerivativeName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getEnergyTermParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumComputedValues;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumEnergyParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumEnergyTerms;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumGlobalParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumPerParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setComputedValueParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setEnergyTermParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setExclusionParticles;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setFunctionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setGlobalParameterDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setGlobalParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setPerParticleParameterName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_usesPeriodicBoundaryConditions;

/**
 * Custom GB Force.
 */
public class CustomGBForce extends Force {

  /**
   * Create a CustomGBForce.
   */
  public CustomGBForce() {
    super(OpenMM_CustomGBForce_create());
  }

  /**
   * Add a computed value to the force.
   *
   * @param name       The name of the computed value.
   * @param expression The expression for computing the value.
   * @param type       The type of computation to perform.
   * @return The index of the computed value that was added.
   */
  public int addComputedValue(String name, String expression, int type) {
    return OpenMM_CustomGBForce_addComputedValue(pointer, name, expression, type);
  }

  /**
   * Add an energy parameter derivative to the force.
   *
   * @param name The name of the parameter to compute the derivative of the energy with respect to.
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMM_CustomGBForce_addEnergyParameterDerivative(pointer, name);
  }

  /**
   * Add an energy term to the force.
   *
   * @param expression The expression for computing the energy term.
   * @param type       The type of computation to perform.
   * @return The index of the energy term that was added.
   */
  public int addEnergyTerm(String expression, int type) {
    return OpenMM_CustomGBForce_addEnergyTerm(pointer, expression, type);
  }

  /**
   * Add an exclusion to the force.
   *
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   * @return The index of the exclusion that was added.
   */
  public int addExclusion(int particle1, int particle2) {
    return OpenMM_CustomGBForce_addExclusion(pointer, particle1, particle2);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name   The name of the function as it appears in expressions.
   * @param values The tabulated values of the function.
   * @param min    The minimum value of the independent variable for which the function is defined.
   * @param max    The maximum value of the independent variable for which the function is defined.
   * @return The index of the function that was added.
   * @deprecated This method exists only for backward compatibility. Use addTabulatedFunction() instead.
   */
  @Deprecated
  public int addFunction(String name, PointerByReference values, double min, double max) {
    return OpenMM_CustomGBForce_addFunction(pointer, name, values, min, max);
  }

  /**
   * Add a global parameter to the force.
   *
   * @param name         The name of the parameter.
   * @param defaultValue The default value of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMM_CustomGBForce_addGlobalParameter(pointer, name, defaultValue);
  }

  /**
   * Add a particle to the force.
   *
   * @param particleParameters The parameters for the particle.
   * @return The index of the particle that was added.
   */
  public int addParticle(DoubleArray particleParameters) {
    return OpenMM_CustomGBForce_addParticle(pointer, particleParameters.getPointer());
  }

  /**
   * Add a per-particle parameter to the force.
   *
   * @param name The name of the parameter.
   * @return The index of the parameter that was added.
   */
  public int addPerParticleParameter(String name) {
    return OpenMM_CustomGBForce_addPerParticleParameter(pointer, name);
  }

  /**
   * Add a tabulated function that may appear in the energy expression.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomGBForce_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Destroy the force.
   */
  @Override
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomGBForce_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * Get the parameters for a computed value.
   *
   * @param index      The index of the computed value.
   * @param name       The name of the computed value (output).
   * @param expression The expression for computing the value (output).
   * @param type       The type of computation to perform (output).
   */
  public void getComputedValueParameters(int index, PointerByReference name, PointerByReference expression, IntBuffer type) {
    OpenMM_CustomGBForce_getComputedValueParameters(pointer, index, name, expression, type);
  }

  /**
   * Get the parameters for a computed value.
   *
   * @param index      The index of the computed value.
   * @param name       The name of the computed value (output).
   * @param expression The expression for computing the value (output).
   * @param type       The type of computation to perform (output).
   */
  public void getComputedValueParameters(int index, PointerByReference name, PointerByReference expression, IntByReference type) {
    OpenMM_CustomGBForce_getComputedValueParameters(pointer, index, name, expression, type);
  }

  /**
   * Get the cutoff distance.
   *
   * @return The cutoff distance.
   */
  public double getCutoffDistance() {
    return OpenMM_CustomGBForce_getCutoffDistance(pointer);
  }

  /**
   * Get the name of a parameter with respect to which the derivative of the energy should be computed.
   *
   * @param index The index of the parameter derivative.
   * @return The name of the parameter.
   */
  public String getEnergyParameterDerivativeName(int index) {
    Pointer p = OpenMM_CustomGBForce_getEnergyParameterDerivativeName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the parameters for an energy term.
   *
   * @param index      The index of the energy term.
   * @param expression The expression for computing the energy term (output).
   * @param type       The type of computation to perform (output).
   */
  public void getEnergyTermParameters(int index, PointerByReference expression, IntBuffer type) {
    OpenMM_CustomGBForce_getEnergyTermParameters(pointer, index, expression, type);
  }

  /**
   * Get the parameters for an energy term.
   *
   * @param index      The index of the energy term.
   * @param expression The expression for computing the energy term (output).
   * @param type       The type of computation to perform (output).
   */
  public void getEnergyTermParameters(int index, PointerByReference expression, IntByReference type) {
    OpenMM_CustomGBForce_getEnergyTermParameters(pointer, index, expression, type);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntBuffer particle1, IntBuffer particle2) {
    OpenMM_CustomGBForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the particles in an exclusion.
   *
   * @param index     The index of the exclusion.
   * @param particle1 The index of the first particle in the exclusion (output).
   * @param particle2 The index of the second particle in the exclusion (output).
   */
  public void getExclusionParticles(int index, IntByReference particle1, IntByReference particle2) {
    OpenMM_CustomGBForce_getExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function (output).
   * @param values The tabulated values of the function (output).
   * @param min    The minimum value of the independent variable for which the function is defined (output).
   * @param max    The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference values, DoubleBuffer min, DoubleBuffer max) {
    OpenMM_CustomGBForce_getFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Get the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function (output).
   * @param values The tabulated values of the function (output).
   * @param min    The minimum value of the independent variable for which the function is defined (output).
   * @param max    The maximum value of the independent variable for which the function is defined (output).
   */
  public void getFunctionParameters(int index, PointerByReference name, PointerByReference values, DoubleByReference min, DoubleByReference max) {
    OpenMM_CustomGBForce_getFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The default value of the parameter.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMM_CustomGBForce_getGlobalParameterDefaultValue(pointer, index);
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getGlobalParameterName(int index) {
    Pointer p = OpenMM_CustomGBForce_getGlobalParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the nonbonded method.
   *
   * @return The nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMM_CustomGBForce_getNonbondedMethod(pointer);
  }

  /**
   * Get the number of computed values.
   *
   * @return The number of computed values.
   */
  public int getNumComputedValues() {
    return OpenMM_CustomGBForce_getNumComputedValues(pointer);
  }

  /**
   * Get the number of parameters with respect to which the derivative of the energy should be computed.
   *
   * @return The number of parameters.
   */
  public int getNumEnergyParameterDerivatives() {
    return OpenMM_CustomGBForce_getNumEnergyParameterDerivatives(pointer);
  }

  /**
   * Get the number of energy terms.
   *
   * @return The number of energy terms.
   */
  public int getNumEnergyTerms() {
    return OpenMM_CustomGBForce_getNumEnergyTerms(pointer);
  }

  /**
   * Get the number of exclusions.
   *
   * @return The number of exclusions.
   */
  public int getNumExclusions() {
    return OpenMM_CustomGBForce_getNumExclusions(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   * @deprecated This method exists only for backward compatibility. Use getNumTabulatedFunctions() instead.
   */
  @Deprecated
  public int getNumFunctions() {
    return OpenMM_CustomGBForce_getNumFunctions(pointer);
  }

  /**
   * Get the number of global parameters.
   *
   * @return The number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMM_CustomGBForce_getNumGlobalParameters(pointer);
  }

  /**
   * Get the number of particles.
   *
   * @return The number of particles.
   */
  public int getNumParticles() {
    return OpenMM_CustomGBForce_getNumParticles(pointer);
  }

  /**
   * Get the number of per-particle parameters.
   *
   * @return The number of per-particle parameters.
   */
  public int getNumPerParticleParameters() {
    return OpenMM_CustomGBForce_getNumPerParticleParameters(pointer);
  }

  /**
   * Get the number of tabulated functions.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomGBForce_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particleParameters The parameters for the particle (output).
   */
  public void getParticleParameters(int index, DoubleArray particleParameters) {
    OpenMM_CustomGBForce_getParticleParameters(pointer, index, particleParameters.getPointer());
  }

  /**
   * Get the name of a per-particle parameter.
   *
   * @param index The index of the parameter.
   * @return The name of the parameter.
   */
  public String getPerParticleParameterName(int index) {
    Pointer p = OpenMM_CustomGBForce_getPerParticleParameterName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get a reference to a tabulated function.
   *
   * @param index The index of the function.
   * @return A reference to the function.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomGBForce_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function.
   *
   * @param index The index of the function.
   * @return The name of the function.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomGBForce_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the parameters for a computed value.
   *
   * @param index      The index of the computed value.
   * @param name       The name of the computed value.
   * @param expression The expression for computing the value.
   * @param type       The type of computation to perform.
   */
  public void setComputedValueParameters(int index, String name, String expression, int type) {
    OpenMM_CustomGBForce_setComputedValueParameters(pointer, index, name, expression, type);
  }

  /**
   * Set the cutoff distance.
   *
   * @param distance The cutoff distance.
   */
  public void setCutoffDistance(double distance) {
    OpenMM_CustomGBForce_setCutoffDistance(pointer, distance);
  }

  /**
   * Set the parameters for an energy term.
   *
   * @param index      The index of the energy term.
   * @param expression The expression for computing the energy term.
   * @param type       The type of computation to perform.
   */
  public void setEnergyTermParameters(int index, String expression, int type) {
    OpenMM_CustomGBForce_setEnergyTermParameters(pointer, index, expression, type);
  }

  /**
   * Set the particles in an exclusion.
   *
   * @param index     The index of the exclusion.
   * @param particle1 The index of the first particle in the exclusion.
   * @param particle2 The index of the second particle in the exclusion.
   */
  public void setExclusionParticles(int index, int particle1, int particle2) {
    OpenMM_CustomGBForce_setExclusionParticles(pointer, index, particle1, particle2);
  }

  /**
   * Set the parameters for a tabulated function.
   *
   * @param index  The index of the function.
   * @param name   The name of the function.
   * @param values The tabulated values of the function.
   * @param min    The minimum value of the independent variable for which the function is defined.
   * @param max    The maximum value of the independent variable for which the function is defined.
   */
  public void setFunctionParameters(int index, String name, PointerByReference values, double min, double max) {
    OpenMM_CustomGBForce_setFunctionParameters(pointer, index, name, values, min, max);
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index        The index of the parameter.
   * @param defaultValue The default value of the parameter.
   */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMM_CustomGBForce_setGlobalParameterDefaultValue(pointer, index, defaultValue);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMM_CustomGBForce_setGlobalParameterName(pointer, index, name);
  }

  /**
   * Set the nonbonded method.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_CustomGBForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the parameters for a particle.
   *
   * @param index              The index of the particle.
   * @param particleParameters The parameters for the particle.
   */
  public void setParticleParameters(int index, DoubleArray particleParameters) {
    OpenMM_CustomGBForce_setParticleParameters(pointer, index, particleParameters.getPointer());
  }

  /**
   * Set the name of a per-particle parameter.
   *
   * @param index The index of the parameter.
   * @param name  The name of the parameter.
   */
  public void setPerParticleParameterName(int index, String name) {
    OpenMM_CustomGBForce_setPerParticleParameterName(pointer, index, name);
  }

  /**
   * Update the parameters in the context.
   *
   * @param context The context.
   */
  public void updateParametersInContext(Context context) {
    if (context.hasContextPointer()) {
      OpenMM_CustomGBForce_updateParametersInContext(pointer, context.getPointer());
    }
  }

  /**
   * Check if the force uses periodic boundary conditions.
   *
   * @return True if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    int pbc = OpenMM_CustomGBForce_usesPeriodicBoundaryConditions(pointer);
    return pbc == OpenMM_True;
  }
}