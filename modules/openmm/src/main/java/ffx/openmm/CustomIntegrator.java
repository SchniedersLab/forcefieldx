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
import com.sun.jna.ptr.PointerByReference;

import java.nio.IntBuffer;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addComputeGlobal;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addComputePerDof;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addComputeSum;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addGlobalVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addPerDofVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addUpdateContextState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_beginIfBlock;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_beginWhileBlock;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_endBlock;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getComputationStep;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getGlobalVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getGlobalVariableByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getGlobalVariableName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getKineticEnergyExpression;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getNumComputations;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getNumGlobalVariables;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getNumPerDofVariables;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getNumTabulatedFunctions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getPerDofVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getPerDofVariableByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getPerDofVariableName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getTabulatedFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_getTabulatedFunctionName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setGlobalVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setGlobalVariableByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setKineticEnergyExpression;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setPerDofVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setPerDofVariableByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_step;

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
   * Add a computation that computes a global value.
   *
   * @param variable   The name of the variable to store the computed value in.
   * @param expression The expression to evaluate.
   * @return The index of the computation that was added.
   */
  public int addComputeGlobal(String variable, String expression) {
    return OpenMM_CustomIntegrator_addComputeGlobal(pointer, variable, expression);
  }

  /**
   * Add a per-DOF computation to this Integrator.
   *
   * @param variable   The name of the variable to store the computed value in.
   * @param expression The expression to evaluate.
   * @return The index of the computation that was added.
   */
  public int addComputePerDof(String variable, String expression) {
    return OpenMM_CustomIntegrator_addComputePerDof(pointer, variable, expression);
  }

  /**
   * Add a computation that computes a sum over degrees of freedom.
   *
   * @param variable   The name of the variable to store the computed value in.
   * @param expression The expression to evaluate and sum.
   * @return The index of the computation that was added.
   */
  public int addComputeSum(String variable, String expression) {
    return OpenMM_CustomIntegrator_addComputeSum(pointer, variable, expression);
  }

  /**
   * Add a position constraint to this Integrator.
   *
   * @return The index of the computation that was added.
   */
  public int addConstrainPositions() {
    return OpenMM_CustomIntegrator_addConstrainPositions(pointer);
  }

  /**
   * Add a velocity constraint to this Integrator.
   *
   * @return The index of the computation that was added.
   */
  public int addConstrainVelocities() {
    return OpenMM_CustomIntegrator_addConstrainVelocities(pointer);
  }

  /**
   * Add a global variable to this Integrator.
   *
   * @param name         The name of the variable to create.
   * @param initialValue The initial value of the variable.
   * @return The index of the variable that was added.
   */
  public int addGlobalVariable(String name, double initialValue) {
    return OpenMM_CustomIntegrator_addGlobalVariable(pointer, name, initialValue);
  }

  /**
   * Add a per-DOF variable to this Integrator.
   *
   * @param name         The name of the variable to create.
   * @param initialValue The initial value of the variable.
   * @return The index of the variable that was added.
   */
  public int addPerDofVariable(String name, double initialValue) {
    return OpenMM_CustomIntegrator_addPerDofVariable(pointer, name, initialValue);
  }

  /**
   * Add a tabulated function that may appear in expressions.
   *
   * @param name     The name of the function as it appears in expressions.
   * @param function A TabulatedFunction object defining the function.
   * @return The index of the function that was added.
   */
  public int addTabulatedFunction(String name, PointerByReference function) {
    return OpenMM_CustomIntegrator_addTabulatedFunction(pointer, name, function);
  }

  /**
   * Add an update context state to this Integrator.
   *
   * @return The index of the computation that was added.
   */
  public int addUpdateContextState() {
    return OpenMM_CustomIntegrator_addUpdateContextState(pointer);
  }

  /**
   * Begin a new if block.
   *
   * @param condition The condition under which the block should be executed.
   * @return The index of the computation that was added.
   */
  public int beginIfBlock(String condition) {
    return OpenMM_CustomIntegrator_beginIfBlock(pointer, condition);
  }

  /**
   * Begin a new while block.
   *
   * @param condition The condition under which the block should be executed repeatedly.
   * @return The index of the computation that was added.
   */
  public int beginWhileBlock(String condition) {
    return OpenMM_CustomIntegrator_beginWhileBlock(pointer, condition);
  }

  /**
   * Destroy the integrator.
   */
  public void destroy() {
    if (pointer != null) {
      OpenMM_CustomIntegrator_destroy(pointer);
      pointer = null;
    }
  }

  /**
   * End the most recently begun if or while block.
   *
   * @return The index of the computation that was added.
   */
  public int endBlock() {
    return OpenMM_CustomIntegrator_endBlock(pointer);
  }

  /**
   * Get information about a computation in the integrator.
   *
   * @param index      The index of the computation to get.
   * @param type       The type of the computation (output).
   * @param variable   The variable the computation is stored in (output).
   * @param expression The expression to evaluate (output).
   */
  public void getComputationStep(int index, IntByReference type, PointerByReference variable, PointerByReference expression) {
    OpenMM_CustomIntegrator_getComputationStep(pointer, index, type, variable, expression);
  }

  /**
   * Get information about a computation in the integrator.
   *
   * @param index      The index of the computation to get.
   * @param type       The type of the computation (output).
   * @param variable   The variable the computation is stored in (output).
   * @param expression The expression to evaluate (output).
   */
  public void getComputationStep(int index, IntBuffer type, PointerByReference variable, PointerByReference expression) {
    OpenMM_CustomIntegrator_getComputationStep(pointer, index, type, variable, expression);
  }

  /**
   * Get the value of a global variable.
   *
   * @param index The index of the variable to get.
   * @return The value of the variable.
   */
  public double getGlobalVariable(int index) {
    return OpenMM_CustomIntegrator_getGlobalVariable(pointer, index);
  }

  /**
   * Get the value of a global variable, specified by name.
   *
   * @param name The name of the variable to get.
   * @return The value of the variable.
   */
  public double getGlobalVariableByName(String name) {
    return OpenMM_CustomIntegrator_getGlobalVariableByName(pointer, name);
  }

  /**
   * Get the name of a global variable.
   *
   * @param index The index of the variable to get.
   * @return The name of the variable.
   */
  public String getGlobalVariableName(int index) {
    Pointer p = OpenMM_CustomIntegrator_getGlobalVariableName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the expression used to compute the kinetic energy.
   *
   * @return The expression used to compute the kinetic energy.
   */
  public String getKineticEnergyExpression() {
    Pointer p = OpenMM_CustomIntegrator_getKineticEnergyExpression(pointer);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the number of computation steps defined for this integrator.
   *
   * @return The number of computation steps.
   */
  public int getNumComputations() {
    return OpenMM_CustomIntegrator_getNumComputations(pointer);
  }

  /**
   * Get the number of global variables that have been defined.
   *
   * @return The number of global variables.
   */
  public int getNumGlobalVariables() {
    return OpenMM_CustomIntegrator_getNumGlobalVariables(pointer);
  }

  /**
   * Get the number of per-DOF variables that have been defined.
   *
   * @return The number of per-DOF variables.
   */
  public int getNumPerDofVariables() {
    return OpenMM_CustomIntegrator_getNumPerDofVariables(pointer);
  }

  /**
   * Get the number of tabulated functions that have been defined.
   *
   * @return The number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMM_CustomIntegrator_getNumTabulatedFunctions(pointer);
  }

  /**
   * Get the values of a per-DOF variable.
   *
   * @param index    The index of the variable to get.
   * @param variable The values of the variable (output).
   */
  public void getPerDofVariable(int index, PointerByReference variable) {
    OpenMM_CustomIntegrator_getPerDofVariable(pointer, index, variable);
  }

  /**
   * Get the values of a per-DOF variable, specified by name.
   *
   * @param name     The name of the variable to get.
   * @param variable The values of the variable (output).
   */
  public void getPerDofVariableByName(String name, PointerByReference variable) {
    OpenMM_CustomIntegrator_getPerDofVariableByName(pointer, name, variable);
  }

  /**
   * Get the name of a per-DOF variable.
   *
   * @param index The index of the variable to get.
   * @return The name of the variable.
   */
  public String getPerDofVariableName(int index) {
    Pointer p = OpenMM_CustomIntegrator_getPerDofVariableName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Get the random number seed. See setRandomNumberSeed() for details.
   *
   * @return The random number seed.
   */
  public int getRandomNumberSeed() {
    return OpenMM_CustomIntegrator_getRandomNumberSeed(pointer);
  }

  /**
   * Get a reference to a tabulated function that may appear in expressions.
   *
   * @param index The index of the function to get.
   * @return The TabulatedFunction object defining the function.
   */
  public PointerByReference getTabulatedFunction(int index) {
    return OpenMM_CustomIntegrator_getTabulatedFunction(pointer, index);
  }

  /**
   * Get the name of a tabulated function that may appear in expressions.
   *
   * @param index The index of the function to get.
   * @return The name of the function as it appears in expressions.
   */
  public String getTabulatedFunctionName(int index) {
    Pointer p = OpenMM_CustomIntegrator_getTabulatedFunctionName(pointer, index);
    if (p == null) {
      return null;
    }
    return p.getString(0);
  }

  /**
   * Set the value of a global variable.
   *
   * @param index The index of the variable to set.
   * @param value The new value of the variable.
   */
  public void setGlobalVariable(int index, double value) {
    OpenMM_CustomIntegrator_setGlobalVariable(pointer, index, value);
  }

  /**
   * Set the value of a global variable, specified by name.
   *
   * @param name  The name of the variable to set.
   * @param value The new value of the variable.
   */
  public void setGlobalVariableByName(String name, double value) {
    OpenMM_CustomIntegrator_setGlobalVariableByName(pointer, name, value);
  }

  /**
   * Set the expression used to compute the kinetic energy.
   *
   * @param expression The expression used to compute the kinetic energy.
   */
  public void setKineticEnergyExpression(String expression) {
    OpenMM_CustomIntegrator_setKineticEnergyExpression(pointer, expression);
  }

  /**
   * Set the values of a per-DOF variable.
   *
   * @param index    The index of the variable to set.
   * @param variable The new values of the variable.
   */
  public void setPerDofVariable(int index, PointerByReference variable) {
    OpenMM_CustomIntegrator_setPerDofVariable(pointer, index, variable);
  }

  /**
   * Set the values of a per-DOF variable, specified by name.
   *
   * @param name     The name of the variable to set.
   * @param variable The new values of the variable.
   */
  public void setPerDofVariableByName(String name, PointerByReference variable) {
    OpenMM_CustomIntegrator_setPerDofVariableByName(pointer, name, variable);
  }

  /**
   * Set the random number seed. The precise meaning of this parameter is undefined, and is left up
   * to each Platform to interpret in an appropriate way. It is guaranteed that if two simulations
   * are run with different random number seeds, the sequence of random numbers will be different.
   * On the other hand, no guarantees are made about the behavior of simulations that use the same
   * seed. In particular, Platforms are permitted to use non-deterministic algorithms which produce
   * different results on successive runs, even if those runs were initialized identically.
   * <p>
   * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a
   * Context is created from this Force. This is done to ensure that each Context receives unique
   * random seeds without you needing to set them explicitly.
   *
   * @param seed The random number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMM_CustomIntegrator_setRandomNumberSeed(pointer, seed);
  }

  /**
   * Advance a simulation through time by taking a series of time steps.
   *
   * @param steps The number of time steps to take.
   */
  public void step(int steps) {
    OpenMM_CustomIntegrator_step(pointer, steps);
  }
}