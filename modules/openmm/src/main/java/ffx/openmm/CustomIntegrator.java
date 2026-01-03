// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
 * This is an Integrator that can be used to implemented arbitrary, user defined
 * integration algorithms.  It is flexible enough to support a wide range of
 * methods including both deterministic and stochastic integrators, Metropolized
 * integrators, and integrators that must integrate additional quantities along
 * with the particle positions and momenta.
 *
 * <p>To create an integration algorithm, you first define a set of variables the
 * integrator will compute.  Variables come in two types: global variables
 * have a single value, while per-DOF variables have a value for every
 * degree of freedom (x, y, or z coordinate of a particle).  You can define as
 * many variables as you want of each type.  The value of any variable can be
 * computed by the integration algorithm, or set directly by calling a method on
 * the CustomIntegrator.  All variables are persistent between integration
 * steps; once a value is set, it keeps that value until it is changed by the
 * user or recomputed in a later integration step.
 *
 * <p>Next, you define the algorithm as a series of computations.  To execute a
 * time step, the integrator performs the list of computations in order.  Each
 * computation updates the value of one global or per-DOF value.  There are
 * several types of computations that can be done:
 *
 * <ul>
 * <li>Global: You provide a mathematical expression involving only global
 * variables.  It is evaluated and stored into a global variable.</li>
 * <li>Per-DOF: You provide a mathematical expression involving both global and
 * per-DOF variables.  It is evaluated once for every degree of freedom, and
 * the values are stored into a per-DOF variable.</li>
 * <li>Sum: You provide a mathematical expression involving both global and
 * per-DOF variables.  It is evaluated once for every degree of freedom.  All
 * of those values are then added together, and the sum is stored into a global
 * variable.</li>
 * <li>Constrain Positions: The particle positions are updated so that all
 * distance constraints are satisfied.</li>
 * <li>Constrain Velocities: The particle velocities are updated so the net
 * velocity along any constrained distance is 0.</li>
 * </ul>
 *
 * <p>Like all integrators, CustomIntegrator ignores any particle whose mass is 0.
 * It is skipped when doing per-DOF computations, and is not included when
 * computing sums over degrees of freedom.
 *
 * <p>In addition to the variables you define by calling addGlobalVariable() and
 * addPerDofVariable(), the integrator provides the following pre-defined
 * variables:
 *
 * <ul>
 * <li>dt: (global) This is the step size being used by the integrator.</li>
 * <li>energy: (global, read-only) This is the current potential energy of the
 * system.</li>
 * <li>energy0, energy1, energy2, ...: (global, read-only) This is similar to
 * energy, but includes only the contribution from forces in one force group.
 * A single computation step may only depend on a single energy variable
 * (energy, energy0, energy1, etc.).</li>
 * <li>x: (per-DOF) This is the current value of the degree of freedom (the x,
 * y, or z coordinate of a particle).</li>
 * <li>v: (per-DOF) This is the current velocity associated with the degree of
 * freedom (the x, y, or z component of a particle's velocity).</li>
 * <li>f: (per-DOF, read-only) This is the current force acting on the degree of
 * freedom (the x, y, or z component of the force on a particle).</li>
 * <li>f0, f1, f2, ...: (per-DOF, read-only) This is similar to f, but includes
 * only the contribution from forces in one force group.  A single computation
 * step may only depend on a single force variable (f, f0, f1, etc.).</li>
 * <li>m: (per-DOF, read-only) This is the mass of the particle the degree of
 * freedom is associated with.</li>
 * <li>uniform: (either global or per-DOF, read-only) This is a uniformly
 * distributed random number between 0 and 1.  Every time an expression is
 * evaluated, a different value will be used.  When used in a per-DOF
 * expression, a different value will be used for every degree of freedom.
 * Note, however, that if this variable appears multiple times in a single
 * expression, the same value is used everywhere it appears in that
 * expression.</li>
 * <li>gaussian: (either global or per-DOF, read-only) This is a Gaussian
 * distributed random number with mean 0 and variance 1.  Every time an expression
 * is evaluated, a different value will be used.  When used in a per-DOF
 * expression, a different value will be used for every degree of freedom.
 * Note, however, that if this variable appears multiple times in a single
 * expression, the same value is used everywhere it appears in that
 * expression.</li>
 * <li>A global variable is created for every adjustable parameter defined
 * in the integrator's Context.</li>
 * </ul>
 *
 * <p>The following example uses a CustomIntegrator to implement a velocity Verlet
 * integrator:
 *
 * <pre>{@code
 * CustomIntegrator integrator = new CustomIntegrator(0.001);
 * integrator.addComputePerDof("v", "v+0.5*dt*f/m");
 * integrator.addComputePerDof("x", "x+dt*v");
 * integrator.addComputePerDof("v", "v+0.5*dt*f/m");
 * }</pre>
 *
 * <p>The first step updates the velocities based on the current forces.
 * The second step updates the positions based on the new velocities, and the
 * third step updates the velocities again.  Although the first and third steps
 * look identical, the forces used in them are different.  You do not need to
 * tell the integrator that; it will recognize that the positions have changed
 * and know to recompute the forces automatically.
 */
public class CustomIntegrator extends Integrator {

  /**
   * Constructor.
   *
   * @param dt The time step.
   */
  public CustomIntegrator(double dt) {
    super(OpenMM_CustomIntegrator_create(dt));
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
  @Override
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