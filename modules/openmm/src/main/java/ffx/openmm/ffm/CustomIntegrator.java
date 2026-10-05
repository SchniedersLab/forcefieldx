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
package ffx.openmm.ffm;

import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * An integrator whose algorithm is defined by user-supplied expressions, able to represent deterministic,
 * stochastic and Metropolized methods and integrators that also integrate additional quantities.
 *
 * <p>Define global variables ({@link #addGlobalVariable(String, double)}, one value) and per-degree-of-freedom
 * variables ({@link #addPerDofVariable(String, double)}, one value per x, y or z coordinate of each particle), then
 * define the algorithm as an ordered list of computation steps executed once per time step: global ({@link
 * #addComputeGlobal(String, String)}), per-DOF ({@link #addComputePerDof(String, String)}), sum ({@link
 * #addComputeSum(String, String)}), constraint of positions or velocities, a context-state update, and "if"/"while"
 * blocks. Variables persist between steps. Particles with mass 0 are ignored in per-DOF computations and sums.</p>
 *
 * <p>Predefined variables include {@code dt} (global step size), {@code energy} and {@code energy0, energy1, ...}
 * (read-only global energies; a step may depend on only one of them), {@code x}, {@code v}, {@code f} and {@code
 * f0, f1, ...} (per-DOF position, velocity, force; a step may depend on only one force variable), {@code m} (read-only
 * mass), read-only random numbers {@code uniform} and {@code gaussian}, and one global variable per adjustable parameter of
 * the context. Use {@link #addUpdateContextState()} each step so forces such as thermostats and barostats can act, and
 * {@link #setKineticEnergyExpression(String)} when the default {@code m*v*v/2} is not correct (for example for leapfrog
 * methods). The header also describes a {@code deriv()} function for parameter derivatives of the energy.</p>
 *
 * <p>Expressions use the operators + - * / and ^ and the functions sqrt, exp, log, sin, cos, sec, csc, tan, cot,
 * asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta and select
 * (trigonometric functions in radians, log natural); intermediate quantities may follow the main expression after
 * {@code ;}. Per-DOF and sum expressions may also use the vector functions cross, dot, _x, _y, _z and vector. See
 * the OpenMM header for the complete rules.</p>
 *
 * <p>Block conditions are evaluated globally and may use only global variables. Be careful with "while" blocks, since
 * nothing prevents infinite loops. Strings passed to native code are converted to temporary UTF-8 strings.</p>
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.CustomIntegrator} is not. Per-DOF
 * variables use {@link Vec3Array} instead of {@code PointerByReference}; {@link #getPerDofVariable(int)} and {@link
 * #getPerDofVariableByName(String)} return a new caller-owned array. {@link #getTabulatedFunction(int)} returns a
 * borrowed {@link MemorySegment}. {@link #getComputationStep(int)} is implemented by the JNA class but throws {@link
 * UnsupportedOperationException} here because of a C wrapper ABI problem described on that method.</p>
 */
public class CustomIntegrator extends Integrator {

  /**
   * Description of one computation step: its kind, target variable and expression.
   *
   * <p>This record models the data reported by OpenMM's {@code getComputationStep}, but {@link #getComputationStep(int)}
   * currently throws, so instances are not produced by this class. The native {@code type} values are 0 compute global,
   * 1 compute per-DOF, 2 compute sum, 3 constrain positions, 4 constrain velocities, 5 update context state, 6 "if"
   * block start, 7 "while" block start and 8 block end.</p>
   *
   * @param type       native computation-type integer, as listed above.
   * @param variable   variable the step stores its result into; the header says this is an empty string if the step
   *     stores no result.
   * @param expression expression the step evaluates; the header says this is an empty string if the step evaluates no
   *     expression.
   */
  public record ComputationStep(int type, String variable, String expression) {
  }

  /**
   * Create a custom integrator with no variables or computation steps. The FFM runtime is initialized first through
   * {@link OpenMMRuntime#initialize()}.
   *
   * @param stepSize step size with which to integrate the system, in ps; exposed to expressions as {@code dt}.
   */
  public CustomIntegrator(double stepSize) {
    super(create(stepSize));
  }

  /**
   * Add a step that computes a global value each integration step.
   *
   * <p>The expression may involve only global variables.</p>
   *
   * @param variable   name of the global variable that stores the computed value.
   * @param expression mathematical expression involving only global variables.
   * @return index of the computation step that was added.
   */
  public int addComputeGlobal(String variable, String expression) {
    return withPair(variable, expression, OpenMMNative::OpenMM_CustomIntegrator_addComputeGlobal);
  }

  /**
   * Add a step that computes a per-degree-of-freedom value each integration step.
   *
   * <p>The expression may involve global and per-DOF variables and is evaluated for every degree of freedom.</p>
   *
   * @param variable   name of the per-DOF variable that stores the computed value.
   * @param expression mathematical expression involving global and per-DOF variables.
   * @return index of the computation step that was added.
   */
  public int addComputePerDof(String variable, String expression) {
    return withPair(variable, expression, OpenMMNative::OpenMM_CustomIntegrator_addComputePerDof);
  }

  /**
   * Add a step that sums an expression over all degrees of freedom each integration step.
   *
   * <p>The expression may involve global and per-DOF variables; its per-DOF values are added and the sum is stored in a
   * global variable.</p>
   *
   * @param variable   name of the global variable that stores the sum.
   * @param expression mathematical expression involving global and per-DOF variables.
   * @return index of the computation step that was added.
   */
  public int addComputeSum(String variable, String expression) {
    return withPair(variable, expression, OpenMMNative::OpenMM_CustomIntegrator_addComputeSum);
  }

  /**
   * Add a step that updates particle positions so all constraints are satisfied.
   *
   * @return index of the computation step that was added.
   */
  public int addConstrainPositions() {
    return OpenMMNative.OpenMM_CustomIntegrator_addConstrainPositions(getPointer());
  }

  /**
   * Add a step that updates particle velocities so the net velocity along all constrained distances is 0.
   *
   * @return index of the computation step that was added.
   */
  public int addConstrainVelocities() {
    return OpenMMNative.OpenMM_CustomIntegrator_addConstrainVelocities(getPointer());
  }

  /**
   * Define a new global variable.
   *
   * @param name         variable name.
   * @param initialValue initial value of the variable.
   * @return index of the variable that was added.
   */
  public int addGlobalVariable(String name, double initialValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomIntegrator_addGlobalVariable(
            getPointer(), value, initialValue));
  }

  /**
   * Define a new per-degree-of-freedom variable.
   *
   * @param name         variable name.
   * @param initialValue initial value for all degrees of freedom.
   * @return index of the variable that was added.
   */
  public int addPerDofVariable(String name, double initialValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomIntegrator_addPerDofVariable(
            getPointer(), value, initialValue));
  }

  /**
   * Add a tabulated function that may appear in expressions.
   *
   * <p>Native ownership of the function transfers to the integrator, which deletes it when the integrator is deleted;
   * after the call the passed wrapper is invalidated (its handle is cleared without calling the native destructor) and must
   * not be reused or destroyed. The JNA method takes a {@code PointerByReference} and does not invalidate anything.</p>
   *
   * @param name     function name as it appears in expressions.
   * @param function live tabulated function to add; must not be null.
   * @return index of the function that was added.
   * @throws NullPointerException if {@code function} is null.
   * @throws IllegalStateException if the function or this integrator has been destroyed.
   */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomIntegrator_addTabulatedFunction(
            getPointer(), value, function.getPointer()));
    function.invalidate();
    return index;
  }

  /**
   * Add a step that allows forces to update the context state (for example an Andersen thermostat randomizing
   * velocities or a Monte Carlo barostat scaling positions).
   *
   * @return index of the computation step that was added.
   */
  public int addUpdateContextState() {
    return OpenMMNative.OpenMM_CustomIntegrator_addUpdateContextState(getPointer());
  }

  /**
   * Begin an "if" block. Steps up to the matching {@link #endBlock()} execute only if the condition is true.
   *
   * @param condition expression with a comparison operator ({@code =}, {@code <}, {@code >}, {@code !=}, {@code <=} or
   *     {@code >=}) involving only global variables.
   * @return index of the computation step that was added.
   */
  public int beginIfBlock(String condition) {
    return withString(condition, OpenMMNative::OpenMM_CustomIntegrator_beginIfBlock);
  }

  /**
   * Begin a "while" block. Steps up to the matching {@link #endBlock()} execute repeatedly while the condition remains
   * true; take care to avoid infinite loops.
   *
   * @param condition expression with a comparison operator involving only global variables.
   * @return index of the computation step that was added.
   */
  public int beginWhileBlock(String condition) {
    return withString(condition, OpenMMNative::OpenMM_CustomIntegrator_beginWhileBlock);
  }

  /**
   * Destroy the native integrator.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. If the
   * integrator was passed to a {@link Context}, that context's destruction also destroys the integrator, so do not call
   * this for an integrator that a live context owns.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomIntegrator_destroy);
  }

  /**
   * End the most recently begun "if" or "while" block.
   *
   * @return index of the computation step that was added.
   */
  public int endBlock() {
    return OpenMMNative.OpenMM_CustomIntegrator_endBlock(getPointer());
  }

  /**
   * Read the data of a computation step. This method is not supported and always throws.
   *
   * <p>The installed C wrapper declares the string outputs as {@code char**} but writes constructed C++ {@code std::string}
   * objects to those addresses, so it cannot safely be called through FFM; the wrapper must be fixed before this output
   * operation can be exposed. The JNA class provides {@code getComputationStep(int, IntByReference, PointerByReference,
   * PointerByReference)} and an {@code IntBuffer} variant. This is a documented wrapper/header mismatch: the header
   * returns type, variable and expression through output parameters.</p>
   *
   * @param index index of the computation step.
   * @return never returns normally.
   * @throws UnsupportedOperationException always.
   */
  public ComputationStep getComputationStep(int index) {
    throw new UnsupportedOperationException(
        "OpenMM_CustomIntegrator_getComputationStep has an incompatible C wrapper ABI.");
  }

  /**
   * Get the current value of a global variable.
   *
   * @param index index of the variable.
   * @return current value.
   */
  public double getGlobalVariable(int index) {
    return OpenMMNative.OpenMM_CustomIntegrator_getGlobalVariable(getPointer(), index);
  }

  /**
   * Get the current value of a global variable by name.
   *
   * @param name variable name.
   * @return current value.
   */
  public double getGlobalVariableByName(String name) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomIntegrator_getGlobalVariableByName(
            getPointer(), value));
  }

  /**
   * Get the name of a global variable.
   *
   * @param index index of the variable.
   * @return variable name, copied into a Java string.
   */
  public String getGlobalVariableName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomIntegrator_getGlobalVariableName(
        getPointer(), index));
  }

  /**
   * Get the expression used to compute the kinetic energy. It is evaluated for every degree of freedom (excluding
   * mass-0 particles) and the values are summed. The default is {@code m*v*v/2}.
   *
   * @return kinetic-energy expression, copied into a Java string.
   */
  public String getKineticEnergyExpression() {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomIntegrator_getKineticEnergyExpression(getPointer()));
  }

  /**
   * Get the number of computation steps that have been added, including constraint, context-update and block steps.
   *
   * @return number of computation steps.
   */
  public int getNumComputations() {
    return OpenMMNative.OpenMM_CustomIntegrator_getNumComputations(getPointer());
  }

  /**
   * Get the number of global variables that have been defined.
   *
   * @return number of global variables.
   */
  public int getNumGlobalVariables() {
    return OpenMMNative.OpenMM_CustomIntegrator_getNumGlobalVariables(getPointer());
  }

  /**
   * Get the number of per-degree-of-freedom variables that have been defined.
   *
   * @return number of per-DOF variables.
   */
  public int getNumPerDofVariables() {
    return OpenMMNative.OpenMM_CustomIntegrator_getNumPerDofVariables(getPointer());
  }

  /**
   * Get the number of tabulated functions that have been defined.
   *
   * @return number of tabulated functions.
   */
  public int getNumTabulatedFunctions() {
    return OpenMMNative.OpenMM_CustomIntegrator_getNumTabulatedFunctions(getPointer());
  }

  /**
   * Get the values of a per-degree-of-freedom variable.
   *
   * <p>The values are copied into a new native vector array, one {@link Vec3} per particle whose components are the
   * three degrees of freedom. The caller owns the returned array and must close or destroy it. The JNA method fills a
   * {@code PointerByReference} instead.</p>
   *
   * @param index index of the variable.
   * @return new caller-owned array of the variable values.
   */
  public Vec3Array getPerDofVariable(int index) {
    Vec3Array values = new Vec3Array(0);
    OpenMMNative.OpenMM_CustomIntegrator_getPerDofVariable(
        getPointer(), index, values.getPointer());
    return values;
  }

  /**
   * Get the values of a per-degree-of-freedom variable by name.
   *
   * <p>See {@link #getPerDofVariable(int)} for the layout and ownership of the returned array.</p>
   *
   * @param name variable name.
   * @return new caller-owned array of the variable values.
   */
  public Vec3Array getPerDofVariableByName(String name) {
    Vec3Array values = new Vec3Array(0);
    OpenMMStrings.withUtf8String(name, value -> {
      OpenMMNative.OpenMM_CustomIntegrator_getPerDofVariableByName(
          getPointer(), value, values.getPointer());
    });
    return values;
  }

  /**
   * Get the name of a per-degree-of-freedom variable.
   *
   * @param index index of the variable.
   * @return variable name, copied into a Java string.
   */
  public String getPerDofVariableName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomIntegrator_getPerDofVariableName(getPointer(), index));
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)}.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen at context creation.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_CustomIntegrator_getRandomNumberSeed(getPointer());
  }

  /**
   * Get a tabulated function.
   *
   * <p>The returned segment is a borrowed, non-owning native handle valid only while this integrator is alive; it is not
   * wrapped in a Java class and must not be destroyed by the caller. The JNA method returns a {@code PointerByReference}.</p>
   *
   * @param index index of the function.
   * @return borrowed native function handle.
   */
  public MemorySegment getTabulatedFunction(int index) {
    return OpenMMNative.OpenMM_CustomIntegrator_getTabulatedFunction(getPointer(), index);
  }

  /**
   * Get the name of a tabulated function as it appears in expressions.
   *
   * @param index index of the function.
   * @return function name, copied into a Java string.
   */
  public String getTabulatedFunctionName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomIntegrator_getTabulatedFunctionName(getPointer(), index));
  }

  /**
   * Set the value of a global variable.
   *
   * @param index index of the variable.
   * @param value new value.
   */
  public void setGlobalVariable(int index, double value) {
    OpenMMNative.OpenMM_CustomIntegrator_setGlobalVariable(getPointer(), index, value);
  }

  /**
   * Set the value of a global variable by name.
   *
   * @param name  variable name.
   * @param value new value.
   */
  public void setGlobalVariableByName(String name, double value) {
    OpenMMStrings.withUtf8String(name,
        nativeName -> OpenMMNative.OpenMM_CustomIntegrator_setGlobalVariableByName(
            getPointer(), nativeName, value));
  }

  /**
   * Set the expression used to compute the kinetic energy. It is evaluated for every degree of freedom (excluding
   * mass-0 particles) and the values are summed. It may depend only on {@code x}, {@code v}, {@code f}, {@code m} and {@code dt},
   * user-defined variables and context parameters, not on potential energy, single-group forces or random numbers.
   *
   * @param expression kinetic-energy expression.
   */
  public void setKineticEnergyExpression(String expression) {
    OpenMMStrings.withUtf8String(expression,
        nativeExpression -> OpenMMNative.OpenMM_CustomIntegrator_setKineticEnergyExpression(
            getPointer(), nativeExpression));
  }

  /**
   * Set the values of a per-degree-of-freedom variable.
   *
   * <p>The values are copied by OpenMM; the caller keeps ownership of {@code values}. The JNA method takes a {@code
   * PointerByReference}.</p>
   *
   * @param index  index of the variable.
   * @param values new values, one {@link Vec3} per particle (one entry per degree-of-freedom triple); must not be null.
   */
  public void setPerDofVariable(int index, Vec3Array values) {
    OpenMMNative.OpenMM_CustomIntegrator_setPerDofVariable(
        getPointer(), index, values.getPointer());
  }

  /**
   * Set the values of a per-degree-of-freedom variable by name.
   *
   * <p>See {@link #setPerDofVariable(int, Vec3Array)}; the caller keeps ownership of {@code values}.</p>
   *
   * @param name   variable name.
   * @param values new values; must not be null.
   */
  public void setPerDofVariableByName(String name, Vec3Array values) {
    OpenMMStrings.withUtf8String(name,
        nativeName -> OpenMMNative.OpenMM_CustomIntegrator_setPerDofVariableByName(
            getPointer(), nativeName, values.getPointer()));
  }

  /**
   * Set the random-number seed. The meaning is left to each OpenMM platform. Different seeds are guaranteed to give
   * different random-number sequences, but equal seeds carry no guarantee because platforms may use non-deterministic
   * algorithms. A seed of 0 (the default) causes a unique seed to be chosen when a context is created.
   *
   * @param seed random-number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_CustomIntegrator_setRandomNumberSeed(getPointer(), seed);
  }

  /**
   * Advance the simulation by executing the defined algorithm for each time step. This overrides {@link
   * Integrator#step(int)} with the CustomIntegrator native step call.
   *
   * @param steps number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_CustomIntegrator_step(getPointer(), steps);
  }

  /**
   * Create the native integrator after loading the FFM runtime.
   *
   * @param stepSize step size, in ps.
   * @return native integrator handle owned by the new wrapper.
   */
  private static MemorySegment create(double stepSize) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_CustomIntegrator_create(stepSize);
  }

  private int withPair(String variable, String expression, PairCall call) {
    return OpenMMStrings.withUtf8StringResult(variable, nativeVariable ->
        OpenMMStrings.withUtf8StringResult(expression, nativeExpression ->
            call.invoke(getPointer(), nativeVariable, nativeExpression)));
  }

  private int withString(String value, StringCall call) {
    return OpenMMStrings.withUtf8StringResult(value,
        nativeValue -> call.invoke(getPointer(), nativeValue));
  }

  @FunctionalInterface
  private interface PairCall {
    int invoke(MemorySegment target, MemorySegment variable, MemorySegment expression);
  }

  @FunctionalInterface
  private interface StringCall {
    int invoke(MemorySegment target, MemorySegment value);
  }
}
