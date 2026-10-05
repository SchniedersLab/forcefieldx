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
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Custom energy function of collective variables produced by child forces.
 *
 * <p>Each variable name is available in the expression as a scalar whose value is computed by
 * its child {@link Force}. Energy expressions use OpenMM's custom-expression language; energy is
 * in kJ/mol, while the child force determines each variable's physical units. Global parameters
 * and tabulated functions may also be referenced.
 *
 * <p>Adding a wrapped child force or {@link TabulatedFunction} transfers native ownership to this
 * force and invalidates the supplied wrapper. Raw {@link MemorySegment} overloads also transfer
 * native ownership of child/function handles, but cannot invalidate associated Java wrappers;
 * names are borrowed only for the call. Do not destroy transferred handles independently.
 * Returned child, inner-context, and function
 * handles are borrowed and must not be destroyed independently. Context updates do not recompile
 * expressions or replace child forces/declarations; global defaults only affect new contexts.
 * Java strings and scalar arrays returned by getters are copied.
 */
public class CustomCVForce extends Force {

  /** Create a CV force from a custom expression over names of subsequently added variables.
   * @param energy expression evaluated in terms of collective-variable names */
  public CustomCVForce(String energy) {
    super(create(energy));
  }

  /**
   * Add a named child force as a collective variable and transfer its ownership.
   *
   * @param name variable name referenced by the energy expression
   * @param force child force whose native ownership transfers; this Java wrapper is invalidated
   * @return variable index
   */
  public int addCollectiveVariable(String name, Force force) {
    int index = withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCVForce_addCollectiveVariable(
            getPointer(), value, force.getPointer()));
    force.invalidate();
    return index;
  }

  /**
   * Add a variable from a caller-owned name segment and native child-force handle. Native force
   * ownership transfers to this force, but the raw overload cannot invalidate an associated Java
   * wrapper.
   *
   * @param name NUL-terminated UTF-8 variable-name segment
   * @param force native child-force handle
   * @return variable index
   */
  public int addCollectiveVariable(MemorySegment name, MemorySegment force) {
    return OpenMMNative.OpenMM_CustomCVForce_addCollectiveVariable(
        getPointer(), name, force);
  }

  /** Request an energy derivative for an already declared global parameter.
   * @param name exact global parameter name */
  public void addEnergyParameterDerivative(String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCVForce_addEnergyParameterDerivative(getPointer(), value));
  }

  /** Request a derivative using a caller-owned NUL-terminated UTF-8 name segment.
   * @param name native name segment borrowed for this call */
  public void addEnergyParameterDerivative(MemorySegment name) {
    OpenMMNative.OpenMM_CustomCVForce_addEnergyParameterDerivative(getPointer(), name);
  }

  /** Declare a global parameter and default for newly created contexts.
   * @param name expression parameter name
   * @param defaultValue initial value in expression-consistent units
   * @return declaration index */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCVForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /** Declare a global parameter from a caller-owned NUL-terminated UTF-8 name.
   * @param name native name segment borrowed during this call
   * @param defaultValue initial value
   * @return declaration index */
  public int addGlobalParameter(MemorySegment name, double defaultValue) {
    return OpenMMNative.OpenMM_CustomCVForce_addGlobalParameter(
        getPointer(), name, defaultValue);
  }

  /** Add a tabulated function and transfer its native ownership to this force.
   * @param name expression-visible function name
   * @param function wrapper invalidated after ownership transfer
   * @return function index */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCVForce_addTabulatedFunction(
            getPointer(), value, function.getPointer()));
    function.invalidate();
    return index;
  }

  /** Add a function from a caller-owned name segment and wrapped function. Native ownership
   * transfers and the wrapper is invalidated.
   * @param name NUL-terminated UTF-8 segment borrowed during the call
   * @param function function wrapper invalidated after transfer
   * @return function index */
  public int addTabulatedFunction(MemorySegment name, TabulatedFunction function) {
    int index = OpenMMNative.OpenMM_CustomCVForce_addTabulatedFunction(
        getPointer(), name, function.getPointer());
    function.invalidate();
    return index;
  }

  /** Add a function using a caller-owned name segment and native function handle. Native ownership
   * transfers to this force, but an associated Java wrapper is not invalidated automatically.
   * @param name NUL-terminated UTF-8 name segment, borrowed during the call
   * @param function native function handle
   * @return function index */
  public int addTabulatedFunction(MemorySegment name, MemorySegment function) {
    return OpenMMNative.OpenMM_CustomCVForce_addTabulatedFunction(
        getPointer(), name, function);
  }

  /** Destroy this force. */
  /** Destroy the native force, including forces/functions whose ownership was transferred to it. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomCVForce_destroy);
  }

  /** @param index collective-variable index
   *  @return borrowed child-force handle owned by this force; do not destroy it */
  public MemorySegment getCollectiveVariable(int index) {
    return OpenMMNative.OpenMM_CustomCVForce_getCollectiveVariable(getPointer(), index);
  }

  /** @param index collective-variable index
   *  @return copied variable name */
  public String getCollectiveVariableName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCVForce_getCollectiveVariableName(getPointer(), index));
  }

  /** Compute all collective-variable values in the supplied context and return a Java copy.
   * @param context initialized context associated with this force
   * @return values in collective-variable declaration order */
  public double[] getCollectiveVariableValues(Context context) {
    try (DoubleArray values = new DoubleArray(0)) {
      OpenMMNative.OpenMM_CustomCVForce_getCollectiveVariableValues(
          getPointer(), context.getPointer(), values.getPointer());
      return copy(values);
    }
  }

  /** Copy collective-variable values into caller-owned native storage.
   * @param context initialized context associated with this force
   * @param values caller-owned output array with capacity for every variable */
  public void getCollectiveVariableValues(Context context, DoubleArray values) {
    OpenMMNative.OpenMM_CustomCVForce_getCollectiveVariableValues(
        getPointer(), context.getPointer(), values.getPointer());
  }

  /** Copy values into a caller-owned native segment with capacity for every variable.
   * @param context initialized context associated with this force
   * @param values output array segment borrowed during the call */
  public void getCollectiveVariableValues(Context context, MemorySegment values) {
    OpenMMNative.OpenMM_CustomCVForce_getCollectiveVariableValues(
        getPointer(), context.getPointer(), values);
  }

  /** @return copied current energy expression */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomCVForce_getEnergyFunction(getPointer()));
  }

  /** @param index derivative index
   *  @return copied differentiated global parameter name */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCVForce_getEnergyParameterDerivativeName(
            getPointer(), index));
  }

  /** @param index global parameter index
   *  @return default used for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomCVForce_getGlobalParameterDefaultValue(
        getPointer(), index);
  }

  /** @param index global parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCVForce_getGlobalParameterName(getPointer(), index));
  }

  /** Get the borrowed inner context associated with an initialized outer context. Do not destroy
   * the returned handle.
   * @param context outer context containing this force
   * @return borrowed inner context handle */
  public MemorySegment getInnerContext(Context context) {
    return OpenMMNative.OpenMM_CustomCVForce_getInnerContext(
        getPointer(), context.getPointer());
  }

  /** @return number of child forces/collective variables */
  public int getNumCollectiveVariables() {
    return OpenMMNative.OpenMM_CustomCVForce_getNumCollectiveVariables(getPointer());
  }

  /** @return number of requested global-parameter derivatives */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_CustomCVForce_getNumEnergyParameterDerivatives(getPointer());
  }

  /** @return number of declared global parameters */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomCVForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of registered tabulated functions */
  public int getNumTabulatedFunctions() {
    return OpenMMNative.OpenMM_CustomCVForce_getNumTabulatedFunctions(getPointer());
  }

  /** @param index function index
   *  @return borrowed force-owned function handle; do not destroy it */
  public MemorySegment getTabulatedFunction(int index) {
    return OpenMMNative.OpenMM_CustomCVForce_getTabulatedFunction(getPointer(), index);
  }

  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCVForce_getTabulatedFunctionName(getPointer(), index));
  }

  /** Replace the energy expression; existing contexts must be recreated to use it.
   * @param energy custom expression over declared CV names */
  public void setEnergyFunction(String energy) {
    withString(energy, value ->
        OpenMMNative.OpenMM_CustomCVForce_setEnergyFunction(getPointer(), value));
  }

  /** Replace the expression from a live caller-owned NUL-terminated UTF-8 segment.
   * @param energy native expression segment borrowed for this call */
  public void setEnergyFunction(MemorySegment energy) {
    OpenMMNative.OpenMM_CustomCVForce_setEnergyFunction(getPointer(), energy);
  }

  /** Change the default for future contexts, not the value in existing contexts.
   * @param index global parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_CustomCVForce_setGlobalParameterDefaultValue(
        getPointer(), index, value);
  }

  /** Rename a declared global parameter; existing contexts are not reconfigured.
   * @param index global parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCVForce_setGlobalParameterName(getPointer(), index, value));
  }

  /** Rename a global parameter from a caller-owned NUL-terminated UTF-8 segment.
   * @param index global parameter index
   * @param name native name segment borrowed for this call */
  public void setGlobalParameterName(int index, MemorySegment name) {
    OpenMMNative.OpenMM_CustomCVForce_setGlobalParameterName(getPointer(), index, name);
  }

  /** Apply only tabulated-function value changes to a context; function dimensions and
   * domain/range must remain unchanged.
   * @param context context containing this force; no-op without a native handle. This cannot
   *     replace child forces, expressions, declarations, or global defaults. */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_CustomCVForce_updateParametersInContext(
                getPointer(), pointer));
  }

  /** @return whether this force reports use of periodic boundary conditions. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomCVForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(String energy) {
    return withStringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomCVForce_create(value);
    });
  }

  private static double[] copy(DoubleArray values) {
    double[] result = new double[values.getSize()];
    for (int index = 0; index < result.length; index++) {
      result[index] = values.get(index);
    }
    return result;
  }

  private static void withString(String value, Consumer<MemorySegment> action) {
    OpenMMStrings.withUtf8String(value, action);
  }

  private static <T> T withStringResult(String value, Function<MemorySegment, T> action) {
    return OpenMMStrings.withUtf8StringResult(value, action);
  }
}
