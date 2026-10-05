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
 * Custom algebraic interaction for bonds containing a fixed number of particles.
 *
 * <p>The constructor's {@code numParticles} fixes the number and ordering of particle indices in
 * every bond. The expression may use OpenMM geometric functions of the ordered positions (for
 * example distances, angles, and dihedrals), global/per-bond parameters, and custom functions.
 * Distances are in nm, angles in radians, and energy in kJ/mol; arbitrary parameter values are
 * not unit-converted. Per-bond vectors follow declaration order.
 *
 * <p>Java arrays and record snapshots are copied. Native array inputs are borrowed for the call
 * only. Adding a {@link TabulatedFunction} transfers native ownership and invalidates that
 * wrapper; raw native-handle overloads transfer ownership too, but cannot invalidate any
 * associated Java wrapper. Returned function handles are borrowed from the force. Context
 * updates propagate only per-bond and tabulated-function values and are a no-op for a context
 * with no native pointer; expressions, declarations, topology, and global defaults require context
 * recreation.
 *
 * <p>The native header deprecates legacy continuous-function APIs in favor of tabulated
 * functions; these Java compatibility overloads remain. The legacy function-parameter getter is
 * unsupported: its C wrapper declares a {@code char**} output while writing a C++
 * {@code std::string}, an incompatible ABI. Do not treat that output as usable.
 */
public class CustomCompoundBondForce extends Force {

  /**
   * Snapshot of a bond's ordered particle indices and per-bond values.
   *
   * @param particles copied particle indices in expression order
   * @param parameters copied values in per-bond declaration order
   */
  public record BondParameters(int[] particles, double[] parameters) {}
  /**
   * Legacy continuous-function definition.
   *
   * @param name function name referenced by an expression
   * @param values copied tabulated values
   * @param min lower bound of the function's input domain
   * @param max upper bound of the function's input domain
   */
  public record FunctionParameters(String name, double[] values, double min, double max) {}

  /**
   * Create a custom force for bonds with a fixed particle count.
   *
   * @param numParticles number and ordering length of particles in each bond
   * @param energy custom energy expression
   */
  public CustomCompoundBondForce(int numParticles, String energy) {
    super(create(numParticles, energy));
  }

  /**
   * Add a bond using caller-owned native arrays.
   *
   * @param particles native array containing exactly {@code numParticles} ordered particle indices
   * @param parameters native array containing one value per declared per-bond parameter
   * @return index assigned to the bond
   */
  public int addBond(IntArray particles, DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_addBond(
        getPointer(), particles.getPointer(), parameters.getPointer());
  }

  /**
   * Add a bond from Java arrays; both arrays are copied to temporary native storage.
   *
   * @param particles ordered particle indices; length must equal the constructor count
   * @param parameters per-bond values in declaration order
   * @return index assigned to the bond
   */
  public int addBond(int[] particles, double[] parameters) {
    try (IntArray nativeParticles = toNative(particles);
         DoubleArray nativeParameters = toNative(parameters)) {
      return addBond(nativeParticles, nativeParameters);
    }
  }

  /** Request energy differentiation with respect to an already declared global parameter.
   * @param name exact global parameter name */
  public void addEnergyParameterDerivative(String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addEnergyParameterDerivative(
            getPointer(), value));
  }

  /** Add a legacy continuous function; prefer {@link #addTabulatedFunction(String, TabulatedFunction)}.
   * The supplied native array is borrowed synchronously.
   * @param name function name
   * @param values native function values
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  @Deprecated
  public int addFunction(String name, DoubleArray values, double min, double max) {
    return withStringResult(name, nativeName ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addFunction(
            getPointer(), nativeName, values.getPointer(), min, max));
  }

  /** Add a legacy continuous function from a caller-owned native value-array segment.
   * @param name function name
   * @param values caller-owned native value array, borrowed for the call
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  @Deprecated
  public int addFunction(String name, MemorySegment values, double min, double max) {
    return withStringResult(name, nativeName ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addFunction(
            getPointer(), nativeName, values, min, max));
  }

  /** Add a legacy continuous function using caller-owned NUL-terminated UTF-8 name and value-array
   * segments borrowed during the call.
   * @param name native function-name segment
   * @param values native value-array segment
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  @Deprecated
  public int addFunction(
      MemorySegment name, MemorySegment values, double min, double max) {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_addFunction(
        getPointer(), name, values, min, max);
  }

  /** Declare an expression-wide parameter with a default for newly created contexts.
   * @param name expression parameter name
   * @param defaultValue initial value, using expression-consistent units
   * @return declaration index */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /** Declare a per-bond parameter; every bond supplies one value in declaration order.
   * @param name parameter name
   * @return declaration index */
  public int addPerBondParameter(String name) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addPerBondParameter(
            getPointer(), value));
  }

  /** Add a tabulated function and transfer its native ownership to this force.
   * @param name expression-visible function name
   * @param function wrapper invalidated after ownership transfer
   * @return function index */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addTabulatedFunction(
            getPointer(), value, function.getPointer()));
    function.invalidate();
    return index;
  }

  /** Add a function from a Java name and native handle; native ownership transfers to this force,
   * but this overload does not invalidate an associated Java wrapper.
   * @param name function name
   * @param function native function handle */
  public int addTabulatedFunction(String name, MemorySegment function) {
    return withStringResult(name, nativeName ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_addTabulatedFunction(
            getPointer(), nativeName, function));
  }

  /** Add a function from caller-owned handles; native function ownership transfers to this force,
   * but an associated Java wrapper is not invalidated automatically.
   * @param name NUL-terminated UTF-8 name segment, borrowed synchronously
   * @param function native function handle */
  public int addTabulatedFunction(MemorySegment name, MemorySegment function) {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_addTabulatedFunction(
        getPointer(), name, function);
  }

  /** Destroy this force. */
  /** Destroy the native force and its owned child resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomCompoundBondForce_destroy);
  }

  /** @param index bond index
   *  @return independent copies of ordered particle indices and per-bond values */
  public BondParameters getBondParameters(int index) {
    try (IntArray particles = new IntArray(0); DoubleArray parameters = new DoubleArray(0)) {
      OpenMMNative.OpenMM_CustomCompoundBondForce_getBondParameters(
          getPointer(), index, particles.getPointer(), parameters.getPointer());
      return new BondParameters(copy(particles), copy(parameters));
    }
  }

  /** Copy bond values to caller-owned native arrays with sufficient capacities.
   * @param index bond index
   * @param particles output array for ordered particle indices
   * @param parameters output array for per-bond values */
  public void getBondParameters(int index, IntArray particles, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_getBondParameters(
        getPointer(), index, particles.getPointer(), parameters.getPointer());
  }

  /** @return copied current energy expression */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCompoundBondForce_getEnergyFunction(getPointer()));
  }

  /** @param index derivative index
   *  @return copied differentiated global parameter name */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCompoundBondForce_getEnergyParameterDerivativeName(
            getPointer(), index));
  }

  /**
   * This output is unsupported: the C wrapper uses {@code char**}, but its implementation writes
   * a C++ {@code std::string}; that ABI cannot safely return the function name.
   *
   * @param index legacy function index
   * @return never returns
   * @throws UnsupportedOperationException always, because of the incompatible output ABI
   */
  public FunctionParameters getFunctionParameters(int index) {
    throw new UnsupportedOperationException(
        "OpenMM_CustomCompoundBondForce_getFunctionParameters has an incompatible char** output ABI.");
  }

  /** @param index global parameter index
   *  @return default value used for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getGlobalParameterDefaultValue(
        getPointer(), index);
  }

  /** @param index global parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCompoundBondForce_getGlobalParameterName(getPointer(), index));
  }

  /** @return number of stored compound bonds */
  public int getNumBonds() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumBonds(getPointer());
  }

  /** @return number of requested global-parameter derivatives */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumEnergyParameterDerivatives(
        getPointer());
  }

  /** @return legacy function count (deprecated alias for tabulated functions) */
  @Deprecated
  public int getNumFunctions() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumFunctions(getPointer());
  }

  /** @return number of global parameter declarations */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumGlobalParameters(getPointer());
  }

  /** @return fixed number of ordered particles in each bond */
  public int getNumParticlesPerBond() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumParticlesPerBond(getPointer());
  }

  /** @return number of per-bond declarations and values required in each bond */
  public int getNumPerBondParameters() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumPerBondParameters(getPointer());
  }

  /** @return number of registered tabulated functions */
  public int getNumTabulatedFunctions() {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getNumTabulatedFunctions(getPointer());
  }

  /** @param index per-bond parameter index
   *  @return copied parameter name */
  public String getPerBondParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCompoundBondForce_getPerBondParameterName(
            getPointer(), index));
  }

  /** @param index function index
   *  @return borrowed force-owned function handle; do not destroy it */
  public MemorySegment getTabulatedFunction(int index) {
    return OpenMMNative.OpenMM_CustomCompoundBondForce_getTabulatedFunction(getPointer(), index);
  }

  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCompoundBondForce_getTabulatedFunctionName(
            getPointer(), index));
  }

  /** Replace one bond's ordered particle indices and all per-bond values; native arrays are
   * borrowed synchronously.
   * @param index stored bond index
   * @param particles replacement ordered indices
   * @param parameters replacement values in declaration order */
  public void setBondParameters(int index, IntArray particles, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_setBondParameters(
        getPointer(), index, particles.getPointer(), parameters.getPointer());
  }

  /** Replace bond values from caller-owned segments borrowed during the call.
   * @param index stored bond index
   * @param particles native ordered-index array
   * @param parameters native per-bond value array */
  public void setBondParameters(
      int index, MemorySegment particles, MemorySegment parameters) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_setBondParameters(
        getPointer(), index, particles, parameters);
  }

  /** Replace the expression; recreate contexts to use the new expression.
   * @param energy OpenMM custom energy expression */
  public void setEnergyFunction(String energy) {
    withString(energy, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_setEnergyFunction(getPointer(), value));
  }

  /** Replace a legacy continuous-function definition; prefer tabulated-function APIs.
   * @param index legacy function index
   * @param name function name
   * @param values native values, borrowed during the call
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  @Deprecated
  public void setFunctionParameters(
      int index, String name, DoubleArray values, double min, double max) {
    withString(name, nativeName ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_setFunctionParameters(
            getPointer(), index, nativeName, values.getPointer(), min, max));
  }

  /** Replace a legacy function from a caller-owned native value-array segment borrowed for this call.
   * @param index legacy function index
   * @param name function name
   * @param values native value-array segment
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  @Deprecated
  public void setFunctionParameters(
      int index, String name, MemorySegment values, double min, double max) {
    withString(name, nativeName ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_setFunctionParameters(
            getPointer(), index, nativeName, values, min, max));
  }

  /** Replace a legacy function using caller-owned NUL-terminated name and value-array segments.
   * @param index legacy function index
   * @param name native name segment
   * @param values native values segment
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  @Deprecated
  public void setFunctionParameters(
      int index, MemorySegment name, MemorySegment values, double min, double max) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_setFunctionParameters(
        getPointer(), index, name, values, min, max);
  }

  /** Set the default for new contexts, not the current value in live contexts.
   * @param index global parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_setGlobalParameterDefaultValue(
        getPointer(), index, value);
  }

  /** Rename a global parameter; existing contexts are not reconfigured.
   * @param index global parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /** Rename a per-bond parameter; existing contexts are not reconfigured.
   * @param index per-bond parameter index
   * @param name new expression name */
  public void setPerBondParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCompoundBondForce_setPerBondParameterName(
            getPointer(), index, value));
  }

  /** Enable periodic minimum-image displacements for the bond geometry.
   * @param periodic {@code true} enables periodic displacements */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /** Set the native integer boolean flag (zero disables; nonzero enables periodic displacements).
   * @param periodic native boolean integer */
  public void setUsesPeriodicBoundaryConditions(int periodic) {
    OpenMMNative.OpenMM_CustomCompoundBondForce_setUsesPeriodicBoundaryConditions(
        getPointer(), periodic);
  }

  /** Apply only per-bond parameter values and tabulated values in an existing context. Particle
   * membership/order cannot change and bonds cannot be added; tabulated-function dimensions and
   * domain/range must remain unchanged.
   * @param context context containing this force; no-op if no native context pointer is present.
   *     Expression, topology, declarations, and global defaults are not propagated. */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_CustomCompoundBondForce_updateParametersInContext(
                getPointer(), pointer));
  }

  /** @return {@code true} if periodic boundary conditions are enabled for bond geometry. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomCompoundBondForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(int numParticles, String energy) {
    return withStringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomCompoundBondForce_create(numParticles, value);
    });
  }

  private static IntArray toNative(int[] values) {
    IntArray result = new IntArray(values.length);
    for (int index = 0; index < values.length; index++) {
      result.set(index, values[index]);
    }
    return result;
  }

  private static DoubleArray toNative(double[] values) {
    DoubleArray result = new DoubleArray(values.length);
    for (int index = 0; index < values.length; index++) {
      result.set(index, values[index]);
    }
    return result;
  }

  private static int[] copy(IntArray values) {
    int[] result = new int[values.getSize()];
    for (int index = 0; index < result.length; index++) {
      result[index] = values.get(index);
    }
    return result;
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
