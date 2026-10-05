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
 * Custom interactions between particle groups, using each group's weighted centroid as a site.
 *
 * <p>The constructor's positive {@code numGroups} fixes the number of group indices in every
 * bond. The energy expression can use OpenMM custom geometric variables such as distances,
 * angles, and dihedrals between group centroids, global and per-bond parameters, and custom
 * functions. Geometric distances are in nm, angles are in radians, and energy is in kJ/mol.
 * Parameter vectors follow declaration order; each group has particle indices and corresponding
 * weights, which define its centroid.
 *
 * <p>The header's {@code getNumFunctions()} compatibility alias is deprecated in favor of
 * {@link #getNumTabulatedFunctions()}.
 *
 * <p>Java primitive arrays and the record-returning getters are copied. Native-array inputs are
 * borrowed only for the call. Adding a tabulated function transfers native ownership to this
 * force; the typed overload invalidates its {@link TabulatedFunction} wrapper, whereas raw-handle
 * overloads cannot invalidate an associated wrapper. Do not destroy a transferred function
 * independently. Getter handles are borrowed and must not be destroyed by the caller.
 * Updating a context propagates only per-bond parameter and tabulated-function values, not
 * expressions, declarations,
 * structure, or global defaults. A context without a native handle is ignored.
 */
public class CustomCentroidBondForce extends Force {

  /**
   * Snapshot of a bond's group indices and per-bond values.
   *
   * @param groups copied group indices in bond order
   * @param parameters copied values in per-bond declaration order
   */
  public record BondParameters(int[] groups, double[] parameters) {}
  /**
   * Snapshot of a group's centroid definition.
   *
   * @param particles copied particle indices
   * @param weights copied corresponding weights; element {@code i} weights particle {@code i}
   */
  public record GroupParameters(int[] particles, double[] weights) {}

  /**
   * Create a centroid-bond force.
   *
   * @param numGroups fixed number of group centroids participating in each bond
   * @param energy OpenMM custom expression evaluated for each bond
   */
  public CustomCentroidBondForce(int numGroups, String energy) {
    super(create(numGroups, energy));
  }

  /**
   * Add a bond using caller-owned native arrays.
   *
   * @param groups native integer array of exactly {@code numGroups} group indices
   * @param parameters native double array with one value per declared per-bond parameter
   * @return index assigned to the bond
   */
  public int addBond(IntArray groups, DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_addBond(
        getPointer(), groups.getPointer(), parameters.getPointer());
  }

  /**
   * Add a bond from Java arrays, which are copied to temporary native storage.
   *
   * @param groups group indices in bond order; length must equal the constructor's group count
   * @param parameters per-bond values in declaration order
   * @return index assigned to the bond
   */
  public int addBond(int[] groups, double[] parameters) {
    try (IntArray nativeGroups = toNative(groups);
         DoubleArray nativeParameters = toNative(parameters)) {
      return addBond(nativeGroups, nativeParameters);
    }
  }

  /** Request an energy derivative for a previously declared global parameter.
   * @param name exact global parameter name */
  public void addEnergyParameterDerivative(String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_addEnergyParameterDerivative(
            getPointer(), value));
  }

  /** Declare an expression-wide parameter and the default used by newly created contexts.
   * @param name expression parameter name
   * @param defaultValue initial value, in expression-consistent units
   * @return declaration index */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /**
   * Add a centroid group from native arrays.
   *
   * @param particles native particle-index array
   * @param weights native array of corresponding weights; use an empty array for default weights
   * @return index assigned to the group
   */
  public int addGroup(IntArray particles, DoubleArray weights) {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_addGroup(
        getPointer(), particles.getPointer(), weights.getPointer());
  }

  /**
   * Add a centroid group; Java arrays are copied and must have matching lengths when weights are
   * explicitly supplied.
   *
   * @param particles particle indices defining the group
   * @param weights corresponding weights, or an empty array for default weights
   * @return index assigned to the group
   */
  public int addGroup(int[] particles, double[] weights) {
    try (IntArray nativeParticles = toNative(particles);
         DoubleArray nativeWeights = toNative(weights)) {
      return addGroup(nativeParticles, nativeWeights);
    }
  }

  /** Declare a per-bond parameter; each bond supplies one value in declaration order.
   * @param name parameter name
   * @return declaration index */
  public int addPerBondParameter(String name) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_addPerBondParameter(
            getPointer(), value));
  }

  /**
   * Add a tabulated function referenced by the energy expression.
   *
   * @param name expression-visible function name
   * @param function function wrapper whose native ownership transfers; wrapper is invalidated
   * @return function index
   */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_addTabulatedFunction(
            getPointer(), value, function.getPointer()));
    function.invalidate();
    return index;
  }

  /** Add a function from a native handle; native ownership transfers to this force, but any
   * associated Java wrapper is not invalidated automatically.
   * @param name function name, encoded from a temporary Java string
   * @param function native function handle */
  public int addTabulatedFunction(String name, MemorySegment function) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_addTabulatedFunction(
            getPointer(), value, function));
  }

  /** Add a function from caller-owned handles. The function's native ownership transfers to this
   * force, but an associated Java wrapper is not invalidated automatically.
   * @param name NUL-terminated UTF-8 function-name segment, borrowed during the call
   * @param function native function handle */
  public int addTabulatedFunction(MemorySegment name, MemorySegment function) {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_addTabulatedFunction(
        getPointer(), name, function);
  }

  /** Destroy this force. */
  /** Destroy the native force and its owned child resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomCentroidBondForce_destroy);
  }

  /** @param index bond index
   *  @return copied group indices and per-bond values in declaration order */
  public BondParameters getBondParameters(int index) {
    try (IntArray groups = new IntArray(0); DoubleArray parameters = new DoubleArray(0)) {
      OpenMMNative.OpenMM_CustomCentroidBondForce_getBondParameters(
          getPointer(), index, groups.getPointer(), parameters.getPointer());
      return new BondParameters(copy(groups), copy(parameters));
    }
  }

  /**
   * Copy a bond's group indices and values to native arrays.
   *
   * @param index bond index
   * @param groups caller-owned native integer array; caller must provide sufficient capacity
   * @param parameters caller-owned native double array; caller must provide sufficient capacity
   */
  public void getBondParameters(int index, IntArray groups, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_getBondParameters(
        getPointer(), index, groups.getPointer(), parameters.getPointer());
  }

  /** @return copied energy-expression string */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCentroidBondForce_getEnergyFunction(getPointer()));
  }

  /** @param index requested derivative index
   *  @return copied global parameter name */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCentroidBondForce_getEnergyParameterDerivativeName(
            getPointer(), index));
  }

  /** @param index global parameter index
   *  @return default used for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getGlobalParameterDefaultValue(
        getPointer(), index);
  }

  /** @param index global parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCentroidBondForce_getGlobalParameterName(getPointer(), index));
  }

  /** @param index group index
   *  @return copied particle indices and corresponding centroid weights */
  public GroupParameters getGroupParameters(int index) {
    try (IntArray particles = new IntArray(0); DoubleArray weights = new DoubleArray(0)) {
      OpenMMNative.OpenMM_CustomCentroidBondForce_getGroupParameters(
          getPointer(), index, particles.getPointer(), weights.getPointer());
      return new GroupParameters(copy(particles), copy(weights));
    }
  }

  /** Copy a group's particle indices and weights to caller-owned native arrays.
   * @param index group index
   * @param particles caller-owned integer array with sufficient capacity
   * @param weights caller-owned double array with sufficient capacity */
  public void getGroupParameters(int index, IntArray particles, DoubleArray weights) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_getGroupParameters(
        getPointer(), index, particles.getPointer(), weights.getPointer());
  }

  /** @return number of stored bonds */
  public int getNumBonds() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumBonds(getPointer());
  }

  /** @return number of requested global-parameter derivatives */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumEnergyParameterDerivatives(
        getPointer());
  }

  /** @return number of registered functions (deprecated alias for tabulated functions) */
  @Deprecated
  public int getNumFunctions() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumFunctions(getPointer());
  }

  /** @return number of global parameter declarations */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of registered centroid groups */
  public int getNumGroups() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumGroups(getPointer());
  }

  /** @return fixed group count per bond, supplied to the constructor */
  public int getNumGroupsPerBond() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumGroupsPerBond(getPointer());
  }

  /** @return number of per-bond declarations and required values per bond */
  public int getNumPerBondParameters() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumPerBondParameters(getPointer());
  }

  /** @return number of registered tabulated functions */
  public int getNumTabulatedFunctions() {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getNumTabulatedFunctions(getPointer());
  }

  /** @param index per-bond parameter index
   *  @return copied parameter name */
  public String getPerBondParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCentroidBondForce_getPerBondParameterName(getPointer(), index));
  }

  /** @param index function index
   *  @return borrowed handle owned by this force; do not destroy it */
  public MemorySegment getTabulatedFunction(int index) {
    return OpenMMNative.OpenMM_CustomCentroidBondForce_getTabulatedFunction(getPointer(), index);
  }

  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomCentroidBondForce_getTabulatedFunctionName(getPointer(), index));
  }

  /** Replace bond group indices and per-bond values; native arrays are borrowed synchronously.
   * @param index stored bond index
   * @param groups native group-index array of the fixed group count
   * @param parameters per-bond values in declaration order */
  public void setBondParameters(int index, IntArray groups, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_setBondParameters(
        getPointer(), index, groups.getPointer(), parameters.getPointer());
  }

  /** Replace bond values from caller-owned native segments borrowed for this call.
   * @param index stored bond index
   * @param groups native group-index array
   * @param parameters native per-bond value array */
  public void setBondParameters(
      int index, MemorySegment groups, MemorySegment parameters) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_setBondParameters(
        getPointer(), index, groups, parameters);
  }

  /** Replace the energy expression; existing contexts must be recreated to use it.
   * @param energy OpenMM custom energy expression */
  public void setEnergyFunction(String energy) {
    withString(energy, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_setEnergyFunction(getPointer(), value));
  }

  /** Set the default for future contexts; does not set the value in live contexts.
   * @param index global parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_setGlobalParameterDefaultValue(
        getPointer(), index, value);
  }

  /** Rename a global parameter; existing contexts are not reconfigured.
   * @param index global parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /** Replace a group's centroid definition using caller-owned native arrays.
   * @param index group index
   * @param particles particle-index array
   * @param weights corresponding weights */
  public void setGroupParameters(int index, IntArray particles, DoubleArray weights) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_setGroupParameters(
        getPointer(), index, particles.getPointer(), weights.getPointer());
  }

  /** Replace a group's centroid definition using caller-owned native segments borrowed for this call.
   * @param index group index
   * @param particles native particle-index array
   * @param weights native corresponding weight array */
  public void setGroupParameters(int index, MemorySegment particles, MemorySegment weights) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_setGroupParameters(
        getPointer(), index, particles, weights);
  }

  /** Rename a per-bond parameter; existing contexts are not reconfigured.
   * @param index per-bond parameter index
   * @param name new expression name */
  public void setPerBondParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomCentroidBondForce_setPerBondParameterName(
            getPointer(), index, value));
  }

  /** Set the native integer periodic flag (zero disables; nonzero enables periodic displacements).
   * @param periodic native boolean integer flag */
  public void setUsesPeriodicBoundaryConditions(int periodic) {
    OpenMMNative.OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(
        getPointer(), periodic);
  }

  /** Set whether centroid distances use periodic minimum-image displacements.
   * @param periodic {@code true} enables periodic displacements */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    setUsesPeriodicBoundaryConditions(OpenMMBooleans.toNative(periodic) != 0 ? 1 : 0);
  }

  /**
   * Apply only per-bond parameter and tabulated-value changes in an existing context. Centroid
   * group definitions and the groups used by a bond cannot be changed, and bonds cannot be added.
   * For a tabulated function, dimensions and domain/range must remain unchanged.
   *
   * @param context context containing this force; no-op if it has no native context handle.
   *     Expression, topology, group definitions, declarations, and global defaults are not
   *     propagated.
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_CustomCentroidBondForce_updateParametersInContext(
                getPointer(), pointer));
  }

  /** @return {@code true} if periodic displacements are enabled. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomCentroidBondForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(int numGroups, String energy) {
    return withStringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomCentroidBondForce_create(numGroups, value);
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
