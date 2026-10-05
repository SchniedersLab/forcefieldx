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
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;

/**
 * Custom pairwise nonbonded interaction with global and per-particle parameters.
 *
 * <p>The energy expression is evaluated for each eligible particle pair and may use {@code r},
 * distance in nm, and parameter names with {@code 1}/{@code 2} suffixes for the first and second
 * particles. OpenMM custom computed values may define reusable expressions; interaction groups
 * restrict pair selection, and exclusions suppress listed pairs. Energy is in kJ/mol; parameter
 * units are not converted. Nonbonded methods are enum 0 {@code NoCutoff}, 1
 * {@code CutoffNonPeriodic}, and 2 {@code CutoffPeriodic}; cutoffs and switching distances are nm.
 * Long-range correction and switching are separate boolean options.
 *
 * <p>The native header deprecates legacy continuous-function methods and {@code getNumFunctions}
 * in favor of tabulated-function APIs; the Java overloads remain for compatibility.
 *
 * <p>Java arrays are copied into temporary native storage. {@link DoubleArray} and set inputs are
 * borrowed for the synchronous call; returned arrays are independent Java copies. Adding a
 * tabulated function transfers native ownership: the typed overload invalidates the wrapper,
 * while the raw-handle overload cannot invalidate an associated wrapper, which must not be
 * destroyed independently. Getter handles are borrowed. The function/computed-value output getters are
 * unsupported because the wrapper uses {@code char**} outputs where its implementation writes
 * C++ {@code std::string}; do not call them expecting usable results. Context updates propagate
 * only per-particle parameter and tabulated-function values, not
 * expression/declaration changes, exclusions, groups, or global defaults; no-op if the context
 * lacks a native handle.
 */
public class CustomNonbondedForce extends Force {

  /** @param particle1 first particle index
   *  @param particle2 second particle index */
  public record Exclusion(int particle1, int particle2) {}

  /** Create a pairwise custom force.
   * @param energy expression; {@code r} is pair distance in nm and per-particle parameters use
   *     particle suffixes */
  public CustomNonbondedForce(String energy) {
    super(create(energy));
  }

  /** Request energy differentiation with respect to an already declared global parameter.
   * @param name exact global parameter name */
  public void addEnergyParameterDerivative(String name) {
    withString(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(getPointer(), s));
  }

  /** Declare an expression-wide parameter and default for new contexts.
   * @param name expression name
   * @param value default value
   * @return declaration index */
  public int addGlobalParameter(String name, double value) {
    return withStringResult(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_addGlobalParameter(getPointer(), s, value));
  }

  /** Add a tabulated function and transfer ownership to this force.
   * @param name expression-visible name
   * @param function wrapper invalidated after native ownership transfer
   * @return function index */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = withStringResult(name,
        s -> OpenMMNative.OpenMM_CustomNonbondedForce_addTabulatedFunction(
            getPointer(), s, function.getPointer()));
    function.invalidate();
    return index;
  }

  /** Add a tabulated function from caller-owned handles borrowed during the call.
   * @param name NUL-terminated UTF-8 name segment
   * @param function native function handle */
  public int addTabulatedFunction(MemorySegment name, MemorySegment function) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addTabulatedFunction(
        getPointer(), name, function);
  }

  /** Declare one value for every particle in the system.
   * @param name expression parameter name
   * @return declaration index */
  public int addPerParticleParameter(String name) {
    return withStringResult(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_addPerParticleParameter(getPointer(), s));
  }

  /** Add a particle's per-particle values in declaration order; the Java array is copied.
   * @param parameters one value per per-particle declaration
   * @return particle index */
  public int addParticle(double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      return addParticle(values);
    }
  }

  /** Add particle values from a native array borrowed during the call.
   * @param parameters values in declaration order
   * @return particle index */
  public int addParticle(DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addParticle(getPointer(), parameters.getPointer());
  }

  /** Add values from a caller-owned native segment borrowed during the call.
   * @param parameters native per-particle values
   * @return particle index */
  public int addParticle(MemorySegment parameters) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addParticle(getPointer(), parameters);
  }

  /** Exclude one particle pair from the interaction.
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @return exclusion index */
  public int addExclusion(int particle1, int particle2) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addExclusion(getPointer(), particle1, particle2);
  }

  /** Add a legacy continuous one-dimensional function (deprecated by the native API); prefer
   * {@link #addTabulatedFunction(String, TabulatedFunction)}.
   * @param name expression-visible function name
   * @param values function samples
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  public int addFunction(String name, double[] values, double min, double max) {
    try (DoubleArray nativeValues = CustomForceParameters.toNative(values)) {
      return withStringResult(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_addFunction(
          getPointer(), s, nativeValues.getPointer(), min, max));
    }
  }

  /** Add a function from caller-owned name and native function handles. Native function ownership
   * transfers to this force, but an associated Java wrapper is not invalidated automatically.
   * @param name native function name
   * @param values native function samples
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  public int addFunction(
      MemorySegment name, MemorySegment values, double min, double max) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addFunction(
        getPointer(), name, values, min, max);
  }

  /** Declare a named pairwise computed expression available to the energy expression.
   * @param name computed-value name
   * @param expression expression defining the value
   * @return computed-value index */
  public int addComputedValue(String name, String expression) {
    return withStringResult(name, n -> withStringResult(expression, e ->
        OpenMMNative.OpenMM_CustomNonbondedForce_addComputedValue(getPointer(), n, e)));
  }

  /** Restrict evaluation to pairs with one member in each supplied particle set.
   * @param set1 first native set
   * @param set2 second native set
   * @return interaction-group index */
  public int addInteractionGroup(IntSet set1, IntSet set2) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addInteractionGroup(
        getPointer(), set1.getPointer(), set2.getPointer());
  }

  /** Add a group from caller-owned native set handles borrowed during this call.
   * @param set1 first set handle
   * @param set2 second set handle
   * @return interaction-group index */
  public int addInteractionGroup(MemorySegment set1, MemorySegment set2) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_addInteractionGroup(
        getPointer(), set1, set2);
  }

  /** Add exclusions for pairs connected by at most {@code bondCutoff} bonds.
   * @param bonds native bond list
   * @param bondCutoff maximum graph distance in bonds */
  public void createExclusionsFromBonds(BondArray bonds, int bondCutoff) {
    OpenMMNative.OpenMM_CustomNonbondedForce_createExclusionsFromBonds(
        getPointer(), bonds.getPointer(), bondCutoff);
  }

  /** Destroy the native force and its owned tabulated functions. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomNonbondedForce_destroy);
  }

  /** @return copied current energy expression */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomNonbondedForce_getEnergyFunction(getPointer()));
  }

  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomNonbondedForce_getGlobalParameterName(getPointer(), index));
  }

  /** @param index global-parameter index
   *  @return default for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_getGlobalParameterDefaultValue(getPointer(), index);
  }

  /** @param index per-particle declaration index
   *  @return copied parameter name */
  public String getPerParticleParameterName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomNonbondedForce_getPerParticleParameterName(getPointer(), index));
  }

  /** @param index derivative index
   *  @return copied differentiated global parameter name */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomNonbondedForce_getEnergyParameterDerivativeName(getPointer(), index));
  }

  /** @param index particle index
   *  @return copied values in declaration order */
  public double[] getParticleParameters(int index) {
    try (DoubleArray values = new DoubleArray(0)) {
      OpenMMNative.OpenMM_CustomNonbondedForce_getParticleParameters(getPointer(), index, values.getPointer());
      return CustomForceParameters.copy(values);
    }
  }

  /** @param index exclusion index
   *  @return excluded particle indices */
  public Exclusion getExclusionParticles(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment p1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment p2 = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomNonbondedForce_getExclusionParticles(getPointer(), index, p1, p2);
      return new Exclusion(p1.get(ValueLayout.JAVA_INT, 0), p2.get(ValueLayout.JAVA_INT, 0));
    }
  }

  /** @param index function index
   *  @return borrowed force-owned handle; do not destroy it */
  public MemorySegment getTabulatedFunction(int index) {
    return OpenMMNative.OpenMM_CustomNonbondedForce_getTabulatedFunction(getPointer(), index);
  }

  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomNonbondedForce_getTabulatedFunctionName(getPointer(), index));
  }

  /** Copy an interaction group's two sets into caller-owned set objects.
   * @param index interaction-group index
   * @param set1 output set for the first particle group
   * @param set2 output set for the second particle group */
  public void getInteractionGroupParameters(int index, IntSet set1, IntSet set2) {
    OpenMMNative.OpenMM_CustomNonbondedForce_getInteractionGroupParameters(
        getPointer(), index, set1.getPointer(), set2.getPointer());
  }

  /** Copy an interaction group into caller-owned set handles borrowed for this call.
   * @param index interaction-group index
   * @param set1 output native set handle
   * @param set2 output native set handle */
  public void getInteractionGroupParameters(
      int index, MemorySegment set1, MemorySegment set2) {
    OpenMMNative.OpenMM_CustomNonbondedForce_getInteractionGroupParameters(
        getPointer(), index, set1, set2);
  }

  /**
   * Unsupported: the wrapper's {@code char**} name output is implemented by writing a C++
   * {@code std::string}, an incompatible ABI that cannot safely provide this record.
   *
   * @param index continuous-function index
   * @return never returns
   * @throws UnsupportedOperationException always due to the incompatible output ABI
   */
  public FunctionParameters getFunctionParameters(int index) {
    throw invalidStringOutput("getFunctionParameters");
  }

  /**
   * Unsupported: the wrapper's {@code char**} name/expression outputs are implemented by writing
   * C++ strings and are ABI-incompatible.
   *
   * @param index computed-value index
   * @return never returns
   * @throws UnsupportedOperationException always due to the incompatible output ABI
   */
  public ComputedValueParameters getComputedValueParameters(int index) {
    throw invalidStringOutput("getComputedValueParameters");
  }

  /** @return number of particles */
  public int getNumParticles() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumParticles(getPointer()); }
  /** @return number of excluded pairs */
  public int getNumExclusions() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumExclusions(getPointer()); }
  /** @return per-particle declarations and required parameter-vector length */
  public int getNumPerParticleParameters() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumPerParticleParameters(getPointer()); }
  /** @return global parameter declarations */
  public int getNumGlobalParameters() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumGlobalParameters(getPointer()); }
  /** @return registered tabulated functions */
  public int getNumTabulatedFunctions() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumTabulatedFunctions(getPointer()); }
  /** @return legacy function count; prefer {@link #getNumTabulatedFunctions()} */
  public int getNumFunctions() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumFunctions(getPointer()); }
  /** @return computed-value declarations */
  public int getNumComputedValues() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumComputedValues(getPointer()); }
  /** @return interaction groups */
  public int getNumInteractionGroups() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumInteractionGroups(getPointer()); }
  /** @return requested global-parameter energy derivatives */
  public int getNumEnergyParameterDerivatives() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNumEnergyParameterDerivatives(getPointer()); }
  /** @return enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2 {@code CutoffPeriodic} */
  public int getNonbondedMethod() { return OpenMMNative.OpenMM_CustomNonbondedForce_getNonbondedMethod(getPointer()); }
  /** Set enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2 {@code CutoffPeriodic}.
   * @param method native enum value */
  public void setNonbondedMethod(int method) { OpenMMNative.OpenMM_CustomNonbondedForce_setNonbondedMethod(getPointer(), method); }
  /** @return cutoff distance in nm */
  public double getCutoffDistance() { return OpenMMNative.OpenMM_CustomNonbondedForce_getCutoffDistance(getPointer()); }
  /** Set the cutoff distance.
   * @param distance cutoff in nm */
  public void setCutoffDistance(double distance) { OpenMMNative.OpenMM_CustomNonbondedForce_setCutoffDistance(getPointer(), distance); }
  /** @return whether the switching function is enabled */
  public boolean getUseSwitchingFunction() { return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_CustomNonbondedForce_getUseSwitchingFunction(getPointer())); }
  /** Enable or disable switching of the interaction energy at the switching distance.
   * @param use {@code true} to enable switching */
  public void setUseSwitchingFunction(boolean use) { OpenMMNative.OpenMM_CustomNonbondedForce_setUseSwitchingFunction(getPointer(), OpenMMBooleans.toNative(use)); }
  /** @return switching distance in nm */
  public double getSwitchingDistance() { return OpenMMNative.OpenMM_CustomNonbondedForce_getSwitchingDistance(getPointer()); }
  /** Set the distance at which energy switching begins.
   * @param distance switching distance in nm */
  public void setSwitchingDistance(double distance) { OpenMMNative.OpenMM_CustomNonbondedForce_setSwitchingDistance(getPointer(), distance); }
  /** @return whether long-range correction is enabled */
  public boolean getUseLongRangeCorrection() { return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_CustomNonbondedForce_getUseLongRangeCorrection(getPointer())); }
  /** Enable or disable long-range correction.
   * @param use {@code true} to enable correction */
  public void setUseLongRangeCorrection(boolean use) { OpenMMNative.OpenMM_CustomNonbondedForce_setUseLongRangeCorrection(getPointer(), OpenMMBooleans.toNative(use)); }

  /** Replace a particle's values; Java values are copied.
   * @param index particle index
   * @param parameters one value per declaration, in declaration order */
  public void setParticleParameters(int index, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) { setParticleParameters(index, values); }
  }
  /** Replace values from a native array borrowed during the call.
   * @param index particle index
   * @param parameters values in declaration order */
  public void setParticleParameters(int index, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomNonbondedForce_setParticleParameters(getPointer(), index, parameters.getPointer());
  }
  /** Replace values from a caller-owned native segment borrowed during the call.
   * @param index particle index
   * @param parameters native values in declaration order */
  public void setParticleParameters(int index, MemorySegment parameters) {
    OpenMMNative.OpenMM_CustomNonbondedForce_setParticleParameters(
        getPointer(), index, parameters);
  }
  /** Replace an excluded particle pair.
   * @param index exclusion index
   * @param particle1 first particle
   * @param particle2 second particle */
  public void setExclusionParticles(int index, int particle1, int particle2) {
    OpenMMNative.OpenMM_CustomNonbondedForce_setExclusionParticles(getPointer(), index, particle1, particle2);
  }
  /** Replace energy expression; recreate contexts to apply the expression change.
   * @param energy custom pairwise expression; {@code r} is in nm */
  public void setEnergyFunction(String energy) { withString(energy, s -> OpenMMNative.OpenMM_CustomNonbondedForce_setEnergyFunction(getPointer(), s)); }
  /** Rename a global parameter; existing contexts are not reconfigured.
   * @param index parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) { withString(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_setGlobalParameterName(getPointer(), index, s)); }
  /** Set the default for future contexts; live context values are unchanged.
   * @param index parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) { OpenMMNative.OpenMM_CustomNonbondedForce_setGlobalParameterDefaultValue(getPointer(), index, value); }
  /** Rename a per-particle declaration; existing contexts are not reconfigured.
   * @param index declaration index
   * @param name new expression name */
  public void setPerParticleParameterName(int index, String name) { withString(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_setPerParticleParameterName(getPointer(), index, s)); }
  /** Replace a legacy continuous-function definition, deprecated by the native API; prefer the
   * tabulated-function parameter API.
   * @param index function index
   * @param name function name
   * @param values sampled values
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  public void setFunctionParameters(int index, String name, double[] values, double min, double max) {
    try (DoubleArray nativeValues = CustomForceParameters.toNative(values)) {
      withString(name, s -> OpenMMNative.OpenMM_CustomNonbondedForce_setFunctionParameters(
          getPointer(), index, s, nativeValues.getPointer(), min, max));
    }
  }
  /** Replace a function from caller-owned name and sample segments borrowed during the call.
   * @param index function index
   * @param name NUL-terminated UTF-8 name
   * @param values native samples
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  public void setFunctionParameters(
      int index, MemorySegment name, MemorySegment values, double min, double max) {
    OpenMMNative.OpenMM_CustomNonbondedForce_setFunctionParameters(
        getPointer(), index, name, values, min, max);
  }
  /** Replace a computed expression; contexts must be recreated to observe expression changes.
   * @param index computed-value index
   * @param name value name
   * @param expression new expression */
  public void setComputedValueParameters(int index, String name, String expression) {
    withString(name, n -> withString(expression, e ->
        OpenMMNative.OpenMM_CustomNonbondedForce_setComputedValueParameters(
            getPointer(), index, n, e)));
  }
  /** Replace a group using caller-owned set wrappers.
   * @param index interaction-group index
   * @param set1 replacement first particle set
   * @param set2 replacement second particle set */
  public void setInteractionGroupParameters(int index, IntSet set1, IntSet set2) {
    OpenMMNative.OpenMM_CustomNonbondedForce_setInteractionGroupParameters(getPointer(), index, set1.getPointer(), set2.getPointer());
  }
  /** Replace a group using caller-owned set handles borrowed during the call.
   * @param index interaction-group index
   * @param set1 replacement first set handle
   * @param set2 replacement second set handle */
  public void setInteractionGroupParameters(
      int index, MemorySegment set1, MemorySegment set2) {
    OpenMMNative.OpenMM_CustomNonbondedForce_setInteractionGroupParameters(
        getPointer(), index, set1, set2);
  }
  /** Apply only per-particle parameter and tabulated-function value changes to a context. New
   * particles cannot be added and function dimensions/domain/range cannot change.
   * @param context context containing this force; no-op without a native context handle. Does not
   *     propagate expressions, nonbonded method, cutoff, declarations, exclusions, interaction
   *     groups, or global defaults. */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(context, c -> OpenMMNative.OpenMM_CustomNonbondedForce_updateParametersInContext(getPointer(), c));
  }
  /** @return whether the current nonbonded method uses periodic boundary conditions. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_CustomNonbondedForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /** Continuous function definition snapshot.
   * @param name function name
   * @param values copied samples
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  public record FunctionParameters(String name, double[] values, double min, double max) {}
  /** Computed-value expression snapshot.
   * @param name computed-value name
   * @param expression expression defining its value */
  public record ComputedValueParameters(String name, String expression) {}
  private static MemorySegment create(String energy) {
    return withStringResult(energy, s -> { OpenMMRuntime.initialize(); return OpenMMNative.OpenMM_CustomNonbondedForce_create(s); });
  }
  private static void withString(String value, java.util.function.Consumer<MemorySegment> call) {
    OpenMMStrings.withUtf8String(value, call);
  }
  private static <T> T withStringResult(String value, java.util.function.Function<MemorySegment,T> call) {
    return OpenMMStrings.withUtf8StringResult(value, call);
  }
  private static UnsupportedOperationException invalidStringOutput(String method) {
    return new UnsupportedOperationException("OpenMM_CustomNonbondedForce_" + method +
        " has an incompatible char** output ABI in this OpenMM wrapper.");
  }
}
