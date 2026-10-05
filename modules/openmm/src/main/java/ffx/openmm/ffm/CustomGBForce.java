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
 * General-purpose custom generalized-Born force assembled from computed values and energy terms.
 *
 * <p>Expressions use OpenMM's custom-expression language and may reference geometry (including
 * {@code r}), global/per-particle parameters, and computed values in their defined dependency
 * order. A {@code SingleParticle} computation is evaluated independently per particle;
 * {@code ParticlePair} sums over non-excluded pairs; {@code ParticlePairNoExclusions} sums over
 * all pairs including excluded pairs. Integer values for these header enums are 0, 1, and 2,
 * respectively. Nonbonded method values are 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic},
 * and 2 {@code CutoffPeriodic}. Pair cutoff distances are in nm; energy is in kJ/mol, while
 * custom parameter/computed-value units are determined by the expressions.
 *
 * <p>The native header deprecates its legacy continuous-function API in favor of tabulated
 * functions; the matching Java overloads remain for compatibility.
 *
 * <p>Java parameter arrays are copied; native array inputs are borrowed for the call. Tabulated
 * functions added through either overload transfer native ownership to this force. The typed
 * overload invalidates its wrapper; the raw-handle overload cannot invalidate an associated Java
 * wrapper, so do not destroy that handle independently. Returned handles are borrowed. The computed-value,
 * energy-term, and continuous
 * function getters are explicitly unsupported: their C wrappers expose {@code char**} outputs
 * that are implemented by writing C++ {@code std::string} objects, so that output ABI is
 * incompatible and unsafe. Context update skips contexts without a native pointer and only
 * propagates OpenMM-supported parameter edits, not expression structure or global defaults.
 */
public class CustomGBForce extends Force {

  /** Snapshot of a computed-value definition.
   * @param name expression-visible value name
   * @param expression expression used to compute the value
   * @param type OpenMM computation-type enum: 0 single particle, 1 pair sum excluding exclusions,
   *     2 pair sum including exclusions */
  public record ComputedValueParameters(String name, String expression, int type) {}
  /** Snapshot of an energy-term definition.
   * @param expression energy expression
   * @param type OpenMM computation-type enum: 0 single particle, 1 pair sum excluding exclusions,
   *     2 pair sum including exclusions */
  public record EnergyTermParameters(String expression, int type) {}
  /** @param particle1 first excluded particle index
   *  @param particle2 second excluded particle index */
  public record Exclusion(int particle1, int particle2) {}
  /** Snapshot of a continuous one-dimensional function.
   * @param name function name
   * @param values copied tabulated values
   * @param min lower bound of the function's input domain
   * @param max upper bound of the function's input domain */
  public record FunctionParameters(String name, double[] values, double min, double max) {}

  /** Create an empty force; computed values and energy terms may be added afterward. */
  public CustomGBForce() { super(create()); }

  /** Add a computed value with its expression and evaluation mode.
   * @param name value name usable in later expressions
   * @param expression expression evaluated according to {@code type}
   * @param type 0 single-particle, 1 sum over non-excluded pairs, or 2 sum over all pairs
   * @return computed-value index */
  public int addComputedValue(String name, String expression, int type) {
    return strings(name, expression, (n,e) -> OpenMMNative.OpenMM_CustomGBForce_addComputedValue(getPointer(), n, e, type));
  }
  /** Add an energy expression evaluated using the selected computation mode.
   * @param expression energy expression
   * @param type 0 single-particle, 1 pair sum excluding exclusions, or 2 pair sum including them
   * @return energy-term index */
  public int addEnergyTerm(String expression, int type) {
    return withStringResult(expression, e -> OpenMMNative.OpenMM_CustomGBForce_addEnergyTerm(getPointer(), e, type));
  }
  /** Request the energy derivative with respect to a declared global parameter.
   * @param name exact global parameter name */
  public void addEnergyParameterDerivative(String name) {
    withString(name, n -> OpenMMNative.OpenMM_CustomGBForce_addEnergyParameterDerivative(getPointer(), n));
  }
  /** Declare a parameter for which every particle supplies one value in declaration order.
   * @param name expression name
   * @return declaration index */
  public int addPerParticleParameter(String name) {
    return withStringResult(name, n -> OpenMMNative.OpenMM_CustomGBForce_addPerParticleParameter(getPointer(), n));
  }
  /** Declare a global parameter and its default for newly created contexts.
   * @param name expression name
   * @param defaultValue initial value
   * @return declaration index */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, n -> OpenMMNative.OpenMM_CustomGBForce_addGlobalParameter(getPointer(), n, defaultValue));
  }
  /** Add one particle's values in per-particle declaration order; the Java array is copied.
   * @param parameters one value per declared per-particle parameter
   * @return particle index */
  public int addParticle(double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) { return addParticle(values); }
  }
  /** Add one particle using a caller-owned native array borrowed for this call.
   * @param parameters native per-particle values in declaration order
   * @return particle index */
  public int addParticle(DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomGBForce_addParticle(getPointer(), parameters.getPointer());
  }
  /** Add one particle from a caller-owned native segment borrowed for this call.
   * @param parameters native per-particle values in declaration order
   * @return particle index */
  public int addParticle(MemorySegment parameters) {
    return OpenMMNative.OpenMM_CustomGBForce_addParticle(getPointer(), parameters);
  }
  /** Exclude a particle pair from computations using {@code ParticlePair}.
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @return exclusion index */
  public int addExclusion(int particle1, int particle2) {
    return OpenMMNative.OpenMM_CustomGBForce_addExclusion(getPointer(), particle1, particle2);
  }
  /** Add a tabulated function and transfer native ownership to this force.
   * @param name expression-visible function name
   * @param function wrapper invalidated after transfer
   * @return function index */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = withStringResult(name, n -> OpenMMNative.OpenMM_CustomGBForce_addTabulatedFunction(getPointer(), n, function.getPointer()));
    function.invalidate();
    return index;
  }
  /** Add a tabulated function using a caller-owned name segment and native function handle.
   * Native function ownership transfers to this force; an associated wrapper is not invalidated.
   * @param name NUL-terminated UTF-8 function-name segment, borrowed during the call
   * @param function native function handle */
  public int addTabulatedFunction(MemorySegment name, MemorySegment function) {
    return OpenMMNative.OpenMM_CustomGBForce_addTabulatedFunction(
        getPointer(), name, function);
  }
  /** Add a legacy continuous one-dimensional function (deprecated by the native API); prefer
   * {@link #addTabulatedFunction(String, TabulatedFunction)}.
   * @param name expression-visible function name
   * @param values sampled function values
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  public int addFunction(String name, double[] values, double min, double max) {
    try (DoubleArray nativeValues = CustomForceParameters.toNative(values)) {
      return withStringResult(name, n -> OpenMMNative.OpenMM_CustomGBForce_addFunction(getPointer(), n, nativeValues.getPointer(), min, max));
    }
  }
  /** Add a continuous function from caller-owned name and value-array segments borrowed for this call.
   * @param name NUL-terminated UTF-8 function name
   * @param values native sampled values
   * @param min lower input-domain bound
   * @param max upper input-domain bound
   * @return function index */
  public int addFunction(
      MemorySegment name, MemorySegment values, double min, double max) {
    return OpenMMNative.OpenMM_CustomGBForce_addFunction(
        getPointer(), name, values, min, max);
  }
  /** Destroy the force and release native functions owned by it. */
  @Override public void destroy() { destroy(OpenMMNative::OpenMM_CustomGBForce_destroy); }

  /**
   * Unsupported: the wrapper's {@code char**} output is implemented by writing C++ strings and
   * therefore has an incompatible, unsafe ABI.
   *
   * @param index computed-value index
   * @return never returns
   * @throws UnsupportedOperationException always because of the wrapper ABI mismatch
   */
  public ComputedValueParameters getComputedValueParameters(int index) { throw invalidStringOutput("getComputedValueParameters"); }
  /**
   * Unsupported: the {@code char**} output is implemented as a C++ {@code std::string} output,
   * which is ABI-incompatible.
   *
   * @param index energy-term index
   * @return never returns
   * @throws UnsupportedOperationException always because of the wrapper ABI mismatch
   */
  public EnergyTermParameters getEnergyTermParameters(int index) { throw invalidStringOutput("getEnergyTermParameters"); }
  /**
   * Unsupported: the function-name {@code char**} output is written as a C++ string and cannot
   * safely be read through this wrapper ABI.
   *
   * @param index continuous-function index
   * @return never returns
   * @throws UnsupportedOperationException always because of the wrapper ABI mismatch
   */
  public FunctionParameters getFunctionParameters(int index) { throw invalidStringOutput("getFunctionParameters"); }
  /** @param index exclusion index
   *  @return both particle indices in the excluded pair */
  public Exclusion getExclusionParticles(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment p1=arena.allocate(ValueLayout.JAVA_INT),p2=arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomGBForce_getExclusionParticles(getPointer(),index,p1,p2);
      return new Exclusion(p1.get(ValueLayout.JAVA_INT,0),p2.get(ValueLayout.JAVA_INT,0));
    }
  }
  /** @param index particle index
   *  @return copied per-particle values in declaration order */
  public double[] getParticleParameters(int index) {
    try (DoubleArray values = new DoubleArray(0)) {
      OpenMMNative.OpenMM_CustomGBForce_getParticleParameters(getPointer(), index, values.getPointer());
      return CustomForceParameters.copy(values);
    }
  }
  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomGBForce_getGlobalParameterName(getPointer(),index)); }
  /** @param index global-parameter index
   *  @return default value for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) { return OpenMMNative.OpenMM_CustomGBForce_getGlobalParameterDefaultValue(getPointer(),index); }
  /** @param index per-particle declaration index
   *  @return copied parameter name */
  public String getPerParticleParameterName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomGBForce_getPerParticleParameterName(getPointer(),index)); }
  /** @param index derivative index
   *  @return copied differentiated parameter name */
  public String getEnergyParameterDerivativeName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomGBForce_getEnergyParameterDerivativeName(getPointer(),index)); }
  /** @param index function index
   *  @return borrowed force-owned handle; do not destroy it */
  public MemorySegment getTabulatedFunction(int index) { return OpenMMNative.OpenMM_CustomGBForce_getTabulatedFunction(getPointer(),index); }
  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomGBForce_getTabulatedFunctionName(getPointer(),index)); }

  /** @return number of particles */ public int getNumParticles(){return OpenMMNative.OpenMM_CustomGBForce_getNumParticles(getPointer());}
  /** @return number of exclusions */ public int getNumExclusions(){return OpenMMNative.OpenMM_CustomGBForce_getNumExclusions(getPointer());}
  /** @return per-particle parameter declarations */ public int getNumPerParticleParameters(){return OpenMMNative.OpenMM_CustomGBForce_getNumPerParticleParameters(getPointer());}
  /** @return global parameter declarations */ public int getNumGlobalParameters(){return OpenMMNative.OpenMM_CustomGBForce_getNumGlobalParameters(getPointer());}
  /** @return requested energy derivatives */ public int getNumEnergyParameterDerivatives(){return OpenMMNative.OpenMM_CustomGBForce_getNumEnergyParameterDerivatives(getPointer());}
  /** @return registered tabulated functions */ public int getNumTabulatedFunctions(){return OpenMMNative.OpenMM_CustomGBForce_getNumTabulatedFunctions(getPointer());}
  /** @return legacy function count; prefer {@link #getNumTabulatedFunctions()} */
  public int getNumFunctions(){return OpenMMNative.OpenMM_CustomGBForce_getNumFunctions(getPointer());}
  /** @return computed-value definitions */ public int getNumComputedValues(){return OpenMMNative.OpenMM_CustomGBForce_getNumComputedValues(getPointer());}
  /** @return energy-term definitions */ public int getNumEnergyTerms(){return OpenMMNative.OpenMM_CustomGBForce_getNumEnergyTerms(getPointer());}
  /**
   * @return OpenMM enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2
   *     {@code CutoffPeriodic}
   */
  public int getNonbondedMethod(){return OpenMMNative.OpenMM_CustomGBForce_getNonbondedMethod(getPointer());}
  /** Set the OpenMM nonbonded method: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2
   * {@code CutoffPeriodic}.
   * @param method native enum value */
  public void setNonbondedMethod(int method){OpenMMNative.OpenMM_CustomGBForce_setNonbondedMethod(getPointer(),method);}
  /** @return interaction cutoff in nm */ public double getCutoffDistance(){return OpenMMNative.OpenMM_CustomGBForce_getCutoffDistance(getPointer());}
  /** Set the interaction cutoff.
   * @param distance cutoff distance in nm */
  public void setCutoffDistance(double distance){OpenMMNative.OpenMM_CustomGBForce_setCutoffDistance(getPointer(),distance);}

  /** Replace a computed-value definition.
   * @param index computed-value index
   * @param name expression-visible name
   * @param expression replacement expression
   * @param type 0 single particle, 1 pair sum excluding exclusions, 2 pair sum including exclusions */
  public void setComputedValueParameters(int index,String name,String expression,int type) {
    strings(name,expression,(n,e)->{OpenMMNative.OpenMM_CustomGBForce_setComputedValueParameters(getPointer(),index,n,e,type);return null;});
  }
  /** Replace an energy-term expression and evaluation mode.
   * @param index energy-term index
   * @param expression replacement expression
   * @param type 0 single particle, 1 pair sum excluding exclusions, 2 pair sum including exclusions */
  public void setEnergyTermParameters(int index,String expression,int type) {
    withString(expression,e->OpenMMNative.OpenMM_CustomGBForce_setEnergyTermParameters(getPointer(),index,e,type));
  }
  /** Replace all per-particle values; Java values are copied.
   * @param index particle index
   * @param parameters one value per declared per-particle parameter, in declaration order */
  public void setParticleParameters(int index,double[] parameters) {
    try(DoubleArray values=CustomForceParameters.toNative(parameters)){setParticleParameters(index,values);}
  }
  /** Replace per-particle values from a native array borrowed during the call.
   * @param index particle index
   * @param parameters values in declaration order */
  public void setParticleParameters(int index,DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomGBForce_setParticleParameters(getPointer(),index,parameters.getPointer());
  }
  /** Replace per-particle values from a native segment borrowed during the call.
   * @param index particle index
   * @param parameters caller-owned values in declaration order */
  public void setParticleParameters(int index, MemorySegment parameters) {
    OpenMMNative.OpenMM_CustomGBForce_setParticleParameters(
        getPointer(), index, parameters);
  }
  /** Replace an exclusion pair.
   * @param index exclusion index
   * @param particle1 first particle
   * @param particle2 second particle */
  public void setExclusionParticles(int index,int particle1,int particle2){OpenMMNative.OpenMM_CustomGBForce_setExclusionParticles(getPointer(),index,particle1,particle2);}
  /** Replace a legacy continuous-function definition, deprecated by the native API; prefer
   * setting parameters on its {@link TabulatedFunction}.
   * @param index function index
   * @param name expression-visible function name
   * @param values sampled values
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  public void setFunctionParameters(int index,String name,double[] values,double min,double max){
    try(DoubleArray nativeValues=CustomForceParameters.toNative(values)){
      withString(name,n->OpenMMNative.OpenMM_CustomGBForce_setFunctionParameters(getPointer(),index,n,nativeValues.getPointer(),min,max));
    }
  }
  /** Replace a function from caller-owned name and value-array segments borrowed during the call.
   * @param index function index
   * @param name NUL-terminated UTF-8 name
   * @param values native samples
   * @param min lower input-domain bound
   * @param max upper input-domain bound */
  public void setFunctionParameters(
      int index, MemorySegment name, MemorySegment values, double min, double max) {
    OpenMMNative.OpenMM_CustomGBForce_setFunctionParameters(
        getPointer(), index, name, values, min, max);
  }
  /** Rename a global parameter; live contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new name */
  public void setGlobalParameterName(int index,String name){withString(name,n->OpenMMNative.OpenMM_CustomGBForce_setGlobalParameterName(getPointer(),index,n));}
  /** Set the default for future contexts, not live context values.
   * @param index global-parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index,double value){OpenMMNative.OpenMM_CustomGBForce_setGlobalParameterDefaultValue(getPointer(),index,value);}
  /** Rename a per-particle declaration; live contexts are not reconfigured.
   * @param index declaration index
   * @param name new name */
  public void setPerParticleParameterName(int index,String name){withString(name,n->OpenMMNative.OpenMM_CustomGBForce_setPerParticleParameterName(getPointer(),index,n));}
  /** Propagate only per-particle values and tabulated-function values; new particles cannot be
   * added and function dimensions/domain/range cannot change.
   * @param context context containing this force; no-op without a native pointer. Structural
   *     definitions and global defaults are not updated */
  public void updateParametersInContext(Context context){CustomForceParameters.updateContext(context,c->OpenMMNative.OpenMM_CustomGBForce_updateParametersInContext(getPointer(),c));}
  /** @return whether this force uses periodic boundary conditions. */
  @Override public boolean usesPeriodicBoundaryConditions(){return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_CustomGBForce_usesPeriodicBoundaryConditions(getPointer()));}

  private static MemorySegment create(){OpenMMRuntime.initialize();return OpenMMNative.OpenMM_CustomGBForce_create.makeInvoker().apply();}
  private static UnsupportedOperationException invalidStringOutput(String method){return new UnsupportedOperationException("OpenMM_CustomGBForce_"+method+" has an incompatible char** output ABI in this OpenMM wrapper.");}
  private static void withString(String value,java.util.function.Consumer<MemorySegment> call){OpenMMStrings.withUtf8String(value,call);}
  private static <T>T withStringResult(String value,java.util.function.Function<MemorySegment,T> call){return OpenMMStrings.withUtf8StringResult(value,call);}
  private static <T>T strings(String first,String second,java.util.function.BiFunction<MemorySegment,MemorySegment,T> call){
    return withStringResult(first,a->withStringResult(second,b->call.apply(a,b)));
  }
}
