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
 * Custom interaction evaluated over selected sets of a fixed number of particles.
 *
 * <p>The expression labels the ordered particles {@code p1}, {@code p2}, and so on, and can
 * reference their coordinates/geometries, global parameters, per-particle parameters, and
 * tabulated functions. Distances are nm, angles/dihedrals radians, and energy kJ/mol. The
 * nonbonded enum values are 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, and 2
 * {@code CutoffPeriodic}; cutoff is in nm. Permutation mode 0 {@code SinglePermutation} evaluates
 * each set once and requires a symmetric expression; mode 1 {@code UniqueCentralParticle}
 * evaluates once for each particle as central {@code p1}, requiring symmetry among the remaining
 * particles. With cutoffs, mode 0 requires all members mutually within cutoff; mode 1 requires
 * each member within cutoff of the central particle.
 *
 * <p>Particle types are integer labels used by per-position type filters; filters constrain which
 * labels may occupy each expression position. Exclusions remove sets containing the excluded pair.
 * Java arrays are copied; native array/set inputs are caller-owned and borrowed synchronously.
 * Adding a function transfers native ownership; the typed overload invalidates its wrapper, while
 * the raw-handle overload cannot invalidate an associated wrapper, which must not be destroyed
 * independently. Function handles returned by getters are borrowed. Context updates affect supported
 * per-particle values only, not expression, types/selection structure, declarations, or global
 * defaults; an uninitialized context is ignored.
 */
public class CustomManyParticleForce extends Force {

  /**
   * Snapshot of a particle's parameters and type label.
   *
   * @param parameters copied values in per-particle declaration order
   * @param type integer type label used by type filters
   */
  public record ParticleParameters(double[] parameters, int type) {}
  /** @param particle1 first excluded particle index
   *  @param particle2 second excluded particle index */
  public record Exclusion(int particle1, int particle2) {}

  /**
   * Create a force evaluated over sets of a fixed size.
   *
   * @param particlesPerSet number of ordered particle positions in each interaction
   * @param energy expression using the set's {@code p1}, {@code p2}, ... variables
   */
  public CustomManyParticleForce(int particlesPerSet, String energy) {
    super(create(particlesPerSet, energy));
  }

  /** Declare a global expression parameter and its default for new contexts.
   * @param name parameter name
   * @param defaultValue initial value
   * @return declaration index */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, s -> OpenMMNative.OpenMM_CustomManyParticleForce_addGlobalParameter(getPointer(), s, defaultValue));
  }
  /** Declare a global parameter from caller-owned NUL-terminated UTF-8 storage.
   * @param name name segment borrowed during the call
   * @param defaultValue initial value
   * @return declaration index */
  public int addGlobalParameter(MemorySegment name, double defaultValue) {
    return OpenMMNative.OpenMM_CustomManyParticleForce_addGlobalParameter(getPointer(), name, defaultValue);
  }
  /** Declare one per-particle value; each particle supplies values in declaration order.
   * @param name expression parameter name
   * @return declaration index */
  public int addPerParticleParameter(String name) {
    return withStringResult(name, s -> OpenMMNative.OpenMM_CustomManyParticleForce_addPerParticleParameter(getPointer(), s));
  }
  /** Declare a value from caller-owned NUL-terminated UTF-8 storage.
   * @param name name segment borrowed during the call
   * @return declaration index */
  public int addPerParticleParameter(MemorySegment name) {
    return OpenMMNative.OpenMM_CustomManyParticleForce_addPerParticleParameter(getPointer(), name);
  }
  /** Add a particle with copied per-particle values and an integer type label.
   * @param parameters one value per per-particle declaration
   * @param type type label used by position filters
   * @return particle index */
  public int addParticle(double[] parameters, int type) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) { return addParticle(values, type); }
  }
  /** Add a particle using a native parameter array borrowed during the call.
   * @param parameters values in declaration order
   * @param type type label
   * @return particle index */
  public int addParticle(DoubleArray parameters, int type) {
    return OpenMMNative.OpenMM_CustomManyParticleForce_addParticle(getPointer(), parameters.getPointer(), type);
  }
  /** Add a particle from a caller-owned native segment borrowed during the call.
   * @param parameters native values in declaration order
   * @param type type label
   * @return particle index */
  public int addParticle(MemorySegment parameters, int type) {
    return OpenMMNative.OpenMM_CustomManyParticleForce_addParticle(getPointer(), parameters, type);
  }
  /** Add a tabulated function and transfer its native ownership to the force.
   * @param name expression-visible function name
   * @param function wrapper invalidated after transfer
   * @return function index */
  public int addTabulatedFunction(String name, TabulatedFunction function) {
    int index = withStringResult(name, s -> OpenMMNative.OpenMM_CustomManyParticleForce_addTabulatedFunction(getPointer(), s, function.getPointer()));
    function.invalidate();
    return index;
  }
  /** Add a function using a caller-owned name segment and native function handle. Native ownership
   * transfers to this force; an associated wrapper is not invalidated automatically.
   * @param name NUL-terminated UTF-8 name, borrowed during the call
   * @param function native function handle
   * @return function index */
  public int addTabulatedFunction(MemorySegment name, MemorySegment function) {
    return OpenMMNative.OpenMM_CustomManyParticleForce_addTabulatedFunction(getPointer(), name, function);
  }
  /** Exclude all interaction sets containing this particle pair.
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @return exclusion index */
  public int addExclusion(int particle1, int particle2) {
    return OpenMMNative.OpenMM_CustomManyParticleForce_addExclusion(getPointer(), particle1, particle2);
  }
  /** Create exclusions between particles separated by at most the specified number of bonds.
   * @param bonds native bond list
   * @param bondCutoff maximum bond-graph separation included in exclusions */
  public void createExclusionsFromBonds(BondArray bonds, int bondCutoff) {
    OpenMMNative.OpenMM_CustomManyParticleForce_createExclusionsFromBonds(getPointer(), bonds.getPointer(), bondCutoff);
  }
  @Override public void destroy() { destroy(OpenMMNative::OpenMM_CustomManyParticleForce_destroy); }

  /** @return copied current energy expression */
  public String getEnergyFunction() { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomManyParticleForce_getEnergyFunction(getPointer())); }
  /** @param index particle index
   *  @return copied values in declaration order and its integer type label */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena arena = Arena.ofConfined(); DoubleArray values = new DoubleArray(0)) {
      MemorySegment type = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomManyParticleForce_getParticleParameters(getPointer(), index, values.getPointer(), type);
      return new ParticleParameters(CustomForceParameters.copy(values), type.get(ValueLayout.JAVA_INT, 0));
    }
  }
  /** @param index exclusion index
   *  @return excluded particle indices */
  public Exclusion getExclusionParticles(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment p1 = arena.allocate(ValueLayout.JAVA_INT), p2 = arena.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomManyParticleForce_getExclusionParticles(getPointer(), index, p1, p2);
      return new Exclusion(p1.get(ValueLayout.JAVA_INT,0), p2.get(ValueLayout.JAVA_INT,0));
    }
  }
  /** @param index function index
   *  @return borrowed force-owned native handle; do not destroy it */
  public MemorySegment getTabulatedFunction(int index) { return OpenMMNative.OpenMM_CustomManyParticleForce_getTabulatedFunction(getPointer(), index); }
  /** @param index function index
   *  @return copied function name */
  public String getTabulatedFunctionName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomManyParticleForce_getTabulatedFunctionName(getPointer(), index)); }
  /** Copy the allowed integer type labels for one particle position to a caller-owned set.
   * @param index position in the set, from zero to particles-per-set minus one
   * @param types output set with caller-owned lifetime */
  public void getTypeFilter(int index, IntSet types) { OpenMMNative.OpenMM_CustomManyParticleForce_getTypeFilter(getPointer(), index, types.getPointer()); }
  /** @param index global-parameter index
   *  @return copied name */
  public String getGlobalParameterName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomManyParticleForce_getGlobalParameterName(getPointer(), index)); }
  /** @param index global-parameter index
   *  @return default for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) { return OpenMMNative.OpenMM_CustomManyParticleForce_getGlobalParameterDefaultValue(getPointer(), index); }
  /** @param index per-particle declaration index
   *  @return copied name */
  public String getPerParticleParameterName(int index) { return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomManyParticleForce_getPerParticleParameterName(getPointer(), index)); }
  /** @return fixed number of particles in every evaluated set */
  public int getNumParticlesPerSet() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNumParticlesPerSet(getPointer()); }
  /** @return number of particles added to this force */
  public int getNumParticles() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNumParticles(getPointer()); }
  /** @return number of excluded particle pairs */
  public int getNumExclusions() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNumExclusions(getPointer()); }
  /** @return per-particle parameter declarations and required vector length */
  public int getNumPerParticleParameters() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNumPerParticleParameters(getPointer()); }
  /** @return global parameter declarations */
  public int getNumGlobalParameters() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNumGlobalParameters(getPointer()); }
  /** @return registered tabulated functions */
  public int getNumTabulatedFunctions() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNumTabulatedFunctions(getPointer()); }
  /** @return enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2 {@code CutoffPeriodic} */
  public int getNonbondedMethod() { return OpenMMNative.OpenMM_CustomManyParticleForce_getNonbondedMethod(getPointer()); }
  /** Set enum: 0 {@code NoCutoff}, 1 {@code CutoffNonPeriodic}, or 2 {@code CutoffPeriodic}.
   * @param method native enum value */
  public void setNonbondedMethod(int method) { OpenMMNative.OpenMM_CustomManyParticleForce_setNonbondedMethod(getPointer(), method); }
  /** @return enum: 0 {@code SinglePermutation}, 1 {@code UniqueCentralParticle} */
  public int getPermutationMode() { return OpenMMNative.OpenMM_CustomManyParticleForce_getPermutationMode(getPointer()); }
  /** Set the permutation policy.
   * @param mode 0 {@code SinglePermutation} (symmetric in all particles) or 1
   *     {@code UniqueCentralParticle} (symmetric among all except central {@code p1}) */
  public void setPermutationMode(int mode) { OpenMMNative.OpenMM_CustomManyParticleForce_setPermutationMode(getPointer(), mode); }
  /** @return cutoff distance in nm */
  public double getCutoffDistance() { return OpenMMNative.OpenMM_CustomManyParticleForce_getCutoffDistance(getPointer()); }
  /** Set cutoff distance.
   * @param distance cutoff in nm */
  public void setCutoffDistance(double distance) { OpenMMNative.OpenMM_CustomManyParticleForce_setCutoffDistance(getPointer(), distance); }
  /** Replace expression; existing contexts must be recreated to use it.
   * @param energy custom expression over ordered particle positions */
  public void setEnergyFunction(String energy) { withString(energy, s -> OpenMMNative.OpenMM_CustomManyParticleForce_setEnergyFunction(getPointer(), s)); }
  /** Replace expression from a caller-owned NUL-terminated UTF-8 segment borrowed for the call.
   * @param energy native expression segment */
  public void setEnergyFunction(MemorySegment energy) { OpenMMNative.OpenMM_CustomManyParticleForce_setEnergyFunction(getPointer(), energy); }
  /** Rename a global parameter; live contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) { withString(name, s -> OpenMMNative.OpenMM_CustomManyParticleForce_setGlobalParameterName(getPointer(), index, s)); }
  /** Rename using a caller-owned NUL-terminated UTF-8 segment.
   * @param index global-parameter index
   * @param name native name segment borrowed during the call */
  public void setGlobalParameterName(int index, MemorySegment name) { OpenMMNative.OpenMM_CustomManyParticleForce_setGlobalParameterName(getPointer(), index, name); }
  /** Set the default for future contexts, not existing context values.
   * @param index global-parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) { OpenMMNative.OpenMM_CustomManyParticleForce_setGlobalParameterDefaultValue(getPointer(), index, value); }
  /** Rename a per-particle declaration; live contexts are not reconfigured.
   * @param index declaration index
   * @param name new expression name */
  public void setPerParticleParameterName(int index, String name) { withString(name, s -> OpenMMNative.OpenMM_CustomManyParticleForce_setPerParticleParameterName(getPointer(), index, s)); }
  /** Rename from caller-owned NUL-terminated UTF-8 storage.
   * @param index declaration index
   * @param name name segment borrowed during the call */
  public void setPerParticleParameterName(int index, MemorySegment name) { OpenMMNative.OpenMM_CustomManyParticleForce_setPerParticleParameterName(getPointer(), index, name); }
  /** Replace a particle's copied values and type label.
   * @param index particle index
   * @param parameters one value per declaration, in declaration order
   * @param type integer type label */
  public void setParticleParameters(int index, double[] parameters, int type) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) { setParticleParameters(index, values, type); }
  }
  /** Replace particle values from a native array borrowed during the call.
   * @param index particle index
   * @param parameters values in declaration order
   * @param type integer type label */
  public void setParticleParameters(int index, DoubleArray parameters, int type) {
    OpenMMNative.OpenMM_CustomManyParticleForce_setParticleParameters(getPointer(), index, parameters.getPointer(), type);
  }
  /** Replace particle values from a caller-owned segment borrowed during the call.
   * @param index particle index
   * @param parameters native values in declaration order
   * @param type integer type label */
  public void setParticleParameters(int index, MemorySegment parameters, int type) {
    OpenMMNative.OpenMM_CustomManyParticleForce_setParticleParameters(getPointer(), index, parameters, type);
  }
  /** Replace an excluded particle pair.
   * @param index exclusion index
   * @param particle1 first particle
   * @param particle2 second particle */
  public void setExclusionParticles(int index, int particle1, int particle2) { OpenMMNative.OpenMM_CustomManyParticleForce_setExclusionParticles(getPointer(), index, particle1, particle2); }
  /** Set the allowed type labels for one ordered particle position.
   * @param index position in the set, zero-based
   * @param types caller-owned native set borrowed synchronously */
  public void setTypeFilter(int index, IntSet types) { OpenMMNative.OpenMM_CustomManyParticleForce_setTypeFilter(getPointer(), index, types.getPointer()); }
  /** Propagate only per-particle values and tabulated-function values to an existing context.
   * New particles cannot be added and function dimensions/domain/range cannot change.
   * @param context context containing this force; no-op without a native handle. Does not update
   *     expression, nonbonded method, cutoff, type filters, particle topology, declarations, or
   *     global defaults. */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(context, c -> OpenMMNative.OpenMM_CustomManyParticleForce_updateParametersInContext(getPointer(), c));
  }
  /** @return whether the selected nonbonded method uses periodic boundary conditions. */
  @Override public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_CustomManyParticleForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(int particlesPerSet, String energy) {
    return withStringResult(energy, s -> { OpenMMRuntime.initialize(); return OpenMMNative.OpenMM_CustomManyParticleForce_create(particlesPerSet, s); });
  }
  private static void withString(String value, java.util.function.Consumer<MemorySegment> call) { OpenMMStrings.withUtf8String(value, call); }
  private static <T> T withStringResult(String value, java.util.function.Function<MemorySegment,T> call) { return OpenMMStrings.withUtf8StringResult(value, call); }
}
