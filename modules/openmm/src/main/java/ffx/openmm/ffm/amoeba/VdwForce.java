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
package ffx.openmm.ffm.amoeba;

import ffx.openmm.ffm.Context;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.IntArray;
import ffx.openmm.ffm.OpenMMBooleans;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.OpenMMStrings;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.Objects;

/**
 * Models AMOEBA van der Waals interactions using buffered 14-7 or Lennard-Jones 12-6
 * potentials.
 *
 * <p>Native defaults are NoCutoff, Buffered147, alchemical method None, dispersion correction
 * enabled, per-particle rather than type-based parameters, sigma combining rule CUBIC-MEAN,
 * epsilon combining rule HHG, softcore power 5, softcore alpha 0.7, and context parameter name
 * {@code AmoebaVdwLambda}. Parameters may be specified per particle, or by particle type with
 * optional type-pair overrides. Interaction sites can be displaced from their parent particles
 * with a reduction factor. Alchemical interactions use the force's context parameter, whose name
 * is returned by {@link #getLambda()}.</p>
 *
 * <p>The native API accepts an alchemical flag for each particle and uses the configured
 * alchemical method: None (0) leaves interactions unchanged, Decouple (1) preserves full-strength
 * interactions between two alchemical particles, and Annihilate (2) makes interactions involving
 * alchemical particles lambda-dependent, turning off alchemical-pair interactions at lambda zero.
 * The lambda is supplied as a Context parameter in the range [0, 1].
 * Particle parameter changes are copied into an existing Context with
 * {@link #updateParametersInContext(Context)}; non-particle settings require Context
 * reinitialization. Lengths (sigma and cutoff) are in nanometers; epsilon values are in
 * kilojoules per mole.</p>
 */
public class VdwForce extends Force {

  /**
   * Create an AMOEBA van der Waals force.
   *
   * @throws IllegalStateException if the native force cannot be created.
   */
  public VdwForce() {
    super(create());
  }

  /**
   * Add per-particle van der Waals parameters.
   *
   * @param parentIndex parent particle index.
   * @param sigma van der Waals sigma in nanometers.
   * @param epsilon van der Waals well depth in kilojoules per mole.
   * @param reductionFactor fraction of the parent-to-particle distance for the interaction site.
   * @param isAlchemical whether this particle is alchemical.
   * @param scaleFactor dimensionless scale applied to interactions involving this particle.
   * @return index of the added particle.
   */
  public int addParticle(int parentIndex, double sigma, double epsilon, double reductionFactor,
                         boolean isAlchemical, double scaleFactor) {
    return OpenMMNative.OpenMM_AmoebaVdwForce_addParticle(getPointer(), parentIndex, sigma, epsilon,
        reductionFactor, OpenMMBooleans.toNative(isAlchemical), scaleFactor);
  }

  /**
   * Add per-particle van der Waals parameters using the legacy integer boolean representation.
   *
   * @param parentIndex parent particle index.
   * @param sigma van der Waals sigma in nanometers.
   * @param epsilon van der Waals well depth in kilojoules per mole.
   * @param reductionFactor fraction of the parent-to-particle distance for the interaction site.
   * @param isAlchemical nonzero when this particle is alchemical.
   * @param scaleFactor dimensionless scale applied to interactions involving this particle.
   * @return index of the added particle.
   */
  public int addParticle(int parentIndex, double sigma, double epsilon, double reductionFactor,
                         int isAlchemical, double scaleFactor) {
    return addParticle(parentIndex, sigma, epsilon, reductionFactor, isAlchemical != 0, scaleFactor);
  }

  /**
   * Add a particle using type-based parameters.
   *
   * @param parentIndex parent particle index.
   * @param typeIndex particle type index.
   * @param reductionFactor dimensionless fraction of the parent-to-particle distance for the interaction site.
   * @param isAlchemical whether this particle is alchemical.
   * @param scaleFactor dimensionless scale applied to interactions involving this particle.
   * @return index of the added particle.
   */
  public int addParticle(int parentIndex, int typeIndex, double reductionFactor,
                         boolean isAlchemical, double scaleFactor) {
    return OpenMMNative.OpenMM_AmoebaVdwForce_addParticle_1(getPointer(), parentIndex, typeIndex,
        reductionFactor, OpenMMBooleans.toNative(isAlchemical), scaleFactor);
  }

  /**
   * Add a particle using type-based parameters and the legacy integer boolean representation.
   *
   * @param parentIndex parent particle index.
   * @param typeIndex particle type index.
   * @param reductionFactor dimensionless fraction of the parent-to-particle distance for the interaction site.
   * @param isAlchemical nonzero when this particle is alchemical.
   * @param scaleFactor dimensionless scale applied to interactions involving this particle.
   * @return index of the added particle.
   */
  public int addParticle(int parentIndex, int typeIndex, double reductionFactor,
                         int isAlchemical, double scaleFactor) {
    return addParticle(parentIndex, typeIndex, reductionFactor, isAlchemical != 0, scaleFactor);
  }

  /**
   * Add a particle type.
   *
   * @param sigma sigma for particles of this type, in nanometers.
   * @param epsilon well depth for particles of this type, in kilojoules per mole.
   * @return index of the added type.
   */
  public int addParticleType(double sigma, double epsilon) {
    return OpenMMNative.OpenMM_AmoebaVdwForce_addParticleType(getPointer(), sigma, epsilon);
  }

  /**
   * Add a pair-specific override for two particle types.
   *
   * @param type1 first particle type.
   * @param type2 second particle type.
   * @param sigma pair sigma in nanometers.
   * @param epsilon pair well depth in kilojoules per mole.
   * @return index of the added type pair.
   */
  public int addTypePair(int type1, int type2, double sigma, double epsilon) {
    return OpenMMNative.OpenMM_AmoebaVdwForce_addTypePair(getPointer(), type1, type2, sigma, epsilon);
  }

  /**
   * Destroy this force.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_AmoebaVdwForce_destroy);
  }

  /**
   * Get the alchemical method as its native enum value (None 0, Decouple 1, Annihilate 2).
   *
   * @return configured alchemical method.
   */
  public int getAlchemicalMethod() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getAlchemicalMethod(getPointer());
  }

  /**
   * Get the cutoff distance.
   *
   * @return cutoff distance in nanometers.
   * @deprecated Use {@link #getCutoffDistance()}.
   */
  @Deprecated
  public double getCutoff() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getCutoff(getPointer());
  }

  /**
   * Get the nonbonded cutoff distance.
   *
   * @return cutoff distance in nanometers.
   */
  public double getCutoffDistance() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getCutoffDistance(getPointer());
  }

  /**
   * Get the epsilon combining rule.
   *
   * @return copied native combining-rule string (ARITHMETIC, GEOMETRIC, HARMONIC, W-H, or HHG).
   */
  public String getEpsilonCombiningRule() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_AmoebaVdwForce_getEpsilonCombiningRule(getPointer()));
  }

  /**
   * Get the name of the context parameter controlling alchemical interactions.
   *
   * @return copied lambda context parameter name (initially {@code AmoebaVdwLambda}).
   */
  public String getLambda() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_AmoebaVdwForce_Lambda(getPointer()));
  }

  /**
   * Get the nonbonded method as its native enum value (NoCutoff 0, CutoffPeriodic 1).
   *
   * @return configured nonbonded method.
   */
  public int getNonbondedMethod() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getNonbondedMethod(getPointer());
  }

  /**
   * Get the number of particles.
   *
   * @return particle count.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getNumParticles(getPointer());
  }

  /**
   * Get the number of particle types.
   *
   * @return particle type count.
   */
  public int getNumParticleTypes() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getNumParticleTypes(getPointer());
  }

  /**
   * Get the number of type-pair overrides.
   *
   * @return type-pair count.
   */
  public int getNumTypePairs() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getNumTypePairs(getPointer());
  }

  /**
   * Copy the exclusions for a particle into an owned FFM integer array.
   *
   * @param particleIndex particle index.
   * @return caller-owned exclusions array; call {@link IntArray#destroy()} when finished.
   */
  public IntArray getParticleExclusions(int particleIndex) {
    IntArray exclusions = new IntArray(0);
    try {
      OpenMMNative.OpenMM_AmoebaVdwForce_getParticleExclusions(
          getPointer(), particleIndex, exclusions.getPointer());
      return exclusions;
    } catch (RuntimeException | Error exception) {
      exclusions.destroy();
      throw exception;
    }
  }

  /**
   * Get the parameters for a particle.
   *
   * @param particleIndex particle index.
   * @return copied particle parameters.
   */
  public ParticleParameters getParticleParameters(int particleIndex) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment parent = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment sigma = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment reduction = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment alchemical = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment type = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment scale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_AmoebaVdwForce_getParticleParameters(getPointer(), particleIndex, parent,
          sigma, epsilon, reduction, alchemical, type, scale);
      return new ParticleParameters(parent.get(ValueLayout.JAVA_INT, 0),
          sigma.get(ValueLayout.JAVA_DOUBLE, 0), epsilon.get(ValueLayout.JAVA_DOUBLE, 0),
          reduction.get(ValueLayout.JAVA_DOUBLE, 0),
          OpenMMBooleans.fromNative(alchemical.get(ValueLayout.JAVA_INT, 0)),
          type.get(ValueLayout.JAVA_INT, 0), scale.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the potential function (Buffered147 0 or LennardJones 1).
   *
   * @return configured potential function.
   */
  public int getPotentialFunction() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getPotentialFunction(getPointer());
  }

  /**
   * Get the sigma combining rule.
   *
   * @return copied native combining-rule string (ARITHMETIC, GEOMETRIC, or CUBIC-MEAN).
   */
  public String getSigmaCombiningRule() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_AmoebaVdwForce_getSigmaCombiningRule(getPointer()));
  }

  /**
   * Get the softcore alpha.
   *
   * @return softcore alpha.
   */
  public double getSoftcoreAlpha() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getSoftcoreAlpha(getPointer());
  }

  /**
   * Get the softcore power.
   *
   * @return softcore power.
   */
  public int getSoftcorePower() {
    return OpenMMNative.OpenMM_AmoebaVdwForce_getSoftcorePower(getPointer());
  }

  /**
   * Get the parameters for a particle type.
   *
   * @param typeIndex particle type index.
   * @return copied particle-type parameters.
   */
  public TypeParameters getParticleTypeParameters(int typeIndex) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment sigma = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_AmoebaVdwForce_getParticleTypeParameters(getPointer(), typeIndex,
          sigma, epsilon);
      return new TypeParameters(sigma.get(ValueLayout.JAVA_DOUBLE, 0),
          epsilon.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get parameters for a type-pair override.
   *
   * @param pairIndex type-pair index.
   * @return copied type-pair parameters.
   */
  public TypePairParameters getTypePairParameters(int pairIndex) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment type1 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment type2 = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment sigma = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilon = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_AmoebaVdwForce_getTypePairParameters(getPointer(), pairIndex, type1,
          type2, sigma, epsilon);
      return new TypePairParameters(type1.get(ValueLayout.JAVA_INT, 0),
          type2.get(ValueLayout.JAVA_INT, 0), sigma.get(ValueLayout.JAVA_DOUBLE, 0),
          epsilon.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get whether the long-range dispersion correction is enabled.
   *
   * @return true if the correction is enabled.
   */
  public boolean getUseDispersionCorrection() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_AmoebaVdwForce_getUseDispersionCorrection(getPointer()));
  }

  /**
   * Get whether particle parameters are specified by particle type.
   *
   * @return true when particle types are used.
   */
  public boolean getUseParticleTypes() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_AmoebaVdwForce_getUseParticleTypes(getPointer()));
  }

  /**
   * Set the alchemical method (None 0, Decouple 1, Annihilate 2).
   *
   * @param method native alchemical method value.
   */
  public void setAlchemicalMethod(int method) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setAlchemicalMethod(getPointer(), method);
  }

  /**
   * Set the cutoff distance.
   *
   * @param cutoff cutoff distance in nanometers.
   * @deprecated Use {@link #setCutoffDistance(double)}.
   */
  @Deprecated
  public void setCutoff(double cutoff) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setCutoff(getPointer(), cutoff);
  }

  /**
   * Set the nonbonded cutoff distance.
   *
   * @param distance cutoff distance in nanometers.
   */
  public void setCutoffDistance(double distance) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setCutoffDistance(getPointer(), distance);
  }

  /**
   * Set the epsilon combining rule.
   *
   * @param rule native combining-rule string: ARITHMETIC, GEOMETRIC, HARMONIC, W-H, or HHG.
   */
  public void setEpsilonCombiningRule(String rule) {
    OpenMMStrings.withUtf8String(rule, value ->
        OpenMMNative.OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(getPointer(), value));
  }

  /**
   * Set the name of the context parameter controlling alchemical interactions.
   *
   * @param name lambda parameter name.
   */
  public void setLambdaName(String name) {
    OpenMMStrings.withUtf8String(name, value ->
        OpenMMNative.OpenMM_AmoebaVdwForce_setLambdaName(getPointer(), value));
  }

  /**
   * Set the nonbonded method (NoCutoff 0 or CutoffPeriodic 1).
   *
   * @param method native nonbonded method value.
   */
  public void setNonbondedMethod(int method) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setNonbondedMethod(getPointer(), method);
  }

  /**
   * Replace the exclusions for one particle.
   *
   * @param particleIndex particle index.
   * @param exclusions excluded particle indices.
   * @throws NullPointerException if {@code exclusions} is null.
   */
  public void setParticleExclusions(int particleIndex, IntArray exclusions) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setParticleExclusions(getPointer(), particleIndex,
        Objects.requireNonNull(exclusions, "Exclusions cannot be null.").getPointer());
  }

  /**
   * Replace the parameters for a particle.
   *
   * @param particleIndex particle index.
   * @param parentIndex parent particle index.
   * @param sigma van der Waals sigma in nanometers.
   * @param epsilon van der Waals well depth in kilojoules per mole.
   * @param reductionFactor dimensionless interaction-site reduction factor.
   * @param isAlchemical whether this particle is alchemical.
   * @param typeIndex particle type, or -1 when specified per particle.
   * @param scaleFactor dimensionless scale applied to interactions involving this particle.
   */
  public void setParticleParameters(int particleIndex, int parentIndex, double sigma, double epsilon,
                                    double reductionFactor, boolean isAlchemical, int typeIndex,
                                    double scaleFactor) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setParticleParameters(getPointer(), particleIndex,
        parentIndex, sigma, epsilon, reductionFactor, OpenMMBooleans.toNative(isAlchemical),
        typeIndex, scaleFactor);
  }

  /**
   * Replace particle parameters using the legacy integer boolean representation.
   *
   * @param particleIndex particle index.
   * @param parentIndex parent particle index.
   * @param sigma van der Waals sigma in nanometers.
   * @param epsilon van der Waals well depth in kilojoules per mole.
   * @param reductionFactor dimensionless interaction-site reduction factor.
   * @param isAlchemical nonzero if this particle is alchemical.
   * @param typeIndex particle type index, or -1 when parameters are specified per particle.
   * @param scaleFactor dimensionless scale applied to interactions involving this particle.
   */
  public void setParticleParameters(int particleIndex, int parentIndex, double sigma, double epsilon,
                                    double reductionFactor, int isAlchemical, int typeIndex,
                                    double scaleFactor) {
    setParticleParameters(particleIndex, parentIndex, sigma, epsilon, reductionFactor,
        isAlchemical != 0, typeIndex, scaleFactor);
  }

  /**
   * Replace parameters for a particle type.
   *
   * @param typeIndex particle type index.
   * @param sigma type sigma in nanometers.
   * @param epsilon type well depth in kilojoules per mole.
   */
  public void setParticleTypeParameters(int typeIndex, double sigma, double epsilon) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setParticleTypeParameters(getPointer(), typeIndex,
        sigma, epsilon);
  }

  /**
   * Set the potential function (Buffered147 0 or LennardJones 1).
   *
   * @param function native potential-function value.
   */
  public void setPotentialFunction(int function) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setPotentialFunction(getPointer(), function);
  }

  /**
   * Set the sigma combining rule.
   *
   * @param rule native combining-rule string: ARITHMETIC, GEOMETRIC, or CUBIC-MEAN.
   */
  public void setSigmaCombiningRule(String rule) {
    OpenMMStrings.withUtf8String(rule, value ->
        OpenMMNative.OpenMM_AmoebaVdwForce_setSigmaCombiningRule(getPointer(), value));
  }

  /**
   * Set the softcore alpha.
   *
   * @param alpha softcore alpha.
   */
  public void setSoftcoreAlpha(double alpha) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setSoftcoreAlpha(getPointer(), alpha);
  }

  /**
   * Set the softcore power.
   *
   * @param power softcore power.
   */
  public void setSoftcorePower(int power) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setSoftcorePower(getPointer(), power);
  }

  /**
   * Replace parameters for a type-pair override.
   *
   * @param pairIndex type-pair index.
   * @param type1 first particle type.
   * @param type2 second particle type.
   * @param sigma pair sigma in nanometers.
   * @param epsilon pair well depth in kilojoules per mole.
   */
  public void setTypePairParameters(int pairIndex, int type1, int type2, double sigma,
                                    double epsilon) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setTypePairParameters(getPointer(), pairIndex, type1,
        type2, sigma, epsilon);
  }

  /**
   * Set whether to add the long-range dispersion correction.
   *
   * @param useCorrection true to enable the correction.
   */
  public void setUseDispersionCorrection(boolean useCorrection) {
    OpenMMNative.OpenMM_AmoebaVdwForce_setUseDispersionCorrection(getPointer(),
        OpenMMBooleans.toNative(useCorrection));
  }

  /**
   * Copy per-particle parameters into an existing Context. Changes to the nonbonded method,
   * cutoff, potential function, combining rules, alchemical settings, and other force-wide
   * settings require Context reinitialization.
   *
   * @param context context to update.
   * @throws NullPointerException if {@code context} is null.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_AmoebaVdwForce_updateParametersInContext(getPointer(),
        Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Determine whether this force uses periodic boundary conditions.
   *
   * @return true for CutoffPeriodic.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_AmoebaVdwForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Per-particle AMOEBA van der Waals parameters.
   *
   * @param parentIndex parent particle index.
   * @param sigma van der Waals sigma in nanometers.
   * @param epsilon well depth in kilojoules per mole.
   * @param reductionFactor dimensionless fraction locating the interaction site from the parent.
   * @param isAlchemical whether the particle participates in alchemical changes.
   * @param typeIndex particle type index, or -1 when per-particle parameters are used.
   * @param scaleFactor dimensionless scale on interactions involving this particle.
   */
  public record ParticleParameters(int parentIndex, double sigma, double epsilon,
                                   double reductionFactor, boolean isAlchemical, int typeIndex,
                                   double scaleFactor) {
  }

  /**
   * Parameters for one particle type.
   *
   * @param sigma type sigma in nanometers.
   * @param epsilon type well depth in kilojoules per mole.
   */
  public record TypeParameters(double sigma, double epsilon) {
  }

  /**
   * Parameters for one pair-specific type override.
   *
   * @param type1 first particle type index.
   * @param type2 second particle type index.
   * @param sigma pair sigma in nanometers.
   * @param epsilon pair well depth in kilojoules per mole.
   */
  public record TypePairParameters(int type1, int type2, double sigma, double epsilon) {
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_AmoebaVdwForce_create.makeInvoker().apply();
  }
}
