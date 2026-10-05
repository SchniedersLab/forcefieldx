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
import java.util.Objects;

/**
 * An implicit solvation force using the Generalized Born/surface-area (GBSA) OBC model.
 *
 * <p>Create the force, call {@link #addParticle(double, double, double)} once for each particle of the {@link
 * System}, then add it to the system. The number of particles with GBSA parameters must equal the number of
 * system particles or OpenMM throws an exception when a {@link Context} is created. Parameter changes made with
 * {@link #setParticleParameters(int, double, double, double)} do not affect existing contexts until {@link
 * #updateParametersInContext(Context)} is called.</p>
 *
 * <p>The system should also contain a {@link NonbondedForce} with charges identical to those given here; otherwise
 * the results are incorrect. If that force uses a cutoff-nonperiodic or cutoff-periodic method, call {@link
 * NonbondedForce#setReactionFieldDielectric(double)} with 1.0 to turn off the reaction-field approximation, which
 * is not correct when combined with GBSA. The default nonbonded method is {@link NonbondedMethod#NO_CUTOFF}.</p>
 *
 * <p>Context-update scope: {@link #updateParametersInContext(Context)} updates only per-particle parameters. The
 * nonbonded method, cutoff distance, dielectrics, surface-area energy and other settings need a new context, and
 * particles cannot be added.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.GBSAOBCForce} is not. The JNA
 * nonbonded-method getter and setter use raw integers; this class adds the {@link NonbondedMethod} enum and keeps
 * integer variants ({@link #getNonbondedMethodValue()}, {@link #setNonbondedMethod(int)}). {@link
 * #getParticleParameters(int)} returns a record, whereas the JNA version returns void and fills
 * {@code DoubleByReference}/{@code DoubleBuffer} arguments. The OpenMM header documents no defaults or units for
 * the solute and solvent dielectric constants.</p>
 */
public class GBSAOBCForce extends Force {
  /**
   * Long-range nonbonded methods supported by GBSA-OBC. Unlike {@link NonbondedForce.NonbondedMethod}, there are no
   * Ewald, PME or LJ-PME values. Each constant corresponds to the native integer in parentheses.
   *
   * <p>Conversion from the native integer used by {@link GBSAOBCForce#getNonbondedMethod()} throws {@link
   * IllegalStateException} for any value other than 0 through 2.</p>
   */
  public enum NonbondedMethod {
    /**
     * Native value 0. No cutoff is applied and all N^2 interactions are computed exactly, so periodic boundary
     * conditions cannot be used. This is the default.
     */
    NO_CUTOFF(0),
    /**
     * Native value 1 ({@code CutoffNonPeriodic}). Interactions beyond the cutoff distance are ignored.
     */
    CUTOFF_NONPERIODIC(1),
    /**
     * Native value 2 ({@code CutoffPeriodic}). Periodic boundary conditions are used, so each particle interacts
     * only with the nearest periodic copy of each other particle, and interactions beyond the cutoff are ignored.
     */
    CUTOFF_PERIODIC(2);
    private final int nativeValue;

    NonbondedMethod(int nativeValue) {
      this.nativeValue = nativeValue;
    }

    private static NonbondedMethod fromNative(int value) {
      for (NonbondedMethod method : values()) if (method.nativeValue == value) return method;
      throw new IllegalStateException("Unknown OpenMM GBSA-OBC method: " + value);
    }
  }

  /**
   * Create a native force. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * <p>The JNA counterpart calls the native create function directly without this step.</p>
   */
  public GBSAOBCForce() {
    super(create());
  }

  /**
   * Add the GBSA parameters of a particle. Call once for each particle of the system; the i'th call defines the
   * i'th particle.
   *
   * @param charge        particle charge, in units of the proton charge; must match the charge in the
   *     {@link NonbondedForce}.
   * @param radius        GBSA radius of the particle, in nm.
   * @param scalingFactor OBC scaling factor for the particle (dimensionless).
   * @return index of the particle that was added.
   */
  public int addParticle(double charge, double radius, double scalingFactor) {
    return OpenMMNative.OpenMM_GBSAOBCForce_addParticle(getPointer(), charge, radius, scalingFactor);
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_GBSAOBCForce_destroy);
  }

  /**
   * Get the cutoff distance used for nonbonded interactions.
   *
   * @return cutoff distance, in nm; it has no effect for {@link NonbondedMethod#NO_CUTOFF}.
   */
  public double getCutoffDistance() {
    return OpenMMNative.OpenMM_GBSAOBCForce_getCutoffDistance(getPointer());
  }

  /**
   * Get the method used for long-range nonbonded interactions. The JNA counterpart returns the raw native
   * integer; see {@link #getNonbondedMethodValue()} for that representation.
   *
   * @return the method, mapped from OpenMM's native integer.
   * @throws IllegalStateException if the native value is not 0, 1 or 2.
   */
  public NonbondedMethod getNonbondedMethod() {
    return NonbondedMethod.fromNative(OpenMMNative.OpenMM_GBSAOBCForce_getNonbondedMethod(getPointer()));
  }

  /**
   * Get the nonbonded method as OpenMM's native integer (0 NoCutoff, 1 CutoffNonPeriodic, 2 CutoffPeriodic),
   * as the JNA {@code getNonbondedMethod()} does.
   *
   * @return native nonbonded method value, unchecked.
   */
  public int getNonbondedMethodValue() {
    return OpenMMNative.OpenMM_GBSAOBCForce_getNonbondedMethod(getPointer());
  }

  /**
   * Get the number of particles for which GBSA parameters have been defined.
   *
   * @return number of particles.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_GBSAOBCForce_getNumParticles(getPointer());
  }

  /**
   * Get the GBSA parameters of a particle.
   *
   * <p>Native out-parameters are copied into the returned record; the JNA counterpart returns void and fills
   * reference/buffer arguments.</p>
   *
   * @param index index of the particle.
   * @return a {@link ParticleParameters} record.
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment q = a.allocate(ValueLayout.JAVA_DOUBLE), r = a.allocate(ValueLayout.JAVA_DOUBLE), s = a.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_GBSAOBCForce_getParticleParameters(getPointer(), index, q, r, s);
      return new ParticleParameters(q.get(ValueLayout.JAVA_DOUBLE, 0), r.get(ValueLayout.JAVA_DOUBLE, 0), s.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the dielectric constant of the solute.
   *
   * @return solute dielectric constant (dimensionless).
   */
  public double getSoluteDielectric() {
    return OpenMMNative.OpenMM_GBSAOBCForce_getSoluteDielectric(getPointer());
  }

  /**
   * Get the dielectric constant of the solvent.
   *
   * @return solvent dielectric constant (dimensionless).
   */
  public double getSolventDielectric() {
    return OpenMMNative.OpenMM_GBSAOBCForce_getSolventDielectric(getPointer());
  }

  /**
   * Get the energy scale of the surface energy term.
   *
   * @return surface-area energy scale, in kJ/(mol nm^2).
   */
  public double getSurfaceAreaEnergy() {
    return OpenMMNative.OpenMM_GBSAOBCForce_getSurfaceAreaEnergy(getPointer());
  }

  /**
   * Set the cutoff distance used for nonbonded interactions; it has no effect for {@link
   * NonbondedMethod#NO_CUTOFF}. Not updated in existing contexts.
   *
   * @param distance cutoff distance, in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMMNative.OpenMM_GBSAOBCForce_setCutoffDistance(getPointer(), distance);
  }

  /**
   * Set the method used for long-range nonbonded interactions. Not updated in existing contexts.
   *
   * @param method nonbonded method; must not be null.
   * @throws NullPointerException if {@code method} is null.
   */
  public void setNonbondedMethod(NonbondedMethod method) {
    OpenMMNative.OpenMM_GBSAOBCForce_setNonbondedMethod(getPointer(), Objects.requireNonNull(method, "Method cannot be null.").nativeValue);
  }

  /**
   * Set the nonbonded method using OpenMM's native integer (0 NoCutoff, 1 CutoffNonPeriodic, 2 CutoffPeriodic),
   * matching the JNA signature. The value is passed to OpenMM without checking on the Java side.
   *
   * @param method native nonbonded method value.
   */
  public void setNonbondedMethod(int method) {
    OpenMMNative.OpenMM_GBSAOBCForce_setNonbondedMethod(getPointer(), method);
  }

  /**
   * Replace the GBSA parameters of an existing particle. Existing contexts see the change only after {@link
   * #updateParametersInContext(Context)}.
   *
   * @param index         index of the particle.
   * @param charge        particle charge, in units of the proton charge.
   * @param radius        GBSA radius, in nm.
   * @param scalingFactor OBC scaling factor (dimensionless).
   */
  public void setParticleParameters(int index, double charge, double radius, double scalingFactor) {
    OpenMMNative.OpenMM_GBSAOBCForce_setParticleParameters(getPointer(), index, charge, radius, scalingFactor);
  }

  /**
   * Set the dielectric constant of the solute. Not updated in existing contexts.
   *
   * @param dielectric solute dielectric constant (dimensionless).
   */
  public void setSoluteDielectric(double dielectric) {
    OpenMMNative.OpenMM_GBSAOBCForce_setSoluteDielectric(getPointer(), dielectric);
  }

  /**
   * Set the dielectric constant of the solvent. Not updated in existing contexts.
   *
   * @param dielectric solvent dielectric constant (dimensionless).
   */
  public void setSolventDielectric(double dielectric) {
    OpenMMNative.OpenMM_GBSAOBCForce_setSolventDielectric(getPointer(), dielectric);
  }

  /**
   * Set the energy scale of the surface energy term. Not updated in existing contexts.
   *
   * @param energy surface-area energy scale, in kJ/(mol nm^2).
   */
  public void setSurfaceAreaEnergy(double energy) {
    OpenMMNative.OpenMM_GBSAOBCForce_setSurfaceAreaEnergy(getPointer(), energy);
  }

  /**
   * Copy the per-particle parameters stored in this force into an existing {@link Context} without
   * reinitializing it. Only particle parameters are updated; see the class description for what cannot change.
   *
   * @param context live context created from a system containing this force; must not be null.
   * @throws NullPointerException if {@code context} is null.
   * @throws IllegalStateException if this force or the context has been destroyed.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_GBSAOBCForce_updateParametersInContext(getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Determine whether this force uses periodic boundary conditions.
   *
   * <p>Overrides {@link Force#usesPeriodicBoundaryConditions()} with the native GBSA-OBC query; this class has
   * no setter, and the OpenMM header states only that it reports whether the force uses periodic boundary
   * conditions.</p>
   *
   * @return true if the force uses periodic boundary conditions, false otherwise.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_GBSAOBCForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Immutable copy of the GBSA parameters of one particle.
   *
   * @param charge        particle charge, in units of the proton charge.
   * @param radius        GBSA radius, in nm.
   * @param scalingFactor OBC scaling factor (dimensionless).
   */
  public record ParticleParameters(double charge, double radius, double scalingFactor) {
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_GBSAOBCForce_create.makeInvoker().apply();
  }
}
