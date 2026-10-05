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
 * Nonbonded interactions between particles: a Coulomb force for electrostatics and a Lennard-Jones force for
 * van der Waals interactions, with optional cutoffs and periodic methods.
 *
 * <p>Create the force, call {@link #addParticle(double, double, double)} once for each particle in the {@link
 * System}, then add the force to the system. The number of particles with nonbonded parameters must equal the
 * number of particles in the system, or OpenMM throws an exception when a {@link Context} is created. Lennard-Jones
 * pair parameters use the Lorentz-Berthelot rule (arithmetic mean of sigmas, geometric mean of epsilons).</p>
 *
 * <p>Exceptions replace the ordinary interaction for selected pairs, either excluding the pair entirely or giving it
 * modified (for example 1-4) parameters; {@link #createExceptionsFromBonds(BondArray, double, double)} builds them
 * from a bond list. Cutoffs are never applied to exceptions.</p>
 *
 * <p>With a cutoff, Lennard-Jones interactions are sharply truncated by default; enable a switching function with
 * {@link #setUseSwitchingFunction(boolean)} and a {@link #setSwitchingDistance(double) switching distance} smaller
 * than the cutoff to taper them smoothly. A long-range dispersion correction is enabled by default. The default
 * {@link NonbondedMethod} is {@link NonbondedMethod#NO_CUTOFF}, and the PME/LJ-PME separation parameter defaults to 0,
 * meaning that the grid parameters are derived from the Ewald error tolerance.</p>
 *
 * <p>Global parameters ({@link #addGlobalParameter(String, double)}) together with particle and exception parameter
 * offsets modify the effective values as {@code charge = baseCharge + param*chargeScale}, and likewise for sigma and
 * epsilon (and for the charge product, sigma and epsilon of exceptions).</p>
 *
 * <p>Context-update scope: {@link #updateParametersInContext(Context)} copies only per-particle parameters
 * (charge, sigma, epsilon), per-exception charge product, sigma and epsilon, and the charge/sigma/epsilon scales of
 * offsets. The nonbonded method, cutoff, switching, PME and other settings can be changed only by recreating the
 * context, and neither particles/exceptions nor the particle, exception or global parameter a offset refers to
 * can be changed or added.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.NonbondedForce} is not. The JNA class
 * exposes only a subset of the API, with getters filling {@code IntByReference}/{@code DoubleByReference} arguments,
 * integer {@code setNonbondedMethod}/{@code setUseSwitchingFunction}/{@code setUseDispersionCorrection}, and
 * {@code ffx.openmm.BondArray}. This class returns records and the {@link NonbondedMethod} enum, adds boolean
 * overloads (retaining the integer overloads), and adds the nonbonded-method, switching, reaction-field, Ewald, PME
 * and LJ-PME getters, global parameters, offsets, reciprocal-space force group, direct-space, and exception
 * periodic-boundary methods, none of which exist in the JNA class.</p>
 */
public class NonbondedForce extends Force {
  /**
   * The methods OpenMM can use to handle long-range nonbonded interactions. Each constant corresponds to the
   * native integer given in parentheses, which is the value used by {@link #setNonbondedMethod(int)}.
   *
   * <p>{@link #fromNative(int)} rejects any value other than 0 through 5 with an {@link IllegalStateException};
   * that conversion is used by {@link #getNonbondedMethod()}.</p>
   */
  public enum NonbondedMethod {
    /**
     * Native value 0. No cutoff is applied and all N^2 interactions are computed exactly, so periodic
     * boundary conditions cannot be used. This is the OpenMM default.
     */
    NO_CUTOFF(0),
    /**
     * Native value 1 ({@code CutoffNonPeriodic}). Interactions beyond the cutoff distance are ignored and Coulomb
     * interactions closer than the cutoff use the reaction-field method.
     */
    CUTOFF_NONPERIODIC(1),
    /**
     * Native value 2 ({@code CutoffPeriodic}). Periodic boundary conditions are used, each particle interacts only
     * with the nearest periodic copy of each other particle, interactions beyond the cutoff are ignored, and Coulomb
     * interactions inside the cutoff use the reaction-field method.
     */
    CUTOFF_PERIODIC(2),
    /**
     * Native value 3. Periodic boundary conditions are used and Ewald summation computes the Coulomb interaction
     * with all periodic copies.
     */
    EWALD(3),
    /**
     * Native value 4. Periodic boundary conditions are used and Particle-Mesh Ewald summation computes the Coulomb
     * interaction with all periodic copies.
     */
    PME(4),
    /**
     * Native value 5. Periodic boundary conditions are used and Particle-Mesh Ewald summation is used for both the
     * Coulomb and the Lennard-Jones interactions; no switching is used for either.
     */
    LJPME(5);
    private final int nativeValue;

    NonbondedMethod(int nativeValue) {
      this.nativeValue = nativeValue;
    }

    private int nativeValue() {
      return nativeValue;
    }

    private static NonbondedMethod fromNative(int value) {
      for (NonbondedMethod method : values()) if (method.nativeValue == value) return method;
      throw new IllegalStateException("Unknown OpenMM nonbonded method: " + value);
    }
  }

  /**
   * Create a native force. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * <p>The JNA counterpart calls the native create function directly without this step.</p>
   */
  public NonbondedForce() {
    super(create());
  }

  /**
   * Add an exception: a pair interaction calculated with parameters different from the per-particle ones.
   *
   * <p>If the charge product and epsilon are both 0 the interaction is omitted from force and energy calculations.
   * Cutoffs are never applied to exceptions. {@link #createExceptionsFromBonds(BondArray, double, double)} is often
   * more convenient.</p>
   *
   * @param particle1     index of the first particle in the interaction.
   * @param particle2     index of the second particle in the interaction.
   * @param chargeProduct scaled product of the particle charges (strength of the Coulomb interaction), in
   *     proton-charge units squared.
   * @param sigma         Lennard-Jones sigma, in nm.
   * @param epsilon       Lennard-Jones epsilon (well depth), in kJ/mol.
   * @param replace       if true, an existing exception for the same pair is replaced; if false, OpenMM raises an
   *     error in that case. The OpenMM C++ default of {@code false} is not applied by this method; the caller
   *     must supply the value.
   * @return index of the exception that was added.
   */
  public int addException(int particle1, int particle2, double chargeProduct, double sigma,
                          double epsilon, boolean replace) {
    return OpenMMNative.OpenMM_NonbondedForce_addException(getPointer(), particle1, particle2,
        chargeProduct, sigma, epsilon, OpenMMBooleans.toNative(replace));
  }

  /**
   * Add the nonbonded parameters of a particle. Call once for each particle of the system; the i'th call defines
   * the i'th particle.
   *
   * @param charge  particle charge, in units of the proton charge.
   * @param sigma   Lennard-Jones sigma (van der Waals radius), in nm.
   * @param epsilon Lennard-Jones epsilon (well depth), in kJ/mol.
   * @return index of the particle that was added.
   */
  public int addParticle(double charge, double sigma, double epsilon) {
    return OpenMMNative.OpenMM_NonbondedForce_addParticle(getPointer(), charge, sigma, epsilon);
  }

  /**
   * Create exceptions from molecular topology. Pairs separated by one or two bonds are set not to interact, and
   * pairs separated by three bonds (1-4 interactions) have their Coulomb and Lennard-Jones strengths multiplied by
   * the given factors.
   *
   * @param bonds          FFM {@link BondArray} of particle-index pairs that are bonded; must not be null.
   *     The array is only read, and remains owned by the caller.
   * @param coulomb14Scale factor multiplying the Coulomb interaction of 1-4 pairs.
   * @param lj14Scale      factor multiplying the Lennard-Jones interaction of 1-4 pairs.
   * @throws NullPointerException if {@code bonds} is null.
   */
  public void createExceptionsFromBonds(BondArray bonds, double coulomb14Scale, double lj14Scale) {
    OpenMMNative.OpenMM_NonbondedForce_createExceptionsFromBonds(getPointer(),
        Objects.requireNonNull(bonds, "Bonds cannot be null.").getPointer(),
        coulomb14Scale, lj14Scale);
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_NonbondedForce_destroy);
  }

  /**
   * Get the cutoff distance used for nonbonded interactions.
   *
   * @return cutoff distance, in nm; it has no effect for {@link NonbondedMethod#NO_CUTOFF}.
   */
  public double getCutoffDistance() {
    return OpenMMNative.OpenMM_NonbondedForce_getCutoffDistance(getPointer());
  }

  /**
   * Get the parameters of an exception.
   *
   * <p>Native out-parameters are copied into the returned record. The JNA counterpart returns void and fills
   * {@code IntByReference}/{@code DoubleByReference} arguments.</p>
   *
   * @param index index of the exception.
   * @return an {@link ExceptionParameters} record.
   */
  public ExceptionParameters getExceptionParameters(int index) {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment p1 = a.allocate(ValueLayout.JAVA_INT), p2 = a.allocate(ValueLayout.JAVA_INT);
      MemorySegment q = a.allocate(ValueLayout.JAVA_DOUBLE), s = a.allocate(ValueLayout.JAVA_DOUBLE), e = a.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_NonbondedForce_getExceptionParameters(getPointer(), index, p1, p2, q, s, e);
      return new ExceptionParameters(p1.get(ValueLayout.JAVA_INT, 0), p2.get(ValueLayout.JAVA_INT, 0),
          q.get(ValueLayout.JAVA_DOUBLE, 0), s.get(ValueLayout.JAVA_DOUBLE, 0), e.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the Ewald error tolerance, the acceptable fractional force error used to select the reciprocal-space
   * cutoff and separation parameter. It gives no rigorous per-atom guarantee. For PME, it is ignored if the PME
   * alpha has been set to a nonzero value.
   *
   * @return Ewald error tolerance (dimensionless fractional force error).
   */
  public double getEwaldErrorTolerance() {
    return OpenMMNative.OpenMM_NonbondedForce_getEwaldErrorTolerance(getPointer());
  }

  /**
   * Get the number of special interactions calculated differently from the ordinary ones.
   *
   * @return number of exceptions.
   */
  public int getNumExceptions() {
    return OpenMMNative.OpenMM_NonbondedForce_getNumExceptions(getPointer());
  }

  /**
   * Get the number of particles for which nonbonded parameters have been defined.
   *
   * @return number of particles with parameters.
   */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_NonbondedForce_getNumParticles(getPointer());
  }

  /**
   * Get the method used for long-range nonbonded interactions.
   *
   * @return the method, mapped from OpenMM's native integer.
   * @throws IllegalStateException if the native value is not one of 0 through 5.
   */
  public NonbondedMethod getNonbondedMethod() {
    return NonbondedMethod.fromNative(OpenMMNative.OpenMM_NonbondedForce_getNonbondedMethod(getPointer()));
  }

  /**
   * Get the PME parameters configured on this force. If alpha is 0 (the default) the stored grid values are
   * ignored and values are chosen from the Ewald error tolerance. For the values actually used by a context, see
   * {@link #getPMEParametersInContext(Context)}.
   *
   * @return copied {@link PMEParameters}.
   */
  public PMEParameters getPMEParameters() {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment alpha = a.allocate(ValueLayout.JAVA_DOUBLE), nx = a.allocate(ValueLayout.JAVA_INT), ny = a.allocate(ValueLayout.JAVA_INT), nz = a.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_NonbondedForce_getPMEParameters(getPointer(), alpha, nx, ny, nz);
      return new PMEParameters(alpha.get(ValueLayout.JAVA_DOUBLE, 0), nx.get(ValueLayout.JAVA_INT, 0), ny.get(ValueLayout.JAVA_INT, 0), nz.get(ValueLayout.JAVA_INT, 0));
    }
  }

  /**
   * Get the nonbonded parameters of a particle.
   *
   * <p>Native out-parameters are copied into the returned record. The JNA counterpart returns void and fills
   * {@code DoubleByReference} arguments.</p>
   *
   * @param index index of the particle.
   * @return a {@link ParticleParameters} record.
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment q = a.allocate(ValueLayout.JAVA_DOUBLE), s = a.allocate(ValueLayout.JAVA_DOUBLE), e = a.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_NonbondedForce_getParticleParameters(getPointer(), index, q, s, e);
      return new ParticleParameters(q.get(ValueLayout.JAVA_DOUBLE, 0), s.get(ValueLayout.JAVA_DOUBLE, 0), e.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Get the solvent dielectric constant used by the reaction-field approximation.
   *
   * @return reaction-field solvent dielectric constant (dimensionless).
   */
  public double getReactionFieldDielectric() {
    return OpenMMNative.OpenMM_NonbondedForce_getReactionFieldDielectric(getPointer());
  }

  /**
   * Get the distance at which the Lennard-Jones switching function begins to reduce the interaction. It must be
   * less than the cutoff distance.
   *
   * @return switching distance, in nm.
   */
  public double getSwitchingDistance() {
    return OpenMMNative.OpenMM_NonbondedForce_getSwitchingDistance(getPointer());
  }

  /**
   * Get whether a switching function is applied to the Lennard-Jones interaction. The option is ignored when the
   * method is {@link NonbondedMethod#NO_CUTOFF}.
   *
   * @return true if the switching function is enabled.
   */
  public boolean getUseSwitchingFunction() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_NonbondedForce_getUseSwitchingFunction(getPointer()));
  }

  /**
   * Get whether a long-range dispersion correction is added to the energy to approximate Lennard-Jones
   * interactions beyond the cutoff. It depends on the periodic box volume and applies only with periodic boundary
   * conditions. It is enabled by default.
   *
   * @return true if the dispersion correction is enabled.
   */
  public boolean getUseDispersionCorrection() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_NonbondedForce_getUseDispersionCorrection(getPointer()));
  }

  /**
   * Set the cutoff distance used for nonbonded interactions. It has no effect for {@link NonbondedMethod#NO_CUTOFF}.
   * Not updated in existing contexts by {@link #updateParametersInContext(Context)}.
   *
   * @param distance cutoff distance, in nm.
   */
  public void setCutoffDistance(double distance) {
    OpenMMNative.OpenMM_NonbondedForce_setCutoffDistance(getPointer(), distance);
  }

  /**
   * Replace the parameters of an existing exception. Existing contexts see only the charge product, sigma and
   * epsilon changes, after {@link #updateParametersInContext(Context)}; the particle pair cannot change in a context.
   * If the charge product and epsilon are both 0, the interaction is omitted. Cutoffs are never applied to exceptions.
   *
   * @param index          index of the exception.
   * @param particle1      index of the first particle in the interaction.
   * @param particle2      index of the second particle in the interaction.
   * @param chargeProduct scaled product of the particle charges, in proton-charge units squared.
   * @param sigma          Lennard-Jones sigma, in nm.
   * @param epsilon        Lennard-Jones epsilon, in kJ/mol.
   */
  public void setExceptionParameters(int index, int particle1, int particle2, double chargeProduct, double sigma, double epsilon) {
    OpenMMNative.OpenMM_NonbondedForce_setExceptionParameters(getPointer(), index, particle1, particle2, chargeProduct, sigma, epsilon);
  }

  /**
   * Set the Ewald error tolerance. For PME it is ignored if the PME alpha has been set to a nonzero value.
   * Not updated in existing contexts.
   *
   * @param tolerance acceptable fractional force error used to select Ewald/PME parameters (dimensionless).
   */
  public void setEwaldErrorTolerance(double tolerance) {
    OpenMMNative.OpenMM_NonbondedForce_setEwaldErrorTolerance(getPointer(), tolerance);
  }

  /**
   * Set the method used for long-range nonbonded interactions. Not updated in existing contexts.
   *
   * @param method nonbonded method; must not be null.
   * @throws NullPointerException if {@code method} is null.
   */
  public void setNonbondedMethod(NonbondedMethod method) {
    OpenMMNative.OpenMM_NonbondedForce_setNonbondedMethod(getPointer(), Objects.requireNonNull(method, "Method cannot be null.").nativeValue());
  }

  /**
   * Set the nonbonded method by its OpenMM integer value (0 NoCutoff, 1 CutoffNonPeriodic, 2 CutoffPeriodic,
   * 3 Ewald, 4 PME, 5 LJPME), matching the JNA signature. The value is passed to OpenMM without checking on the
   * Java side.
   *
   * @param method native nonbonded method value.
   */
  public void setNonbondedMethod(int method) {
    OpenMMNative.OpenMM_NonbondedForce_setNonbondedMethod(getPointer(), method);
  }

  /**
   * Set the PME separation parameter and grid dimensions. If alpha is 0 (the default) these values are ignored and
   * chosen from the Ewald error tolerance. Not updated in existing contexts.
   *
   * @param alpha separation parameter.
   * @param nx    number of grid points along the X axis.
   * @param ny    number of grid points along the Y axis.
   * @param nz    number of grid points along the Z axis.
   */
  public void setPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMMNative.OpenMM_NonbondedForce_setPMEParameters(getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Set the separation parameter and grid dimensions of the Lennard-Jones dispersion term of LJ-PME. If alpha is 0
   * (the default) these values are ignored and chosen from the Ewald error tolerance. Not updated in existing contexts.
   *
   * @param alpha separation parameter.
   * @param nx    number of dispersion grid points along the X axis.
   * @param ny    number of dispersion grid points along the Y axis.
   * @param nz    number of dispersion grid points along the Z axis.
   */
  public void setLJPMEParameters(double alpha, int nx, int ny, int nz) {
    OpenMMNative.OpenMM_NonbondedForce_setLJPMEParameters(getPointer(), alpha, nx, ny, nz);
  }

  /**
   * Replace the nonbonded parameters of an existing particle. Existing contexts see the change only after
   * {@link #updateParametersInContext(Context)}.
   *
   * @param index   index of the particle.
   * @param charge  particle charge, in units of the proton charge.
   * @param sigma   Lennard-Jones sigma, in nm.
   * @param epsilon Lennard-Jones epsilon, in kJ/mol.
   */
  public void setParticleParameters(int index, double charge, double sigma, double epsilon) {
    OpenMMNative.OpenMM_NonbondedForce_setParticleParameters(getPointer(), index, charge, sigma, epsilon);
  }

  /**
   * Set the solvent dielectric constant used by the reaction-field approximation. Not updated in existing
   * contexts.
   *
   * @param dielectric solvent dielectric constant (dimensionless).
   */
  public void setReactionFieldDielectric(double dielectric) {
    OpenMMNative.OpenMM_NonbondedForce_setReactionFieldDielectric(getPointer(), dielectric);
  }

  /**
   * Set the distance at which the Lennard-Jones switching function begins to reduce the interaction; it must be
   * less than the cutoff distance. Not updated in existing contexts.
   *
   * @param distance switching distance, in nm.
   */
  public void setSwitchingDistance(double distance) {
    OpenMMNative.OpenMM_NonbondedForce_setSwitchingDistance(getPointer(), distance);
  }

  /**
   * Set whether a switching function is applied to the Lennard-Jones interaction (ignored for
   * {@link NonbondedMethod#NO_CUTOFF}). Not updated in existing contexts.
   *
   * @param use true to enable the switching function.
   */
  public void setUseSwitchingFunction(boolean use) {
    OpenMMNative.OpenMM_NonbondedForce_setUseSwitchingFunction(getPointer(), OpenMMBooleans.toNative(use));
  }

  /**
   * Set the switching-function flag using OpenMM's integer boolean (1 true, 0 false), matching the JNA
   * signature. The value is passed through unchanged.
   *
   * @param use native boolean value.
   */
  public void setUseSwitchingFunction(int use) {
    OpenMMNative.OpenMM_NonbondedForce_setUseSwitchingFunction(getPointer(), use);
  }

  /**
   * Set whether the long-range dispersion correction is applied. Not updated in existing contexts.
   *
   * @param use true to enable the correction (the default).
   */
  public void setUseDispersionCorrection(boolean use) {
    OpenMMNative.OpenMM_NonbondedForce_setUseDispersionCorrection(getPointer(), OpenMMBooleans.toNative(use));
  }

  /**
   * Set the dispersion-correction flag using OpenMM's integer boolean (1 true, 0 false), matching the JNA
   * signature. The value is passed through unchanged.
   *
   * @param use native boolean value.
   */
  public void setUseDispersionCorrection(int use) {
    OpenMMNative.OpenMM_NonbondedForce_setUseDispersionCorrection(getPointer(), use);
  }

  /**
   * Get the number of global parameters that have been added.
   *
   * @return number of global parameters.
   */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_NonbondedForce_getNumGlobalParameters(getPointer());
  }

  /**
   * Get the number of particle parameter offsets that have been added.
   *
   * @return number of particle offsets.
   */
  public int getNumParticleParameterOffsets() {
    return OpenMMNative.OpenMM_NonbondedForce_getNumParticleParameterOffsets(getPointer());
  }

  /**
   * Get the number of exception parameter offsets that have been added.
   *
   * @return number of exception offsets.
   */
  public int getNumExceptionParameterOffsets() {
    return OpenMMNative.OpenMM_NonbondedForce_getNumExceptionParameterOffsets(getPointer());
  }

  /**
   * Add a global parameter that parameter offsets may depend on. Its default value is the initial value in newly
   * created contexts and can later be changed through the context.
   *
   * @param name         parameter name; converted to a temporary native UTF-8 string.
   * @param defaultValue default (initial) value of the parameter.
   * @return index of the added parameter.
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_NonbondedForce_addGlobalParameter(getPointer(), value, defaultValue));
  }

  /**
   * Get the name of a global parameter.
   *
   * @param index index of the parameter.
   * @return parameter name, copied into a Java string.
   */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_NonbondedForce_getGlobalParameterName(getPointer(), index));
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index index of the parameter.
   * @return default value.
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_NonbondedForce_getGlobalParameterDefaultValue(getPointer(), index);
  }

  /**
   * Set the name of a global parameter.
   *
   * @param index index of the parameter.
   * @param name  new parameter name.
   */
  public void setGlobalParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_NonbondedForce_setGlobalParameterName(getPointer(), index, value));
  }

  /**
   * Set the default value of a global parameter.
   *
   * @param index index of the parameter.
   * @param value new default value.
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_NonbondedForce_setGlobalParameterDefaultValue(getPointer(), index, value);
  }

  /**
   * Get the LJ-PME dispersion parameters configured on this force. If alpha is 0 (the default) the stored grid
   * values are ignored and chosen from the Ewald error tolerance.
   *
   * @return copied {@link PMEParameters}.
   */
  public PMEParameters getLJPMEParameters() {
    try (Arena a = Arena.ofConfined()) {
      MemorySegment alpha = a.allocate(ValueLayout.JAVA_DOUBLE), nx = a.allocate(ValueLayout.JAVA_INT), ny = a.allocate(ValueLayout.JAVA_INT), nz = a.allocate(ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_NonbondedForce_getLJPMEParameters(getPointer(), alpha, nx, ny, nz);
      return new PMEParameters(alpha.get(ValueLayout.JAVA_DOUBLE, 0), nx.get(ValueLayout.JAVA_INT, 0), ny.get(ValueLayout.JAVA_INT, 0), nz.get(ValueLayout.JAVA_INT, 0));
    }
  }

  /**
   * Get the PME parameters actually used in a context. Platforms may restrict grid sizes, so these can differ
   * slightly from the configured or Ewald-tolerance-derived values.
   *
   * @param context live context containing this force; must not be null.
   * @return copied {@link PMEParameters}.
   * @throws NullPointerException if {@code context} is null.
   */
  public PMEParameters getPMEParametersInContext(Context context) {
    return getGridParameters(context, false);
  }

  /**
   * Get the LJ-PME dispersion parameters actually used in a context. Platforms may restrict grid sizes, so
   * these can differ slightly from the configured values.
   *
   * @param context live context containing this force; must not be null.
   * @return copied {@link PMEParameters}.
   * @throws NullPointerException if {@code context} is null.
   */
  public PMEParameters getLJPMEParametersInContext(Context context) {
    return getGridParameters(context, true);
  }

  /**
   * Add an offset to the per-particle parameters of a particle, driven by a global parameter.
   *
   * @param parameter     name of a global parameter that was already added with {@link #addGlobalParameter(String, double)}.
   * @param particleIndex index of the affected particle.
   * @param chargeScale   scale; this times the parameter value is added to the particle's charge.
   * @param sigmaScale    scale; this times the parameter value is added to the particle's sigma.
   * @param epsilonScale  scale; this times the parameter value is added to the particle's epsilon.
   * @return index of the added offset.
   */
  public int addParticleParameterOffset(String parameter, int particleIndex, double chargeScale,
                                        double sigmaScale, double epsilonScale) {
    return OpenMMStrings.withUtf8StringResult(parameter, name ->
        OpenMMNative.OpenMM_NonbondedForce_addParticleParameterOffset(
            getPointer(), name, particleIndex, chargeScale, sigmaScale, epsilonScale));
  }

  /**
   * Get a particle parameter offset. The native parameter-name string is copied into the returned record.
   *
   * @param index index of the offset, as returned by {@link #addParticleParameterOffset}.
   * @return a {@link ParticleParameterOffset} record.
   */
  public ParticleParameterOffset getParticleParameterOffset(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment parameter = arena.allocate(ValueLayout.ADDRESS);
      MemorySegment particleIndex = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment chargeScale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment sigmaScale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilonScale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_NonbondedForce_getParticleParameterOffset(
          getPointer(), index, parameter, particleIndex, chargeScale, sigmaScale, epsilonScale);
      return new ParticleParameterOffset(
          OpenMMStrings.copy(parameter.get(ValueLayout.ADDRESS, 0)),
          particleIndex.get(ValueLayout.JAVA_INT, 0),
          chargeScale.get(ValueLayout.JAVA_DOUBLE, 0), sigmaScale.get(ValueLayout.JAVA_DOUBLE, 0),
          epsilonScale.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace a particle parameter offset. Existing contexts see changed scales after {@link
   * #updateParametersInContext(Context)}, but not a changed particle or global parameter.
   *
   * @param index         index of the offset.
   * @param parameter     name of an existing global parameter.
   * @param particleIndex index of the affected particle.
   * @param chargeScale   scale applied to the parameter value and added to charge.
   * @param sigmaScale    scale applied to the parameter value and added to sigma.
   * @param epsilonScale  scale applied to the parameter value and added to epsilon.
   */
  public void setParticleParameterOffset(int index, String parameter, int particleIndex,
                                         double chargeScale, double sigmaScale, double epsilonScale) {
    OpenMMStrings.withUtf8String(parameter, name ->
        OpenMMNative.OpenMM_NonbondedForce_setParticleParameterOffset(
            getPointer(), index, name, particleIndex, chargeScale, sigmaScale, epsilonScale));
  }

  /**
   * Add an offset to the parameters of an exception, driven by a global parameter.
   *
   * @param parameter         name of a global parameter that was already added.
   * @param exceptionIndex    index of the affected exception.
   * @param chargeProductScale scale; this times the parameter value is added to the exception's charge product.
   * @param sigmaScale        scale; this times the parameter value is added to the exception's sigma.
   * @param epsilonScale      scale; this times the parameter value is added to the exception's epsilon.
   * @return index of the added offset.
   */
  public int addExceptionParameterOffset(String parameter, int exceptionIndex,
                                         double chargeProductScale, double sigmaScale, double epsilonScale) {
    return OpenMMStrings.withUtf8StringResult(parameter, name ->
        OpenMMNative.OpenMM_NonbondedForce_addExceptionParameterOffset(
            getPointer(), name, exceptionIndex, chargeProductScale, sigmaScale, epsilonScale));
  }

  /**
   * Get an exception parameter offset. The native parameter-name string is copied into the returned record.
   *
   * @param index index of the offset, as returned by {@link #addExceptionParameterOffset}.
   * @return an {@link ExceptionParameterOffset} record.
   */
  public ExceptionParameterOffset getExceptionParameterOffset(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment parameter = arena.allocate(ValueLayout.ADDRESS);
      MemorySegment exceptionIndex = arena.allocate(ValueLayout.JAVA_INT);
      MemorySegment chargeProductScale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment sigmaScale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment epsilonScale = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_NonbondedForce_getExceptionParameterOffset(
          getPointer(), index, parameter, exceptionIndex, chargeProductScale, sigmaScale, epsilonScale);
      return new ExceptionParameterOffset(
          OpenMMStrings.copy(parameter.get(ValueLayout.ADDRESS, 0)),
          exceptionIndex.get(ValueLayout.JAVA_INT, 0),
          chargeProductScale.get(ValueLayout.JAVA_DOUBLE, 0),
          sigmaScale.get(ValueLayout.JAVA_DOUBLE, 0), epsilonScale.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace an exception parameter offset. Existing contexts see changed scales after {@link
   * #updateParametersInContext(Context)}, but not a changed exception or global parameter.
   *
   * @param index              index of the offset.
   * @param parameter          name of an existing global parameter.
   * @param exceptionIndex     index of the affected exception.
   * @param chargeProductScale scale applied to the parameter value and added to the charge product.
   * @param sigmaScale         scale applied to the parameter value and added to sigma.
   * @param epsilonScale       scale applied to the parameter value and added to epsilon.
   */
  public void setExceptionParameterOffset(int index, String parameter, int exceptionIndex,
                                          double chargeProductScale, double sigmaScale, double epsilonScale) {
    OpenMMStrings.withUtf8String(parameter, name ->
        OpenMMNative.OpenMM_NonbondedForce_setExceptionParameterOffset(
            getPointer(), index, name, exceptionIndex, chargeProductScale, sigmaScale, epsilonScale));
  }

  /**
   * Get the force group that reciprocal-space Ewald/PME interactions belong to. Together with {@link
   * #getForceGroup()} (direct space) it lets multiple-time-step integrators evaluate them at different intervals.
   *
   * @return group index from 0 through 31, or -1 (the default) to use the same group as direct space.
   */
  public int getReciprocalSpaceForceGroup() {
    return OpenMMNative.OpenMM_NonbondedForce_getReciprocalSpaceForceGroup(getPointer());
  }

  /**
   * Set the force group for reciprocal-space Ewald/PME interactions.
   *
   * @param group group index from 0 through 31, or -1 to use the same group as direct space.
   */
  public void setReciprocalSpaceForceGroup(int group) {
    OpenMMNative.OpenMM_NonbondedForce_setReciprocalSpaceForceGroup(getPointer(), group);
  }

  /**
   * Get whether direct-space interactions are included in forces and energies. Excluding them is useful when
   * a different force replaces the direct-space calculation while this force supplies reciprocal space.
   *
   * @return true if direct-space interactions are included.
   */
  public boolean getIncludeDirectSpace() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_NonbondedForce_getIncludeDirectSpace(getPointer()));
  }

  /**
   * Set whether direct-space interactions are included in forces and energies.
   *
   * @param include true to include direct-space interactions.
   */
  public void setIncludeDirectSpace(boolean include) {
    OpenMMNative.OpenMM_NonbondedForce_setIncludeDirectSpace(getPointer(), OpenMMBooleans.toNative(include));
  }

  /**
   * Set direct-space inclusion using OpenMM's integer boolean (1 true, 0 false); the value is passed unchanged.
   *
   * @param include native boolean value.
   */
  public void setIncludeDirectSpace(int include) {
    OpenMMNative.OpenMM_NonbondedForce_setIncludeDirectSpace(getPointer(), include);
  }

  /**
   * Get whether periodic boundary conditions are applied to exceptions. This is usually inappropriate because
   * exceptions normally represent bonded pairs. It matters only when periodic conditions are applied to other
   * interactions; it is ignored for {@link NonbondedMethod#NO_CUTOFF} and {@link NonbondedMethod#CUTOFF_NONPERIODIC}.
   * Cutoffs are never applied to exceptions.
   *
   * @return true if exceptions use periodic boundary conditions.
   */
  public boolean getExceptionsUsePeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_NonbondedForce_getExceptionsUsePeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Set whether periodic boundary conditions are applied to exceptions (see {@link
   * #getExceptionsUsePeriodicBoundaryConditions()}).
   *
   * @param periodic true to apply periodic boundary conditions to exceptions.
   */
  public void setExceptionsUsePeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_NonbondedForce_setExceptionsUsePeriodicBoundaryConditions(getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Set exception periodic-boundary handling using OpenMM's integer boolean (1 true, 0 false); the value is
   * passed unchanged.
   *
   * @param periodic native boolean value.
   */
  public void setExceptionsUsePeriodicBoundaryConditions(int periodic) {
    OpenMMNative.OpenMM_NonbondedForce_setExceptionsUsePeriodicBoundaryConditions(getPointer(), periodic);
  }

  /**
   * Copy parameters stored in this force into an existing {@link Context} without reinitializing it.
   *
   * <p>Call {@link #setParticleParameters}, {@link #setExceptionParameters}, or the offset setters first. Only
   * particle and exception parameters and the offset scales are updated; see the class description for what cannot
   * be changed. New particles and exceptions cannot be added through this method.</p>
   *
   * @param context live context created from a system containing this force; must not be null.
   * @throws NullPointerException if {@code context} is null.
   * @throws IllegalStateException if this force or the context has been destroyed.
   */
  public void updateParametersInContext(Context context) {
    OpenMMNative.OpenMM_NonbondedForce_updateParametersInContext(getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Determine whether this force uses periodic boundary conditions, as reported by OpenMM for the configured
   * {@link NonbondedMethod}.
   *
   * <p>The JNA implementation compares the native value with {@code OpenMM_True}; this one converts it through
   * {@link OpenMMBooleans#fromNative(int)}.</p>
   *
   * @return true if the force uses periodic boundary conditions.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_NonbondedForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Immutable copy of per-particle nonbonded parameters.
   *
   * @param charge  particle charge, in units of the proton charge.
   * @param sigma   Lennard-Jones sigma, in nm.
   * @param epsilon Lennard-Jones epsilon, in kJ/mol.
   */
  public record ParticleParameters(double charge, double sigma, double epsilon) {
  }

  /**
   * Immutable copy of the parameters of one exception.
   *
   * @param particle1     index of the first particle.
   * @param particle2     index of the second particle.
   * @param chargeProduct scaled product of the particle charges, in proton-charge units squared.
   * @param sigma         Lennard-Jones sigma, in nm.
   * @param epsilon       Lennard-Jones epsilon, in kJ/mol.
   */
  public record ExceptionParameters(int particle1, int particle2, double chargeProduct, double sigma, double epsilon) {
  }

  /**
   * Immutable copy of PME or LJ-PME separation and grid parameters. A zero alpha in configured parameters means
   * the values are derived from the Ewald error tolerance.
   *
   * @param alpha separation parameter.
   * @param nx    number of grid points along the X axis.
   * @param ny    number of grid points along the Y axis.
   * @param nz    number of grid points along the Z axis.
   */
  public record PMEParameters(double alpha, int nx, int ny, int nz) {
  }

  /**
   * Immutable copy of a particle parameter offset. The effective value is the base value plus the global
   * parameter value times the matching scale.
   *
   * @param parameter     name of the global parameter.
   * @param particleIndex index of the affected particle.
   * @param chargeScale   scale added to charge per unit parameter value.
   * @param sigmaScale    scale added to sigma per unit parameter value.
   * @param epsilonScale  scale added to epsilon per unit parameter value.
   */
  public record ParticleParameterOffset(String parameter, int particleIndex, double chargeScale,
                                        double sigmaScale, double epsilonScale) {
  }

  /**
   * Immutable copy of an exception parameter offset. The effective value is the base value plus the global
   * parameter value times the matching scale.
   *
   * @param parameter          name of the global parameter.
   * @param exceptionIndex     index of the affected exception.
   * @param chargeProductScale scale added to the charge product per unit parameter value.
   * @param sigmaScale         scale added to sigma per unit parameter value.
   * @param epsilonScale       scale added to epsilon per unit parameter value.
   */
  public record ExceptionParameterOffset(String parameter, int exceptionIndex,
                                         double chargeProductScale, double sigmaScale, double epsilonScale) {
  }

  private PMEParameters getGridParameters(Context context, boolean ljpme) {
    Objects.requireNonNull(context, "Context cannot be null.");
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment alpha = arena.allocate(ValueLayout.JAVA_DOUBLE), nx = arena.allocate(ValueLayout.JAVA_INT), ny = arena.allocate(ValueLayout.JAVA_INT), nz = arena.allocate(ValueLayout.JAVA_INT);
      if (ljpme) {
        OpenMMNative.OpenMM_NonbondedForce_getLJPMEParametersInContext(
            getPointer(), context.getPointer(), alpha, nx, ny, nz);
      } else {
        OpenMMNative.OpenMM_NonbondedForce_getPMEParametersInContext(
            getPointer(), context.getPointer(), alpha, nx, ny, nz);
      }
      return new PMEParameters(alpha.get(ValueLayout.JAVA_DOUBLE, 0), nx.get(ValueLayout.JAVA_INT, 0),
          ny.get(ValueLayout.JAVA_INT, 0), nz.get(ValueLayout.JAVA_INT, 0));
    }
  }

  private static MemorySegment create() {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_NonbondedForce_create.makeInvoker().apply();
  }
}
