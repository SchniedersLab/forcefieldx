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
import ffx.openmm.ffm.bindings.OpenMM_Vec3;
import java.lang.foreign.Arena;
import java.lang.foreign.MemorySegment;
import java.lang.foreign.ValueLayout;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Implements the Alchemical Transfer Method (ATM) energy from a pair of system energies.
 *
 * <p>The force evaluates child forces for an initial state ({@code u0}) and a target state
 * ({@code u1}), using the per-particle displacement vectors supplied here, and applies an
 * algebraic energy expression to those energies and any global parameters. Add one displacement
 * pair per system particle, in the same order as particles in the system. The target-state
 * energy uses {@code displacement1}; the initial-state energy uses {@code displacement0}. Both
 * vectors are added to the corresponding particle coordinates. Displacements are in nm; energies
 * use kJ/mol. The expression constructor accepts the usual OpenMM algebraic expression syntax and
 * names {@code u0} and {@code u1}. The parameterized constructor selects OpenMM's default
 * soft-core softplus expression and initializes its named parameters.
 *
 * <p>Child forces added through the {@link Force} overload transfer to ATM ownership; the
 * {@link MemorySegment} overload passes a raw native handle and cannot invalidate any Java
 * wrapper for it. Handles returned by {@link #getForce(int)} are borrowed from this force.
 *
 * <p>This FFM wrapper does not expose the native header's static helpers for the default
 * parameter names (such as {@code Lambda1()}, {@code Alpha()}, and {@code Direction()}); use
 * those literal parameter names when configuring the parameterized expression.
 */
public class ATMForce extends Force {

  /**
   * Copied displacement vectors for a particle's target and initial states.
   *
   * @param displacement1 target-state displacement in nm
   * @param displacement0 initial-state displacement in nm
   */
  public record ParticleParameters(Vec3 displacement1, Vec3 displacement0) {}
  /**
   * The two state energies and the resulting ATM expression energy.
   *
   * @param u0 initial-state energy in kJ/mol
   * @param u1 target/displaced-state energy in kJ/mol
   * @param energy value of this ATM force's energy expression in kJ/mol
   */
  public record PerturbationEnergy(double u0, double u1, double energy) {}

  /**
   * Create an ATM force with a custom energy expression.
   *
   * @param energy OpenMM algebraic expression for the ATM energy; it can refer to the child-state
   *     energies {@code u0} and {@code u1} and to global parameters added separately
   */
  public ATMForce(String energy) {
    super(create(energy));
  }

  /**
   * Create an ATM force using OpenMM's default soft-core softplus energy expression.
   *
   * <p>These values initialize the named global parameters in newly created contexts. The
   * dimensionless {@code lambda1} and {@code lambda2} are conventionally between 0 and 1;
   * {@code alpha} is in (kJ/mol)<sup>-1</sup>; {@code uh}, {@code w0}, {@code umax}, and
   * {@code ubcore} are in kJ/mol; {@code acore} is dimensionless; and {@code direction} is
   * dimensionless and should be {@code 1} for forward transfer or {@code -1} for backward
   * transfer. The corresponding native parameter names are {@code Lambda1}, {@code Lambda2},
   * {@code Alpha}, {@code Uh}, {@code W0}, {@code Umax}, {@code Ubcore}, {@code Acore}, and
   * {@code Direction}.
   *
   * @param lambda1 default {@code Lambda1} value, dimensionless
   * @param lambda2 default {@code Lambda2} value, dimensionless
   * @param alpha default {@code Alpha} value in (kJ/mol)<sup>-1</sup>
   * @param uh default {@code Uh} value in kJ/mol
   * @param w0 default {@code W0} value in kJ/mol
   * @param umax default {@code Umax} value in kJ/mol
   * @param ubcore default {@code Ubcore} value in kJ/mol
   * @param acore default {@code Acore} value, dimensionless
   * @param direction default {@code Direction} value, conventionally {@code 1} or {@code -1}
   */
  public ATMForce(
      double lambda1, double lambda2, double alpha, double uh, double w0,
      double umax, double ubcore, double acore, double direction) {
    super(create(lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore, direction));
  }

  /**
   * Request energy-derivative calculation with respect to a previously added global parameter.
   *
   * @param name exact name of the global parameter
   */
  public void addEnergyParameterDerivative(String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_ATMForce_addEnergyParameterDerivative(getPointer(), value));
  }

  /**
   * Request energy-derivative calculation using a caller-provided native string.
   *
   * @param name native null-terminated string naming a previously added global parameter; valid
   *     for the duration of this call
   */
  public void addEnergyParameterDerivative(MemorySegment name) {
    OpenMMNative.OpenMM_ATMForce_addEnergyParameterDerivative(getPointer(), name);
  }

  /**
   * Add a child force whose energy contributes to the ATM state energies.
   *
   * <p>Ownership transfers to this ATM force and the passed Java façade is invalidated; it must
   * not be used or destroyed independently afterward.
   *
   * @param force child force to add
   * @return index of the child force
   */
  public int addForce(Force force) {
    int index = OpenMMNative.OpenMM_ATMForce_addForce(getPointer(), force.getPointer());
    force.invalidate();
    return index;
  }

  /**
   * Add a child force by raw native handle. Native ownership transfers to this ATM force, but
   * this overload cannot invalidate any Java façade that may wrap the handle.
   *
   * @param force native handle for a child force whose ownership can be transferred
   * @return index of the child force
   */
  public int addForce(MemorySegment force) {
    return OpenMMNative.OpenMM_ATMForce_addForce(getPointer(), force);
  }

  /**
   * Add a global parameter and its default value for newly created contexts. A running context's
   * value can be changed through its parameter API.
   *
   * @param name parameter name referenced by the energy expression
   * @param defaultValue initial value for newly created contexts
   * @return index of the added parameter
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_ATMForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /**
   * Add a global parameter using a caller-provided native string.
   *
   * @param name native null-terminated parameter name, valid for the duration of this call
   * @param defaultValue initial value for newly created contexts
   * @return index of the added parameter
   */
  public int addGlobalParameter(MemorySegment name, double defaultValue) {
    return OpenMMNative.OpenMM_ATMForce_addGlobalParameter(
        getPointer(), name, defaultValue);
  }

  /**
   * Add a particle's target-state and initial-state displacement vectors.
   *
   * @param displacement1 target-state displacement in nm
   * @param displacement0 initial-state displacement in nm
   * @return index of the added particle; add particles in system particle order
   */
  public int addParticle(Vec3 displacement1, Vec3 displacement0) {
    try (Arena arena = Arena.ofConfined()) {
      return OpenMMNative.OpenMM_ATMForce_addParticle(
          getPointer(), displacement1.toNative(arena), displacement0.toNative(arena));
    }
  }

  /**
   * Add a particle using caller-owned native vector structs. OpenMM copies their vector values
   * during this call; the caller retains ownership of the structs.
   *
   * @param displacement1 target-state displacement vector in nm
   * @param displacement0 initial-state displacement vector in nm
   * @return index of the added particle; add particles in system particle order
   */
  public int addParticle(MemorySegment displacement1, MemorySegment displacement0) {
    return OpenMMNative.OpenMM_ATMForce_addParticle(
        getPointer(), displacement1, displacement0);
  }

  /** Destroy this force and its owned child forces. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_ATMForce_destroy);
  }

  /**
   * Get the energy expression as a Java string copy.
   *
   * @return configured energy expression
   */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_ATMForce_getEnergyFunction(getPointer()));
  }

  /**
   * Get the name of a requested energy-parameter derivative as a Java string copy.
   *
   * @param index index among requested derivatives
   * @return corresponding global parameter name
   */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_ATMForce_getEnergyParameterDerivativeName(getPointer(), index));
  }

  /**
   * Get a borrowed native handle for a child force owned by this ATM force.
   *
   * <p>The returned handle is not owned by the caller and becomes invalid when this ATM force is
   * destroyed. The raw handle API does not create a child {@link Force} façade.
   *
   * @param index index of the child force
   * @return borrowed child-force handle
   */
  public MemorySegment getForce(int index) {
    return OpenMMNative.OpenMM_ATMForce_getForce(getPointer(), index);
  }

  /**
   * Get the default value of a global parameter.
   *
   * @param index index of the parameter
   * @return default value used for newly created contexts
   */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_ATMForce_getGlobalParameterDefaultValue(getPointer(), index);
  }

  /**
   * Get the name of a global parameter as a Java string copy.
   *
   * @param index index of the parameter
   * @return parameter name
   */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_ATMForce_getGlobalParameterName(getPointer(), index));
  }

  /** @return number of global-parameter energy derivatives requested */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_ATMForce_getNumEnergyParameterDerivatives(getPointer());
  }

  /** @return number of child forces owned by this ATM force */
  public int getNumForces() {
    return OpenMMNative.OpenMM_ATMForce_getNumForces(getPointer());
  }

  /** @return number of global parameters defined for this force */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_ATMForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of particle displacement pairs defined for this force */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_ATMForce_getNumParticles(getPointer());
  }

  /**
   * Get a particle's target-state and initial-state displacements as copied Java values.
   *
   * @param index index of the particle in ATM particle order
   * @return copied target displacement followed by copied initial displacement, both in nm
   */
  public ParticleParameters getParticleParameters(int index) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment displacement1 = arena.allocate(OpenMM_Vec3.layout());
      MemorySegment displacement0 = arena.allocate(OpenMM_Vec3.layout());
      OpenMMNative.OpenMM_ATMForce_getParticleParameters(
          getPointer(), index, displacement1, displacement0);
      return new ParticleParameters(
          Vec3.fromNative(displacement1), Vec3.fromNative(displacement0));
    }
  }

  /**
   * Evaluate and return the two state energies and this force's energy in a context.
   *
   * @param context context whose current state is evaluated
   * @return copied initial-state energy {@code u0}, target/displaced-state energy {@code u1}, and
   *     ATM expression energy, all in kJ/mol
   */
  public PerturbationEnergy getPerturbationEnergy(Context context) {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment u1 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment u0 = arena.allocate(ValueLayout.JAVA_DOUBLE);
      MemorySegment energy = arena.allocate(ValueLayout.JAVA_DOUBLE);
      OpenMMNative.OpenMM_ATMForce_getPerturbationEnergy(
          getPointer(), context.getPointer(), u1, u0, energy);
      return new PerturbationEnergy(
          u0.get(ValueLayout.JAVA_DOUBLE, 0), u1.get(ValueLayout.JAVA_DOUBLE, 0),
          energy.get(ValueLayout.JAVA_DOUBLE, 0));
    }
  }

  /**
   * Replace the ATM energy expression.
   *
   * @param energy new OpenMM algebraic expression using {@code u0}, {@code u1}, and any defined
   *     global parameters
   */
  public void setEnergyFunction(String energy) {
    withString(energy, value ->
        OpenMMNative.OpenMM_ATMForce_setEnergyFunction(getPointer(), value));
  }

  /**
   * Replace the ATM energy expression from a caller-provided native string.
   *
   * @param energy native null-terminated expression string, valid for the duration of this call
   */
  public void setEnergyFunction(MemorySegment energy) {
    OpenMMNative.OpenMM_ATMForce_setEnergyFunction(getPointer(), energy);
  }

  /**
   * Set a global parameter's default value for newly created contexts.
   *
   * @param index index of the parameter
   * @param value new default value
   */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_ATMForce_setGlobalParameterDefaultValue(getPointer(), index, value);
  }

  /**
   * Rename a global parameter.
   *
   * @param index index of the parameter
   * @param name new parameter name
   */
  public void setGlobalParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_ATMForce_setGlobalParameterName(getPointer(), index, value));
  }

  /**
   * Rename a global parameter using a caller-provided native string.
   *
   * @param index index of the parameter
   * @param name native null-terminated parameter name, valid for the duration of this call
   */
  public void setGlobalParameterName(int index, MemorySegment name) {
    OpenMMNative.OpenMM_ATMForce_setGlobalParameterName(getPointer(), index, name);
  }

  /**
   * Replace a particle's target-state and initial-state displacement vectors.
   *
   * @param index index of the particle in ATM particle order
   * @param displacement1 target-state displacement in nm
   * @param displacement0 initial-state displacement in nm
   */
  public void setParticleParameters(
      int index, Vec3 displacement1, Vec3 displacement0) {
    try (Arena arena = Arena.ofConfined()) {
      OpenMMNative.OpenMM_ATMForce_setParticleParameters(
          getPointer(), index, displacement1.toNative(arena), displacement0.toNative(arena));
    }
  }

  /**
   * Replace displacements from caller-owned native vector structs; OpenMM copies their values.
   *
   * @param index index of the particle in ATM particle order
   * @param displacement1 target-state displacement struct in nm
   * @param displacement0 initial-state displacement struct in nm
   */
  public void setParticleParameters(
      int index, MemorySegment displacement1, MemorySegment displacement0) {
    OpenMMNative.OpenMM_ATMForce_setParticleParameters(
        getPointer(), index, displacement1, displacement0);
  }

  /**
   * Copy changed per-particle displacements into an existing context without reinitializing it.
   * The number of particles cannot be changed by this operation.
   *
   * @param context context to update; this FFM wrapper performs no update if it has no native
   *     context handle
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(
        context, pointer ->
            OpenMMNative.OpenMM_ATMForce_updateParametersInContext(getPointer(), pointer));
  }

  /** @return native OpenMM's periodic-boundary-condition status for this force */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_ATMForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(String energy) {
    return withStringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_ATMForce_create(value);
    });
  }

  private static MemorySegment create(
      double lambda1, double lambda2, double alpha, double uh, double w0,
      double umax, double ubcore, double acore, double direction) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_ATMForce_create_2(
        lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore, direction);
  }

  private static void withString(String value, Consumer<MemorySegment> action) {
    OpenMMStrings.withUtf8String(value, action);
  }

  private static <T> T withStringResult(String value, Function<MemorySegment, T> action) {
    return OpenMMStrings.withUtf8StringResult(value, action);
  }
}
