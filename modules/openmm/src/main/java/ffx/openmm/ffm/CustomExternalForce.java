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

/**
 * Custom algebraic external potential applied independently to selected particles.
 *
 * <p>The energy expression can use particle coordinates {@code x}, {@code y}, and {@code z}
 * (nm), global parameters, per-particle parameters, and OpenMM's custom-expression operators,
 * functions, and tabulated functions. Energy is in kJ/mol; parameter values have no automatic
 * unit conversion. Each particle's parameter vector contains one value per declared per-particle
 * parameter in declaration order.
 *
 * <p>When using periodic systems, the expression receives the particle's coordinates; it is not
 * automatically made translationally periodic. Stored per-particle value edits affect existing
 * contexts only after {@link #updateParametersInContext(Context)}. That method cannot update the
 * expression, topology, declarations, or global defaults; recreate contexts after such structural
 * changes. Global parameter values in a live context use the context parameter API.
 *
 * <p>Java arrays are copied through temporary native arrays. {@link DoubleArray} arguments are
 * borrowed synchronously and stay caller-owned; returned record arrays are copied into Java.
 * Java strings are encoded as temporary UTF-8 data, while returned strings are copied.
 */
public class CustomExternalForce extends Force {

  /**
   * Snapshot of a selected particle and its per-particle values.
   *
   * @param particle particle index
   * @param parameters copied values in declaration order
   */
  public record ParticleParameters(int particle, double[] parameters) {}

  /**
   * Create an external force from an OpenMM custom expression.
   *
   * @param energy expression evaluated per particle; coordinates {@code x}, {@code y}, and
   *     {@code z} are in nm
   */
  public CustomExternalForce(String energy) {
    super(create(energy));
  }

  /**
   * Declare an expression-wide parameter and its default for newly created contexts.
   *
   * @param name expression parameter name
   * @param defaultValue initial value, with units chosen consistently with the expression
   * @return index assigned to the declaration
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomExternalForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /**
   * Add a selected particle and all its per-particle values.
   *
   * @param particle particle index
   * @param parameters values in declaration order, copied during this call; length must equal
   *     {@link #getNumPerParticleParameters()}
   * @return index assigned to the selected particle
   */
  public int addParticle(int particle, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      return addParticle(particle, values);
    }
  }

  /**
   * Add a selected particle using a caller-owned native parameter array borrowed for this call.
   *
   * @param particle particle index
   * @param parameters values in declaration order, one for every declared per-particle parameter
   * @return index assigned to the selected particle
   */
  public int addParticle(int particle, DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomExternalForce_addParticle(
        getPointer(), particle, parameters.getPointer());
  }

  /**
   * Declare a parameter with a separate value for each selected particle.
   *
   * @param name expression parameter name
   * @return index assigned to the declaration
   */
  public int addPerParticleParameter(String name) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomExternalForce_addPerParticleParameter(
            getPointer(), value));
  }

  /** Destroy the native force and release its owned resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomExternalForce_destroy);
  }

  /** @return Java copy of the current energy expression. */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomExternalForce_getEnergyFunction(getPointer()));
  }

  /** @param index global-parameter index
   *  @return default value for contexts created after this setting */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomExternalForce_getGlobalParameterDefaultValue(getPointer(), index);
  }

  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomExternalForce_getGlobalParameterName(getPointer(), index));
  }

  /** @return number of declared global parameters */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomExternalForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of particles selected for this external force */
  public int getNumParticles() {
    return OpenMMNative.OpenMM_CustomExternalForce_getNumParticles(getPointer());
  }

  /** @return number of declarations and required values per selected particle */
  public int getNumPerParticleParameters() {
    return OpenMMNative.OpenMM_CustomExternalForce_getNumPerParticleParameters(getPointer());
  }

  /**
   * Get one selected particle's index and a copy of its values.
   *
   * @param index index into the force's selected-particle list
   * @return snapshot with values in per-particle declaration order
   */
  public ParticleParameters getParticleParameters(int index) {
    try (java.lang.foreign.Arena arena = java.lang.foreign.Arena.ofConfined();
         DoubleArray parameters = new DoubleArray(0)) {
      var particle = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomExternalForce_getParticleParameters(
          getPointer(), index, particle, parameters.getPointer());
      return new ParticleParameters(particle.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          CustomForceParameters.copy(parameters));
    }
  }

  /** @param index per-particle parameter index
   *  @return copied parameter name */
  public String getPerParticleParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomExternalForce_getPerParticleParameterName(getPointer(), index));
  }

  /** Replace the expression; contexts must be recreated to use the new expression.
   * @param energy expression in which particle coordinates are in nm */
  public void setEnergyFunction(String energy) {
    OpenMMStrings.withUtf8String(energy,
        value -> OpenMMNative.OpenMM_CustomExternalForce_setEnergyFunction(getPointer(), value));
  }

  /** Set the default for future contexts, not the current value in existing contexts.
   * @param index global-parameter index
   * @param defaultValue new default value */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMMNative.OpenMM_CustomExternalForce_setGlobalParameterDefaultValue(
        getPointer(), index, defaultValue);
  }

  /** Rename a declared global parameter; existing contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomExternalForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /**
   * Replace a selected particle's index and all per-particle values.
   *
   * @param index selected-particle list index
   * @param particle replacement particle index
   * @param parameters replacement values in declaration order, copied during the call
   */
  public void setParticleParameters(int index, int particle, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      setParticleParameters(index, particle, values);
    }
  }

  /**
   * Replace one selected particle from a caller-owned native array borrowed for this call.
   *
   * @param index selected-particle list index
   * @param particle replacement particle index
   * @param parameters replacement values in declaration order
   */
  public void setParticleParameters(int index, int particle, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomExternalForce_setParticleParameters(
        getPointer(), index, particle, parameters.getPointer());
  }

  /** Rename a declared per-particle parameter; existing contexts are not reconfigured.
   * @param index per-particle parameter index
   * @param name new expression name */
  public void setPerParticleParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomExternalForce_setPerParticleParameterName(
            getPointer(), index, value));
  }

  /**
   * Apply supported per-particle value changes in a context.
   *
   * @param context context containing this force; no-op if it has no native context pointer.
   *     Changes to expression, selection topology, declarations, and global defaults are not
   *     propagated by this operation.
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(context,
        value -> OpenMMNative.OpenMM_CustomExternalForce_updateParametersInContext(
            getPointer(), value));
  }

  /** @return whether this force reports use of periodic boundary conditions. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomExternalForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static java.lang.foreign.MemorySegment create(String energy) {
    return OpenMMStrings.withUtf8StringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomExternalForce_create(value);
    });
  }
}
