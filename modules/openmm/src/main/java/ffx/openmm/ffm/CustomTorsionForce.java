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

/**
 * Custom algebraic interactions for ordered quadruples of particles.
 *
 * <p>The expression is evaluated per torsion and may use {@code theta}, the torsion angle in
 * radians, global and per-torsion parameters, and OpenMM custom-expression operators, functions,
 * and tabulated functions. Energy is in kJ/mol; parameter values are not converted and must use
 * units consistent with the expression. A torsion's values are ordered according to the
 * per-torsion declarations.
 *
 * <p>Changing stored torsion values does not modify an existing context until
 * {@link #updateParametersInContext(Context)} is called. Only supported torsion parameter values
 * are propagated; changing the expression, declarations, topology, or global defaults requires
 * context recreation. The context parameter API changes live global parameter values.
 *
 * <p>Java-array values are copied to temporary native arrays; {@link DoubleArray} and
 * {@link MemorySegment} inputs are borrowed synchronously and remain caller-owned. The
 * {@code MemorySegment} string overloads additionally require a live, NUL-terminated UTF-8
 * segment for the duration of the call. Returned record arrays and strings are Java copies.
 */
public class CustomTorsionForce extends Force {

  /**
   * Snapshot of a torsion's ordered particle indices and per-torsion values.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param particle3 third particle index
   * @param particle4 fourth particle index
   * @param parameters copied values in declaration order
   */
  public record TorsionParameters(
      int particle1, int particle2, int particle3, int particle4, double[] parameters) {}

  /**
   * Create a custom torsion force.
   *
   * @param energy expression evaluated for each torsion; {@code theta} is in radians
   */
  public CustomTorsionForce(String energy) {
    super(create(energy));
  }

  /**
   * Request the energy derivative with respect to a previously declared global parameter.
   *
   * @param name exact global parameter name
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_addEnergyParameterDerivative(
            getPointer(), value));
  }

  /**
   * Request a derivative using a caller-owned NUL-terminated UTF-8 name.
   *
   * @param name live native name segment, borrowed for this call
   */
  public void addEnergyParameterDerivative(MemorySegment name) {
    OpenMMNative.OpenMM_CustomTorsionForce_addEnergyParameterDerivative(getPointer(), name);
  }

  /**
   * Declare an expression-wide parameter with a default for newly created contexts.
   *
   * @param name expression parameter name
   * @param defaultValue initial value; units are determined by the expression
   * @return index assigned to the declaration
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /**
   * Declare a global parameter from a caller-owned NUL-terminated UTF-8 name.
   *
   * @param name live native name segment, borrowed for this call
   * @param defaultValue initial value
   * @return index assigned to the declaration
   */
  public int addGlobalParameter(MemorySegment name, double defaultValue) {
    return OpenMMNative.OpenMM_CustomTorsionForce_addGlobalParameter(
        getPointer(), name, defaultValue);
  }

  /**
   * Declare a parameter whose value can differ for each torsion.
   *
   * @param name expression parameter name
   * @return index assigned to the declaration
   */
  public int addPerTorsionParameter(String name) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_addPerTorsionParameter(
            getPointer(), value));
  }

  /**
   * Declare a per-torsion parameter from a caller-owned NUL-terminated UTF-8 name.
   *
   * @param name live native name segment, borrowed for this call
   * @return index assigned to the declaration
   */
  public int addPerTorsionParameter(MemorySegment name) {
    return OpenMMNative.OpenMM_CustomTorsionForce_addPerTorsionParameter(getPointer(), name);
  }

  /**
   * Add an ordered four-particle torsion and its per-torsion values.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param particle3 third particle index
   * @param particle4 fourth particle index
   * @param parameters values in declaration order, copied for this call; one per declaration
   * @return index assigned to the torsion
   */
  public int addTorsion(
      int particle1, int particle2, int particle3, int particle4, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      return addTorsion(particle1, particle2, particle3, particle4, values);
    }
  }

  /**
   * Add a torsion using a caller-owned native array borrowed for this call.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param particle3 third particle index
   * @param particle4 fourth particle index
   * @param parameters per-torsion values in declaration order
   * @return index assigned to the torsion
   */
  public int addTorsion(
      int particle1, int particle2, int particle3, int particle4, DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomTorsionForce_addTorsion(
        getPointer(), particle1, particle2, particle3, particle4, parameters.getPointer());
  }

  /** Destroy the native force and release its owned resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomTorsionForce_destroy);
  }

  /** @return Java copy of the current energy expression. */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomTorsionForce_getEnergyFunction(getPointer()));
  }

  /** @param index derivative index, from zero to count minus one
   *  @return copied name of the differentiated global parameter */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomTorsionForce_getEnergyParameterDerivativeName(
            getPointer(), index));
  }

  /** @param index global-parameter index
   *  @return default value for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomTorsionForce_getGlobalParameterDefaultValue(
        getPointer(), index);
  }

  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomTorsionForce_getGlobalParameterName(getPointer(), index));
  }

  /** @return number of requested global-parameter energy derivatives */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_CustomTorsionForce_getNumEnergyParameterDerivatives(getPointer());
  }

  /** @return number of declared global parameters */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomTorsionForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of declared per-torsion parameters and values required per torsion */
  public int getNumPerTorsionParameters() {
    return OpenMMNative.OpenMM_CustomTorsionForce_getNumPerTorsionParameters(getPointer());
  }

  /** @return number of torsions stored in this force */
  public int getNumTorsions() {
    return OpenMMNative.OpenMM_CustomTorsionForce_getNumTorsions(getPointer());
  }

  /** @param index per-torsion parameter index
   *  @return copied parameter name */
  public String getPerTorsionParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomTorsionForce_getPerTorsionParameterName(getPointer(), index));
  }

  /**
   * Return one torsion and copies of all its parameter values.
   *
   * @param index torsion index
   * @return ordered particle indices and values in declaration order
   */
  public TorsionParameters getTorsionParameters(int index) {
    try (java.lang.foreign.Arena arena = java.lang.foreign.Arena.ofConfined();
         DoubleArray parameters = new DoubleArray(0)) {
      var particle1 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      var particle2 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      var particle3 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      var particle4 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomTorsionForce_getTorsionParameters(
          getPointer(), index, particle1, particle2, particle3, particle4, parameters.getPointer());
      return new TorsionParameters(particle1.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          particle2.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          particle3.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          particle4.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          CustomForceParameters.copy(parameters));
    }
  }

  /** Replace the expression; existing contexts must be recreated to use it.
   * @param energy OpenMM custom expression; {@code theta} is in radians */
  public void setEnergyFunction(String energy) {
    OpenMMStrings.withUtf8String(energy,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_setEnergyFunction(getPointer(), value));
  }

  /** Replace the expression from a caller-owned NUL-terminated UTF-8 segment.
   * @param energy live native expression segment, borrowed for this call */
  public void setEnergyFunction(MemorySegment energy) {
    OpenMMNative.OpenMM_CustomTorsionForce_setEnergyFunction(getPointer(), energy);
  }

  /** Set the default used in future contexts, not the current value of existing contexts.
   * @param index global-parameter index
   * @param defaultValue new default value */
  public void setGlobalParameterDefaultValue(int index, double defaultValue) {
    OpenMMNative.OpenMM_CustomTorsionForce_setGlobalParameterDefaultValue(
        getPointer(), index, defaultValue);
  }

  /** Rename a global parameter; existing contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /** Rename a global parameter using a live NUL-terminated UTF-8 segment.
   * @param index global-parameter index
   * @param name caller-owned segment borrowed for this call */
  public void setGlobalParameterName(int index, MemorySegment name) {
    OpenMMNative.OpenMM_CustomTorsionForce_setGlobalParameterName(getPointer(), index, name);
  }

  /** Rename a per-torsion parameter; existing contexts are not reconfigured.
   * @param index per-torsion parameter index
   * @param name new expression name */
  public void setPerTorsionParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_setPerTorsionParameterName(
            getPointer(), index, value));
  }

  /** Rename a per-torsion parameter using a live NUL-terminated UTF-8 segment.
   * @param index per-torsion parameter index
   * @param name caller-owned segment borrowed for this call */
  public void setPerTorsionParameterName(int index, MemorySegment name) {
    OpenMMNative.OpenMM_CustomTorsionForce_setPerTorsionParameterName(getPointer(), index, name);
  }

  /**
   * Replace one torsion's ordered indices and all its parameter values.
   *
   * @param index stored torsion index
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param particle3 third particle index
   * @param particle4 fourth particle index
   * @param parameters replacement values in declaration order, copied during this call
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3,
      int particle4, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      setTorsionParameters(index, particle1, particle2, particle3, particle4, values);
    }
  }

  /**
   * Replace a torsion from a caller-owned native parameter array borrowed for this call.
   *
   * @param index stored torsion index
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param particle3 third particle index
   * @param particle4 fourth particle index
   * @param parameters replacement values in declaration order
   */
  public void setTorsionParameters(int index, int particle1, int particle2, int particle3,
      int particle4, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomTorsionForce_setTorsionParameters(
        getPointer(), index, particle1, particle2, particle3, particle4, parameters.getPointer());
  }

  /** Set whether torsion geometry uses periodic minimum-image displacements.
   * @param periodic {@code true} to use periodic displacements */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_CustomTorsionForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Apply supported per-torsion value changes to an existing context.
   *
   * @param context context containing this force; no-op if it has no native context pointer.
   *     Does not propagate expression, topology, declaration, or global-default changes.
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(context,
        value -> OpenMMNative.OpenMM_CustomTorsionForce_updateParametersInContext(
            getPointer(), value));
  }

  /** @return {@code true} if periodic boundary conditions are enabled. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomTorsionForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static java.lang.foreign.MemorySegment create(String energy) {
    return OpenMMStrings.withUtf8StringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomTorsionForce_create(value);
    });
  }
}
