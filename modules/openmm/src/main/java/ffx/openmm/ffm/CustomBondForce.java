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
 * Custom algebraic interactions between pairs of particles.
 *
 * <p>The expression is evaluated once per bond and may use {@code r}, the interparticle distance
 * in nm, global parameters, per-bond parameters, and OpenMM's custom-expression operators,
 * functions, and tabulated functions. The energy unit is kJ/mol. Parameter values are supplied
 * without conversion and must be consistent with the expression. A bond's parameter vector must
 * contain one value per declared per-bond parameter, in declaration order.
 *
 * <p>Edits to stored bond values do not affect an existing context until
 * {@link #updateParametersInContext(Context)} is called. That operation propagates only parameter
 * values that OpenMM permits to be updated; changing the expression, parameter declarations,
 * topology, or global defaults requires a new context. Global values in a live context are
 * controlled through the context parameter API.
 *
 * <p>Java primitive arrays are copied to temporary native arrays; {@link DoubleArray} arguments
 * are borrowed synchronously and remain caller-owned. Returned record arrays are independent Java
 * copies. Java {@link String} values are passed as temporary UTF-8 native strings and returned
 * strings are copied.
 */
public class CustomBondForce extends Force {

  /**
   * Snapshot of a bond's particle indices and per-bond values.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param parameters copied values in per-bond declaration order
   */
  public record BondParameters(int particle1, int particle2, double[] parameters) {}

  /**
   * Create a custom bond force.
   *
   * @param energy expression evaluated for every bond; {@code r} is distance in nm
   */
  public CustomBondForce(String energy) {
    super(create(energy));
  }

  /**
   * Add a bond and all its per-bond values.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param parameters values in declaration order, copied during the call; length must equal
   *     {@link #getNumPerBondParameters()}
   * @return index assigned to the bond
   */
  public int addBond(int particle1, int particle2, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      return addBond(particle1, particle2, values);
    }
  }

  /**
   * Add a bond using a caller-owned native array borrowed for this call.
   *
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param parameters per-bond values in declaration order, one for each declaration
   * @return index assigned to the bond
   */
  public int addBond(int particle1, int particle2, DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomBondForce_addBond(
        getPointer(), particle1, particle2, parameters.getPointer());
  }

  /**
   * Request an energy derivative with respect to a declared global parameter.
   *
   * @param name exact global parameter name
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomBondForce_addEnergyParameterDerivative(getPointer(), value));
  }

  /**
   * Declare a global parameter and the value used in newly created contexts.
   *
   * @param name expression parameter name
   * @param defaultValue initial value; units are determined by the expression
   * @return index assigned to the parameter
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomBondForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /**
   * Declare a parameter whose value can differ for each bond.
   *
   * @param name expression parameter name
   * @return index assigned to the declaration
   */
  public int addPerBondParameter(String name) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomBondForce_addPerBondParameter(getPointer(), value));
  }

  /** Destroy the native force and release its owned resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomBondForce_destroy);
  }

  /**
   * Get one bond's indices and independent copies of its parameter values.
   *
   * @param index bond index
   * @return snapshot with values in per-bond declaration order
   */
  public BondParameters getBondParameters(int index) {
    try (java.lang.foreign.Arena arena = java.lang.foreign.Arena.ofConfined();
         DoubleArray parameters = new DoubleArray(0)) {
      var particle1 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      var particle2 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomBondForce_getBondParameters(
          getPointer(), index, particle1, particle2, parameters.getPointer());
      return new BondParameters(particle1.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          particle2.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          CustomForceParameters.copy(parameters));
    }
  }

  /** @return Java copy of the current energy expression. */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomBondForce_getEnergyFunction(getPointer()));
  }

  /** @param index derivative index, from zero to count minus one
   *  @return name of the differentiated global parameter */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomBondForce_getEnergyParameterDerivativeName(getPointer(), index));
  }

  /** @param index global-parameter index
   *  @return default for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomBondForce_getGlobalParameterDefaultValue(getPointer(), index);
  }

  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomBondForce_getGlobalParameterName(
        getPointer(), index));
  }

  /** @return number of bonds stored in this force */
  public int getNumBonds() {
    return OpenMMNative.OpenMM_CustomBondForce_getNumBonds(getPointer());
  }

  /** @return number of requested global-parameter energy derivatives */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_CustomBondForce_getNumEnergyParameterDerivatives(getPointer());
  }

  /** @return number of declared global parameters */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomBondForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of declarations and required values in each bond parameter vector */
  public int getNumPerBondParameters() {
    return OpenMMNative.OpenMM_CustomBondForce_getNumPerBondParameters(getPointer());
  }

  /** @param index per-bond parameter index
   *  @return copied parameter name */
  public String getPerBondParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomBondForce_getPerBondParameterName(getPointer(), index));
  }

  /**
   * Replace one bond's particle indices and all per-bond values.
   *
   * @param index stored bond index
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param parameters replacement values in declaration order, copied for this call, with one
   *     value per declared per-bond parameter
   */
  public void setBondParameters(int index, int particle1, int particle2, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      setBondParameters(index, particle1, particle2, values);
    }
  }

  /**
   * Replace bond values from a caller-owned native array borrowed for this call.
   *
   * @param index stored bond index
   * @param particle1 first particle index
   * @param particle2 second particle index
   * @param parameters replacement per-bond values in declaration order
   */
  public void setBondParameters(int index, int particle1, int particle2, DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomBondForce_setBondParameters(
        getPointer(), index, particle1, particle2, parameters.getPointer());
  }

  /** Replace the expression; existing contexts must be recreated to use the new expression.
   * @param energy OpenMM custom expression with {@code r} measured in nm */
  public void setEnergyFunction(String energy) {
    OpenMMStrings.withUtf8String(energy,
        value -> OpenMMNative.OpenMM_CustomBondForce_setEnergyFunction(getPointer(), value));
  }

  /** Set the default for future contexts; does not change a live context's current value.
   * @param index global-parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_CustomBondForce_setGlobalParameterDefaultValue(getPointer(), index, value);
  }

  /** Rename a declared global parameter; existing contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new name */
  public void setGlobalParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomBondForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /** Rename a declared per-bond parameter; existing contexts are not reconfigured.
   * @param index per-bond parameter index
   * @param name new name */
  public void setPerBondParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomBondForce_setPerBondParameterName(
            getPointer(), index, value));
  }

  /**
   * Set whether bond distances use minimum-image periodic displacements.
   *
   * @param periodic {@code true} to use periodic displacements
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_CustomBondForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Apply supported bond-value changes to a context.
   *
   * @param context context containing this force; the helper is a no-op if it has no native
   *     context pointer. This does not propagate expression, declaration, topology, or global
   *     default changes.
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(context,
        value -> OpenMMNative.OpenMM_CustomBondForce_updateParametersInContext(getPointer(), value));
  }

  /** @return {@code true} if periodic boundary conditions are enabled for bond distances. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomBondForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static java.lang.foreign.MemorySegment create(String energy) {
    return OpenMMStrings.withUtf8StringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomBondForce_create(value);
    });
  }
}
