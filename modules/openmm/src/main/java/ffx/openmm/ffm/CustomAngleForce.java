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
 * Custom algebraic interactions for ordered triples of particles.
 *
 * <p>The constructor expression is evaluated once for each angle. It may use {@code theta} (the
 * angle in radians), global parameters, per-angle parameters, and supported OpenMM custom
 * expression operators, functions, and tabulated functions. OpenMM reports energy in kJ/mol;
 * parameter values are unitless at this API boundary and must use units consistent with the
 * expression. Per-angle values correspond positionally to parameters declared by
 * {@link #addPerAngleParameter(String)}.
 *
 * <p>Changing a force definition after creating a {@link Context} does not change that context
 * automatically. {@link #updateParametersInContext(Context)} propagates supported per-angle
 * parameter changes only; it does not recompile changed expressions or update parameter names,
 * global-parameter defaults, or other structural settings. Global parameter values in a running
 * context are changed through the context parameter API.
 *
 * <p>Java-array overloads copy values into temporary native arrays for the duration of the call.
 * {@link DoubleArray} overloads borrow the supplied native array synchronously; callers retain
 * ownership and must keep it alive for the call. Getters returning records contain independent
 * Java-array copies. Strings passed as {@link String} are encoded as temporary UTF-8 strings;
 * returned strings are copied into Java strings.
 */
public class CustomAngleForce extends Force {

  /**
   * Snapshot of one angle's particle indices and per-angle values.
   *
   * @param particle1 index of the first particle
   * @param particle2 index of the central particle
   * @param particle3 index of the third particle
   * @param parameters copied values in declaration order of the per-angle parameters
   */
  public record AngleParameters(int particle1, int particle2, int particle3, double[] parameters) {}

  /**
   * Create a force with the supplied OpenMM custom energy expression.
   *
   * @param energy expression evaluated for each angle; {@code theta} denotes its angle in radians
   */
  public CustomAngleForce(String energy) {
    super(create(energy));
  }

  /**
   * Add an angle with values for all declared per-angle parameters.
   *
   * @param particle1 first particle index
   * @param particle2 central particle index
   * @param particle3 third particle index
   * @param parameters values in the order parameters were declared; the array is copied during
   *     this call and its length must equal {@link #getNumPerAngleParameters()}
   * @return index assigned to the angle
   */
  public int addAngle(int particle1, int particle2, int particle3, double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      return addAngle(particle1, particle2, particle3, values);
    }
  }

  /**
   * Add an angle using a borrowed native array of per-angle values in declaration order.
   *
   * @param particle1 first particle index
   * @param particle2 central particle index
   * @param particle3 third particle index
   * @param parameters caller-owned native array, borrowed only for this call, with one value per
   *     declared per-angle parameter
   * @return index assigned to the angle
   */
  public int addAngle(int particle1, int particle2, int particle3, DoubleArray parameters) {
    return OpenMMNative.OpenMM_CustomAngleForce_addAngle(
        getPointer(), particle1, particle2, particle3, parameters.getPointer());
  }

  /**
   * Request computation of the energy derivative for a previously declared global parameter.
   *
   * @param name exact name passed to {@link #addGlobalParameter(String, double)}
   */
  public void addEnergyParameterDerivative(String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomAngleForce_addEnergyParameterDerivative(
            getPointer(), value));
  }

  /**
   * Declare a global expression parameter and its default for newly created contexts.
   *
   * @param name parameter name referenced by the expression
   * @param defaultValue default parameter value (in units chosen consistently with the expression)
   * @return index assigned to the parameter
   */
  public int addGlobalParameter(String name, double defaultValue) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomAngleForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /**
   * Declare a per-angle expression parameter. Each angle subsequently added must provide exactly
   * one value for each declaration, in declaration order.
   *
   * @param name parameter name referenced by the expression
   * @return index assigned to the declaration
   */
  public int addPerAngleParameter(String name) {
    return OpenMMStrings.withUtf8StringResult(name,
        value -> OpenMMNative.OpenMM_CustomAngleForce_addPerAngleParameter(getPointer(), value));
  }

  /** Destroy the native force and release its owned resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomAngleForce_destroy);
  }

  /**
   * Return a snapshot of one angle and its copied per-angle values.
   *
   * @param index angle index
   * @return particle indices and values in per-angle declaration order
   */
  public AngleParameters getAngleParameters(int index) {
    try (java.lang.foreign.Arena arena = java.lang.foreign.Arena.ofConfined();
         DoubleArray parameters = new DoubleArray(0)) {
      var particle1 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      var particle2 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      var particle3 = arena.allocate(java.lang.foreign.ValueLayout.JAVA_INT);
      OpenMMNative.OpenMM_CustomAngleForce_getAngleParameters(
          getPointer(), index, particle1, particle2, particle3, parameters.getPointer());
      return new AngleParameters(particle1.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          particle2.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          particle3.get(java.lang.foreign.ValueLayout.JAVA_INT, 0),
          CustomForceParameters.copy(parameters));
    }
  }

  /** @return a Java copy of the current energy expression. */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomAngleForce_getEnergyFunction(getPointer()));
  }

  /** @param index derivative index; valid indices are from zero to count minus one
   *  @return name of the global parameter whose energy derivative was requested */
  public String getEnergyParameterDerivativeName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomAngleForce_getEnergyParameterDerivativeName(getPointer(), index));
  }

  /** @param index global-parameter index
   *  @return default value used when a context is created */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomAngleForce_getGlobalParameterDefaultValue(getPointer(), index);
  }

  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomAngleForce_getGlobalParameterName(
        getPointer(), index));
  }

  /** @return number of angle terms currently stored in this force */
  public int getNumAngles() {
    return OpenMMNative.OpenMM_CustomAngleForce_getNumAngles(getPointer());
  }

  /** @return number of global-parameter energy derivatives requested */
  public int getNumEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_CustomAngleForce_getNumEnergyParameterDerivatives(getPointer());
  }

  /** @return number of declared global parameters */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomAngleForce_getNumGlobalParameters(getPointer());
  }

  /** @return number of declared per-angle parameters (the required parameter-vector length) */
  public int getNumPerAngleParameters() {
    return OpenMMNative.OpenMM_CustomAngleForce_getNumPerAngleParameters(getPointer());
  }

  /** @param index per-angle parameter index
   *  @return copied parameter name */
  public String getPerAngleParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomAngleForce_getPerAngleParameterName(getPointer(), index));
  }

  /**
   * Replace an angle's indices and all per-angle values.
   *
   * @param index stored angle index
   * @param particle1 first particle index
   * @param particle2 central particle index
   * @param particle3 third particle index
   * @param parameters replacement values in declaration order; copied for this call and required
   *     to contain one value per declared per-angle parameter
   */
  public void setAngleParameters(int index, int particle1, int particle2, int particle3,
      double[] parameters) {
    try (DoubleArray values = CustomForceParameters.toNative(parameters)) {
      setAngleParameters(index, particle1, particle2, particle3, values);
    }
  }

  /**
   * Replace an angle using a caller-owned native array borrowed for this call.
   *
   * @param index stored angle index
   * @param particle1 first particle index
   * @param particle2 central particle index
   * @param particle3 third particle index
   * @param parameters replacement per-angle values in declaration order
   */
  public void setAngleParameters(int index, int particle1, int particle2, int particle3,
      DoubleArray parameters) {
    OpenMMNative.OpenMM_CustomAngleForce_setAngleParameters(
        getPointer(), index, particle1, particle2, particle3, parameters.getPointer());
  }

  /** Replace the energy expression; existing contexts require recreation to observe expression changes.
   * @param energy OpenMM custom expression, with {@code theta} in radians */
  public void setEnergyFunction(String energy) {
    OpenMMStrings.withUtf8String(energy,
        value -> OpenMMNative.OpenMM_CustomAngleForce_setEnergyFunction(getPointer(), value));
  }

  /** Set the default used for subsequently created contexts; this does not change an existing
   * context's current parameter value.
   * @param index global-parameter index
   * @param value new default value */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_CustomAngleForce_setGlobalParameterDefaultValue(getPointer(), index, value);
  }

  /** Rename a declared global parameter; existing contexts are not reconfigured by this setter.
   * @param index global-parameter index
   * @param name new expression parameter name */
  public void setGlobalParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomAngleForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /** Rename a declared per-angle parameter; existing contexts are not reconfigured by this setter.
   * @param index per-angle parameter index
   * @param name new expression parameter name */
  public void setPerAngleParameterName(int index, String name) {
    OpenMMStrings.withUtf8String(name,
        value -> OpenMMNative.OpenMM_CustomAngleForce_setPerAngleParameterName(
            getPointer(), index, value));
  }

  /**
   * Set whether the angle's particle displacements use periodic boundary conditions.
   *
   * @param periodic {@code true} to use the context periodic box, otherwise use unwrapped
   *     coordinate displacements
   */
  public void setUsesPeriodicBoundaryConditions(boolean periodic) {
    OpenMMNative.OpenMM_CustomAngleForce_setUsesPeriodicBoundaryConditions(
        getPointer(), OpenMMBooleans.toNative(periodic));
  }

  /**
   * Apply supported per-angle parameter changes to an existing context.
   *
   * @param context context containing this force; if it has no native context handle this helper
   *     intentionally does nothing. OpenMM restricts this update to modifiable angle parameters,
   *     not changes to the expression, declarations, or topology.
   */
  public void updateParametersInContext(Context context) {
    CustomForceParameters.updateContext(context,
        value -> OpenMMNative.OpenMM_CustomAngleForce_updateParametersInContext(getPointer(), value));
  }

  /** @return {@code true} if this force uses periodic displacements. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomAngleForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static java.lang.foreign.MemorySegment create(String energy) {
    return OpenMMStrings.withUtf8StringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomAngleForce_create(value);
    });
  }
}
