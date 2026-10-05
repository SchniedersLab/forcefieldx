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
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Custom energy depending on periodic box volume and box-vector components.
 *
 * <p>The expression may reference {@code v} (box volume in nm^3), {@code ax}, {@code bx},
 * {@code by}, {@code cx}, {@code cy}, and {@code cz} (periodic box-vector components in nm), and
 * global parameters. The first box vector is {@code (ax,0,0)} and the second is {@code (bx,by,0)}.
 * The force contributes energy (kJ/mol) but applies no particle forces because it has no
 * particle-position dependence; it is intended for volume/box-dependent terms such as pressure
 * matching. Parameter units are not converted.
 *
 * <p>This force always reports periodic-boundary-condition use. Java strings are passed through
 * temporary UTF-8 storage and returned strings are copied. {@link MemorySegment} string overloads
 * borrow a caller-owned NUL-terminated UTF-8 segment only for the call.
 */
public class CustomVolumeForce extends Force {

  /**
   * Create a force from an expression of box volume, box-vector components, and global parameters.
   *
   * @param energy custom expression; {@code v} is nm^3 and box components are nm
   */
  public CustomVolumeForce(String energy) {
    super(create(energy));
  }

  /** Declare an expression-wide parameter with a default for newly created contexts.
   * @param name expression parameter name
   * @param defaultValue initial value, in expression-consistent units
   * @return index assigned to the parameter */
  public int addGlobalParameter(String name, double defaultValue) {
    return withStringResult(name, value ->
        OpenMMNative.OpenMM_CustomVolumeForce_addGlobalParameter(
            getPointer(), value, defaultValue));
  }

  /** Declare a global parameter from a caller-owned NUL-terminated UTF-8 name segment.
   * @param name name segment borrowed during the call
   * @param defaultValue initial value
   * @return index assigned to the parameter */
  public int addGlobalParameter(MemorySegment name, double defaultValue) {
    return OpenMMNative.OpenMM_CustomVolumeForce_addGlobalParameter(
        getPointer(), name, defaultValue);
  }

  /** Destroy this force. */
  /** Destroy the native force and release its owned resources. */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CustomVolumeForce_destroy);
  }

  /** @return copied current energy expression */
  public String getEnergyFunction() {
    return OpenMMStrings.copy(OpenMMNative.OpenMM_CustomVolumeForce_getEnergyFunction(getPointer()));
  }

  /** @param index global-parameter index
   *  @return default for newly created contexts */
  public double getGlobalParameterDefaultValue(int index) {
    return OpenMMNative.OpenMM_CustomVolumeForce_getGlobalParameterDefaultValue(
        getPointer(), index);
  }

  /** @param index global-parameter index
   *  @return copied parameter name */
  public String getGlobalParameterName(int index) {
    return OpenMMStrings.copy(
        OpenMMNative.OpenMM_CustomVolumeForce_getGlobalParameterName(getPointer(), index));
  }

  /** @return number of declared global parameters */
  public int getNumGlobalParameters() {
    return OpenMMNative.OpenMM_CustomVolumeForce_getNumGlobalParameters(getPointer());
  }

  /** Replace the expression; recreate contexts to use it.
   * @param energy custom expression in terms of box variables and global parameters */
  public void setEnergyFunction(String energy) {
    withString(energy, value ->
        OpenMMNative.OpenMM_CustomVolumeForce_setEnergyFunction(getPointer(), value));
  }

  /** Replace the expression from a caller-owned NUL-terminated UTF-8 segment borrowed for the call.
   * @param energy live native expression segment */
  public void setEnergyFunction(MemorySegment energy) {
    OpenMMNative.OpenMM_CustomVolumeForce_setEnergyFunction(getPointer(), energy);
  }

  /** Set the default for future contexts; existing context values are unchanged.
   * @param index global-parameter index
   * @param value new default */
  public void setGlobalParameterDefaultValue(int index, double value) {
    OpenMMNative.OpenMM_CustomVolumeForce_setGlobalParameterDefaultValue(
        getPointer(), index, value);
  }

  /** Rename a global parameter; existing contexts are not reconfigured.
   * @param index global-parameter index
   * @param name new expression name */
  public void setGlobalParameterName(int index, String name) {
    withString(name, value ->
        OpenMMNative.OpenMM_CustomVolumeForce_setGlobalParameterName(
            getPointer(), index, value));
  }

  /** Rename a parameter from a caller-owned NUL-terminated UTF-8 segment.
   * @param index global-parameter index
   * @param name segment borrowed during the call */
  public void setGlobalParameterName(int index, MemorySegment name) {
    OpenMMNative.OpenMM_CustomVolumeForce_setGlobalParameterName(getPointer(), index, name);
  }

  /** @return always {@code true}; this force depends on the periodic box. */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CustomVolumeForce_usesPeriodicBoundaryConditions(getPointer()));
  }

  private static MemorySegment create(String energy) {
    return withStringResult(energy, value -> {
      OpenMMRuntime.initialize();
      return OpenMMNative.OpenMM_CustomVolumeForce_create(value);
    });
  }

  private static void withString(String value, Consumer<MemorySegment> action) {
    OpenMMStrings.withUtf8String(value, action);
  }

  private static <T> T withStringResult(String value, Function<MemorySegment, T> action) {
    return OpenMMStrings.withUtf8StringResult(value, action);
  }
}
