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

/**
 * Immutable snapshot of selected simulation data from an OpenMM {@link Context}.
 *
 * <p>Instances own the native state handle and must be closed or destroyed. Vector data is copied
 * into packed {@code x,y,z} arrays, so the returned arrays remain valid after this state is
 * destroyed. Only data requested when the state was created is available; querying omitted
 * positions, velocities, forces, energies, parameters, or parameter derivatives is an OpenMM
 * error. Vector results are independent Java arrays, while parameter-map accessors expose borrowed
 * opaque native handles rather than Java maps.</p>
 */
public class State extends OpenMMHandle {

  /**
   * Wrap a native state handle and assume responsibility for releasing it.
   *
   * @param pointer native state handle whose lifetime is transferred to this wrapper.
   */
  public State(MemorySegment pointer) {
    super(pointer);
    OpenMMRuntime.initialize();
  }

  /**
   * Release this native state snapshot. Repeated calls have no effect.
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_State_destroy);
  }

  /**
   * Get the sum of OpenMM data-type flags identifying the data included in this state.
   *
   * @return bitwise combination of OpenMM state-data flags.
   */
  public int getDataTypes() {
    return OpenMMNative.OpenMM_State_getDataTypes(getPointer());
  }

  /**
   * Get the force acting on each particle as packed {@code x,y,z} components in
   * kJ/(mol nm).
   *
   * <p>Forces must have been requested when the state was created.</p>
   *
   * @return independent array with three force components per particle.
   */
  public double[] getForces() {
    return copyVec3Array(OpenMMNative.OpenMM_State_getForces(getPointer()));
  }

  /**
   * Get total kinetic energy in kJ/mol.
   *
   * <p>Energy data must have been requested when the state was created. OpenMM computes this using
   * the associated integrator's velocity convention.</p>
   *
   * @return kinetic energy in kJ/mol.
   */
  public double getKineticEnergy() {
    return OpenMMNative.OpenMM_State_getKineticEnergy(getPointer());
  }

  /**
   * Get derivatives of potential energy with respect to adjustable context parameters.
   *
   * <p>This returns the borrowed opaque C-wrapper handle for OpenMM's parameter map, not a Java
   * {@code Map}. The state must have been created with parameter derivatives requested, and the
   * handle becomes invalid when this state is destroyed.</p>
   *
   * @return borrowed {@code OpenMM_ParameterArray} handle, valid only while this state remains
   *         alive.
   */
  public MemorySegment getEnergyParameterDerivatives() {
    return OpenMMNative.OpenMM_State_getEnergyParameterDerivatives(getPointer());
  }

  /**
   * Get the adjustable context parameter values recorded in this state.
   *
   * <p>This returns the borrowed opaque C-wrapper handle for OpenMM's parameter map, not a Java
   * {@code Map}. The state must have been created with parameter values requested, and the handle
   * becomes invalid when this state is destroyed.</p>
   *
   * @return borrowed {@code OpenMM_ParameterArray} handle, valid only while this state remains
   *         alive.
   */
  public MemorySegment getParameters() {
    return OpenMMNative.OpenMM_State_getParameters(getPointer());
  }

  /**
   * Get the periodic box vectors, copied into Java values.
   *
   * @return copied periodic box vectors in nanometers.
   */
  public System.PeriodicBoxVectors getPeriodicBoxVectors() {
    try (Arena arena = Arena.ofConfined()) {
      MemorySegment a = new Vec3(0.0, 0.0, 0.0).toNative(arena);
      MemorySegment b = new Vec3(0.0, 0.0, 0.0).toNative(arena);
      MemorySegment c = new Vec3(0.0, 0.0, 0.0).toNative(arena);
      OpenMMNative.OpenMM_State_getPeriodicBoxVectors(getPointer(), a, b, c);
      return new System.PeriodicBoxVectors(
          Vec3.fromNative(a), Vec3.fromNative(b), Vec3.fromNative(c));
    }
  }

  /**
   * Get the periodic box volume in cubic nanometers.
   *
   * @return periodic box volume in cubic nanometers.
   */
  public double getPeriodicBoxVolume() {
    return OpenMMNative.OpenMM_State_getPeriodicBoxVolume(getPointer());
  }

  /**
   * Get each particle's position as packed {@code x,y,z} components in nanometers.
   *
   * <p>Positions must have been requested when the state was created.</p>
   *
   * @return independent array with three position components per particle.
   */
  public double[] getPositions() {
    return copyVec3Array(OpenMMNative.OpenMM_State_getPositions(getPointer()));
  }

  /**
   * Get total potential energy in kJ/mol.
   *
   * <p>Energy data must have been requested when the state was created.</p>
   *
   * @return potential energy in kJ/mol.
   */
  public double getPotentialEnergy() {
    return OpenMMNative.OpenMM_State_getPotentialEnergy(getPointer());
  }

  /**
   * Get the integration step count recorded in this state.
   *
   * @return number of completed integration steps.
   */
  public long getStepCount() {
    return OpenMMNative.OpenMM_State_getStepCount(getPointer());
  }

  /**
   * Get simulation time in picoseconds.
   *
   * @return simulation time in picoseconds.
   */
  public double getTime() {
    return OpenMMNative.OpenMM_State_getTime(getPointer());
  }

  /**
   * Get each particle's velocity as packed {@code x,y,z} components in nanometers per picosecond.
   *
   * <p>Velocities must have been requested when the state was created.</p>
   *
   * @return independent array with three velocity components per particle.
   */
  public double[] getVelocities() {
    return copyVec3Array(OpenMMNative.OpenMM_State_getVelocities(getPointer()));
  }

  /**
   * Copy a borrowed native vector array.
   *
   * @param array native array owned by the state.
   * @return packed Java vector components.
   */
  private static double[] copyVec3Array(MemorySegment array) {
    int size = OpenMMNative.OpenMM_Vec3Array_getSize(array);
    double[] values = new double[size * 3];
    for (int index = 0; index < size; index++) {
      Vec3 vector = Vec3.fromNative(OpenMMNative.OpenMM_Vec3Array_get(array, index));
      int offset = index * 3;
      values[offset] = vector.x();
      values[offset + 1] = vector.y();
      values[offset + 2] = vector.z();
    }
    return values;
  }
}
