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
 * Prevents the center of mass of a {@link System} from drifting. At each time step it calculates the
 * center-of-mass momentum and adjusts individual particle velocities to make it zero.
 *
 * <p>The frequency (in time steps) at which removal occurs is set at construction and can be changed with {@link
 * #setFrequency(int)}. The OpenMM header does not state whether a frequency change affects an existing {@link
 * Context}, and this class provides no {@code updateParametersInContext}.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.CMMotionRemover} is not. The OpenMM
 * C++ constructor defaults the frequency to 1; this class, like the JNA one, always requires the value.</p>
 */
public class CMMotionRemover extends Force {
  /**
   * Create a center-of-mass motion remover. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}, which the JNA counterpart does not do.
   *
   * @param frequency frequency, in time steps, at which center-of-mass motion is removed.
   */
  public CMMotionRemover(int frequency) {
    super(create(frequency));
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_CMMotionRemover_destroy);
  }

  /**
   * Get the frequency at which center-of-mass motion is removed.
   *
   * @return frequency, in time steps.
   */
  public int getFrequency() {
    return OpenMMNative.OpenMM_CMMotionRemover_getFrequency(getPointer());
  }

  /**
   * Set the frequency at which center-of-mass motion is removed.
   *
   * @param frequency frequency, in time steps.
   */
  public void setFrequency(int frequency) {
    OpenMMNative.OpenMM_CMMotionRemover_setFrequency(getPointer(), frequency);
  }

  /**
   * Determine whether this force uses periodic boundary conditions.
   *
   * <p>The OpenMM header's inline implementation of this query returns false, so this force does not use periodic
   * boundary conditions. The value is converted through {@link OpenMMBooleans#fromNative(int)}.</p>
   *
   * @return false for this force, as reported by OpenMM.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_CMMotionRemover_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @param frequency frequency, in time steps, at which center-of-mass motion is removed.
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create(int frequency) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_CMMotionRemover_create(frequency);
  }
}
