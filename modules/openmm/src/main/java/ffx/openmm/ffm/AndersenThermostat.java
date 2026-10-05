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
 * Uses the Andersen method to maintain constant temperature, with randomly occurring collisions that
 * reassign particle velocities from a heat bath.
 *
 * <p>The heat-bath temperature and collision frequency set here are defaults for contexts created afterwards;
 * the OpenMM header states that changing them affects new contexts but not ones that already exist. This class
 * does not expose the native context-parameter names ({@code Temperature()} and {@code CollisionFrequency()}) that the
 * OpenMM header defines for changing these values in a live context, and it provides no {@code
 * updateParametersInContext}.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.AndersenThermostat} is not.</p>
 */
public class AndersenThermostat extends Force {
  /**
   * Create an Andersen thermostat. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}, which the JNA counterpart does not do.
   *
   * @param defaultTemperature        default temperature of the heat bath, in K.
   * @param defaultCollisionFrequency default collision frequency, in 1/ps.
   */
  public AndersenThermostat(double defaultTemperature, double defaultCollisionFrequency) {
    super(create(defaultTemperature, defaultCollisionFrequency));
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_AndersenThermostat_destroy);
  }

  /**
   * Get the default collision frequency.
   *
   * @return default collision frequency, in 1/ps.
   */
  public double getDefaultCollisionFrequency() {
    return OpenMMNative.OpenMM_AndersenThermostat_getDefaultCollisionFrequency(getPointer());
  }

  /**
   * Get the default temperature of the heat bath.
   *
   * @return default heat-bath temperature, in K.
   */
  public double getDefaultTemperature() {
    return OpenMMNative.OpenMM_AndersenThermostat_getDefaultTemperature(getPointer());
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)} for its interpretation.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen when a context is created.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_AndersenThermostat_getRandomNumberSeed(getPointer());
  }

  /**
   * Set the default collision frequency. This affects contexts created afterwards, not ones that already exist.
   *
   * @param frequency default collision frequency, in 1/ps.
   */
  public void setDefaultCollisionFrequency(double frequency) {
    OpenMMNative.OpenMM_AndersenThermostat_setDefaultCollisionFrequency(getPointer(), frequency);
  }

  /**
   * Set the default temperature of the heat bath. This affects contexts created afterwards, not ones that
   * already exist.
   *
   * @param temperature default heat-bath temperature, in K.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMMNative.OpenMM_AndersenThermostat_setDefaultTemperature(getPointer(), temperature);
  }

  /**
   * Set the random-number seed.
   *
   * <p>The precise meaning is left to each OpenMM platform. Different seeds are guaranteed to give different
   * collision sequences, but no guarantee is made for equal seeds: platforms may use non-deterministic algorithms.
   * A seed of 0 (the default) causes a unique seed to be chosen when a context is created from this force.</p>
   *
   * @param seed random-number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_AndersenThermostat_setRandomNumberSeed(getPointer(), seed);
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
        OpenMMNative.OpenMM_AndersenThermostat_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @param temperature default heat-bath temperature, in K.
   * @param frequency   default collision frequency, in 1/ps.
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create(double temperature, double frequency) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_AndersenThermostat_create(temperature, frequency);
  }
}
