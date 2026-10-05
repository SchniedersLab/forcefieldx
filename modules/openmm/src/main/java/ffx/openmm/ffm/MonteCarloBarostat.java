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
import java.util.Objects;

/**
 * Adjusts the size of the periodic box with Monte Carlo moves to simulate constant pressure, scaling the
 * whole box isotropically.
 *
 * <p>The barostat assumes the simulation is run at constant temperature and needs that temperature because it
 * affects the Monte Carlo acceptance probability, but it does not regulate temperature itself; use another
 * mechanism such as a Langevin integrator or an {@link AndersenThermostat}.</p>
 *
 * <p>Context-update scope: the default pressure and temperature (and any surface tension) apply to contexts
 * created afterwards, not to ones that already exist, as stated by the OpenMM header. The header has no {@code
 * updateParametersInContext} for this force, and does not state whether other changes (frequency, seed, modes)
 * reach an existing context. Native context-parameter names ({@code Pressure()}, {@code Temperature()} and
 * similar) defined by the header are not exposed by this class.</p>
 *
 * <p>Ownership: this wrapper owns the native force until it is added to a {@link System}. OpenMM then
 * assumes native ownership, and the system invalidates this wrapper when the system is destroyed, so a
 * force added to a system should not also be destroyed directly.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.MonteCarloBarostat} is not. The C++ constructor defaults the
 * frequency to 25 time steps; this class always requires a value, as the JNA class does.</p>
 */
public class MonteCarloBarostat extends Force {
  /**
   * Create a Monte Carlo barostat. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * @param pressure    default pressure acting on the system, in bar.
   * @param temperature default temperature at which the system is maintained, in K.
   * @param frequency   frequency, in time steps, at which Monte Carlo pressure changes are attempted; 0 disables
   *     the barostat.
   */
  public MonteCarloBarostat(double pressure, double temperature, int frequency) {
    super(create(pressure, temperature, frequency));
  }

  /**
   * Compute the instantaneous pressure of a system to which this barostat is applied.
   *
   * <p>The pressure is computed from the molecular virial using a finite difference to obtain the derivative of
   * potential energy with respect to volume. For systems in equilibrium, its time average should equal the applied
   * pressure, but fluctuations can be very large.</p>
   *
   * @param context live context of the system; must not be null.
   * @return instantaneous pressure, in bar (the header states no unit for this return value; bar is the unit of the
   *     barostat pressure).
   * @throws NullPointerException if {@code context} is null.
   */
  public double computeCurrentPressure(Context context) {
    return OpenMMNative.OpenMM_MonteCarloBarostat_computeCurrentPressure(
        getPointer(), Objects.requireNonNull(context, "Context cannot be null.").getPointer());
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_MonteCarloBarostat_destroy);
  }

  /**
   * Get the default pressure acting on the system.
   *
   * @return default pressure, in bar.
   */
  public double getDefaultPressure() {
    return OpenMMNative.OpenMM_MonteCarloBarostat_getDefaultPressure(getPointer());
  }

  /**
   * Get the default temperature at which the system is maintained.
   *
   * @return default temperature, in K.
   */
  public double getDefaultTemperature() {
    return OpenMMNative.OpenMM_MonteCarloBarostat_getDefaultTemperature(getPointer());
  }

  /**
   * Get the frequency at which Monte Carlo pressure changes are attempted.
   *
   * @return frequency in time steps; 0 means the barostat is disabled.
   */
  public int getFrequency() {
    return OpenMMNative.OpenMM_MonteCarloBarostat_getFrequency(getPointer());
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)}.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen at context creation.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_MonteCarloBarostat_getRandomNumberSeed(getPointer());
  }

  /**
   * Set the default pressure acting on the system. This affects contexts created afterwards, not ones that
   * already exist.
   *
   * @param pressure default pressure, in bar.
   */
  public void setDefaultPressure(double pressure) {
    OpenMMNative.OpenMM_MonteCarloBarostat_setDefaultPressure(getPointer(), pressure);
  }

  /**
   * Set the default temperature at which the system is maintained. This affects contexts created afterwards,
   * not ones that already exist.
   *
   * @param temperature default temperature, in K.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMMNative.OpenMM_MonteCarloBarostat_setDefaultTemperature(getPointer(), temperature);
  }

  /**
   * Set the frequency at which Monte Carlo pressure changes are attempted.
   *
   * @param frequency frequency in time steps; 0 disables the barostat.
   */
  public void setFrequency(int frequency) {
    OpenMMNative.OpenMM_MonteCarloBarostat_setFrequency(getPointer(), frequency);
  }

  /**
   * Set the random-number seed. Different seeds are guaranteed to give different Monte Carlo move sequences,
   * but equal seeds carry no guarantee because platforms may use non-deterministic algorithms. A seed of 0 (the
   * default) causes a unique seed to be chosen when a context is created from this force.
   *
   * @param seed random-number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_MonteCarloBarostat_setRandomNumberSeed(getPointer(), seed);
  }

  /**
   * Determine whether this force uses periodic boundary conditions, as reported by OpenMM.
   *
   * <p><b>Header mismatch:</b> although a barostat needs a periodic box, the OpenMM header's inline implementation
   * returns false, and this method converts that native value through {@link OpenMMBooleans#fromNative(int)}. It
   * therefore returns false, contrary to the earlier statement in this class that it returns true.</p>
   *
   * @return false, as reported by OpenMM for this force.
   */
  @Override
  public boolean usesPeriodicBoundaryConditions() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_MonteCarloBarostat_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @param pressure    default pressure, in bar.
   * @param temperature default temperature, in K.
   * @param frequency   frequency in time steps; 0 disables the barostat.
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create(double pressure, double temperature, int frequency) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_MonteCarloBarostat_create(pressure, temperature, frequency);
  }
}
