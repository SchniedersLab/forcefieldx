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
import java.util.Objects;

/**
 * Adjusts the periodic box with Monte Carlo moves to simulate constant pressure. Unlike {@link
 * MonteCarloBarostat}, each move scales only one axis, so the box shape as well as its size may change.
 * A different pressure can be given for each axis, and individual axes can be kept fixed.
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
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.MonteCarloAnisotropicBarostat} is not. The JNA class
 * uses {@code OpenMM_Vec3} for the pressure and integers for the scale flags and getters; this class uses {@link
 * Vec3} and booleans. The C++ constructor defaults the three scale flags to true and the frequency to 25; this class
 * always requires them. The header's {@code computeCurrentPressure(Context)} (a per-axis pressure vector) is not
 * mapped by this class, nor by the JNA class.</p>
 */
public class MonteCarloAnisotropicBarostat extends Force {
  /**
   * Create an anisotropic Monte Carlo barostat. The pressure vector is converted to temporary native memory, so
   * the caller keeps ownership of {@code pressure}.
   *
   * @param pressure    default pressure acting along the X, Y and Z axes, in bar; must not be null.
   * @param temperature default temperature at which the system is maintained, in K.
   * @param scaleX      whether the X dimension of the periodic box may change size.
   * @param scaleY      whether the Y dimension of the periodic box may change size.
   * @param scaleZ      whether the Z dimension of the periodic box may change size.
   * @param frequency   frequency, in time steps, at which Monte Carlo pressure changes are attempted; 0 disables
   *     the barostat.
   * @throws NullPointerException if {@code pressure} is null.
   */
  public MonteCarloAnisotropicBarostat(
      Vec3 pressure, double temperature, boolean scaleX, boolean scaleY, boolean scaleZ,
      int frequency) {
    super(create(pressure, temperature, scaleX, scaleY, scaleZ, frequency));
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_MonteCarloAnisotropicBarostat_destroy);
  }

  /**
   * Get the default pressure along each axis.
   *
   * @return a copied {@link Vec3} of the default X, Y and Z pressures, in bar.
   */
  public Vec3 getDefaultPressure() {
    return Vec3.fromNative(
        OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getDefaultPressure(getPointer()));
  }

  /**
   * Get the default temperature at which the system is maintained.
   *
   * @return default temperature, in K.
   */
  public double getDefaultTemperature() {
    return OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getDefaultTemperature(getPointer());
  }

  /**
   * Get the frequency at which Monte Carlo pressure changes are attempted.
   *
   * @return frequency in time steps; 0 means the barostat is disabled.
   */
  public int getFrequency() {
    return OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getFrequency(getPointer());
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)}.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen at context creation.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getRandomNumberSeed(getPointer());
  }

  /**
   * Get whether the X dimension of the periodic box may change size.
   *
   * @return true if X may change; the JNA getter returns an int.
   */
  public boolean getScaleX() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getScaleX(getPointer()));
  }

  /**
   * Get whether the Y dimension of the periodic box may change size.
   *
   * @return true if Y may change; the JNA getter returns an int.
   */
  public boolean getScaleY() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getScaleY(getPointer()));
  }

  /**
   * Get whether the Z dimension of the periodic box may change size.
   *
   * @return true if Z may change; the JNA getter returns an int.
   */
  public boolean getScaleZ() {
    return OpenMMBooleans.fromNative(
        OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_getScaleZ(getPointer()));
  }

  /**
   * Set the default pressure along each axis. This affects contexts created afterwards, not ones that already
   * exist.
   *
   * @param pressure X, Y and Z pressures, in bar; must not be null.
   * @throws NullPointerException if {@code pressure} is null.
   */
  public void setDefaultPressure(Vec3 pressure) {
    try (Arena arena = Arena.ofConfined()) {
      OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_setDefaultPressure(
          getPointer(), Objects.requireNonNull(pressure, "Pressure cannot be null.").toNative(arena));
    }
  }

  /**
   * Set the default temperature at which the system is maintained. This affects contexts created afterwards,
   * not ones that already exist.
   *
   * @param temperature default temperature, in K.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_setDefaultTemperature(
        getPointer(), temperature);
  }

  /**
   * Set the frequency at which Monte Carlo pressure changes are attempted.
   *
   * @param frequency frequency in time steps; 0 disables the barostat.
   */
  public void setFrequency(int frequency) {
    OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_setFrequency(getPointer(), frequency);
  }

  /**
   * Set the random-number seed. Different seeds are guaranteed to give different Monte Carlo move sequences,
   * but equal seeds carry no guarantee because platforms may use non-deterministic algorithms. A seed of 0 (the
   * default) causes a unique seed to be chosen when a context is created from this force.
   *
   * @param seed random-number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_setRandomNumberSeed(getPointer(), seed);
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
        OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_usesPeriodicBoundaryConditions(
            getPointer()));
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @param pressure    default X, Y and Z pressures, in bar.
   * @param temperature default temperature, in K.
   * @param scaleX      whether X may change.
   * @param scaleY      whether Y may change.
   * @param scaleZ      whether Z may change.
   * @param frequency   frequency in time steps; 0 disables the barostat.
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create(
      Vec3 pressure, double temperature, boolean scaleX, boolean scaleY, boolean scaleZ,
      int frequency) {
    OpenMMRuntime.initialize();
    try (Arena arena = Arena.ofConfined()) {
      return OpenMMNative.OpenMM_MonteCarloAnisotropicBarostat_create(
          Objects.requireNonNull(pressure, "Pressure cannot be null.").toNative(arena),
          temperature, OpenMMBooleans.toNative(scaleX), OpenMMBooleans.toNative(scaleY),
          OpenMMBooleans.toNative(scaleZ), frequency);
    }
  }
}
