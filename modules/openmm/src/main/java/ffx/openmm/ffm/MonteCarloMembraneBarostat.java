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
 * A Monte Carlo barostat designed for membrane simulations, assuming the membrane lies in the XY plane. The
 * acceptance criterion includes an isotropic-pressure term depending on box volume and a surface-tension term
 * depending on the XY cross-sectional area. Pressure and surface tension have opposite senses: a larger pressure
 * tends to shrink the box, a larger surface tension tends to enlarge it.
 *
 * <p>{@link XYMode} chooses whether X and Y scale together or independently, and {@link ZMode} whether Z varies
 * freely, is fixed, or varies inversely to X and Y so the volume is constant. In {@link ZMode#CONSTANT_VOLUME}
 * pressure has no effect, only surface tension.</p>
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
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.MonteCarloMembraneBarostat} is not. The JNA class
 * uses integer mode arguments and getters; this class uses the {@link XYMode} and {@link ZMode} enums, translating
 * them through native constant accessors. The C++ constructor defaults the frequency to 25; this class always
 * requires it. The header's {@code computeCurrentPressure(Context)} is not mapped by this class or the JNA class.
 * The header documents surface tension in bar*nm, but the parameter description of {@code setDefaultSurfaceTension}
 * says "bar"; this class uses bar*nm.</p>
 */
public class MonteCarloMembraneBarostat extends Force {
  /**
   * How the X and Y box axes may change. Constants map to OpenMM's {@code XYIsotropic} (0) and {@code
   * XYAnisotropic} (1).
   */
  public enum XYMode {ISOTROPIC, ANISOTROPIC}

  /**
   * How the Z box axis may change. Constants map to OpenMM's {@code ZFree} (0), {@code ZFixed} (1) and {@code
   * ConstantVolume} (2). Native values are obtained through native constant accessors, and an unrecognized native value
   * read by {@link #getZMode()} gives an {@link IllegalStateException}.
   */
  public enum ZMode {FREE, FIXED, CONSTANT_VOLUME}

  /**
   * Create a membrane Monte Carlo barostat. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   *
   * @param pressure       default pressure acting on the system, in bar.
   * @param surfaceTension default surface tension acting on the system, in bar*nm.
   * @param temperature    default temperature at which the system is maintained, in K.
   * @param xyMode         behavior of the X and Y axes; must not be null.
   * @param zMode          behavior of the Z axis; must not be null.
   * @param frequency      frequency, in time steps, at which Monte Carlo volume changes are attempted; 0 disables
   *     the barostat.
   * @throws NullPointerException if {@code xyMode} or {@code zMode} is null.
   */
  public MonteCarloMembraneBarostat(double pressure, double surfaceTension, double temperature,
                                    XYMode xyMode, ZMode zMode, int frequency) {
    super(create(pressure, surfaceTension, temperature, xyMode, zMode, frequency));
  }

  /**
   * Destroy the native force.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect.
   * Do not call this for a force whose ownership was transferred to a {@link System}.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_MonteCarloMembraneBarostat_destroy);
  }

  /**
   * Get the default pressure acting on the system.
   *
   * @return default pressure, in bar.
   */
  public double getDefaultPressure() {
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getDefaultPressure(getPointer());
  }

  /**
   * Get the default surface tension acting on the system.
   *
   * @return default surface tension, in bar*nm.
   */
  public double getDefaultSurfaceTension() {
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getDefaultSurfaceTension(getPointer());
  }

  /**
   * Get the default temperature at which the system is maintained.
   *
   * @return default temperature, in K.
   */
  public double getDefaultTemperature() {
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getDefaultTemperature(getPointer());
  }

  /**
   * Get the frequency at which Monte Carlo volume changes are attempted.
   *
   * @return frequency in time steps; 0 means the barostat is disabled.
   */
  public int getFrequency() {
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getFrequency(getPointer());
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)}.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen at context creation.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getRandomNumberSeed(getPointer());
  }

  /**
   * Get the behavior of the X and Y axes.
   *
   * @return {@link XYMode#ANISOTROPIC} if the native value equals the anisotropic constant, otherwise {@link
   *     XYMode#ISOTROPIC} (no unknown-value check is made); the JNA getter returns the native integer.
   */
  public XYMode getXYMode() {
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getXYMode(getPointer()) == OpenMMNative.OpenMM_MonteCarloMembraneBarostat_XYAnisotropic() ? XYMode.ANISOTROPIC : XYMode.ISOTROPIC;
  }

  /**
   * Get the behavior of the Z axis.
   *
   * @return the Z mode; the JNA getter returns the native integer.
   * @throws IllegalStateException if the native value is not a known Z mode.
   */
  public ZMode getZMode() {
    int mode = OpenMMNative.OpenMM_MonteCarloMembraneBarostat_getZMode(getPointer());
    if (mode == OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ZFree()) return ZMode.FREE;
    if (mode == OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ZFixed()) return ZMode.FIXED;
    if (mode == OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ConstantVolume()) return ZMode.CONSTANT_VOLUME;
    throw new IllegalStateException("Unknown OpenMM membrane barostat Z mode: " + mode);
  }

  /**
   * Set the default pressure acting on the system. This affects contexts created afterwards, not ones that
   * already exist.
   *
   * @param pressure default pressure, in bar.
   */
  public void setDefaultPressure(double pressure) {
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setDefaultPressure(getPointer(), pressure);
  }

  /**
   * Set the default surface tension acting on the system. This affects contexts created afterwards, not ones
   * that already exist.
   *
   * @param surfaceTension default surface tension, in bar*nm.
   */
  public void setDefaultSurfaceTension(double surfaceTension) {
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setDefaultSurfaceTension(getPointer(), surfaceTension);
  }

  /**
   * Set the default temperature at which the system is maintained. This affects contexts created afterwards,
   * not ones that already exist.
   *
   * @param temperature default temperature, in K.
   */
  public void setDefaultTemperature(double temperature) {
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setDefaultTemperature(getPointer(), temperature);
  }

  /**
   * Set the frequency at which Monte Carlo volume changes are attempted.
   *
   * @param frequency frequency in time steps; 0 disables the barostat.
   */
  public void setFrequency(int frequency) {
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setFrequency(getPointer(), frequency);
  }

  /**
   * Set the random-number seed. Different seeds are guaranteed to give different Monte Carlo move sequences,
   * but equal seeds carry no guarantee because platforms may use non-deterministic algorithms. A seed of 0 (the
   * default) causes a unique seed to be chosen when a context is created from this force.
   *
   * @param seed random-number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setRandomNumberSeed(getPointer(), seed);
  }

  /**
   * Set the behavior of the X and Y axes. The header does not state whether this reaches an existing context.
   *
   * @param mode XY mode; must not be null.
   * @throws NullPointerException if {@code mode} is null.
   */
  public void setXYMode(XYMode mode) {
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setXYMode(getPointer(), Objects.requireNonNull(mode, "XY mode cannot be null.") == XYMode.ANISOTROPIC ? OpenMMNative.OpenMM_MonteCarloMembraneBarostat_XYAnisotropic() : OpenMMNative.OpenMM_MonteCarloMembraneBarostat_XYIsotropic());
  }

  /**
   * Set the behavior of the Z axis. The header does not state whether this reaches an existing context.
   *
   * @param mode Z mode; must not be null.
   * @throws NullPointerException if {@code mode} is null.
   */
  public void setZMode(ZMode mode) {
    int nativeMode = switch (Objects.requireNonNull(mode, "Z mode cannot be null.")) {
      case FREE -> OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ZFree();
      case FIXED -> OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ZFixed();
      case CONSTANT_VOLUME -> OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ConstantVolume();
    };
    OpenMMNative.OpenMM_MonteCarloMembraneBarostat_setZMode(getPointer(), nativeMode);
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
    return OpenMMBooleans.fromNative(OpenMMNative.OpenMM_MonteCarloMembraneBarostat_usesPeriodicBoundaryConditions(getPointer()));
  }

  /**
   * Create the native force after loading the FFM runtime.
   *
   * @param p  default pressure, in bar.
   * @param s  default surface tension, in bar*nm.
   * @param t  default temperature, in K.
   * @param xy XY mode; must not be null.
   * @param z  Z mode; must not be null.
   * @param f  frequency in time steps; 0 disables the barostat.
   * @return native force handle owned by the new wrapper.
   */
  private static MemorySegment create(double p, double s, double t, XYMode xy, ZMode z, int f) {
    OpenMMRuntime.initialize();
    int xyValue = Objects.requireNonNull(xy, "XY mode cannot be null.") == XYMode.ANISOTROPIC
        ? OpenMMNative.OpenMM_MonteCarloMembraneBarostat_XYAnisotropic()
        : OpenMMNative.OpenMM_MonteCarloMembraneBarostat_XYIsotropic();
    int zValue = switch (Objects.requireNonNull(z, "Z mode cannot be null.")) {
      case FREE -> OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ZFree();
      case FIXED -> OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ZFixed();
      case CONSTANT_VOLUME -> OpenMMNative.OpenMM_MonteCarloMembraneBarostat_ConstantVolume();
    };
    return OpenMMNative.OpenMM_MonteCarloMembraneBarostat_create(p, s, t, xyValue, zValue, f);
  }
}
