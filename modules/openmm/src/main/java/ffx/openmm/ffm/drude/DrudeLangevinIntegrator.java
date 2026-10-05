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
package ffx.openmm.ffm.drude;

import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * Langevin integrator for systems containing Drude particles, with separate thermostats for
 * ordinary and Drude internal coordinates.
 *
 * <p>The ordinary-particle thermostat also acts on the center-of-mass motion of each Drude pair;
 * the second thermostat acts on the pairs' relative internal motion and is typically set to a
 * lower temperature. The optional hard wall limits the distance between a Drude particle and its
 * parent; it defaults to 0.02 nm and is disabled by setting the distance to zero. The System must
 * contain a DrudeForce so the integrator can identify Drude particles.</p>
 */
public class DrudeLangevinIntegrator extends DrudeIntegrator {

  /**
   * Create a Drude Langevin integrator.
   *
   * @param stepSize         Integration time step in picoseconds.
   * @param temperature      Target temperature of the main heat bath in kelvin.
   * @param friction         Friction coefficient coupling ordinary degrees of freedom to the main bath,
   *                         in ps<sup>-1</sup>.
   * @param drudeTemperature Target temperature of the heat bath for Drude internal coordinates,
   *                         in kelvin.
   * @param drudeFriction    Friction coefficient coupling Drude internal coordinates to their heat
   *                         bath, in ps<sup>-1</sup>.
   */
  public DrudeLangevinIntegrator(
      double stepSize, double temperature, double friction,
      double drudeTemperature, double drudeFriction) {
    super(create(stepSize, temperature, friction, drudeTemperature, drudeFriction));
  }

  /**
   * Release the native integrator handle.
   *
   * <p>The integrator must not be used after it is destroyed.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_DrudeLangevinIntegrator_destroy);
  }

  /**
   * Compute the instantaneous temperature of Drude internal coordinates.
   *
   * <p>This is calculated from the kinetic energy of the relative internal motion of Drude pairs
   * and should remain close, on average, to the configured Drude temperature.</p>
   *
   * @return Instantaneous Drude internal-coordinate temperature in kelvin.
   */
  public double computeDrudeTemperature() {
    return OpenMMNative.OpenMM_DrudeLangevinIntegrator_computeDrudeTemperature(getPointer());
  }

  /**
   * Compute the instantaneous temperature of the ordinary system degrees of freedom.
   *
   * <p>This includes ordinary particles and the center-of-mass motion of Drude pairs, but excludes
   * their relative internal motion. On average it should be approximately the temperature returned
   * by {@link #getTemperature()}.</p>
   *
   * @return Instantaneous system temperature in kelvin.
   */
  public double computeSystemTemperature() {
    return OpenMMNative.OpenMM_DrudeLangevinIntegrator_computeSystemTemperature(getPointer());
  }

  /**
   * Get the friction coefficient coupling Drude internal coordinates to their heat bath.
   *
   * @return Drude friction coefficient in ps<sup>-1</sup>.
   */
  public double getDrudeFriction() {
    return OpenMMNative.OpenMM_DrudeLangevinIntegrator_getDrudeFriction(getPointer());
  }

  /**
   * Get the friction coefficient coupling ordinary degrees of freedom to the main heat bath.
   *
   * @return Main-bath friction coefficient in ps<sup>-1</sup>.
   */
  public double getFriction() {
    return OpenMMNative.OpenMM_DrudeLangevinIntegrator_getFriction(getPointer());
  }

  /**
   * Get the target temperature of the main heat bath.
   *
   * @return Main-bath temperature in kelvin.
   */
  public double getTemperature() {
    return OpenMMNative.OpenMM_DrudeLangevinIntegrator_getTemperature(getPointer());
  }

  /**
   * Set the friction coefficient coupling Drude internal coordinates to their heat bath.
   *
   * @param friction Drude friction coefficient in ps<sup>-1</sup>.
   */
  public void setDrudeFriction(double friction) {
    OpenMMNative.OpenMM_DrudeLangevinIntegrator_setDrudeFriction(getPointer(), friction);
  }

  /**
   * Set the friction coefficient coupling ordinary degrees of freedom to the main heat bath.
   *
   * @param friction Main-bath friction coefficient in ps<sup>-1</sup>.
   */
  public void setFriction(double friction) {
    OpenMMNative.OpenMM_DrudeLangevinIntegrator_setFriction(getPointer(), friction);
  }

  /**
   * Set the target temperature of the main heat bath.
   *
   * @param temperature Main-bath temperature in kelvin.
   */
  public void setTemperature(double temperature) {
    OpenMMNative.OpenMM_DrudeLangevinIntegrator_setTemperature(getPointer(), temperature);
  }

  /**
   * Advance the simulation by the requested number of time steps.
   *
   * @param steps Number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_DrudeLangevinIntegrator_step(getPointer(), steps);
  }

  private static MemorySegment create(
      double stepSize, double temperature, double friction,
      double drudeTemperature, double drudeFriction) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_DrudeLangevinIntegrator_create(
        temperature, friction, drudeTemperature, drudeFriction, stepSize);
  }
}
