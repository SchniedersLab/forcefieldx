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

import ffx.openmm.ffm.Integrator;
import ffx.openmm.ffm.OpenMMRuntime;
import ffx.openmm.ffm.bindings.OpenMMNative;

import java.lang.foreign.MemorySegment;

/**
 * Base class for integrators that simulate systems containing Drude particles.
 *
 * <p>The maximum Drude-parent separation is enforced as a hard-wall constraint. Its documented
 * default is 0.02 nm; setting it to zero omits the constraint. Drude integrators use a
 * {@link ffx.openmm.ffm.System} containing a DrudeForce to identify Drude particles.</p>
 */
public class DrudeIntegrator extends Integrator {

  /**
   * Wrap an existing native Drude integrator handle and take ownership of it.
   *
   * <p>Calling {@link #destroy()} releases the native integrator; do not use this wrapper or the
   * handle after destruction.</p>
   *
   * @param pointer Native Drude integrator handle to own.
   */
  public DrudeIntegrator(MemorySegment pointer) {
    super(pointer);
  }

  /**
   * Create a base Drude integrator with the specified time step.
   *
   * @param stepSize Integration time step in picoseconds.
   */
  public DrudeIntegrator(double stepSize) {
    super(create(stepSize));
  }

  /**
   * Release the native integrator handle.
   *
   * <p>The integrator must not be used after it is destroyed.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_DrudeIntegrator_destroy);
  }

  /**
   * Get the heat-bath temperature assigned to Drude internal coordinates.
   *
   * @return Drude internal-coordinate heat-bath temperature in kelvin.
   */
  public double getDrudeTemperature() {
    return OpenMMNative.OpenMM_DrudeIntegrator_getDrudeTemperature(getPointer());
  }

  /**
   * Get the hard-wall limit on the distance between each Drude particle and its parent.
   *
   * @return Maximum Drude-parent distance in nm; zero means that the hard-wall constraint is
   * omitted. The documented default is 0.02 nm.
   */
  public double getMaxDrudeDistance() {
    return OpenMMNative.OpenMM_DrudeIntegrator_getMaxDrudeDistance(getPointer());
  }

  /**
   * Get the random-number seed configured for this integrator.
   *
   * @return Configured random-number seed; see {@link #setRandomNumberSeed(int)}.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_DrudeIntegrator_getRandomNumberSeed(getPointer());
  }

  /**
   * Set the heat-bath temperature assigned to Drude internal coordinates.
   *
   * @param temperature Drude internal-coordinate heat-bath temperature in kelvin.
   */
  public void setDrudeTemperature(double temperature) {
    OpenMMNative.OpenMM_DrudeIntegrator_setDrudeTemperature(getPointer(), temperature);
  }

  /**
   * Set the hard-wall limit on the distance between each Drude particle and its parent.
   *
   * <p>Setting the distance to zero omits the hard-wall constraint. The documented default is
   * 0.02 nm.</p>
   *
   * @param distance Maximum Drude-parent distance in nm, or zero to omit the constraint.
   */
  public void setMaxDrudeDistance(double distance) {
    OpenMMNative.OpenMM_DrudeIntegrator_setMaxDrudeDistance(getPointer(), distance);
  }

  /**
   * Set the seed used to generate random forces.
   *
   * <p>The precise interpretation is platform-dependent. Different seeds guarantee different
   * sequences of random forces, but using the same seed does not guarantee identical results:
   * platforms may use nondeterministic algorithms. A seed of zero (the default) causes a unique
   * seed to be selected when a Context is created.</p>
   *
   * @param seed Random-number seed; zero requests a unique seed per Context.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_DrudeIntegrator_setRandomNumberSeed(getPointer(), seed);
  }

  /**
   * Advance the simulation by the requested number of time steps.
   *
   * @param steps Number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_DrudeIntegrator_step(getPointer(), steps);
  }

  private static MemorySegment create(double stepSize) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_DrudeIntegrator_create(stepSize);
  }
}
