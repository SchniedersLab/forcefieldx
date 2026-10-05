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
 * Simulates a {@link System} using Brownian dynamics with a fixed step size.
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>FFM/JNA differences: this class is final, whereas the JNA {@link ffx.openmm.BrownianIntegrator} is not.</p>
 */
public class BrownianIntegrator extends Integrator {

  /**
   * Create a Brownian integrator. The FFM runtime is initialized first through {@link OpenMMRuntime#initialize()}.
   *
   * @param temperature   temperature of the heat bath, in K.
   * @param frictionCoeff friction coefficient coupling the system to the heat bath, in 1/ps.
   * @param stepSize      step size with which to integrate the system, in ps.
   */
  public BrownianIntegrator(double temperature, double frictionCoeff, double stepSize) {
    super(create(temperature, frictionCoeff, stepSize));
  }

  /**
   * Destroy the native integrator.
   *
   * <p>The native handle is released once and this wrapper is invalidated; repeated calls have no effect. If the
   * integrator was passed to a {@link Context}, that context's destruction also destroys the integrator, so do not call
   * this for an integrator that a live context owns.</p>
   */
  @Override
  public void destroy() {
    destroy(OpenMMNative::OpenMM_BrownianIntegrator_destroy);
  }

  /**
   * Get the friction coefficient that determines how strongly the system is coupled to the heat bath.
   *
   * @return friction coefficient, in 1/ps.
   */
  public double getFriction() {
    return OpenMMNative.OpenMM_BrownianIntegrator_getFriction(getPointer());
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)}.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen at context creation.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_BrownianIntegrator_getRandomNumberSeed(getPointer());
  }

  /**
   * Get the temperature of the heat bath.
   *
   * @return heat-bath temperature, in K.
   */
  public double getTemperature() {
    return OpenMMNative.OpenMM_BrownianIntegrator_getTemperature(getPointer());
  }

  /**
   * Set the friction coefficient that determines how strongly the system is coupled to the heat bath.
   *
   * @param coeff friction coefficient, in 1/ps.
   */
  public void setFriction(double coeff) {
    OpenMMNative.OpenMM_BrownianIntegrator_setFriction(getPointer(), coeff);
  }

  /**
   * Set the random-number seed. The meaning is left to each OpenMM platform. Different seeds are guaranteed to
   * give different random-force sequences, but equal seeds carry no guarantee because platforms may use
   * non-deterministic algorithms. A seed of 0 (the default) causes a unique seed to be chosen when a context is
   * created from this integrator.
   *
   * @param seed random-number seed.
   */
  public void setRandomNumberSeed(int seed) {
    OpenMMNative.OpenMM_BrownianIntegrator_setRandomNumberSeed(getPointer(), seed);
  }

  /**
   * Set the temperature of the heat bath.
   *
   * @param temp heat-bath temperature, in K.
   */
  public void setTemperature(double temp) {
    OpenMMNative.OpenMM_BrownianIntegrator_setTemperature(getPointer(), temp);
  }

  /**
   * Advance the simulation by a series of fixed time steps. This overrides {@link Integrator#step(int)} with the
   * Brownian native step call.
   *
   * @param steps number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_BrownianIntegrator_step(getPointer(), steps);
  }

  /**
   * Create the native integrator after loading the FFM runtime.
   *
   * @param temperature   heat-bath temperature, in K.
   * @param frictionCoeff friction coefficient, in 1/ps.
   * @param stepSize      step size, in ps.
   * @return native integrator handle owned by the new wrapper.
   */
  private static MemorySegment create(double temperature, double frictionCoeff, double stepSize) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, stepSize);
  }
}
