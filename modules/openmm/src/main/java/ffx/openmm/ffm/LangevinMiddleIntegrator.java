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
 * Simulates a {@link System} using Langevin dynamics with the LFMiddle discretization (J. Phys. Chem. A 2019,
 * 123, 28, 6056-6079), which tends to give more accurate configurational sampling than other discretizations. The
 * algorithm is closely related to BAOAB (Proc. R. Soc. A 472: 20160138): both give identical trajectories, but
 * LFMiddle returns half-step (leapfrog) velocities while BAOAB returns on-step velocities, and the former give a more
 * accurate sampling of the thermal ensemble.
 *
 * <p>Ownership: an integrator is bound to one {@link Context}, created by passing the integrator to a {@link
 * Context} constructor. Following the existing {@link Context} documentation, destroying that context also destroys
 * its integrator, so an integrator bound to a context should not also be destroyed directly. Until then this wrapper
 * owns the native integrator.</p>
 *
 * <p>Unlike most integrator wrappers in this package, this class is not final because {@link LangevinIntegrator}
 * extends it. The OpenMM C++ and C constructors take (temperature, frictionCoeff, stepSize); this Java constructor,
 * like the JNA class, takes (dt, temp, gamma), and reorders them when calling the native function. The JNA {@link
 * ffx.openmm.LangevinMiddleIntegrator} is also not final.</p>
 */
public class LangevinMiddleIntegrator extends Integrator {

  /**
   * Create an LFMiddle Langevin integrator. The FFM runtime is initialized first through {@link
   * OpenMMRuntime#initialize()}.
   *
   * <p>Note the argument order differs from the native constructor (temperature, friction, step size).</p>
   *
   * @param dt    step size with which to integrate the system, in ps.
   * @param temp  temperature of the heat bath, in K.
   * @param gamma friction coefficient coupling the system to the heat bath, in 1/ps.
   */
  public LangevinMiddleIntegrator(double dt, double temp, double gamma) {
    super(create(dt, temp, gamma));
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
    destroy(OpenMMNative::OpenMM_LangevinMiddleIntegrator_destroy);
  }

  /**
   * Get the friction coefficient that determines how strongly the system is coupled to the heat bath.
   *
   * @return friction coefficient, in 1/ps.
   */
  public double getFriction() {
    return OpenMMNative.OpenMM_LangevinMiddleIntegrator_getFriction(getPointer());
  }

  /**
   * Get the random-number seed. See {@link #setRandomNumberSeed(int)}.
   *
   * @return random-number seed; 0 (the default) means a unique seed is chosen at context creation.
   */
  public int getRandomNumberSeed() {
    return OpenMMNative.OpenMM_LangevinMiddleIntegrator_getRandomNumberSeed(getPointer());
  }

  /**
   * Get the temperature of the heat bath.
   *
   * @return heat-bath temperature, in K.
   */
  public double getTemperature() {
    return OpenMMNative.OpenMM_LangevinMiddleIntegrator_getTemperature(getPointer());
  }

  /**
   * Set the friction coefficient that determines how strongly the system is coupled to the heat bath.
   *
   * @param gamma friction coefficient, in 1/ps.
   */
  public void setFriction(double gamma) {
    OpenMMNative.OpenMM_LangevinMiddleIntegrator_setFriction(getPointer(), gamma);
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
    OpenMMNative.OpenMM_LangevinMiddleIntegrator_setRandomNumberSeed(getPointer(), seed);
  }

  /**
   * Set the temperature of the heat bath.
   *
   * @param temp heat-bath temperature, in K.
   */
  public void setTemperature(double temp) {
    OpenMMNative.OpenMM_LangevinMiddleIntegrator_setTemperature(getPointer(), temp);
  }

  /**
   * Advance the simulation by a series of fixed time steps. This overrides {@link Integrator#step(int)} with the
   * LangevinMiddle native step call.
   *
   * @param steps number of time steps to take.
   */
  @Override
  public void step(int steps) {
    OpenMMNative.OpenMM_LangevinMiddleIntegrator_step(getPointer(), steps);
  }

  /**
   * Create the native integrator after loading the FFM runtime.
   *
   * @param dt          step size, in ps.
   * @param temperature heat-bath temperature, in K.
   * @param friction    friction coefficient, in 1/ps.
   * @return native integrator handle owned by the new wrapper.
   *
   * <p>The native create function is called with (temperature, friction, dt).</p>
   */
  private static MemorySegment create(double dt, double temperature, double friction) {
    OpenMMRuntime.initialize();
    return OpenMMNative.OpenMM_LangevinMiddleIntegrator_create(temperature, friction, dt);
  }
}
